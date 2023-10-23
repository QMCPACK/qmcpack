#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 4 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/shared_main.hpp"

#include <boost/multi_array.hpp>

#include<iostream>
#include<vector>
#include<cstdlib>

void mySum(void *in_, void *out_, int *len, MPI_Datatype *dtype) {
  double* in = reinterpret_cast<double*>(in_);
  double* out = reinterpret_cast<double*>(out_);
  MPI_Aint sz,lb;
  MPI_Type_get_extent(*dtype,&lb,&sz);
  for(int i=0, iend=sz/sizeof(double); i<iend; i++)
    out[i] += in[i];
}

namespace mpi3 = boost::mpi3;

using std::cout;
using std::endl;
using boost::extents;

// run with 4 procs
int mpi3::main(int argc,char* argv[], mpi3::shared_communicator& node){

//  MPI_Init(&argc,&argv);;

//  int rank, size;
//  MPI_Comm_size(MPI_COMM_WORLD, &size);
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nw=16;
  int ni=2*624*624;
  double t1,t2;

  boost::multi_array<std::complex<double>,2> A(extents[ni][nw*node.size()]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> B(extents[ni][nw]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> C(extents[ni][nw]);  // 2x2 matrix

  std::vector<int> nr{10000,100000,250000,500000,ni};

  for(std::size_t n=0; n != nr.size(); n++) {
	ni = nr[n];
	mpi3::type tA = mpi3::double_.vector(ni, 2*nw, 2*A.strides()[0]);
//  MPI_Datatype typeA;
//  MPI_Type_vector(ni, 2*nw, 2*A.strides()[0], MPI_DOUBLE, &typeA);
	tA.commit();
//  MPI_Type_commit(&typeA);
//  MPI_Op myOpSum;
//  MPI_Op_create(mySum,true,&myOpSum);
//  mpi3::operation my_sum = 
	mpi3::commutative_operation my_sum(&mySum);
//  MPI_Allreduce(MPI_IN_PLACE,A.data(),2*ni*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	node.all_reduce_in_place_n(A.data(), A.num_elements(), std::plus<>{});
  
	int nt = 4;
	node.barrier();
//  MPI_Barrier(MPI_COMM_WORLD);
	{
		boost::timer::auto_cpu_timer t("All reduce");
//  t1=getTime();
		for(int i = 0; i != nt; i++)
			node.all_reduce_in_place_n(A.data(), A.num_elements(), std::plus<>{});
///			MPI_Allreduce(MPI_IN_PLACE,A.data(),2*ni*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	}
//  t2=getTime(); 
//  double tv1=(t2-t1)/nt;
	{
		for(int r=0, p=0; r != node.size(); r++, p+=nw) 
//		if(rank==r) node.reduce_in_place_n(A.data() + p, 1, my_sum, 
			MPI_Reduce(MPI_IN_PLACE, A.data()+p, 1, tA.impl_, myOpSum, r, node.impl_);
		else
			MPI_Reduce(A.data()+p, NULL, 1, tA.impl_, myOpSum, r, MPI_COMM_WORLD);
	}
	node.barrier();
	{
		boost::timer::auto_cpu_timer t;
		for(int i=0; i != nt; i++)
			for(int r=0, p=0; r != size; r++, p+=nw)
				if(rank==r)
					MPI_Reduce(MPI_IN_PLACE, A.data()+p, 1, A.impl_, myOpSum, r, node.impl_);
				else
					MPI_Reduce(A.data()+p, NULL, 1, typeA, myOpSum, r, node.impl_);
	}
	std::vector<mpi3::request> req(node.size()); 
//	std::vector<mpi3::status>  st(node.size()); 
	for(int r=0; r != node.size(); r++)
		if(node.rank()==r)
			MPI_Ireduce(C.data(), B.data(),2*ni*nw,MPI_DOUBLE,MPI_SUM,r, node.impl_, &req[r].impl_);
		else
			MPI_Ireduce(C.data(), B.data(),2*ni*nw,MPI_DOUBLE,MPI_SUM,r, node.impl_, &req[r].impl_);
	mpi3::wait_all(req);
	node.barrier();

	{
		boost::timer::auto_cpu_timer t("multiple Reduce");
		for(int i=0; i!=nt; i++){
			for(int r=0; r != node.size(); r++)
				if(node.rank() == r)
					MPI_Ireduce(C.data(), B.data(), 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_, &req[r].impl_);
				else
					MPI_Ireduce(C.data(), B.data(), 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_, &req[r].impl_);
				mpi3::wait_all(req);
		}
	}
	
	for(int r=0; r != node.size(); r++)
		if(node.rank()==r)
			MPI_Reduce(MPI_IN_PLACE, B.data(), 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_);
		else
			MPI_Reduce(B.data(), NULL, 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_);

	node.barrier();

	{
		boost::timer::auto_cpu_timer t("multiple Ireduce");
		for(int i=0; i != nt; i++) 
			for(int r=0; r<size; r++)
				if(rank==r)
					MPI_Reduce(MPI_IN_PLACE, B.data(), 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_);
				else
					MPI_Reduce(B.data(), NULL, 2*ni*nw, MPI_DOUBLE, MPI_SUM, r, node.impl_);
	}
  if(rank==0) cout<<" n, times (Allreduce, multiple Reduce, multiple IReduce): " <<ni <<" " <<tv1 <<"  " <<tv2 <<" " <<tv3 <<" " <<tv4  <<endl;
  MPI_Type_free(&typeA);

  }

  MPI_Finalize();

  return 0;
}

#include<iostream>
#include<vector>
#include<cstdlib>
#include<mpi.h>
#include <sys/time.h>
#include <ctime>
#include <boost/multi_array.hpp>

double getTime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
}

void mySum(void *in_, void *out_, int *len, MPI_Datatype *dtype) {
  double* in = reinterpret_cast<double*>(in_);
  double* out = reinterpret_cast<double*>(out_);
  MPI_Aint sz,lb;
  MPI_Type_get_extent(*dtype,&lb,&sz);
  for(int i=0, iend=sz/sizeof(double); i<iend; i++)
    out[i] += in[i];
}

using std::cout;
using std::cin;
using std::endl;
using boost::extents;

// run with 4 procs
int main(int argc,char* argv[])
{

  MPI_Init(&argc,&argv);;

  int global_rank, global_size;
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);

  MPI_Info info;
  MPI_Comm MPI_COMM_NODE_LOCAL;
  MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,global_rank,info,&MPI_COMM_NODE_LOCAL);

  int local_rank,local_size;
  MPI_Comm_size(MPI_COMM_NODE_LOCAL, &local_size);
  MPI_Comm_rank(MPI_COMM_NODE_LOCAL, &local_rank);

  MPI_Comm MPI_COMM_GROUP;
  MPI_Comm_split(MPI_COMM_WORLD,local_rank,global_rank,&MPI_COMM_GROUP);

  int rank,size;
  MPI_Comm_size(MPI_COMM_GROUP, &size);
  MPI_Comm_rank(MPI_COMM_GROUP, &rank);

  if(size*local_size != global_size) {
    std::cerr<<" Proc split error. \n";
    exit(1);
  }

  int nw=16;
  int ni=10000/local_size;
  double t1,t2;
  int blk=ni;
  if(argc>1) blk=atoi(argv[1]);
  std::vector<int> nr{1,10,25,50,75};

  boost::multi_array<std::complex<double>,2> A(extents[ni*75][nw*size]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> A_(extents[ni*75][nw*size]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> B(extents[ni*75][nw]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> C(extents[ni*75][nw]);  // 2x2 matrix


  for(int n=0; n<nr.size(); n++) {

  int ni_ = ni*nr[n];
  if(blk > ni_) continue; 
  int nblk = ni_/blk;
  ni_ = nblk*blk;

  MPI_Allreduce(A_.data(),A.data(),2*ni_*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_GROUP);

  int nt=4;
  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++)
    MPI_Allreduce(A_.data(),A.data(),2*ni_*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_GROUP);
  MPI_Barrier(MPI_COMM_WORLD);
  t2=getTime();
  double tv1=(t2-t1)/nt;

  int nmmax=4096;
  std::vector<MPI_Request> req(nmmax);
  std::vector<MPI_Status> st(nmmax);  
  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++) {
    int cnt=0;
    for(int b=0; b<nblk; b++) {
      for(int r=0; r<size; r++) {
        if(rank==r)
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_GROUP,&req[cnt++]);
        else
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_GROUP,&req[cnt++]);
      }
      if(cnt==nmmax) {
        MPI_Waitall(cnt,req.data(),st.data());
        cnt=0;
      }
    }
    if(cnt>0) {
      MPI_Waitall(cnt,req.data(),st.data());
      cnt=0;
    }  
  }
  MPI_Barrier(MPI_COMM_WORLD);
  t2=getTime();
  double tv5=(t2-t1)/nt;

  if(global_rank==0) cout<<" n, times (Allreduce, multiple IReduce): " <<ni_*local_size <<" " <<ni_ <<" " <<tv1 <<" " <<tv5 <<endl;

  }

  MPI_Finalize();

  return 0;
}

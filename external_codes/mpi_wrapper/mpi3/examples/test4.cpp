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

  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nw=16;
  int ni=2*624*624;
  double t1,t2;
  int blk=10000;
  if(argc>1) blk=atoi(argv[1]);

  boost::multi_array<std::complex<double>,2> A(extents[ni][nw*size]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> B(extents[ni][nw]);  // 2x2 matrix
  boost::multi_array<std::complex<double>,2> C(extents[ni][nw]);  // 2x2 matrix

  std::vector<int> nr{10000,100000,250000,500000,ni};

  for(int n=0; n<nr.size(); n++) {

  ni = nr[n];
  if(blk > ni) continue; 
  int nblk = ni/blk;
  ni = nblk*blk;

  MPI_Allreduce(MPI_IN_PLACE,A.data(),2*ni*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  int nt=4;
  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++)
    MPI_Allreduce(MPI_IN_PLACE,A.data(),2*ni*nw*size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  t2=getTime();
  double tv1=(t2-t1)/nt;

  std::vector<MPI_Request> req(size*nblk); 
  std::vector<MPI_Status> st(size*nblk); 
  for(int r=0; r<size; r++) {
    if(rank==r)
      for(int b=0; b<nblk; b++)
        MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
    else
      for(int b=0; b<nblk; b++)
        MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
  }
  MPI_Waitall(size*nblk,req.data(),st.data());

  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++) {
    for(int r=0; r<size; r++) {
      if(rank==r)
        for(int b=0; b<nblk; b++)
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
      else
        for(int b=0; b<nblk; b++)
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
    }
    MPI_Waitall(size*nblk,req.data(),st.data());
  }
  t2=getTime();
  double tv3=(t2-t1)/nt;


  for(int b=0; b<nblk; b++)
    for(int r=0; r<size; r++) {
      if(rank==r)
        MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
      else
        MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
    }
  MPI_Waitall(size*nblk,req.data(),st.data());

  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++) {
    for(int b=0; b<nblk; b++)
      for(int r=0; r<size; r++) {
        if(rank==r)
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
        else
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[r*nblk+b]);
      }
    MPI_Waitall(size*nblk,req.data(),st.data());
  }
  t2=getTime();
  double tv4=(t2-t1)/nt;

  int nmmax=4096;
  req.resize(nmmax);
  st.resize(nmmax);
  MPI_Barrier(MPI_COMM_WORLD);
  t1=getTime();
  for(int i=0; i<nt; i++) {
    int cnt=0;
    for(int b=0; b<nblk; b++) {
      for(int r=0; r<size; r++) {
        if(rank==r)
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[cnt++]);
        else
          MPI_Ireduce(C.data()+b*blk,B.data()+b*blk,2*blk*nw,MPI_DOUBLE,MPI_SUM,r,MPI_COMM_WORLD,&req[cnt++]);
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
  t2=getTime();
  double tv5=(t2-t1)/nt;

  if(rank==0) cout<<" n, times (Allreduce, multiple IReduce 1, multiple IReduce 2): " <<ni <<" " <<tv1 <<" " <<tv3 <<" " <<tv4 <<" " <<tv5 <<endl;

  }

  MPI_Finalize();

  return 0;
}

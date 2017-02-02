#include<iostream>
#include<cstdlib>
#include <mpi.h>

using namespace std;

#include "sys/sysinfo.h"

inline size_t freemem()
{
  struct sysinfo si;
  sysinfo(&si);
  si.freeram+=si.bufferram;
  return si.freeram>>20;
}

#include "dv.h"

int main(int argc, char* argv[])
{
  int rank, nproc;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  qmcplusplus::SMDenseVector<double> buff; 

  if(rank==0) std::cout<<"Memory: " <<freemem() <<std::endl;
  buff.setup(rank==0,std::string("Buffer_0"),nproc,MPI_COMM_WORLD);
  buff.resize(1000000);
  if(rank==0) std::cout<<"Memory: " <<freemem() <<std::endl;
  std::cout<<"rank: " <<rank <<" " <<nproc <<" " <<buff.size() <<std::endl;

  for(int i=0; i<5; i++)
   *(buff.values()+rank*5+i) = rank; 

  if(rank==0) 
   for(int i=0; i<10; i++)
    std::cout<<"  " <<i <<" " <<*(buff.values()+i) <<std::endl; 

  buff.resize(100000000);
  if(rank==0) std::cout<<"Memory: " <<freemem() <<std::endl;
  std::cout<<"rank: " <<rank <<" " <<nproc <<" " <<buff.size() <<std::endl;

  for(int i=0; i<10; i++)
   *(buff.values()+rank*10+i) = rank+10;

  if(rank==0)
   for(int i=0; i<20; i++)
    std::cout<<"  " <<i <<" " <<*(buff.values()+i) <<std::endl;

  buff.resize(200000000);
  if(rank==0) std::cout<<"Memory: " <<freemem() <<std::endl;

  buff.resize(20,true);
  if(rank==0) std::cout<<"Memory: " <<freemem() <<std::endl;


  MPI_Finalize();

  return 0;
}

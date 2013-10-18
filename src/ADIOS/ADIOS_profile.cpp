
#include "ADIOS/ADIOS_profile.h"

namespace ADIOS_PROFILE
{

void profile_adios_size(Communicate* myComm, OUTPUT_T op, uint64_t adios_groupsize, uint64_t adios_totalsize)
{
#ifdef IO_PROFILE
  for(int i=0; i<myComm->size(); i++)
  {
    MPI_Barrier(myComm->getMPI());
    if(myComm->rank()==i)
    {
      if(op==TRACES)
      {
        cout<<myComm->rank()<<" profile adios trace write size is "<<adios_groupsize<<" "<<adios_totalsize<<endl;
      }
      else if(op == CKPOINT)
      {
        cout<<myComm->rank()<<" profile adios checkpoint write size is "<<adios_groupsize<<" "<<adios_totalsize<<endl;
      }
    }
  }
#endif
}

void profile_adios_init(int nBlock)
{
#ifdef IO_PROFILE
  int size = nBlock*sizeof(double);
  comp_times = (double *)malloc(size);
  trace_times= (double *)malloc(size);
  checkpoint_times= (double *)malloc(size);
  bzero(comp_times, size);
  bzero(trace_times, size);
  bzero(checkpoint_times, size);
  start = MPI_Wtime();
#endif
}

void profile_adios_finalize(Communicate* myComm, int nBlocks)
{
#ifdef IO_PROFILE
  double end_time = MPI_Wtime() - start;
  for(int i=0; i<myComm->size(); i++)
  {
    MPI_Barrier(myComm->getMPI());
    if(i==myComm->rank())
    {
      cout<<myComm->rank()<<" profile"<<endl;
      double total_comp_time = 0.0;
      double total_trace_time = 0.0;
      double total_checkpoint_time = 0.0;
      for(int j=0; j<nBlocks; j++)
      {
        printf("comp time %f trace time %f checkpoint_time %f\n", comp_times[j], trace_times[j], checkpoint_times[j]);
        //cout<<"comp time "<<comp_times[j]<<" trace time "<<trace_times[j]<<" checkpoint_times "<<checkpoint_times[j]<<endl;
        total_comp_time += comp_times[j];
        total_trace_time += trace_times[j];
        total_checkpoint_time += checkpoint_times[j];
      }
      printf("total time is %f comp time %f trace time %f checkpoint time %f rank %d\n", end_time, total_comp_time, total_trace_time, total_checkpoint_time, myComm->rank());
      //cout<<"total time is "<<end_time<<" comp time "<<total_comp_time<<" trace time "<<total_trace_time<<" checkpoint time "<<total_checkpoint_time<<" rank "<<myComm->rank()<<endl;
    }
    MPI_Barrier(myComm->getMPI());
  }
  free(comp_times);
  free(trace_times);
  free(checkpoint_times);
#endif
}

void profile_adios_start_comp(int block)
{
#ifdef IO_PROFILE
  comp_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_comp(int block)
{
#ifdef IO_PROFILE
  double tmp = comp_times[block];
  comp_times[block]=MPI_Wtime() - tmp;
#endif
}

void profile_adios_start_trace(int block)
{
#ifdef IO_PROFILE
  trace_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_trace(int block)
{
#ifdef IO_PROFILE
  double tmp = trace_times[block];
  trace_times[block]=MPI_Wtime() - tmp;
#endif
}

void profile_adios_start_checkpoint(int block)
{
#ifdef IO_PROFILE
  checkpoint_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_checkpoint(int block)
{
#ifdef IO_PROFILE
  double tmp = checkpoint_times[block];
  checkpoint_times[block]=MPI_Wtime() - tmp;
#endif
}

};

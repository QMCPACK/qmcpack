
#include "ADIOS/ADIOS_profile.h"

namespace ADIOS_PROFILE
{

void profile_adios_size(Communicate* myComm, OUTPUT_T op, uint64_t adios_groupsize, uint64_t adios_totalsize)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
	if(op == TRACES)
	{
		trace_data_grp[trace_index]=adios_groupsize;
		trace_data_total[trace_index]=adios_totalsize;
		trace_index++;
	}
	else if(op==CKPOINT)
	{
		ckp_data_grp[ckp_index]=adios_groupsize;
		ckp_data_total[ckp_index]=adios_totalsize;
		ckp_index++;
	}
#endif
}

void profile_adios_init(int nBlock)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  int size = nBlock*sizeof(double);
  comp_times = (double *)malloc(size);
  trace_times= (double *)malloc(size);
  checkpoint_times= (double *)malloc(size);
  bzero(comp_times, size);
  bzero(trace_times, size);
  bzero(checkpoint_times, size);
	int data_size = nBlock*sizeof(int);
	ckp_data_grp = (int *)malloc(data_size);
	trace_data_grp = (int *)malloc(data_size);
	ckp_data_total = (int *)malloc(data_size);
	trace_data_total = (int *)malloc(data_size);
	bzero(ckp_data_grp, data_size);
	bzero(trace_data_grp, data_size);
	bzero(ckp_data_total, data_size);
	bzero(trace_data_total, data_size);
	trace_index = 0;
	ckp_index = 0;
  start = MPI_Wtime();
#endif
}

void profile_adios_finalize(Communicate* myComm, int nBlocks)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
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
				if(j<trace_index)
					cout<<myComm->rank()<<" profile adios trace write size is "<<trace_data_grp[j]<<" "<<trace_data_total[j]<<endl;
				if(j<ckp_index)
					cout<<myComm->rank()<<" profile adios checkpoint write size is "<<ckp_data_grp[j]<<" "<<ckp_data_total[j]<<endl;
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
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  comp_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_comp(int block)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  double tmp = comp_times[block];
  comp_times[block]=MPI_Wtime() - tmp;
#endif
}

void profile_adios_start_trace(int block)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  trace_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_trace(int block)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  double tmp = trace_times[block];
  trace_times[block]=MPI_Wtime() - tmp;
#endif
}

void profile_adios_start_checkpoint(int block)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  checkpoint_times[block]=MPI_Wtime();
#endif
}

void profile_adios_end_checkpoint(int block)
{
#if defined(HAVE_ADIOS) && defined(IO_PROFILE)
  double tmp = checkpoint_times[block];
  checkpoint_times[block]=MPI_Wtime() - tmp;
#endif
}

};



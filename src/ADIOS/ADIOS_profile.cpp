//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, Oak Ridge National Laboratory
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "ADIOS/ADIOS_profile.h"

#if (defined HAVE_ADIOS) && (defined IO_PROFILE)
namespace ADIOS_PROFILE
{

void profile_adios_size(Communicate* myComm, OUTPUT_T op, uint64_t adios_groupsize, uint64_t adios_totalsize)
{
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
}

void profile_adios_init(int nBlock)
{
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
}

void profile_adios_finalize(Communicate* myComm, int nBlocks)
{
  double end_time = MPI_Wtime() - start;
  for(int i=0; i<myComm->size(); i++)
  {
    MPI_Barrier(myComm->getMPI());
    if(i==myComm->rank())
    {
      std::cout <<myComm->rank()<<" profile"<< std::endl;
      double total_comp_time = 0.0;
      double total_trace_time = 0.0;
      double total_checkpoint_time = 0.0;
      for(int j=0; j<nBlocks; j++)
      {
        printf("comp time %f trace time %f checkpoint_time %f\n", comp_times[j], trace_times[j], checkpoint_times[j]);
        //cout<<"comp time "<<comp_times[j]<<" trace time "<<trace_times[j]<<" checkpoint_times "<<checkpoint_times[j]<< std::endl;
        total_comp_time += comp_times[j];
        total_trace_time += trace_times[j];
        total_checkpoint_time += checkpoint_times[j];
				if(j<trace_index)
					cout<<myComm->rank()<<" profile adios trace write size is "<<trace_data_grp[j]<<" "<<trace_data_total[j]<< std::endl;
				if(j<ckp_index)
					cout<<myComm->rank()<<" profile adios checkpoint write size is "<<ckp_data_grp[j]<<" "<<ckp_data_total[j]<< std::endl;
      }
      printf("total time is %f comp time %f trace time %f checkpoint time %f rank %d\n", end_time, total_comp_time, total_trace_time, total_checkpoint_time, myComm->rank());
      //cout<<"total time is "<<end_time<<" comp time "<<total_comp_time<<" trace time "<<total_trace_time<<" checkpoint time "<<total_checkpoint_time<<" rank "<<myComm->rank()<< std::endl;
    }
    MPI_Barrier(myComm->getMPI());
  }
  free(comp_times);
  free(trace_times);
  free(checkpoint_times);
}

void profile_adios_start_comp(int block)
{
  comp_times[block]=MPI_Wtime();
}

void profile_adios_end_comp(int block)
{
  double tmp = comp_times[block];
  comp_times[block]=MPI_Wtime() - tmp;
}

void profile_adios_start_trace(int block)
{
  trace_times[block]=MPI_Wtime();
}

void profile_adios_end_trace(int block)
{
  double tmp = trace_times[block];
  trace_times[block]=MPI_Wtime() - tmp;
}

void profile_adios_start_checkpoint(int block)
{
  checkpoint_times[block]=MPI_Wtime();
}

void profile_adios_end_checkpoint(int block)
{
  double tmp = checkpoint_times[block];
  checkpoint_times[block]=MPI_Wtime() - tmp;
}

void updateBlock(int block){
  ADIOS_PROFILE::block = block;
}


void profile_init(int rank){
  myrank = rank;
  comp_start = 0.0;
  comp_end = 0.0;
  comm_start = 0.0;
  comm_end = 0.0;
  io_open_start = 0.0;
  io_open_end = 0.0;
  io_group_start = 0.0;
  io_group_end = 0.0;
  io_write_start = 0.0;
  io_write_end = 0.0;
  io_close_start = 0.0;
  io_close_end = 0.0;
  io_start = 0.0;
  io_end = 0.0;
  comp_total = 0.0;
  comm_total = 0.0;
  io_total = 0.0;
  times.clear();
}

void profile_final(){
  //char buf1[20];
  //sprintf(buf1, "output%d", myrank);

  //FILE * fp;
  //fp = fopen (buf1, "a");
  //for(int i=0; i<times.size(); i++){
    //char buf[200];
    //bzero(buf, 200);
    //sprintf(buf, "%f %d  %d\n", times[i].time, times[i].t_attr, times[i].block);
    //fputs(buf, fp);
  //}
  //fflush(fp);
  //fclose(fp);
                 

  //ofstream myfile;
  //myfile.open(buf);
  //for(int i=0; i<times.size(); i++){
    //myfile <<times[i].time<<"\t"<<times[i].t_attr<<"\t"<<times[i].block<<"\t"<<times[i].step<< std::endl;
  //}
  //myfile.close();
}

void comp_s(){
  comp_start = MPI_Wtime();
}

void comp_e(){
  comp_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" comp "<<comp_end-comp_start<< std::endl;
  //TIME_INFO t;
  //t.time = comp_end - comp_start;
  //t.t_attr = COMP;
  //t.block = block;
  //times.push_back(t);
  //comp_total += t.time;
}

void comm_s(){
  comm_start = MPI_Wtime();
}

void comm_e(){
  comm_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" comm "<<comm_end-comm_start<< std::endl;
  //TIME_INFO t;
  //t.time = comm_end - comm_start;
  //t.t_attr = COMM;
  //t.block = block;
  //times.push_back(t);
  //comm_total += t.time;
}

void io_open_s(){
  io_open_start = MPI_Wtime();
}

void io_open_e(){
  io_open_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" open "<<io_open_end-io_open_start<< std::endl;
  //TIME_INFO t;
  //t.time = io_open_end - io_open_start;
  //t.t_attr = IO_OPEN;
  //t.block = block;
  //times.push_back(t);
}

void io_group_s(){
  io_group_start = MPI_Wtime();
}

void io_group_e(){
  io_group_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" group "<<io_group_end-io_group_start<< std::endl;
  //TIME_INFO t;
  //t.time = io_group_end - io_group_start;
  //t.t_attr = IO_GROUP;
  //t.block = block;
  //times.push_back(t);
}

void io_write_s(){
  io_write_start = MPI_Wtime();
}

void io_write_e(){
  io_write_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" write "<<io_write_end-io_write_start<< std::endl;
  //TIME_INFO t;
  //t.time = io_write_end - io_write_start;
  //t.t_attr = IO_WRITE;
  //t.block = block;
  //times.push_back(t);
}

void io_close_s(){
  io_close_start = MPI_Wtime();
}

void io_close_e(){
  io_close_end = MPI_Wtime();
  qmcplusplus::app_log()<<block<<" close "<<io_close_end-io_close_start<< std::endl;
  //TIME_INFO t;
  //t.time = io_close_end - io_close_start;
  //t.t_attr = IO_CLOSE;
  //t.block = block;
  //times.push_back(t);
}

void io_s(){
  io_start = MPI_Wtime();
}

void io_e(){
  io_end = MPI_Wtime();
  //TIME_INFO t;
  //t.time = io_end - io_start;
  //t.t_attr = IO;
  //t.block = block;
  //times.push_back(t);
  //io_total += t.time;
}
};
#endif



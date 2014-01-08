
#ifndef ADIOS_ADIOS_PROFILE_H
#define ADIOS_ADIOS_PROFILE_H

#include <string>
#include <Configuration.h>
#include <stdint.h>
#include <strings.h>

namespace ADIOS_PROFILE
{

typedef enum
{
  TRACES,
  CKPOINT
}	OUTPUT_T;

static double start = 0.0;
static double *comp_times;
static double *trace_times;
static double *checkpoint_times;
static int *ckp_data_grp;
static int *trace_data_grp;
static int *ckp_data_total;
static int *trace_data_total;
static int ckp_index;
static int trace_index;

#if (defined HAVE_ADIOS) && (defined IO_PROFILE)
void profile_adios_size(Communicate* myComm, OUTPUT_T op, uint64_t adios_groupsize, uint64_t adios_totalsize);
void profile_adios_init(int nBlock);
void profile_adios_finalize(Communicate* myComm, int nBlock);
void profile_adios_start_comp(int block);
void profile_adios_start_trace(int block);
void profile_adios_start_checkpoint(int block);
void profile_adios_end_comp(int block);
void profile_adios_end_trace(int block);
void profile_adios_end_checkpoint(int block);
#else
inline void profile_adios_size(Communicate* myComm, OUTPUT_T op, uint64_t adios_groupsize, uint64_t adios_totalsize){}
inline void profile_adios_init(int nBlock){}
inline void profile_adios_finalize(Communicate* myComm, int nBlock){}
inline void profile_adios_start_comp(int block){}
inline void profile_adios_start_trace(int block){}
inline void profile_adios_start_checkpoint(int block){}
inline void profile_adios_end_comp(int block){}
inline void profile_adios_end_trace(int block){}
inline void profile_adios_end_checkpoint(int block){}
#endif
}

#endif

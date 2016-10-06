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
    
    


#if defined(ADIOS_VERIFY) && defined(HAVE_ADIOS)
#include "ADIOS/ADIOS_verify.h"

namespace IO_VERIFY
{
void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, unsigned long origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  int size;
  vi = adios_inq_var(fp, name);
  if (vi->ndim > 0)
  {
    std::cout <<name<<" verification not passed, not a scalar"<< std::endl;
    return;
  }
  size = count*adios_type_size(vi->type, vi->value);
  unsigned long* mem= (unsigned long* )malloc(size);
  adios_schedule_read(fp, NULL, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  if(mem[0] == origin)
  {
    std::cout <<name<<" verification passed"<<mem[0]<< std::endl;
  }
  else
  {
    std::cout <<name<<" verification not passed, readin: "<<mem[0]<<" writeout: "<<origin<< std::endl;
  }
  adios_free_varinfo (vi);
  free(mem);
}

void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, int* origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  int size;
  vi = adios_inq_var(fp, name);
  if (vi->ndim > 0)
  {
    std::cout <<name<<" verification not passed, not a scalar"<< std::endl;
    return;
  }
  size = count*adios_type_size(vi->type, vi->value);
  int * mem= (int * )malloc(size);
  ADIOS_SELECTION *sel = adios_selection_writeblock(OHMMS::Controller->rank());
  adios_schedule_read(fp, sel, name, 0, 1, mem);
  //adios_schedule_read(fp, NULL, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  if(mem[0] == *origin)
  {
    std::cout <<name<<" verification passed "<<mem[0]<< std::endl;
  }
  else
  {
    std::cout <<name<<" verification not passed, readin: "<<mem[0]<<" writeout: "<<*origin<< std::endl;
  }
  adios_free_varinfo (vi);
  adios_selection_delete(sel);
  free(mem);
}


void adios_checkpoint_verify_intarray_variables(ADIOS_FILE* fp, const char* name, int* origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  int size;
  vi = adios_inq_var(fp, name);
  adios_inq_var_blockinfo (fp, vi);
  if (vi->ndim > 0)
  {
    count*=vi->dims[0];
    for (int j = 1; j < vi->ndim; j++)
    {
      count *= vi->dims[j];
    }
  }
  size = count*adios_type_size(vi->type, vi->value);
  int *mem= (int *)malloc(size);
  adios_schedule_read(fp, NULL, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  for(int i=0; i<count; i++)
  {
    if(mem[i] == origin[i])
    {
      std::cout <<name<<"["<<i<<"]verification passed "<<mem[i]<< std::endl;
    }
    else
    {
      std::cout <<name<<"["<<i<<"]verification not passed, readin: "<<mem[i]<<" writeout: "<<origin[i]<< std::endl;
    }
  }
  adios_free_varinfo (vi);
  free(mem);
}


void adios_checkpoint_verify_variables(ADIOS_FILE* fp, const char* name, RealType* origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  int size;
  vi = adios_inq_var(fp, name);
  adios_inq_var_blockinfo (fp, vi);
  if (vi->ndim > 0)
  {
    count*=vi->dims[0];
    for (int j = 1; j < vi->ndim; j++)
    {
      count *= vi->dims[j];
    }
  }
  size = count*adios_type_size(vi->type, vi->value);
  RealType *mem= (RealType *)malloc(size);
  adios_schedule_read(fp, NULL, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  for(int i=0; i<count; i++)
  {
    if(mem[i] == origin[i])
    {
      std::cout <<name<<"["<<i<<"]verification passed "<<mem[i]<< std::endl;
    }
    else
    {
      std::cout <<name<<"["<<i<<"]verification not passed, readin: "<<mem[i]<<" writeout: "<<origin[i]<< std::endl;
    }
  }
  adios_free_varinfo (vi);
  free(mem);
}

void adios_checkpoint_verify_random_variables(ADIOS_FILE* fp, const char* name, uint_type* origin)
{
  ADIOS_VARINFO *vi;
  int count_int = 1;
  int size;
  uint64_t *start;
  uint64_t *count;
  vi = adios_inq_var(fp, name);
  adios_inq_var_blockinfo (fp, vi);
  if (vi->ndim > 0)
  {
    start = (uint64_t *)malloc(vi->ndim * sizeof(uint64_t));
    count = (uint64_t *)malloc(vi->ndim * sizeof(uint64_t));
  }
  for (int j=0; j<vi->nblocks[0]; j++)
  {
    if(j == OHMMS::Controller->rank())
    {
      for (int k=0; k<vi->ndim; k++)
      {
        start[k] = vi->blockinfo[j].start[k];
        count[k] = vi->blockinfo[j].count[k];
        count_int *= count[k];
        //cout<<OHMMS::Controller->rank()<<" count "<<start[k]<<" "<<count[k]<< std::endl;
      }
    }
  }
  size = count_int*adios_type_size(vi->type, vi->value);
  uint_type *mem= (uint_type*)malloc(size);
  ADIOS_SELECTION *sel = adios_selection_boundingbox(vi->ndim, start, count);
  adios_schedule_read(fp, sel, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  int flag = 0;
  for(int i=0; i<count_int; i++)
  {
    if(mem[i] == origin[i])
    {
      //cout<<name<<"["<<i<<"]verification passed "<<mem[i]<< std::endl;
    }
    else
    {
      flag = 1;
      std::cout <<name<<"["<<i<<"]verification not passed, readin: "<<mem[i]<<" writeout: "<<origin[i]<< std::endl;
    }
  }
  if (flag == 0) std::cout <<name<<" verification passed "<< std::endl;
  else std::cout <<name<<" verification not passed "<< std::endl;
  adios_free_varinfo (vi);
  adios_selection_delete(sel);
  free(start);
  free(count);
  free(mem);
}


void adios_checkpoint_verify_local_variables(ADIOS_FILE* fp, const char* name, OHMMS_PRECISION* origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  unsigned long size = 1;
  vi = adios_inq_var(fp, name);
  adios_inq_var_blockinfo(fp, vi);
  for (int j=0; j<vi->nblocks[0]; j++)
  {
    if(OHMMS::Controller->rank() == j)
    {
      for (int k=0; k<vi->ndim; k++)
      {
        count *= vi->blockinfo[j].count[k];
      }
    }
  }
  size = count * adios_type_size(vi->type, vi->value);
  OHMMS_PRECISION* mem = (OHMMS_PRECISION*)malloc(size);
  ADIOS_SELECTION *sel = adios_selection_writeblock(OHMMS::Controller->rank());
  adios_schedule_read(fp, sel, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  int flag = 0;
  for(int i=0; i<count; i++)
  {
    if(mem[i] == origin[i])
    {
    }
    else
    {
      flag = 1;
      std::cout <<OHMMS::Controller->rank()<<" "<<name<<"["<<i<<"]verification not passed, readin: "<<mem[i]<<" writeout: "<<origin[i]<< std::endl;
    }
  }
  if(flag == 0) std::cout <<OHMMS::Controller->rank()<<" "<<name<<" verification passed"<< std::endl;
  adios_selection_delete(sel);
  adios_free_varinfo (vi);
  free(mem);
}

void adios_trace_verify_local_variables(ADIOS_FILE* fp, const char* name, double* origin)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  unsigned long size = 1;
  vi = adios_inq_var(fp, name);
  adios_inq_var_blockinfo(fp, vi);
  for (int j=0; j<vi->nblocks[0]; j++)
  {
    if(OHMMS::Controller->rank() == j)
    {
      for (int k=0; k<vi->ndim; k++)
      {
        count *= vi->blockinfo[j].count[k];
      }
    }
  }
  size = count * adios_type_size(vi->type, vi->value);
  double* mem = (double*)malloc(size);
  ADIOS_SELECTION *sel = adios_selection_writeblock(OHMMS::Controller->rank());
  adios_schedule_read(fp, sel, name, 0, 1, mem);
  adios_perform_reads(fp, 1);
  int flag = 0;
  for(int i=0; i<count; i++)
  {
    if(mem[i] == origin[i])
    {
    }
    else
    {
      flag = 1;
      std::cout <<OHMMS::Controller->rank()<<" "<<name<<"["<<i<<"]verification not passed, readin: "<<mem[i]<<" writeout: "<<origin[i]<< std::endl;
    }
  }
  if(flag == 0) std::cout <<OHMMS::Controller->rank()<<" "<<name<<" verification passed"<< std::endl;
  adios_selection_delete(sel);
  adios_free_varinfo (vi);
  free(mem);
}

};
#endif

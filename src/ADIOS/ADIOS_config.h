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
    
    


#ifndef ADIOS_ADIOS_CONFIG_H
#define ADIOS_ADIOS_CONFIG_H

#include <Configuration.h>
#include <stdint.h>
#include <cstring>
#ifdef HAVE_ADIOS
#include "Utilities/RandomGenerator.h"
#include <adios_read.h>
extern "C" {
#include <adios_error.h>
}

namespace ADIOS
{
void initialize(bool use_hdf5, bool use_adios);
void initialze(std::string &xml_filename, bool use_hdf5, bool use_adios);
void initialize(const char*value);
bool useADIOS();
bool useHDF5();
const std::string& get_adios_xml();

static std::string adiosname;
static ADIOS_FILE * openfp;

void set_adios_init(bool b);
bool get_adios_init();

bool getRdADIOS();
bool getRdHDF5();
bool getFirstOpen();
void setFirstOpen(bool a);
string getTraceFileName();
typedef qmcplusplus::RandomGenerator_t::uint_type uint_type;

inline int open(const std::string &fname, MPI_Comm comm)
{
  int ret = 0;
  adiosname = fname;
  if(fname.find("config.bp")>= fname.size())
    adiosname.append(".config.bp");
  openfp = adios_read_open_file(adiosname.c_str(),
                                        ADIOS_READ_METHOD_BP,
                                        comm);
  if(openfp == NULL)
  {
    qmcplusplus::app_error() <<"Fail to open adios file "<<adiosname<<" Abort."<< std::endl;
  }
  else
  {
    qmcplusplus::app_log()<<adiosname<<" is open "<< std::endl;
    ret = 1;
  }
  return ret;
}

template <class T>
void read(T data, const std::string& aname)
{
  ADIOS_VARINFO *vi;
  int count = 1;
  int size;

  char * name = new char [aname.length()+1];
  std::strcpy (name, aname.c_str());
	if ( openfp == NULL)
	{
		qmcplusplus::app_error()<<"openfp is null "<< std::endl;
	}
  vi = adios_inq_var(openfp, name);
  adios_inq_var_blockinfo (openfp, vi);
  if (vi->ndim > 0)
  {
    count*=vi->dims[0];
    for (int j = 1; j < vi->ndim; j++)
    {
      count *= vi->dims[j];
    }
  }
  size = count*adios_type_size(vi->type, vi->value);
  adios_schedule_read(openfp, NULL, name, 0, 1, data);
  adios_perform_reads(openfp, 1);
  adios_free_varinfo (vi);
}

template<class T, class C>
void read_scalar(T& data, const std::string& aname, C& index)
{
  ADIOS_VARINFO *vi;
  int size;

  char *name = new char[aname.length()+1];
  std::strcpy(name, aname.c_str());
  if (openfp == NULL)
  {
    qmcplusplus::app_error()<<"openfp is null "<< std::endl;
  }
  vi = adios_inq_var(openfp, name);
  adios_inq_var_blockinfo(openfp, vi);
  index = vi->nblocks[0];
  size = index*adios_type_size(vi->type, vi->value);
  data = (int *)malloc(size);

  ADIOS_SELECTION *sel;

  for(int i=0; i<index; i++)
  {
    sel=adios_selection_writeblock(i);
    adios_schedule_read(openfp, sel, name, 0, 1, &data[i]);
  }
  adios_perform_reads(openfp, 1);
  adios_free_varinfo(vi);
  adios_selection_delete(sel);
}

template<class T>
void read_walkers(T& data, const std::string& aname)
{
  ADIOS_VARINFO *vi;
  int size;

  char *name = new char[aname.length()+1];
  std::strcpy(name, aname.c_str());
  if (openfp == NULL)
  {
    qmcplusplus::app_error()<<"openfp is null "<< std::endl;
  }
  vi = adios_inq_var(openfp, name);
  adios_inq_var_blockinfo(openfp, vi);
  ADIOS_SELECTION *sel;
  int index = 0;
  for(int i = 0; i<vi->nblocks[0];i++)
  {
   int start = vi->blockinfo[i].start[0];
   int count = vi->blockinfo[i].count[0];
   sel=adios_selection_writeblock(i);
   adios_schedule_read(openfp, sel, name, 0, 1, &data[index]);
   index += count;
  }
  adios_perform_reads(openfp, 1);
  adios_free_varinfo(vi);
  adios_selection_delete(sel);
}

template <class T, class S>
void read_random(T& data, S& shape, const std::string& aname)
{
  ADIOS_VARINFO *vi;
  int size;

  char * name = new char [aname.length()+1];
  std::strcpy (name, aname.c_str());
  if ( openfp == NULL)
  {
    qmcplusplus::app_error()<<"openfp is null "<< std::endl;
  }

	uint64_t DIM = 2;
  uint64_t * start= (uint64_t *)malloc(sizeof(uint64_t)*DIM);
  uint64_t * count= (uint64_t *)malloc(sizeof(uint64_t)*DIM);
	start[0] = 0;
	start[1] = 0;
	count[0] = 1;

  vi = adios_inq_var(openfp, name);
  if (vi->ndim != DIM)
  {
 		qmcplusplus::app_error() <<"random number dimension is not 2."<< std::endl; 
	}  
	
	ADIOS_SELECTION *sel;
	uint_type * data1 = (uint_type * ) malloc(sizeof(uint_type)*vi->dims[0]*vi->dims[1]);

	shape[0] = vi->dims[0];
	shape[1] = vi->dims[1];
	
	count[1] = vi->dims[1];
  sel = adios_selection_boundingbox(DIM, start, count);	
  adios_schedule_read(openfp, sel, name, 0, 1, data1);
  adios_perform_reads(openfp, 1);
	data.assign(data1, data1+vi->dims[0]*vi->dims[1]);
	free(data1);
  adios_free_varinfo (vi);
	adios_selection_delete(sel);
}

inline void close()
{
	int ret = adios_read_close(openfp);
	if(ret)
		qmcplusplus::app_error() <<"Fail to close adios file "<<adiosname<<" Abort."<< std::endl;
  else
    qmcplusplus::app_log()<<adiosname<<" is closed "<< std::endl;
}

};

#else
//empty inline functions to avoid compiler macros
namespace ADIOS
{
  inline bool useADIOS() { return false;}
  inline bool useHDF5() { return true;}
}
#endif

#endif

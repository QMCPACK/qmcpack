#ifndef ADIOS_ADIOS_CONFIG_H
#define ADIOS_ADIOS_CONFIG_H

#include <Configuration.h>
#include <stdint.h>
#include <cstring>
#ifdef HAVE_ADIOS
#include "Utilities/RandomGenerator.h"
#include <adios_read.h>

namespace ADIOS
{
void initialize(bool use_hdf5, bool use_adios);
void initialze(std::string &xml_filename, bool use_hdf5, bool use_adios);
void initialize(const char*value);
bool useADIOS();
bool useHDF5();
const std::string& get_adios_xml();

static string adiosname;
static ADIOS_FILE * openfp;

bool getRdADIOS();
bool getRdHDF5();
typedef qmcplusplus::RandomGenerator_t::uint_type uint_type;

inline int open(const string &fname, MPI_Comm comm)
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
    qmcplusplus::app_error() <<"Fail to open adios file "<<adiosname<<" Abort."<<endl;
  }
  else
  {
    qmcplusplus::app_log()<<adiosname<<" is open "<<endl;
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
		qmcplusplus::app_error()<<"openfp is null "<<endl;
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

template <class T, class S>
void read_random(T& data, S& shape, const std::string& aname)
{
  ADIOS_VARINFO *vi;
  int size;

  char * name = new char [aname.length()+1];
  std::strcpy (name, aname.c_str());
  if ( openfp == NULL)
  {
    qmcplusplus::app_error()<<"openfp is null "<<endl;
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
 		qmcplusplus::app_error() <<"random number dimension is not 2."<<endl; 
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
		qmcplusplus::app_error() <<"Fail to close adios file "<<adiosname<<" Abort."<<endl;
  else
    qmcplusplus::app_log()<<adiosname<<" is closed "<<endl;
}

};

#endif

#endif

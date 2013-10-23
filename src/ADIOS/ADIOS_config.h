
#ifndef ADIOS_ADIOS_CONFIG_H
#define ADIOS_ADIOS_CONFIG_H

#include <string>
#include <Configuration.h>
#include <cstring>
#ifdef HAVE_ADIOS
#include <adios_read.h>

namespace ADIOS
{
void initialize(bool use_hdf5, bool use_adios);
void initialze(std::string &xml_filename, bool use_hdf5, bool use_adios);
bool useADIOS();
bool useHDF5();
const std::string& get_adios_xml();

static string adiosname;
static ADIOS_FILE * openfp;

//void open(const string &fname, MPI_Comm comm);

inline void open(const string &fname, MPI_Comm comm)
{
  adiosname = fname;
  if(fname.find("config.bp")>= fname.size())
    adiosname.append(".config.bp");
  openfp = adios_read_open_file(adiosname.c_str(),
                                        ADIOS_READ_METHOD_BP,
                                        comm);
  if(openfp == NULL)
    qmcplusplus::app_error() <<"Fail to open adios file "<<adiosname<<" Abort."<<endl;
  else
    qmcplusplus::app_log()<<adiosname<<" is open "<<endl;
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

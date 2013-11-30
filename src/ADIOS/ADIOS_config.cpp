#include "ADIOS/ADIOS_config.h"

static bool UseHDF5 = false;
static bool UseADIOS = false;
static std::string adios_xml_filename;
const static std::string empty("");
static bool rdADIOS = false;
static bool rdHDF5 = false;

namespace ADIOS
{
void initialize(bool use_hdf5, bool use_adios)
{
  adios_xml_filename = empty;
  UseHDF5 = use_hdf5;
  UseADIOS = use_adios;
}

void initialze(std::string &xml_filename, bool use_hdf5, bool use_adios)
{
  adios_xml_filename = xml_filename;
  UseHDF5 = use_hdf5;
  UseADIOS = use_adios;
}

void initialize(const char *value)
 {
    if(!strcmp(value, "adios"))
    {
      rdADIOS=true;
      rdHDF5=false;
      qmcplusplus::app_log() << "Checkpoint restart from a .bp file" << std::endl;
    }
    else if(!strcmp(value, "hdf"))
    {
      rdADIOS=false;
      rdHDF5=true;
      qmcplusplus::app_log() << "Checkpoint restart from a .h5 file" << std::endl;
    }
    else
    {
      rdADIOS=false;
      rdHDF5=false;
      qmcplusplus::app_warning() << "Checkpoint restart method "<<value<<" is not supported."<< std::endl;
    }
}

bool getRdADIOS()
{
  return rdADIOS;
}

bool getRdHDF5()
{
  return rdHDF5;
}

#ifdef HAVE_ADIOS
bool useADIOS()
{
  return UseADIOS;
}

bool useHDF5()
{
  return UseHDF5;
}
#endif

const std::string& get_adios_xml()
{
  if (!useADIOS() && !useHDF5())
  {
    //qmcplusplus::app_warning() << "Attempted to retrieve adios xml filename before initializing." << std::endl;
    return empty;
  }
  else
    return adios_xml_filename;
}
};

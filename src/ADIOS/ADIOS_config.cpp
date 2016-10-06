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
    
    


#include "ADIOS/ADIOS_config.h"

static bool UseHDF5 = false;
static bool UseADIOS = false;
static std::string adios_xml_filename;
const static std::string empty("");
static bool rdADIOS = false;
static bool rdHDF5 = false;
static bool adios_first_open = true;
static std::string trace_file = "trace.bp";
static bool adios_init_flag = false;

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

void set_adios_init(bool b){
  adios_init_flag = b;
}

bool get_adios_init(){
  return adios_init_flag;
}

bool getFirstOpen(){
  return adios_first_open;
}

void setFirstOpen(bool a){
  adios_first_open = a;
}

std::string getTraceFileName(){
  return trace_file;
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

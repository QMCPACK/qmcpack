//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file OhmmsInfo.cpp
 * @brief Definition of OhmmsInfo class.
 */
#include "config.h"
#include <cstdlib>
#include <cstdio>
#include "Utilities/OhmmsInfo.h"

bool OhmmsInfo::Writeable = false;
OhmmsInform* OhmmsInfo::Debug = 0;
OhmmsInform* OhmmsInfo::Warn = 0;
OhmmsInform* OhmmsInfo::Error = 0;
OhmmsInform* OhmmsInfo::Log = 0;

OhmmsInfo::OhmmsInfo(int argc, char** argv, int master)
{
  if(argc)
  {
    initialize(argv[0],master);
  }
  else
  {
    initialize("ohmms",master);
  }
}

OhmmsInfo::~OhmmsInfo()
{
}

OhmmsInfo::OhmmsInfo(const std::string& fin_name, int rank, int gid, int num_groups)
{
  Writeable = (rank == 0);
  if(Log==0)//check if this is the first time
  {
    Warn = new OhmmsInform("WARNING",false,Writeable);
    Error = new OhmmsInform("ERROR",false,Writeable);
    Log = new OhmmsInform(" ",false,Writeable);
  }
  if(Writeable && num_groups>1)
  {
    char fn[128];
    sprintf(fn,"%s.g%03d.qmc",fin_name.c_str(),gid);
    Log->set(fn);
    Warn->set(*Log,"WARNING");
    Error->set(*Log,"ERROR");
  }
#if defined(PRINT_DEBUG)
  if(Debug==0)
  {
    char fn[128];
    sprintf(fn,"%s.p%03d.debug",fin_name.c_str(),rank);
    Debug = new OhmmsInform("DEBUG",false,true);
    Debug->set(fn);
    Debug->getStream().setf(std::ios::scientific, std::ios::floatfield);
    Debug->getStream().precision(6);
  }
#endif
}

void OhmmsInfo::initialize(const char* froot, int master)
{
  //if(master > 0)
  //  Writeable = false;
  //else
  //  Writeable = true;
  // initialize the estimator to record data
  Warn = new OhmmsInform("WARNING",false,Writeable);
  Error = new OhmmsInform("ERROR",false,Writeable);
  Log = new OhmmsInform(" ",false,Writeable);
  Error->setStdError();
#ifdef PRINT_DEBUG
  Debug = new OhmmsInform("DEBUG",false,Writeable);
#endif
  Log->getStream().setf(std::ios::scientific, std::ios::floatfield);
  Log->getStream().precision(6);
  Log->getStream() << "<?xml version=\"1.0\"?>\n";
//    bool useone = true; //always share the std
//    if(useone)
//    {
//      // share a std::ostream of Warrning
//      Warn->set(froot);
//      Error->set(*Warn);
//      Log->set(*Warn);
//  #ifdef PRINT_DEBUG
//      Debug->set(*Warn);
//  #endif
//    }
//    else
//    {
//      char fname[128];
//      sprintf(fname, "%s.warn",froot);
//      Warn->set(fname);
//      sprintf(fname, "%s.error",froot);
//      Error->set(fname);
//      sprintf(fname, "%s.log",froot);
//      Log->set(fname);
//  #ifdef PRINT_DEBUG
//      sprintf(fname, "%s.debug",froot);
//      Debug->set(fname);
//  #endif
//    }
}

void OhmmsInfo::die(const char* msg)
{
  Warn->getStream() << msg << std::endl;
  Warn->getStream() << "Stop the execution." << std::endl;
  exit(1);
}

/** flush the buffer
 */
void OhmmsInfo::flush()
{
  //if(!Writeable) {
  Log->getStream().flush();
  Error->getStream().flush();
  Warn->getStream().flush();
#ifdef PRINT_DEBUG
  Debug->getStream().flush();
#endif
  //}
}
//std::ostream& app_log(){ return OhmmsInfo::Log->getStream();}
//std::ostream& app_error(){ return OhmmsInfo::Error->getStream();}
//std::ostream& app_warning(){ return OhmmsInfo::Warn->getStream();}
//std::ostream& app_debug(){ return OhmmsInfo::Debug->getStream();}

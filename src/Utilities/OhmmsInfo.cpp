//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OhmmsInfo.cpp
 * @brief Definition of OhmmsInfo class.
 */
#include "Utilities/OhmmsInfo.h"

bool OhmmsInfo::Writeable = false;
OhmmsInform* OhmmsInfo::Debug = 0;
OhmmsInform* OhmmsInfo::Warn = 0;
OhmmsInform* OhmmsInfo::Error = 0;
OhmmsInform* OhmmsInfo::Log = 0;

OhmmsInfo::OhmmsInfo(int argc, char** argv, int master)
{
  if(argc) {
    initialize(argv[0],master);
  } else {
    initialize("ohmms",master);
  }
}

OhmmsInfo::~OhmmsInfo() {

}

void OhmmsInfo::initialize(const char* froot, int master){

  if(master > 0) 
    Writeable = false;
  else 
    Writeable = true;

  // initialize the estimator to record data
  Warn = new OhmmsInform("WARNING",false,Writeable);
  Error = new OhmmsInform("ERROR",false,Writeable);
  Log = new OhmmsInform(" ",false,Writeable);
#ifdef PRINT_DEBUG
  Debug = new OhmmsInform("DEBUG",false,Writeable);
#endif

  Log->getStream().setf(std::ios::scientific, std::ios::floatfield);
  Log->getStream().precision(6);

//    bool useone = true; //always share the std
//    if(useone) 
//    {
//      // share a ostream of Warrning
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

void OhmmsInfo::die(const char* msg) {
  Warn->getStream() << msg << std::endl;
  Warn->getStream() << "Stop the execution." << std::endl;
  exit(1);
}

/** flush the buffer
 */
void OhmmsInfo::flush() {

  if(!Writeable) {
    Log->getStream().flush();
    Error->getStream().flush();
    Warn->getStream().flush();
#ifdef PRINT_DEBUG
    Debug->getStream().flush();
#endif
  }
}
//std::ostream& app_log(){ return OhmmsInfo::Log->getStream();}
//std::ostream& app_error(){ return OhmmsInfo::Error->getStream();}
//std::ostream& app_warning(){ return OhmmsInfo::Warn->getStream();}
//std::ostream& app_debug(){ return OhmmsInfo::Debug->getStream();}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

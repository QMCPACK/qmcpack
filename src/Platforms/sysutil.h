//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_SYSTEM_UTILITIES_H
#define OHMMS_SYSTEM_UTILITIES_H

/*!\file sysutil.h
 * Function declarations to get system information.
 */
#include <string>

//!< return the user name
std::string getUserName();

//!< return the host name
std::string getHostName();

//!< return the date and time
std::string getDateAndTime();

/** get the time and date with a format
 */
std::string getDateAndTime(const char* format);

#ifdef __linux__
#include "sys/sysinfo.h"

inline size_t freemem()
{
  struct sysinfo si;
  sysinfo(&si);
  si.freeram+=si.bufferram;
  return si.freeram>>20;
  //return (si.freeram + si.bufferram);
}
#else

inline size_t freemem()
{
  return 0;
}

#endif

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

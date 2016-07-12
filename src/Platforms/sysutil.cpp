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
/*!\file sysutil.cpp
 * Implement function to get system information in Unix environment.
 */

#include <string>
#include <sstream>
#include <iostream>
using std::string;
#include <time.h>


// Dummy version of getHostName, in case its needed
#if 0
#if defined(_CRAYMPI) || defined(XT_CATAMOUNT)
string getHostName()
{
  return "jaguar";
}
#endif
#endif

#include <sys/utsname.h>

string getHostName()
{
  utsname mysys;
  uname(&mysys);
  return std::string(mysys.nodename);
}

string getDateAndTime()
{
  time_t now;
  time(&now);
  return ctime(&now);
}

string getDateAndTime(const char* format)
{
  time_t now;
  time(&now);
  tm* now_c = localtime(&now);
  char d[32];
  strftime(d,32,format,now_c);
  return std::string(d);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

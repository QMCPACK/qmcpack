//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/*!\file sysutil.cpp
 * Implement function to get system information in Unix environment.
 */

#include <string>
using std::string;

#include <time.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <pwd.h>

string getUserName() {
  struct passwd *who;
  if((who = getpwuid(getuid())) != NULL) {
    return who->pw_name;
  }
  return "auser";
}

string getDateAndTime() {
  time_t now;
  time(&now);
  return ctime(&now);
}

string getHostName() {
  utsname mysys;
  uname(&mysys);
  return string(mysys.nodename);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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


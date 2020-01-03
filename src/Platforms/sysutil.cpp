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


#include "Platforms/sysutil.h"
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
  strftime(d, 32, format, now_c);
  return std::string(d);
}

#ifdef __linux__
#include <sys/sysinfo.h>
#include <sys/resource.h>
#endif

size_t freemem()
{
#ifdef __linux__
  struct sysinfo si;
  sysinfo(&si);
  si.freeram += si.bufferram;
  return si.freeram >> 20;
#else
  return 0;
#endif
}

/* returns heap memory usage in MB */
size_t memusage()
{
#ifdef __linux__
  struct rusage RU; /* heap memory usage */
  getrusage(RUSAGE_SELF, &RU);
  return RU.ru_maxrss >> 10;
#else
  return 0;
#endif
}

void print_mem(const char* title, std::ostream& log)
{
  char msg[256];
  sprintf(msg, "===== %s =====\n", title);
  log << msg;
  sprintf(msg, "Available memory on node 0, free + buffers : %zu MB\n", freemem());
  log << msg;
  sprintf(msg, "Memory footprint by rank 0 on node 0       : %zu MB\n", memusage());
  log << msg;
  sprintf(msg, "==================================================\n");
  log << msg;
}

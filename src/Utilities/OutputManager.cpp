//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////



#include <Utilities/OutputManager.h>


InfoStream infoSummary(&std::cout);
InfoStream infoLog(&std::cout);
InfoStream infoError(&std::cerr);
InfoStream infoDebug(&std::cout);

// default verbosity is high for now.  After converting to app_summary, the default will be low.
OutputManagerClass outputManager(Verbosity::HIGH);


void
OutputManagerClass::setVerbosity(Verbosity level)
{
  global_verbosity_level = level;
  if (isActive(Verbosity::DEBUG)) {
    infoSummary.resume();
    infoLog.resume();
    infoDebug.resume();
  } else if (isActive(Verbosity::HIGH)) {
    infoSummary.resume();
    infoLog.resume();
    infoDebug.pause();
  } else if (isActive(Verbosity::LOW)) {
    infoSummary.resume();
    infoLog.pause();
    infoDebug.pause();
  }
}

bool OutputManagerClass::isActive(Verbosity level)
{
  return level <= global_verbosity_level;
}

void OutputManagerClass::pause()
{
  infoSummary.pause();
  infoLog.pause();
}

void OutputManagerClass::resume()
{
  infoSummary.resume();
  infoLog.resume();
}

void OutputManagerClass::shutOff()
{
  infoSummary.shutOff();
  infoLog.shutOff();
  infoError.shutOff();
  infoDebug.shutOff();
}


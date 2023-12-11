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


/** @file OutputManager.h
 * @brief Declaration of OutputManager class.
 */
#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include "InfoStream.h"


enum class Verbosity
{
  LOW,
  HIGH,
  DEBUG
};

extern InfoStream infoSummary;
extern InfoStream infoLog;
extern InfoStream infoError;
extern InfoStream infoDebug;

class OutputManagerClass
{
  Verbosity global_verbosity_level;

public:
  OutputManagerClass(Verbosity level = Verbosity::LOW) { setVerbosity(level); }

  void setVerbosity(Verbosity level);

  bool isActive(Verbosity level) const;

  bool isDebugActive() const { return isActive(Verbosity::DEBUG); }

  bool isHighActive() const { return isActive(Verbosity::HIGH); }

  /// Pause the summary and log streams
  void pause();

  /// Resume the summary and log streams
  void resume();

  /// Permanently shut off all streams
  void shutOff();
};

extern OutputManagerClass outputManager;

namespace qmcplusplus
{
inline std::ostream& app_summary() { return infoSummary.getStream(); }

inline std::ostream& app_log() { return infoLog.getStream(); }

inline std::ostream& app_error() { return infoError.getStream() << "ERROR "; }

inline std::ostream& app_warning() { return infoLog.getStream() << "WARNING "; }

inline std::ostream& app_debug_stream() { return infoDebug.getStream(); }

// From https://stackoverflow.com/questions/11826554/standard-no-op-output-stream
// If debugging is not active, this skips evaluation of the arguments
#define app_debug                        \
  if (!outputManager.isDebugActive()) {} \
  else                                   \
    app_debug_stream


// Keep these macros temporarily until all output uses streams
#define LOGMSG(msg)                             \
  {                                             \
    qmcplusplus::app_log() << msg << std::endl; \
  }
#define ERRORMSG(msg)                \
  {                                  \
    app_error() << msg << std::endl; \
  }
#define WARNMSG(msg)                   \
  {                                    \
    app_warning() << msg << std::endl; \
  }
#ifdef PRINT_DEBUG
#define DEBUGMSG(msg)                \
  {                                  \
    app_debug() << msg << std::endl; \
  }
#else
#define DEBUGMSG(msg)
#endif

#define XMLReport(msg)

}; // namespace qmcplusplus

#endif

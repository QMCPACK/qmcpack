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

#include <Utilities/InfoStream.h>


enum class Verbosity {LOW, HIGH, DEBUG};
enum class LogType {SUMMARY, APP, ERROR, DEBUG};

extern InfoStream infoSummary;
extern InfoStream infoLog;
extern InfoStream infoError;
extern InfoStream infoDebug;

class OutputManagerClass
{
  Verbosity global_verbosity_level;

public:
  OutputManagerClass(Verbosity level=Verbosity::LOW) { setVerbosity(level); }

  void setVerbosity(Verbosity level);

  bool isActive(Verbosity level);

  bool isDebugActive()
  {
    return isActive(Verbosity::DEBUG);
  }

  bool isHighActive()
  {
    return isActive(Verbosity::HIGH);
  }

  std::ostream& getStream(LogType log) {
    switch (log) {
    case LogType::SUMMARY:
      return infoSummary.getStream();
    case LogType::APP:
      return infoLog.getStream();
    case LogType::ERROR:
      return infoError.getStream();
    case LogType::DEBUG:
      return infoDebug.getStream();
    }
  }

  /// Pause the summary and log streams
  void pause();

  /// Resume the summary and log streams
  void resume();

  /// Permanently shut off all streams
  void shutOff();
};

extern OutputManagerClass outputManager;

namespace qmcplusplus {

inline std::ostream& app_summary()
{
  return outputManager.getStream(LogType::SUMMARY);
}

inline std::ostream& app_log()
{
  return outputManager.getStream(LogType::APP);
}

inline std::ostream& app_error()
{
  outputManager.getStream(LogType::ERROR) << "ERROR ";
  return outputManager.getStream(LogType::ERROR);
}

inline std::ostream& app_warning()
{
  outputManager.getStream(LogType::ERROR) << "WARNING ";
  return outputManager.getStream(LogType::ERROR);
}

inline std::ostream& app_debug_stream()
{
  return outputManager.getStream(LogType::DEBUG);
}

// From https://stackoverflow.com/questions/11826554/standard-no-op-output-stream
// If debugging is not active, this skips evaluation of the arguments
#define app_debug if (!outputManager.isDebugActive()) {} else app_debug_stream



// Keep these macros temporarily until all output uses streams
#define LOGMSG(msg) {qmcplusplus::app_log() << msg << std::endl;}
#define ERRORMSG(msg) {app_error() << msg << std::endl;}
#define WARNMSG(msg) {app_warning() << msg << std::endl;}
#ifdef PRINT_DEBUG
#define DEBUGMSG(msg) {app_debug() << msg << std::endl;}
#else
#define DEBUGMSG(msg)
#endif

#define XMLReport(msg)

};

#endif

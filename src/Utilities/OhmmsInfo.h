//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003 - by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file OhmmsInfo.h
 * @brief Declaration of OhmmsInfo class.
 */
#ifndef OHMMS_OHMMSINFO_H
#define OHMMS_OHMMSINFO_H

#include "Utilities/OhmmsInform.h"
/** Control object for run-time information
 *
 * Similar class to PoomaInfo  of Pooma with very limited functions
 */
class OhmmsInfo
{

public:

  static bool Writeable;
  static OhmmsInform *Debug;
  static OhmmsInform *Warn;
  static OhmmsInform *Error;
  static OhmmsInform *Log;

  static void initialize(const char* froot, int master);
  static void die(const char*);
  static void flush();

  /** constructor using command-line arguments
   * @param argc number of arguments
   * @param argv arguments
   * @param master selects an ostream that does writing
   */
  OhmmsInfo(int argc, char** argv, int master=-1);
  /** constructor
   * @param fin_name input file name
   * @param rank node rank
   * @param group groupd rank
   * @param multi_run true, open one ostream per group
   */
  OhmmsInfo(const std::string& fin_name, int rank=0, int gid=0, int num_groups=1);
  ~OhmmsInfo();
private:
  OhmmsInfo() { }
};

//extern std::ostream& app_log();
//extern std::ostream& app_debug();
//extern std::ostream& app_warn();
//extern std::ostream& app_error();

/**run-time messages
 * - LOGMGS log message
 * - ERRORMSG error message
 * - WARNMSG warning message
 * - DEBUGMSG debug message
 */
#ifdef DONOTUSEOHMMSINFO
#define LOGMSG(msg)
#define ERRORMSG(msg)
#define WARNMSG(msg)
#define DEBUGMSG(msg)
#define XMLReport(msg)
#else
#define LOGMSG(msg) \
 { OhmmsInfo::Log->getStream() << msg << std::endl;}
#define ERRORMSG(msg) \
 { OhmmsInfo::Error->getStream() << "ERROR " << msg << std::endl;}
#define WARNMSG(msg) \
 { OhmmsInfo::Warn->getStream() << "WARN " << msg << std::endl;}
#define XMLReport(msg)
//#define XMLReport(msg) \
//{std::cout<< "XML " << msg << std::endl;}

#ifdef PRINT_DEBUG
#define DEBUGMSG(msg) OhmmsInfo::Debug->getStream() << msg << std::endl
#else
#define DEBUGMSG(msg)
#endif
#endif
#endif//OHMMS_OHMMSINFO_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

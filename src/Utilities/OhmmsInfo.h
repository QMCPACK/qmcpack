//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
/*! \class OhmmsInfo
 *  \brief Similar class to PoomaInfo  of Pooma with very limited functions
 */
#ifndef OHMMS_OHMMSINFO_H
#define OHMMS_OHMMSINFO_H

#include "Utilities/OhmmsInform.h"
class OhmmsInfo {

public:

  static OhmmsInform *Debug;
  static OhmmsInform *Warn;
  static OhmmsInform *Error;
  static OhmmsInform *Log;
  static void initialize(const char* froot, int master);
  static void die(const char*);
  
  OhmmsInfo(int argc, char** argv, int master=-1); 
  ~OhmmsInfo(); 
  OhmmsInfo(){ }

};

/**references to run-time streams
 *@warning NOT utilized in real applications
 */
extern std::ostream& log();
extern std::ostream& error();
extern std::ostream& warning();
extern std::ostream& debug();

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
  { if(OhmmsInfo::Log->open()) OhmmsInfo::Log->getStream() << msg << std::endl;}
#define ERRORMSG(msg) \
  { if(OhmmsInfo::Error->open()) OhmmsInfo::Error->getStream() << msg << std::endl;}
#define WARNMSG(msg) \
 { if(OhmmsInfo::Warn->open()) OhmmsInfo::Warn->getStream() << msg << std::endl;}
#define XMLReport(msg) \
{std::cout<< "XML " << msg << std::endl;}

#ifdef PRINT_DEBUG
#define DEBUGMSG(msg) { OhmmsInfo::Debug->getStream() << msg << std::endl;}
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

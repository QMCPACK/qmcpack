//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Jeongnim Kim
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
/** @file ProgressReportEngine.h
 * @brief declaration of ProgressReportEngine
 */

#ifndef QMCPLUSPLUS_PROGRESSREPORTENGINE_H
#define QMCPLUSPLUS_PROGRESSREPORTENGINE_H

#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"
#include "Utilities/Timer.h"
#include "OhmmsData/OhmmsElementBase.h"

namespace qmcplusplus
{
/**
 *
 * Final class and should not be derived.
 */
class ReportEngine
{
public:

  inline ReportEngine(const std::string& cname, const std::string& fname, int atype=1):
    ReportType(atype),ClassName(cname), FuncName(fname), LogBuffer(*OhmmsInfo::Log)
  {
    LogBuffer << "  " << ClassName << "::" << FuncName << "\n";
    //if(ReportType)
    //  LogBuffer << ("<echo className=\""+ClassName+"\" funcName=\""+FuncName+"\">\n");
    //else
    //  LogBuffer << ("<"+ClassName+">\n");
    LogBuffer.flush();//always flush
  }

  inline ~ReportEngine()
  {
    //if(ReportType)
    //  LogBuffer << "</echo>\n";
    //else
    //  LogBuffer << ("</"+ClassName+">\n");
    LogBuffer.flush();
  }

  inline void flush()
  {
    LogBuffer << '\n';
    LogBuffer.flush();
  }

  inline void startWarning()
  {
    LogBuffer << "<warning>\n";
  }
  inline void endWarning()
  {
    LogBuffer << "</warning>\n";
  }

  inline void warning(const std::string& msg)
  {
    LogBuffer << ("<warning>"+msg+"</warning>\n");
  }

  inline void startError()
  {
    LogBuffer << "<error node=\"" << OHMMS::Controller->rank() << "\">\n";
  }

  inline void endError(bool fatal=false)
  {
    LogBuffer << "</error>\n";
    if(fatal)
      APP_ABORT(ClassName+"::"+FuncName);
  }

  inline void error(const std::string& msg, bool fatal=false)
  {
    LogBuffer << ("<error>"+msg+"</error>\n");
    if(fatal)
      APP_ABORT(ClassName+"::"+FuncName);
  }

  void echo(xmlNodePtr cur, bool recursive=false);

private:
  ///type of report
  int ReportType;
  /** class Name
   */
  std::string ClassName;
  /** name of the current  member function
   */
  std::string FuncName;
  /** arguments
   */
  std::vector<std::string> ArgList;
  /** stream for log message
   */
  OhmmsInform& LogBuffer;
  //disable copy constructor
  ReportEngine(const ReportEngine& a):LogBuffer(*OhmmsInfo::Log) {}
};

// templated version of operator<< for Inform objects
template<class T>
inline
ReportEngine& operator<<(ReportEngine& o, const T& val)
{
  app_log()<< val;
  return o;
}
}
#endif // QMCPLUSPLUS_MPIOBJECTBASE_H
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2468 $   $Date: 2008-02-22 09:27:30 -0500 (Fri, 22 Feb 2008) $
 * $Id: Communicate.h 2468 2008-02-22 14:27:30Z jnkim $
 ***************************************************************************/

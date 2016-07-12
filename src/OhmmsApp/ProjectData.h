//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_PROJECTDATA_H__
#define QMCPLUSPLUS_PROJECTDATA_H__

#include "OhmmsData/OhmmsElementBase.h"
//#include <vector>
//#include <string>
//#include <iostream>
#include "Message/Communicate.h"

namespace qmcplusplus
{

/**class ProjectData
 *\brief Encapsulate data for a project
 *
 * Default: myName = "Project"
 * Should not modify the name, since composite types, such as MDRunData, use the name.
 *
 */
struct ProjectData: public OhmmsElementBase
{

  /// constructor
  ProjectData(const char* aname=0);

  bool get(std::ostream& os) const;
  bool put( std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();

  ///increment a series number and reset m_projectroot
  void advance();

  ///roll-back a series number and reset m_projectroot by one
  void rewind();

  ///set the title
  inline void setTitle(const std::string& atitle)
  {
    m_title=atitle;
    reset();
  }

  void  setCommunicator(Communicate* c);

  ///returns the name of the project
  inline const char* CurrentMainRoot() const
  {
    return m_projectmain.c_str();
  }

  ///returns the name of the project
  inline const char* CurrentRoot() const
  {
    return m_projectroot.c_str();
  }

  ///returns the name of the project
  inline const char* NextRoot() const
  {
    return m_nextroot.c_str();
  }

  /** return the root of the previous sequence
   * @param oldroot is composed by the m_title and m_series
   */
  bool PreviousRoot( std::string& oldroot) const;

  ///title of the project
  std::string m_title;

  ///name of the host where the job is running
  std::string m_host;

  ///date when the job is executed
  std::string m_date;

  ///main root for all the output engines
  std::string m_projectmain;

  ///processor-dependent root for all the output engines
  std::string m_projectroot;

  ///root for the next run
  std::string m_nextroot;

  ///series index
  int m_series;

  ///communicator
  Communicate* myComm;

  ///the xml node for <Project/>
  xmlNodePtr m_cur;
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

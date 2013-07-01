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

#include "Configuration.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
void ReportEngine::echo(xmlNodePtr cur, bool recursive)
{
  if(cur==NULL)
    return;
  app_log()<< "<input node=\""<<(const char*)(cur->name) << "\"";
  xmlAttrPtr att = cur->properties;
  char atext[1024];
  while(att != NULL)
  {
    sprintf(atext,"  %s=\"%s\"",(const char*)(att->name),(const char*)(att->children->content));
    app_log() << atext;
    att = att->next;
  }
  app_log() << "/>\n";
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2468 $   $Date: 2008-02-22 09:27:30 -0500 (Fri, 22 Feb 2008) $
 * $Id: Communicate.h 2468 2008-02-22 14:27:30Z jnkim $
 ***************************************************************************/

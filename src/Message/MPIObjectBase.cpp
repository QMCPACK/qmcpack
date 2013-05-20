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
#include "Message/MPIObjectBase.h"

namespace qmcplusplus
{

MPIObjectBase::MPIObjectBase(Communicate* c): ReportLevel(1),
  myComm(0), ClassName("MPIObjectBase")
{
  initCommunicator(c);
  if(myComm->rank())
    ReportLevel=0;
}

//MPIObjectBase::MPIObjectBase(const MPIObjectBase& a): myComm(0)
//{}

MPIObjectBase::~MPIObjectBase()
{}

void MPIObjectBase::initCommunicator(Communicate* c)
{
  if(myComm && myComm == c)
    return;
  myComm = c ? c:OHMMS::Controller;
}

void MPIObjectBase::setReportLevel(int level)
{
  //demote the level if not the head node
  ReportLevel=(myComm->rank())?0:level;
  //if(ReportLevel)
  //{//inherit the global info streams
  //  LogBuffer.set(*OhmmsInfo::Log,ClassName);
  //  WarnBuffer.set(*OhmmsInfo::Warn);
  //  ErrorBuffer.set(*OhmmsInfo::Error);
  //}
}

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2468 $   $Date: 2008-02-22 09:27:30 -0500 (Fri, 22 Feb 2008) $
 * $Id: Communicate.h 2468 2008-02-22 14:27:30Z jnkim $
 ***************************************************************************/

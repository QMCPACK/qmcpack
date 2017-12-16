//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




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

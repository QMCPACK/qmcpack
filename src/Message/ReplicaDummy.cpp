//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




#include "Message/ReplicaControl.h"

extern "C" {
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
}

// instantiate the Communicator
//Communicate* ReplicaControl::OhmmsComm = new Communicate;

ReplicaControl::ReplicaControl()
{
  ReplicaID = std::vector<int>(1);
  NumReplica = OhmmsComm->getNumNodes();
  ReplicaID[0] = OhmmsComm->getNodeID();
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}

ReplicaControl::ReplicaControl(int argc, char **argv)
{
  ReplicaID = std::vector<int>(1);
  NumReplica = OhmmsComm->getNumNodes();
  ReplicaID[0] = OhmmsComm->getNodeID();
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}

void ReplicaControl::init(int nproc)
{
  if(nproc != ReplicaID.size())
  {
    ReplicaID = std::vector<int>(nproc);
    for(int i=0; i<nproc; i++)
      ReplicaID[i] = i;
  }
} // number of processors to

void ReplicaControl::resize(int n)
{
  if(WallClockTime.size() < n)
  {
    WallClockTime.erase(WallClockTime.begin(), WallClockTime.end());
    MDTime.erase(MDTime.begin(), MDTime.end());
    HMDTime.erase(HMDTime.begin(), HMDTime.end());
    WallClockTime = std::vector<double>(2*n,0.0e0);
    MDTime   = std::vector<double>(2*n,0.0e0);
    HMDTime  = std::vector<double>(2*n,0.0e0);
  }
}

bool ReplicaControl::send(const ParticlePos_t& p, int node, int btag, bool bcast)
{
  return true;
}

bool
ReplicaControl::recv(ParticlePos_t& p, int& node, int btag)
{
  return true;
}


bool ReplicaControl::sendDirections(int& inode)
{
  return true;
}

int ReplicaControl::recvDirections()
{
  return 0;
}

bool ReplicaControl::reportTransition(double drsq)  //slaves
{
  return true;
}

void ReplicaControl::catchTransition(int& remNode, double& wallt, double& rsqmax)  // master
{
}

bool ReplicaControl::sendTransitionTime(double t0)
{
  return true;
}

bool ReplicaControl::recvTransitionTime(double& t0)
{
  return true;
}

bool ReplicaControl::sendExactTime(double wallt, bool havewon)
{
  return true;
}

void ReplicaControl::collectExactTimes(int itrans, double t0)
{
}

void ReplicaControl::reportTimeData(std::ostream& os,int i_trans)
{
  os << std::setw(18) << MDTime[i_trans] << std::setw(18) << TotMD_t << " ";
}



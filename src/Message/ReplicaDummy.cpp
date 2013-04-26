//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//Jeongnim Kim
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
  ReplicaID = vector<int>(1);
  NumReplica = OhmmsComm->getNumNodes();
  ReplicaID[0] = OhmmsComm->getNodeID();
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}

ReplicaControl::ReplicaControl(int argc, char **argv)
{
  ReplicaID = vector<int>(1);
  NumReplica = OhmmsComm->getNumNodes();
  ReplicaID[0] = OhmmsComm->getNodeID();
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}

void ReplicaControl::init(int nproc)
{
  if(nproc != ReplicaID.size())
  {
    ReplicaID = vector<int>(nproc);
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
    WallClockTime = vector<double>(2*n,0.0e0);
    MDTime   = vector<double>(2*n,0.0e0);
    HMDTime  = vector<double>(2*n,0.0e0);
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

void ReplicaControl::reportTimeData(ostream& os,int i_trans)
{
  os << setw(18) << MDTime[i_trans] << setw(18) << TotMD_t << " ";
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * POOMA_VERSION_ID: $Id$
 ***************************************************************************/

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
    
    



#include <mpi.h>
#include "Message/ReplicaControl.h"

extern "C" {
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
}

ReplicaControl::ReplicaControl()
{
  //MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NumReplica);
  MPI_Comm_rank(MPI_COMM_WORLD,&ReplicaID);
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}


ReplicaControl::ReplicaControl(int argc, char **argv)
{
  //MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&NumReplica);
  MPI_Comm_rank(MPI_COMM_WORLD,&ReplicaID);
  TotMD_t = 0.0e0;
  TotHMD_t = 0.0e0;
}

void resize(int n)
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

ReplicaControl::~ReplicaControl()
{
}

bool
ReplicaControl::send(ParticlePos_t& p, int node, int btag, bool bcast)
{
  if(bcast)
    // broadcasting positions to all the nodes
  {
    int size = OHMMS_DIM*newp.size();
    int ierr = MPI_Bcast(&p[0][0], size, MPI_DOUBLE, ReplicaID, btag,
                         MPI_COMM_WORLD);
    return true;
  }
  else
  {
    if(node < 0)
      return true; // do nothing
    MPI_Request request;
    int size = OHMMS_DIM*p.size();
    int errstat = MPI_Isend(&p[0][0], size, MPI_DOUBLE, node, btag,
                            MPI_COMM_WORLD, &request);
    return true;
  }
  return false;
}

bool
ReplicaControl::recv(ParticlePos_t& p, int& node, int btag)
{
  {
    int size = OHMMS_DIM*p.size();
    MPI_Request request;
    MPI_Status status;
    MPI_Irecv(&p[0][0], size, node, btag, MPI_COMM_WORLD, &request);
    bool flag = false;
    while (!flag)
    {
      MPI_Test(&request, &flag, &status);
    }
    return true;
  }
  bool ReplicaControl::sendDirections(int& inode)
  {
    int msg[2];
    if(inode < 0)
    {
      msg[0] = PARMD::finalize_tag;
      MPI_Bcast(&msg[0], 1, MPI_INT, ReplicaID,PARMD::command_tag,MPI_COMM_WORLD);
      return true;
    }
    MPI_Request request;
    MPI_Status status;
    msg[0] = PARMD::winner_tag;
    msg[1] = PARMD::losers_tag;
    bool flag = false;
    for(int i=1; i<NumReplica; i++)
    {
      if(i == inode)
      {
        MPI_Isend(&msg[0],1,MPI_INT,i,PARMD::command_tag,
                  MPI_COMM_WORLD, &request);
      }
      else
      {
        MPI_Isend(&msg[1],1,MPI_INT,i,PARMD::command_tag,
                  MPI_COMM_WORLD, &request);
      }
      flag = false;
      while (!flag)
      {
        MPI_Test(&request, &flag, &status);
      }
    }
    return true;
  }
  int ReplicaControl::recvDirections()
  {
    int whattodo;
    MPI_Request request;
    MPI_Status status;
    MPI_Irecv(&whattodo, 1, MPI_INT, 0, PARMD::command_tag,
              MPI_COMM_WORLD, &request);
    bool flag = false;
    while (!flag)
    {
      MPI_Test(&request, &flag, &status);
    }
    return whatodo;
  }
  bool ReplicaControl::reportTransition(double drsq)  //slaves
  {
    double newdata[2];
    MPI_Request request;
    newdata[0] = WallClockTime[FrameNum];
    newdata[1] = drsq;
    MPI_Isend(newdata,2,MPI_DOUBLE,0,PARMD::transition, MPI_COMM_WORLD, &request);
    return true;
  }
  void ReplicaControl::catchTransition(int& remNode, double& wallt, double& rsqmax)  // master
  {
    double newdata[2];
    static MPI_Request request;
    static MPI_Status status;
    static bool flag = true;
    remNode = -1;
    if(flag)
      // initiate recv
    {
      DEBUGMSG("Open a channel to catch the transition\n");
      MPI_Irecv(newdata, 2, MPI_INT, MPI_ANY_SOURCE,
                PARMD::transition_tag, MPI_COMM_WORLD, &request);
      flag = false;
    }
    MPI_Test(&request, &flag, &status);
    if(flag)
      // a transition is reported
    {
      wallt = newdata[0];
      rsqmax = newdata[1];
      remNode =  status.MPI_SOURCE; // copy the remote node
    }
  }
  bool ReplicaControl::sendTransitionTime(double t0)
  {
    MPI_Bcast(&t0, 1, MPI_DOUBLE, ReplicaID, PARMD::stop_tag, MPI_COMM_WORLD);
    return true;
  }
  bool ReplicaControl::recvTransitionTime(double& t0)
  {
    double trecv = 0;
    MPI_Status status;
    MPI_Recv(&trecv, 1, MPI_DOUBLE, 0, PARMD::stop_tag, MPI_COMM_WORLD, &status);
    t0 = trecv;
    return true;
  }
  bool ReplicaControl::sendExactTime(double wallt, bool havewon)
  {
    int masternode = 0;
    double newtime[3];// std::vector<double> newtime(3);
    if(havewon)
    {
      newtime[0] = MDTime[FrameNum];
      newtime[1] = WallClockTime[FrameNum];
    }
    else
    {
      int iframe=FrameNum;
      while(WallClockTime[iframe] >= wallt && iframe > 0)
      {
        iframe--;
      }
      newtime[0] = MDTime[iframe];
      newtime[1] = WallClockTime[iframe];
    }
    MPI_Request request;
    MPI_Isend(newtime,3,MPI_DOUBLE,0,PARMD::exacttime_tag, MPI_COMM_WORLD, &request);
  }
  void ReplicaControl::collectExactTimes(int itrans, double t0)
  {
    MPI_Status status;
    MPI_Request request;
    int btag = PARMD::exacttime_tag;
    int numrecvd =  NumReplica-1;
    double newtime[3];
    bool flag;
    CumMD_t = t0; // constant time to be added to a new cumm. time
    while(numrecvd > 0)
    {
      MPI_Irecv(newtime, 3, MPI_DOUBLE, MPI_ANY_SOURCE,
                PARMD::exacttime_tag, MPI_COMM_WORLD, &request);
      flag = false;
      while(!false)
      {
        MPI_Test(&request, &flag, &status);
      }
      CumMD_t += newtime[0];
      numrecvd--;
    }
    MDTime[itrans] = CumMD_t;
    TotMD_t += CumMD_t;
  }
  void ReplicaControl::reportTimeData(std::ostream& os,int i_trans)
  {
    os << std::setw(18) << MDTime[i_trans] << std::setw(18) << TotMD_t << " ";
  }
//  bool ReplicaControl::sendNewConfig(ParticlePos_t& newp){
//    int size = OHMMS_DIM*newp.size();
//    MPI_Bcast(&newp[0][0], size, MPI_DOUBLE, ReplicaID, PARMD_newconfig_tag,
//              MPI_COMM_WORLD);
//    return true;
//  }
//  void ReplicaControl::recvNewConfig(ParticlePos_t& newp){
//    int btag = PARMD_newconfig_tag;
//    int size = OHMMS_DIM*newp.size();
//    MPI_Request request;
//    MPI_Status status;
//    int masternode = 0;
//    MPI_Irecv(&newp[0][0], size, masternode, btag, MPI_COMM_WORLD, &request);
//    bool flag = false;
//    while (!flag) {
//      MPI_Test(&request, &flag, &status);
//    }
//  }
//  bool ReplicaControl::sendMinConfig(ParticlePos_t& minp) {
//    MPI_Request request;
//    MPI_Status status;
//    int size = OHMMS_DIM*newp.size();
//    int errstat = MPI_Isend(&minp[0][0], size, MPI_DOUBLE, 0, PARMD_minconfig_tag,
//    		          MPI_COMM_WORLD, &request);
//    return true;
//  }
//  void ReplicaControl::recvMinConfig(ParticlePos_t& minp, int& inode) {
//    int btag = PARMD_minconfig_tag;
//    int size = OHMMS_DIM*newp.size();
//    MPI_Request request;
//    MPI_Status status;
//    MPI_Irecv(&newp[0][0], size, inode, btag, MPI_COMM_WORLD, &request);
//    bool flag = false;
//    while (!flag) {
//      MPI_Test(&request, &flag, &status);
//    }
//  }
//  void ReplicaControl::makeNodeConnection(const char *opt) {
//    if(!strcmp(opt,"chain")) {
//      RightNode =ReplicaID+1;
//      LeftNode = RightNode-2;
//      if(RightNode == NumReplica) RightNode = -1; // not periodic
//    }
//  }

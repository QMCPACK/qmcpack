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
    
    





/******************************************************************* 
// ReplicaControl using MPI
// Responsible for book keeping for parallel replica with or without hmd
// Using timer routines of Parmd7, Art Voter
*******************************************************************/
#ifndef OHMMS_REPLICACONTROL_H
#define OHMMS_REPLICACONTROL_H
#include <vector>

#include "Particle/OhmmsParticle.h"
#include "Message/Communicate.h"

namespace PARMD
{
const int transition_tag = 10;
const int exacttime_tag = 11;
const int newconfig_tag = 12;
const int minconfig_tag = 13;
const int command_tag = 14;
const int stop_tag = 15;
const int pluscycle_tag = 16;
const int minuscycle_tag = 17;
const int winner_tag = 1;
const int losers_tag = 0;
const int noevent_tag = -1;
const int finalize_tag = 99;
}

class ReplicaControl
{

public:

  typedef OHMMS::Particle_t::ParticlePos_t ParticlePos_t;
  enum ReplicaTags { NEWPOS, CUMTIME, TRANSITION};

  //Communicate* OhmmsComm;

  ReplicaControl();
  ReplicaControl(int, char**);
  ~ReplicaControl();

  void init(int nproc); // number of processors to
  inline int getNumNodes() const
  {
    OhmmsComm->getNumNodes();
  }
  inline int getNodeID() const
  {
    OhmmsComm->getNodeID();
  }
  inline bool master() const
  {
    return OhmmsComm->master();
  }

  // Time-related funtions
  inline void saveMDTime(double nowT, double wallt, int iframe)
  {
    FrameNum = iframe; // match the frame index
    WallClockTime[FrameNum] = wallt;
    MDTime[FrameNum] = nowT;
  }

  inline double getWallClockTime()
  {
    return WallClockTime[FrameNum];
  }
  void resetTimer()   // for each transition, reset parameters
  {
    FrameNum = 0;
    CumMD_t = 0.0e0;
    CumHMD_t = 0.0e0;
  }

  bool send(const ParticlePos_t&, int node, int tag, bool bcast=false);
  bool recv(ParticlePos_t&, int& node, int tag);

//    // send/recv of new configurations
//    bool sendNewConfig(ParticlePos_t& );
//    void recvNewConfig(ParticlePos_t& );

//    // send/recv of a mininum configurations
//    bool sendMinConfig(ParticlePos_t& minp);
//    void recvMinConfig(ParticlePos_t& minp, int& inode);
//    bool sendNextNode(ParticlePos_t& , int );
//    bool recvNextNode(ParticlePos_t& , int );

  // send/irecv of the transition time
  bool reportTransition(double);
  void catchTransition(int& inode, double&, double&);

  // send transition time if a transition is detected
  // "the" master sends the first-in transition wall-clock time to slaves
  bool sendTransitionTime(double );
  bool recvTransitionTime(double&);

  // sends MD time at the wall-clock time when the winner found a transition
  // "the" master collect MD times sent from all other nodes
  bool sendExactTime(double, bool havewon = false );
  void collectExactTimes(int itrans, double t0);

  bool sendDirections(int& winner);
  int recvDirections();

  void reportTimeData(std::ostream& ,int itrans);

  void resize(int n);

protected:

  int NumReplica;
  std::vector<int> ReplicaID;

  double CumMD_t, CumHMD_t; // cummulative time over replica
  double TotMD_t, TotHMD_t; // total time
  int            FrameNum;
  int LeftNode, RightNode;

  std::vector<double> WallClockTime;  // wall-clock time
  std::vector<double> MDTime;    // mdtime
  std::vector<double> HMDTime; // hyper time will be added
};

#endif


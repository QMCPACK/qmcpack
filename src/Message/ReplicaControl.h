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
///////////////////////////////////////////////////////////////////////////
// ReplicaControl using MPI 
// Responsible for book keeping for parallel replica with or without hmd
// Using timer routines of Parmd7, Art Voter 
///////////////////////////////////////////////////////////////////////////
#ifndef OHMMS_REPLICACONTROL_H
#define OHMMS_REPLICACONTROL_H
#include <vector>
using std::vector;
#include "Particle/OhmmsParticle.h"
#include "Message/Communicate.h"

namespace PARMD {
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

class ReplicaControl {

public:

  typedef OHMMS::Particle_t::ParticlePos_t ParticlePos_t;  
  enum ReplicaTags { NEWPOS, CUMTIME, TRANSITION};

  //Communicate* OhmmsComm;

  ReplicaControl();
  ReplicaControl(int, char**);
  ~ReplicaControl();

  void init(int nproc); // number of processors to 
  inline int getNumNodes() const {OhmmsComm->getNumNodes();}
  inline int getNodeID() const {OhmmsComm->getNodeID();}
  inline bool master() const { return OhmmsComm->master();}

  // Time-related funtions
  inline void saveMDTime(double nowT, double wallt, int iframe){
    FrameNum = iframe; // match the frame index 
    WallClockTime[FrameNum] = wallt;
    MDTime[FrameNum] = nowT;
  }

  inline double getWallClockTime() { return WallClockTime[FrameNum];}
  void resetTimer() { // for each transition, reset parameters
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
  
  void reportTimeData(ostream& ,int itrans);

  void resize(int n);

protected:  

  int NumReplica;
  vector<int> ReplicaID;

  double CumMD_t, CumHMD_t; // cummulative time over replica
  double TotMD_t, TotHMD_t; // total time
  int            FrameNum;
  int LeftNode, RightNode;

  vector<double> WallClockTime;  // wall-clock time
  vector<double> MDTime;    // mdtime
  vector<double> HMDTime; // hyper time will be added
};

#endif
  
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * POOMA_VERSION_ID: $Id$ 
 ***************************************************************************/

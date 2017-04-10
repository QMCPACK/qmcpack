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
    
    
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_TESTING_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTCLE_TESTING_H
#include "QMCDrivers/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

//class DMCPbyPUpdate;
class DMCUpdateBase;

/** @ingroup QMCDrivers ParticleByParticle
 *@brief Implements the DMC algorithm using particle-by-particle move.
 */
class DMCPbyP: public QMCDriver
{
public:
  /// Constructor.
  DMCPbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
          QMCHamiltonian& h);
  ///destructor
  ~DMCPbyP();

  bool run();
  bool put(xmlNodePtr cur);

private:

  ///Index to determine what to do when node crossing is detected
  IndexType KillNodeCrossing;
  ///Total number of accepted steps per block
  IndexType nAcceptTot;
  ///Total number of rejected steps per block
  IndexType nRejectTot;
  ///Interval between branching
  IndexType BranchInterval;
  ///Interval between branching
  IndexType NonLocalMoveIndex;
  ///hdf5 file name for Branch conditions
  std::string BranchInfo;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  ///input std::string to determine to use reconfiguration
  std::string Reconfiguration;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  ///update engine
  DMCUpdateBase *Mover;
  /// Copy Constructor (disabled)
  DMCPbyP(const DMCPbyP& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMCPbyP& operator=(const DMCPbyP&)
  {
    return *this;
  }

  bool dmcWithReconfiguration();
  bool dmcWithBranching();

};
}

#endif

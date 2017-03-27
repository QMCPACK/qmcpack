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
    
    
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class DMCUpdateBase;
class HamiltonianPool;

/** @ingroup QMCDrivers
 *@brief A dummy QMCDriver for testing
 */
class DMCPbyPOMP: public QMCDriver
{
public:

  /// Constructor.
  DMCPbyPOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  void makeClones(HamiltonianPool& hpool, int np=-1);

  bool run();
  bool put(xmlNodePtr cur);

private:
  ///Index to determine what to do when node crossing is detected
  IndexType KillNodeCrossing;
  ///Column index for Populaton
  IndexType PopIndex;
  ///Column index for E_T
  IndexType EtrialIndex;
  ///Total number of accepted steps per block
  IndexType nAcceptTot;
  ///Total number of rejected steps per block
  IndexType nRejectTot;
  ///hdf5 file name for Branch conditions
  std::string BranchInfo;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  /// Copy Constructor (disabled)
  DMCPbyPOMP(const DMCPbyPOMP& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMCPbyPOMP& operator=(const DMCPbyPOMP&)
  {
    return *this;
  }

  int NumThreads;
  std::string Reconfiguration;
  std::string BenchMarkRun;

  std::vector<DMCUpdateBase*> Movers;
  std::vector<ParticleSet*> wClones;
  std::vector<TrialWaveFunction*> psiClones;
  std::vector<QMCHamiltonian*> hClones;
  std::vector<RandomGenerator_t*> Rng;
  std::vector<BranchEngineType*> branchClones;
  std::vector<int> wPerNode;

  void resetRun();
  void benchMark();
  void dmcWithBranching();
  void dmcWithReconfiguration();
};
}

#endif

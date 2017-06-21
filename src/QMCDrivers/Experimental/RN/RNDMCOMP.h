//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file RNDMCOMP.h
 * @brief Define OpenMP-able DMC Driver.
 */
#ifndef QMCPLUSPLUS_RNDMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_RNDMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief DMC driver using OpenMP paragra
 *
 * This is the main DMC driver with MPI/OpenMP loops over the walkers.
 */
class RNDMCOMP: public QMCDriver, public CloneManager
{
public:

  /// Constructor.
  RNDMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h,
           HamiltonianPool& hpool,WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur);
  void setTau(RealType i);
  void resetComponents(xmlNodePtr cur);

private:
  TrialWaveFunction& Guide;

  ///Interval between branching
  IndexType BranchInterval;
  ///hdf5 file name for Branch conditions
  std::string BranchInfo;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  ///input std::string to use Alternate mover
  std::string useAlternate;
  ///input to control maximum age allowed for walkers.
  IndexType mover_MaxAge;
  IndexType myRNWarmupSteps;


  void resetUpdateEngines();
  /// Copy Constructor (disabled)
//     RNDMCOMP(const RNDMCOMP& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
//     RNDMCOMP& operator=(const RNDMCOMP&) { return *this;}

//     //parameters for released node calculation.
//     int Eindex, nSteps;
};
}

#endif

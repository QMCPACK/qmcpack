//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  string BranchInfo;
  ///input string to determine kill walkers or not
  string KillWalker;
  ///input string to determine swap walkers among mpi processors
  string SwapWalkers;
  ///input string to determine to use nonlocal move
  string NonLocalMove;
  ///input string to use Alternate mover
  string useAlternate;
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
/***************************************************************************
 * $RCSfile: RNDMCOMP.h,v $   $Author: jnkim $
 * $Revision: 1604 $   $Date: 2007-01-06 11:00:24 -0600 (Sat, 06 Jan 2007) $
 * $Id: RNDMCOMP.h 1604 2007-01-06 17:00:24Z jnkim $
 ***************************************************************************/

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
/** @file DMCOMP.h
 * @brief Define OpenMP-able DMC Driver.
 */
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief DMC driver using OpenMP paragra
 *
 * This is the main DMC driver with MPI/OpenMP loops over the walkers.
 */
class DMCOMP: public QMCDriver, public CloneManager
{
public:

  /// Constructor.
  DMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
         HamiltonianPool& hpool,WaveFunctionPool& ppool);

  bool run();
  bool put(xmlNodePtr cur);
  void setTau(RealType i);
  void resetComponents(xmlNodePtr cur);

private:
  ///Index to determine what to do when node crossing is detected
  IndexType KillNodeCrossing;
  ///Interval between branching
  IndexType BranchInterval;
  ///hdf5 file name for Branch conditions
  string BranchInfo;
  ///input string to determine kill walkers or not
  string KillWalker;
  ///input string to determine swap walkers among mpi processors
  string SwapWalkers;
  ///input string to determine to use reconfiguration
  string Reconfiguration;
  ///input string to determine to use nonlocal move
  string NonLocalMove;
  ///input string to benchmark OMP performance
  string BenchMarkRun;
  ///input string to use fast gradient
  string UseFastGrad;
  ///input to control maximum age allowed for walkers.
  IndexType mover_MaxAge;


  void resetUpdateEngines();
  void benchMark();
  /// Copy Constructor (disabled)
  DMCOMP(const DMCOMP& a): QMCDriver(a), CloneManager(a) { }
  /// Copy operator (disabled).
  DMCOMP& operator=(const DMCOMP&)
  {
    return *this;
  }
};
}

#endif
/***************************************************************************
 * $RCSfile: DMCOMP.h,v $   $Author: jnkim $
 * $Revision: 1604 $   $Date: 2007-01-06 11:00:24 -0600 (Sat, 06 Jan 2007) $
 * $Id: DMCOMP.h 1604 2007-01-06 17:00:24Z jnkim $
 ***************************************************************************/

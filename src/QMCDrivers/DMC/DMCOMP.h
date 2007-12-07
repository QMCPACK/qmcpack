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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#define QMCPLUSPLUS_DMC_PARTICLEBYPARTICLE_OPNEMP_H
#include "QMCDrivers/QMCDriver.h" 
#include "QMCDrivers/CloneManager.h" 
namespace qmcplusplus {

  class DMCUpdateBase;

  /** @ingroup QMCDrivers 
   *@brief A dummy QMCDriver for testing
   */
  class DMCOMP: public QMCDriver, public CloneManager {
  public:

    /// Constructor.
    DMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
        HamiltonianPool& hpool);

    bool run();
    bool put(xmlNodePtr cur);
 
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

    void resetUpdateEngines();
    void benchMark();
    /// Copy Constructor (disabled)
    DMCOMP(const DMCOMP& a): QMCDriver(a), CloneManager(a) { }
    /// Copy operator (disabled).
    DMCOMP& operator=(const DMCOMP&) { return *this;}
  };
}

#endif
/***************************************************************************
 * $RCSfile: DMCOMP.h,v $   $Author: jnkim $
 * $Revision: 1604 $   $Date: 2007-01-06 11:00:24 -0600 (Sat, 06 Jan 2007) $
 * $Id: DMCOMP.h 1604 2007-01-06 17:00:24Z jnkim $ 
 ***************************************************************************/

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
#ifndef QMCPLUSPLUS_DMC_MASTERDRIVER_H
#define QMCPLUSPLUS_DMC_MASTERDRIVER_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers
 *@brief A dummy QMCDriver for testing
 */
class DMC: public QMCDriver
{
public:

  /// Constructor.
  DMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

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
  ///update engine
  QMCUpdateBase* Mover;
  ///initialize
  void resetUpdateEngine();
  /// Copy Constructor (disabled)
  DMC(const DMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMC& operator=(const DMC&)
  {
    return *this;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile: DMC.h,v $   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: DMC.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/

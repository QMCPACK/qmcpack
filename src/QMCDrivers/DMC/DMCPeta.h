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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DMCPeta_MASTERDRIVER_PETA_H
#define QMCPLUSPLUS_DMCPeta_MASTERDRIVER_PETA_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

class QMCUpdateBase;

/** @ingroup QMCDrivers
 *@brief A dummy QMCDriver for testing
 */
class DMCPeta: public QMCDriver
{
public:

  /// Constructor.
  DMCPeta(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

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
  DMCPeta(const DMCPeta& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMCPeta& operator=(const DMCPeta&)
  {
    return *this;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile: DMCPeta.h,v $   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: DMCPeta.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/

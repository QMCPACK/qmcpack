//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    


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
  std::string BranchInfo;
  ///input std::string to determine kill walkers or not
  std::string KillWalker;
  ///input std::string to determine swap walkers among mpi processors
  std::string SwapWalkers;
  ///input std::string to determine to use reconfiguration
  std::string Reconfiguration;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  ///input std::string to benchmark OMP performance
  std::string BenchMarkRun;
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

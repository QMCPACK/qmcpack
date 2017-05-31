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
    
    
/**@file DMCMoveAll.h
 * @brief Declaration of DMCMoveAll
 */
#ifndef QMCPLUSPLUS_DMC_MOVEALL_H
#define QMCPLUSPLUS_DMC_MOVEALL_H

#include "QMCDrivers/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

class DMCUpdateBase;

/** @ingroup QMCDrivers  WalkerByWalker
 *@brief implements the DMC algorithm using walker-by-walker move.
 */
class DMCMoveAll: public QMCDriver
{

public:
  /// Constructor.
  DMCMoveAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
  ///destructor
  ~DMCMoveAll();

  bool run();
  bool put(xmlNodePtr q);

private:

  DMCUpdateBase *Mover;
  IndexType KillNodeCrossing;
  IndexType BranchInterval;
  ///Interval between branching
  IndexType NonLocalMoveIndex;
  std::string Reconfiguration;
  std::string KillWalker;
  ///input std::string to determine to use nonlocal move
  std::string NonLocalMove;
  /// Copy Constructor (disabled)
  DMCMoveAll(const DMCMoveAll& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMCMoveAll& operator=(const DMCMoveAll&)
  {
    return *this;
  }

  bool dmcWithBranching();
  bool dmcWithReconfiguration();
};
}
#endif

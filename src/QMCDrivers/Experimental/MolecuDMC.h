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
    
    
/**@file MolecuDMC.h
 * @brief Declaration of MolecuDMC
 */
#ifndef QMCPLUSPLUS_DMC_MOLECU_H
#define QMCPLUSPLUS_DMC_MOLECU_H

#include "QMCDrivers/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers  WalkerByWalker
 *@brief implements the DMC algorithm using walker-by-walker move.
 */
class MolecuDMC: public QMCDriver
{

public:
  /// Constructor.
  MolecuDMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
  ///destructor
  ~MolecuDMC();

  template<class BRANCHER> void advanceKillNodeCrossing(BRANCHER& Branch);
  template<class BRANCHER> void advanceRejectNodeCrossing(BRANCHER& Branch);

  bool run();
  bool put(xmlNodePtr q);

  void setBranchInfo(const std::string& fname);

private:

  IndexType KillNodeCrossing;
  std::string KillWalker;
  ///hdf5 file name for Branch conditions
  std::string BranchInfo;
  /// Copy Constructor (disabled)
  MolecuDMC(const MolecuDMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  MolecuDMC& operator=(const MolecuDMC&)
  {
    return *this;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

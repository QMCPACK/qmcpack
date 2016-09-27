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
    
    


#ifndef QMCPLUSPLUS_DMC_MOLECU_WOS_H
#define QMCPLUSPLUS_DMC_MOLECU_WOS_H

#include "QMC/QMCDriver.h"

namespace qmcplusplus
{

class WOSPotential;

/** Implements the DMC algorithm. */
class DMCWOS: public QMCDriver
{

public:
  /// Constructor.
  DMCWOS(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  template<class BRANCHER>
  void advanceWalkerByWalker(BRANCHER& Branch);

  template<class BRANCHER>
  void advanceAllWalkers(BRANCHER& Branch);

  bool run();
  bool put(xmlNodePtr q);

  void setBranchInfo(const std::string& fname);

private:
  /// Copy Constructor (disabled)
  DMCWOS(const DMCWOS& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DMCWOS& operator=(const DMCWOS&)
  {
    return *this;
  }

  /// Time step to work with WOS
  RealType Tau_var;

  ///pointer to WOSPotential
  WOSPotential *wos_ref;

  ///hdf5 file name for Branch conditions
  std::string BranchInfo;

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

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

  void setBranchInfo(const string& fname);

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

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
/**@file MolecuDMC.h
 * @brief Declaration of MolecuDMC
 */
#ifndef QMCPLUSPLUS_DMC_MOLECU_H
#define QMCPLUSPLUS_DMC_MOLECU_H

#include "QMCDrivers/QMCDriver.h" 
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  /** @ingroup QMCDrivers  WalkerByWalker
   *@brief implements the DMC algorithm using walker-by-walker move. 
   */
  class MolecuDMC: public QMCDriver {

  public:
    /// Constructor.
    MolecuDMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
    ///destructor
    ~MolecuDMC();

    template<class BRANCHER> void advanceKillNodeCrossing(BRANCHER& Branch);
    template<class BRANCHER> void advanceRejectNodeCrossing(BRANCHER& Branch);

    bool run();
    bool put(xmlNodePtr q);

    void setBranchInfo(const string& fname);

  private:

    IndexType KillNodeCrossing;
    string KillWalker;
    ///hdf5 file name for Branch conditions
    std::string BranchInfo;
    /// Copy Constructor (disabled)
    MolecuDMC(const MolecuDMC& a): QMCDriver(a) { }
    /// Copy operator (disabled).
    MolecuDMC& operator=(const MolecuDMC&) { return *this;}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

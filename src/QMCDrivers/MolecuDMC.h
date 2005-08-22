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
#ifndef OHMMS_QMC_DMC_MOLECU_H
#define OHMMS_QMC_DMC_MOLECU_H

#include "QMCDrivers/QMCDriver.h" 
#include "Utilities/OhmmsInfo.h"
//#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/MolecuFixedNodeBranch.h"

namespace ohmmsqmc {

  /** @ingroup QMCDrivers  WalkerByWalker
   *@brief implements the DMC algorithm using walker-by-walker move. 
   */
  class MolecuDMC: public QMCDriver {

  public:
    typedef MolecuFixedNodeBranch<RealType> BranchEngineType;
    //typedef SimpleFixedNodeBranch<RealType> BranchEngineType;
    /// Constructor.
    MolecuDMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
    ///destructor
    ~MolecuDMC();

    template<class BRANCHER>
    void advanceWalkerByWalker(BRANCHER& Branch);

    bool run();
    bool put(xmlNodePtr q);

    void setBranchInfo(const string& fname);

  private:

    ///hdf5 file name for Branch conditions
    std::string BranchInfo;
    ///branch engine
    BranchEngineType *branchEngine;
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

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
#ifndef OHMMS_QMC_MOLECUAPPS_H
#define OHMMS_QMC_MOLECUAPPS_H

#include "QMC/QMCApps.h"

namespace ohmmsqmc {

  /**
   @brief An application to perform QMC simulations of molecules.
   *
   *Molecules consist of fixed centers (ions) and quantum particles
   *(electrons).  The TrialWaveFunction and QMCHamiltonian are 
   *defined accordingly.
   *Use xml input files in Data/HFAtomicSTO
   */
  class MolecuApps: public QMCApps {

  public:

    ///constructor
    MolecuApps(int argc, char** argv);

    ///destructor
    ~MolecuApps();

    ///initialization with a file
    bool init();

  private:
    bool setParticleSets(xmlNodePtr aroot);
    bool setWavefunctions(xmlNodePtr aroot);
    bool setHamiltonian(xmlNodePtr aroot);
    bool setSimulation(xmlNodePtr aroot);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

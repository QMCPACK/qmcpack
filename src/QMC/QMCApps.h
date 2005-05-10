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
#ifndef OHMMS_QMC_QMCAPPS_H
#define OHMMS_QMC_QMCAPPS_H
#include "OhmmsApp/ProjectData.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /** Main class to run the entire process. */
  class QMCApps {

  public:

    ///constructor: initialize the seed for the random number generator
    QMCApps(int argc, char** argv);

    ///destructor
    virtual ~QMCApps();

    ///initialization with a file
    virtual bool init() = 0;

    ///run it
    bool run(xmlNodePtr aroot);

    void saveXml();

  protected:

    xmlDocPtr m_doc;
    xmlNodePtr m_root;
    xmlXPathContextPtr m_context;
    vector<xmlNodePtr> m_walkerset;

    ///project description
    OHMMS::ProjectData myProject;

    ///random number controller
    OHMMS::RandomNumberControl myRandomControl;

    ///name of file that stores configurations (walkers)
    string PrevConfigFile;

    ///the walker ensemble
    MCWalkerConfiguration el;

    ///the ions
    ParticleSet ion;

    ///the wavefunction
    TrialWaveFunction Psi;

    ///the Hamiltonian
    QMCHamiltonian H;

    bool setMCWalkers(xmlNodePtr aroot);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

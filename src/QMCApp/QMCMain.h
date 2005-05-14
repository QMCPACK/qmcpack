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
#ifndef OHMMS_QMC_MAINAPPLICATIONS_H
#define OHMMS_QMC_MAINAPPLICATIONS_H

#include "QMCApp/QMCAppBase.h"

namespace ohmmsqmc {

  class MCWalkerConfiguration;
  class QMCDriver;
  class ParticleSetPool;
  class WaveFunctionPool;
  class HamiltonianPool;

  /** Main application to perform QMC simulations 
   *
   * This is a generalized QMC application which can handle multiple particle sets,
   * walker configurations, hamiltonians and wavefunctions
   */
  class QMCMain: public QMCAppBase {

  public:

    ///constructor
    QMCMain(int argc, char** argv);

    ///destructor
    ~QMCMain();

    bool validateXML();
    bool execute();

  private:

    std::string curMethod;
    QMCDriver *qmcDriver;
    MCWalkerConfiguration *qmcSystem;

    ParticleSetPool* ptclPool;
    WaveFunctionPool* psiPool;
    HamiltonianPool* hamPool;

    ///previous configuration file for next qmc node
    string PrevConfigFile;

    ///name of file that stores configurations (walkers)
    vector<xmlNodePtr> m_walkerset;

    bool runQMC(xmlNodePtr cur);
    bool setMCWalkers(xmlXPathContextPtr cur);
    void processContext(xmlXPathContextPtr context_);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

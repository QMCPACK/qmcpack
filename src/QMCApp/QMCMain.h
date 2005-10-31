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
/**@file QMCMain.h
 * @brief Declaration of QMCMain class
 */
#ifndef QMCPLUSPLUS_MAINAPPLICATIONS_H
#define QMCPLUSPLUS_MAINAPPLICATIONS_H

#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus {

  //forward declaration
  class MCWalkerConfiguration;
  class QMCDriver;
  class ParticleSetPool;
  class WaveFunctionPool;
  class HamiltonianPool;

  /** @ingroup qmcapp
   * @brief Main application to perform QMC simulations 
   *
   * This is a generalized QMC application which can handle multiple ParticleSet,
   * TrialWaveFunction and QMCHamiltonian objects.
   */
  class QMCMain: public QMCAppBase {

  public:

    /*! enum for QMC Run Type */
    enum QMCRunType {DUMMY_RUN, /*!< dummy */
      VMC_RUN, /*!< VMC type: vmc, vmc-ptcl, vmc-multiple, vmc-ptcl-multiple */
      DMC_RUN, /*!< DMC type: dmc, dmc-ptcl*/
      RMC_RUN, /*!< RMC type: rmc, rmc-ptcl */
      OPTIMIZE_RUN /*!< Optimization */
    };

    ///constructor
    QMCMain(int argc, char** argv);

    ///destructor
    ~QMCMain();

    bool validateXML();
    bool execute();

  private:

    ///type of 
    QMCRunType curRunType;

    ///name of the current QMCriver
    std::string curMethod;

    ///flag to indicate that a qmc is the first QMC
    bool FirstQMC;

    /** current QMCDriver
     */
    QMCDriver *qmcDriver;

    /** current MCWalkerConfiguration
     */
    MCWalkerConfiguration *qmcSystem;

    /** ParticleSet Pool
     */
    ParticleSetPool* ptclPool;

    /** TrialWaveFunction Pool
     */
    WaveFunctionPool* psiPool;

    /** QMCHamiltonian Pool
     */
    HamiltonianPool* hamPool;

    ///previous configuration file for next qmc node
    string PrevConfigFile;

    ///xml mcwalkerset elements for output
    vector<xmlNodePtr> m_walkerset;
    ///xml mcwalkerset read-in elements 
    vector<xmlNodePtr> m_walkerset_in;

    ///execute <qmc/> element
    bool runQMC(xmlNodePtr cur);

    ///add <mcwalkerset/> elements to continue a run
    bool setMCWalkers(xmlXPathContextPtr cur);

    /** add unique ParticleSet, TrialWaveFunction and QMCHamiltonian elements to Pool objects
     */
    void processPWH(xmlNodePtr cur);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

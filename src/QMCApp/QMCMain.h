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

#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus {

  /** @ingroup qmcapp
   * @brief Main application to perform QMC simulations 
   *
   * This is a generalized QMC application which can handle multiple ParticleSet,
   * TrialWaveFunction and QMCHamiltonian objects.
   */
  class QMCMain: public QMCDriverFactory, 
                 public QMCAppBase {

  public:

    ///constructor
    QMCMain(Communicate* c);

    ///destructor
    ~QMCMain();

    bool validateXML();
    bool execute();

  private:

    ///flag to indicate that a qmc is the first QMC
    bool FirstQMC;

    ///previous configuration file for next qmc node
    string PrevConfigFile;

    ///xml mcwalkerset elements for output
    vector<xmlNodePtr> m_walkerset;
    ///xml mcwalkerset read-in elements 
    vector<xmlNodePtr> m_walkerset_in;
    ///qmc sections
    vector<pair<xmlNodePtr,bool> > m_qmcaction;
    ///pointer to the last node of the main inputfile
    xmlNodePtr lastInputNode;
    ///execute <qmc/> element
    bool runQMC(xmlNodePtr cur);

    ///add <mcwalkerset/> elements to continue a run
    bool setMCWalkers(xmlXPathContextPtr cur);

    /** add unique ParticleSet, TrialWaveFunction and QMCHamiltonian elements to Pool objects
     * @param cur xml node
     * @return false, if contains qmc actions
     */
    bool processPWH(xmlNodePtr cur);

    /** execute loop **/
    void executeLoop(xmlNodePtr cur);
    /** execute qmc 
     * @param cur qmc xml node
     * @param noloop if true, this qmc section is not in a loop.
     * @return true, if a section is successfully executed.
     */
    bool executeQMCSection(xmlNodePtr cur, bool noloop=true);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

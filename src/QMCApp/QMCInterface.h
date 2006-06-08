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
/**@file QMCInterface.h
 * @brief Declaration of QMCInterface class
 */
#ifndef QMCPLUSPLUS_INTERFACEAPPLICATIONS_H
#define QMCPLUSPLUS_INTERFACEAPPLICATIONS_H

#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/QMCAppBase.h"

namespace qmcplusplus {

  /** @ingroup qmcapp
   * @brief Interface application to perform QMC simulations 
   *
   * This is a generalized QMC application which can handle multiple ParticleSet,
   * TrialWaveFunction and QMCHamiltonian objects.
   */
  class QMCInterface: public QMCDriverFactory, 
                 public QMCAppBase {

  public:

    ///constructor
    QMCInterface(int argc, char** argv);

    ///destructor
    ~QMCInterface();

    bool validateXML();
    bool initialize();
    bool SetVMC(int nblocks);
    bool SetVMCMultiple(int nblocks);
		bool process();
    bool execute();
    void SetPtclPos(int id, double* newR);
    void SetPtclPos(string set, int id, double* newR);
//    vector<double>* GetData();
//		double GetData(string estimator, string tag);
  private:

    /// pointer for xmlNode to be parsed in qmcDriver->process()
    xmlNodePtr runInfoNode;

    ///flag to indicate that a qmc is the first QMC
    bool FirstQMC;

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

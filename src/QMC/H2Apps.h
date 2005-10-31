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
#ifndef QMCPLUSPLUS_H2APP_H
#define QMCPLUSPLUS_H2APP_H
#include "QMC/QMCApps.h"
namespace qmcplusplus {

  class H2Apps: public QMCApps {

  public:

    ///constructor
    H2Apps(int argc, char** argv): QMCApps(argc,argv) { }

    ///initialization with a file
    bool init();

  protected:
    bool setParticleSets(xmlpp::Node*);
    bool setElectrons(xmlpp::Node*);
    bool setIons(xmlpp::Node*);
    bool setWavefunctions(xmlpp::Node*);
    bool setHamiltonian(xmlpp::Node*);  
    bool setSimulation(xmlpp::Node*);  
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

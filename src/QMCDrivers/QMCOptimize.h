//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_GENERIC_OPTIMIZATION_H
#define QMCPLUSPLUS_GENERIC_OPTIMIZATION_H

#include "QMCDrivers/QMCDriver.h" 
#include "QMCDrivers/QMCCostFunction.h"

namespace qmcplusplus {

  /** @ingroup QMCDrivers
   * @brief Implements wave-function optimization
   *
   * Optimization by correlated sampling method with configurations 
   * generated from VMC.
   */

  class QMCOptimize: public QMCDriver
  {
  public:

    ///Constructor.
    QMCOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
	    QMCHamiltonian& h);
    
    ///Destructor
    ~QMCOptimize();

    ///Run the Optimization algorithm.
    bool run();
    ///process xml node
    bool put(xmlNodePtr cur);
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);

    void setWaveFunctionNode(xmlNodePtr cur) { wfNode=cur; }

  private:

    ///index to denote the partition id
    int PartID;
    ///total number of partitions that will share a set of configuratons
    int NumParts;
    ///target cost function to optimize
    QMCCostFunction* optTarget;
    ///xml node to be dumped
    xmlNodePtr wfNode;
    ///xml node for optimizer
    xmlNodePtr optNode;
    ///method for optimization, default conjugate gradient
    string optmethod;
    ///list of files storing configurations  
    vector<string> ConfigFile;

    ///Copy Constructor (disabled).
    QMCOptimize(const QMCOptimize& a): QMCDriver(a) { }  
    ///Copy operator (disabled).
    QMCOptimize& operator=(const QMCOptimize&) { return *this;}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Jordan Vincent
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
#ifndef OHMMS_QMC_VMC_OPT_H
#define OHMMS_QMC_VMC_OPT_H

#include <deque>
#include "Configuration.h"
#include "QMC/QMCDriver.h" 
#include "Optimize/Minimize.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /** Performs optimization of the wavefunction using 
   *the correlated sampling method with configurations 
   *generated from VMC.
   */

  class VMC_OPT: public QMCDriver, public MinimizeFunction 
  {
  public:
    enum {ENERGY, ENERGYSQ};

    ///Constructor.
    VMC_OPT(MCWalkerConfiguration& w, 
	    TrialWaveFunction& psi, 
	    QMCHamiltonian& h, 
	    xmlNodePtr q);
    
    ///Run the Optimization algorithm.
    bool run();
    bool put(xmlNodePtr cur);
    void run_vmc();
    ValueType correlatedSampling();
    ///assign optimization parameter i
    scalar& Params(int i) { return OptParams[i]; }
    ///return optimization parameter i
    scalar Params(int i) const { return OptParams[i]; }
    scalar Cost();

    ///return the number of optimizable parameters
    int NumParams() { return OptParams.size(); }
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);

    void WriteStuff() {}

  private:

    bool resetWaveFunctions();
    bool checkParameters();
    bool putOptParams();  

    RealType evalCost();

    ///storage for previous values of the cost function
    std::deque<scalar> costList;
    ///storage for previous sets of parameters
    std::deque<vector<scalar> > paramList;
    ///parameters to be optimized
    vector<scalar> OptParams;
    ///ID tag for each optimizable parameter
    vector<string> IDtag;  
    ///list of files storing configurations  
    vector<string> ConfigFile;
    ///method for optimization, default conjugate gradient
    string optmethod;
    ///number of times cost function evaluated
    int NumCostCalls;
    ///total number of samples to use in correlated sampling
    int NumSamples;
    ///conjugate gradient tolerance
    RealType cg_tolerance;
    ///conjugate gradient stepsize
    RealType cg_stepsize;
    ///conjugate gradient epsilon
    RealType cg_epsilon;
    ///weights for energy and variance in the cost function
    RealType w_en, w_var;
    ///value of the cost function
    RealType CostValue;
    ///Hamiltonians that depend on the optimization: KE
    QMCHamiltonian H_KE;
    ///Copy Constructor (disabled).
    VMC_OPT(const VMC_OPT& a): QMCDriver(a) { }  
    ///Copy operator (disabled).
    VMC_OPT& operator=(const VMC_OPT&) { return *this;}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

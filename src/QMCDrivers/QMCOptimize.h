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
#ifndef OHMMS_QMC_QMCOptimize_H
#define OHMMS_QMC_QMCOptimize_H

#include <deque>
#include "Configuration.h"
#include "QMCDrivers/QMCDriver.h" 
#include "Optimize/Minimize.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /** @ingroup QMCDrivers
   * @brief Implements wave-function optimization
   *
   * Optimization by correlated sampling method with configurations 
   * generated from VMC.
   */

  class QMCOptimize: public QMCDriver, public MinimizeFunction 
  {
  public:

    enum FieldIndex_OPT {LOGPSI_FIXED=0, LOGPSI_FREE=1, ENERGY_TOT=2, ENERGY_FIXED=3, ENERGY_NEW=4, REWEIGHT=5};
    enum SumIndex_OPT {SUM_E_BARE, SUM_ESQ_BARE, SUM_ABSE_BARE,
      SUM_E_WGT, SUM_ESQ_WGT, SUM_ABSE_WGT, SUM_WGT, SUM_WGTSQ};
    //enum {ENERGY, ENERGYSQ};

    ///Constructor.
    QMCOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, 
	    QMCHamiltonian& h);
    
    ///Destructor
    ~QMCOptimize();

    ///Run the Optimization algorithm.
    bool run();
    bool put(xmlNodePtr cur);
    RealType correlatedSampling();
    ///assign optimization parameter i
    scalar& Params(int i) { return OptParams[i]; }
    ///return optimization parameter i
    scalar Params(int i) const { return OptParams[i]; }
    scalar Cost();

    ///return the number of optimizable parameters
    int NumParams() { return OptParams.size(); }
    ///add a configuration file to the list of files
    void addConfiguration(const string& a);

    void WriteStuff();

  private:

    ///boolean to turn on/off the psi^2/psi^2_old for correlated sampling
    bool UseWeight;
    ///bollean to turn on/off the use of anonymous buffer
    bool needBuffering;
    ///bollean to turn on/off the use of anonymous buffer for the ratio
    bool hamiltonianNeedRatio;
    ///index to denote the partition id
    int PartID;
    ///total number of partitions that will share a set of configuratons
    int NumParts;
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
    RealType w_en, w_var, w_abs;
    ///value of the cost function
    RealType CostValue;
    ///target energy
    RealType Etarget;
    ///Hamiltonians that depend on the optimization: KE
    QMCHamiltonian H_KE;

    /** Sum of energies and weights for averages
     *
     * SumValues[k] where k is one of SumIndex_opt
     */
    TinyVector<RealType,8> SumValue;
    /** Saved properties of all the walkers
     *
     * Records(iw,field_id) returns the field_id value of the iw-th walker
     * field_id is one of FieldIndex_opt
     */
    Matrix<ValueType> Records;
    /** Fixed  Laplacian , \f$\nabla^2\ln\Psi\f$, components */
    ParticleSet::ParticleLaplacian_t dL;

    void checkConfigurations();
    bool resetWaveFunctions();
    bool checkParameters();
    bool putOptParams();  
    RealType evalCost();

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

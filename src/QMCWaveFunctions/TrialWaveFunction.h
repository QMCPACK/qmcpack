//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_QMC_TRIALWAVEFUNCTION_H
#define OHMMS_QMC_TRIALWAVEFUNCTION_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "Optimize/VarList.h"

/**@file TrialWaveFunction.h
 *@brief Declaration of a TrialWaveFunction 
 */
namespace ohmmsqmc {

  /** class for many-body trial wave function 
   *
   *A many-body trial wave function is represented by
   *\f[
   *\Psi({\bf R}) = \prod_i \Theta_i({\bf R}),
   *\f]
   *where each function \f$\Theta_i({\bf R})\f$ is an OrbitalBase.
   *A Composite Pattern is used to handle \f$\prod\f$ operations.
   *Each OrbitalBase should provide proper evaluate functions
   *for the value, gradient and laplacian values.
   *
   *@todo TrialWaveFunction should be a derived class of VarRegistry<RealType>
   */
  class TrialWaveFunction {

  public:

    typedef OrbitalBase::RealType  RealType;
    typedef OrbitalBase::ValueType ValueType;
    typedef OrbitalBase::PosType   PosType;
    typedef OrbitalBase::GradType  GradType;

    TrialWaveFunction();

    ~TrialWaveFunction();

    ValueType evaluate(ParticleSet& P);

    void evaluate(WalkerSetRef& W, OrbitalBase::ValueVectorType& psi);

    void add(OrbitalBase* aterm);

    void reset();

    ///resize the internal storage with number of walkers if necessary
    void resizeByWalkers(int nwalkers);

    /**a list of real variables to be optimized
     * Each builder for a trial wavefuncion is responsible for registering
     * the variables to be optimized.
     */
    VarRegistry<RealType> VarList;

  private:

    ///cannot use copy constructor
    TrialWaveFunction(const TrialWaveFunction&) {}
    
    ///a list of OrbitalBases constituting many-body wave functions
    vector<OrbitalBase*> Z;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


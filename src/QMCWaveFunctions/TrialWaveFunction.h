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

    typedef OrbitalBase::RealType   RealType;
    typedef OrbitalBase::ValueType  ValueType;
    typedef OrbitalBase::PosType    PosType;
    typedef OrbitalBase::GradType   GradType;

    /**a list of real variables to be optimized
     *
     * Each builder for a trial wavefuncion is responsible for registering
     * the variables to be optimized.
     */
    VarRegistry<RealType> VarList;

    TrialWaveFunction();

    ~TrialWaveFunction();

    inline RealType getSign() const { return SignValue;}
    inline ValueType getLogPsi() const { return LogValue;}

    ///Add an OrbitalBase 
    void add(OrbitalBase* aterm);
    ///reset OrbitalBase s during optimization
    void reset();
    /**resize the internal storage
     * @param nwalkers number of walkers 
     *
     * Not used anymore
     */
    void resizeByWalkers(int nwalkers);

    ///Check if aname-ed Single-Particle-Orbital set exists
    bool hasSPOSet(const string& aname);
    ///add a Single-Particle-Orbital set
    void addSPOSet(OhmmsElementBase* spo);
    ///return the aname-ed Single-Particle-Orbital set. 
    OhmmsElementBase* getSPOSet(const string& aname);

    /** evalaute the values of the wavefunction, gradient and laplacian  for a walkers */
    ValueType evaluate(ParticleSet& P);

    /** evalaute the log of the trial wave function */
    ValueType evaluateLog(ParticleSet& P);
    ValueType evaluateLog(ParticleSet& P, bool all);

    /** evalaute the log of the trial wave function */
    ValueType evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& fixedG,
        ParticleSet::ParticleLaplacian_t& fixedL);

    /** functions to handle particle-by-particle update */
    ValueType ratio(ParticleSet& P, int iat);
    void update(ParticleSet& P, int iat);

    ValueType ratio(ParticleSet& P, int iat, 
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

    void restore(int iat);
    void update2(ParticleSet& P, int iat);

    void registerData(ParticleSet& P, PooledData<RealType>& buf);
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf);

    /** evalaute the values of the wavefunction, gradient and laplacian  for all the walkers */
    void evaluate(WalkerSetRef& W, OrbitalBase::ValueVectorType& psi);

  private:

    ///the size of gradient component (QMCTraits::DIM)*the number of particles 
    int TotalDim;

    ///index of the active particle
    int WorkingPtcl;

    ///sign of the trial wave function
    RealType SignValue;

    ///log of the trial wave function
    ValueType LogValue;

    ///a list of OrbitalBases constituting many-body wave functions
    vector<OrbitalBase*> Z;

    ///a list of single-particle-orbital set
    vector<OhmmsElementBase*> SPOSet;

    ///differential gradients
    ParticleSet::ParticleGradient_t delta_G;

    ///differential laplacians
    ParticleSet::ParticleLaplacian_t delta_L;

    ///cannot use copy constructor
    TrialWaveFunction(const TrialWaveFunction&) {}
    
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


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
/**@file TrialWaveFunction.h
 *@brief Declaration of a TrialWaveFunction 
 */
#ifndef OHMMS_QMC_TRIALWAVEFUNCTION_H
#define OHMMS_QMC_TRIALWAVEFUNCTION_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "Optimize/VarList.h"
/**@defgroup MBWfs Many-body wave function group
 * @brief Classes to handle many-body trial wave functions
 */

namespace ohmmsqmc {

  /** @ingroup MBWfs
   * @brief Class to represent a many-body trial wave function 
   *
   *A many-body trial wave function is represented by
   *\f[
   *\Psi({\bf R}) = \prod_i \psi_i({\bf R}),
   *\f]
   *where each function \f$\psi_i({\bf R})\f$ is an OrbitalBase
   (see OrbitalComponent).
   *A Composite Pattern is used to handle \f$\prod\f$ operations.
   *Each OrbitalBase should provide proper evaluate functions
   *for the value, gradient and laplacian values.
   */
  class TrialWaveFunction {

  public:

    typedef OrbitalBase::RealType   RealType;
    typedef OrbitalBase::ValueType  ValueType;
    typedef OrbitalBase::PosType    PosType;
    typedef OrbitalBase::GradType   GradType;
    typedef OrbitalBase::BufferType BufferType;

    /**a list of real variables to be optimized
     *
     * Each builder for a trial wavefuncion is responsible for registering
     * the variables to be optimized.
     */
    VarRegistry<RealType> VarList;

    ///differential gradients
    ParticleSet::ParticleGradient_t G;

    ///differential laplacians
    ParticleSet::ParticleLaplacian_t L;

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

    ValueType registerData(ParticleSet& P, BufferType& buf);
    void copyFromBuffer(ParticleSet& P, BufferType& buf);
    ValueType evaluate(ParticleSet& P, BufferType& buf);

    void dumpToBuffer(ParticleSet& P, BufferType& buf);
    void dumpFromBuffer(ParticleSet& P, BufferType& buf);
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
  /**@}*/
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


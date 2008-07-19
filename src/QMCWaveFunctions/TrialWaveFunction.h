/////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_H

#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Utilities/NewTimer.h"
/**@defgroup MBWfs Many-body wave function group
 * @brief Classes to handle many-body trial wave functions
 */

namespace qmcplusplus {

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
  class TrialWaveFunction: public MPIObjectBase, public OhmmsElementBase 
  {

  public:

    typedef OrbitalBase::RealType           RealType;
    typedef OrbitalBase::ValueType          ValueType;
    typedef OrbitalBase::PosType            PosType;
    typedef OrbitalBase::GradType           GradType;
    typedef OrbitalBase::BufferType         BufferType;
    typedef OrbitalBase::OptimizableSetType OptimizableSetType;

    ///differential gradients
    ParticleSet::ParticleGradient_t G;
    ///differential laplacians
    ParticleSet::ParticleLaplacian_t L;
    /**a list of real variables to be optimized
     *
     * Each builder for a trial wavefuncion is responsible for registering
     * the variables to be optimized.
     */
    OptimizableSetType VarList;

    TrialWaveFunction(Communicate* c);

    ~TrialWaveFunction();

    inline int size() const { return Z.size();}
    inline RealType getPhase() const { return PhaseValue;}
    inline RealType getLogPsi() const { return LogValue;}

    ///Add an OrbitalBase 
    //void addOrbital(OrbitalBase* aterm);
    void addOrbital(OrbitalBase* aterm, const string& aname);

    ///write to ostream
    bool get(std::ostream& ) const;

    ///read from istream
    bool put(std::istream& );

    ///read from xmlNode
    bool put(xmlNodePtr cur);
    ///implement the virtual function
    void reset();
    ///reset member data
    void resetParameters(OptimizableSetType& optVariables);

    /**resize the internal storage
     * @param nwalkers number of walkers 
     *
     * Not used anymore
     */
    void resizeByWalkers(int nwalkers);

    /** recursively change the ParticleSet whose G and L are evaluated */
    void resetTargetParticleSet(ParticleSet& P);

    /////Check if aname-ed Single-Particle-Orbital set exists
    //bool hasSPOSet(const string& aname);
    /////add a Single-Particle-Orbital set
    //void addSPOSet(OhmmsElementBase* spo);
    /////return the aname-ed Single-Particle-Orbital set. 
    //OhmmsElementBase* getSPOSet(const string& aname);

    /** evalaute the values of the wavefunction, gradient and laplacian  for a walkers */
    ValueType evaluate(ParticleSet& P);

    /** evalaute the log of the trial wave function */
    RealType evaluateLog(ParticleSet& P);

    RealType evaluateDeltaLog(ParticleSet& P);

    void evaluateDeltaLog(ParticleSet& P, 
        RealType& logpsi_fixed,
        RealType& logpsi_opt,
        ParticleSet::ParticleGradient_t& fixedG,
        ParticleSet::ParticleLaplacian_t& fixedL);

    /** functions to handle particle-by-particle update */
    RealType ratio(ParticleSet& P, int iat);
    void update(ParticleSet& P, int iat);

    RealType ratio(ParticleSet& P, int iat, 
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL);

//    RealType logRatio(ParticleSet& P, int iat, 
//		    ParticleSet::ParticleGradient_t& dG,
//		    ParticleSet::ParticleLaplacian_t& dL);


    void rejectMove(int iat);
    void acceptMove(ParticleSet& P, int iat);

    RealType registerData(ParticleSet& P, BufferType& buf);
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);
    void copyFromBuffer(ParticleSet& P, BufferType& buf);
    RealType evaluate(ParticleSet& P, BufferType& buf);

    void dumpToBuffer(ParticleSet& P, BufferType& buf);
    void dumpFromBuffer(ParticleSet& P, BufferType& buf);

    /** evalaute the values of the wavefunction, gradient and laplacian  for all the walkers */
    //void evaluate(WalkerSetRef& W, OrbitalBase::ValueVectorType& psi);

    void reverse();

    TrialWaveFunction* makeClone(ParticleSet& tqp) const;

  private:

    ///control how ratio is calculated
    bool Ordered;
    ///the size of ParticleSet
    int NumPtcls;

    ///the size of gradient component (QMCTraits::DIM)*the number of particles 
    int TotalDim;

    ///index of the active particle
    int WorkingPtcl;

    ///sign of the trial wave function
    RealType PhaseValue;

    ///log of the trial wave function
    RealType LogValue;

    ///a list of OrbitalBases constituting many-body wave functions
    vector<OrbitalBase*> Z;

    ///a list of single-particle-orbital set
    //vector<OhmmsElementBase*> SPOSet;

    ///differential gradients
    ParticleSet::ParticleGradient_t delta_G;

    ///differential laplacians
    ParticleSet::ParticleLaplacian_t delta_L;

    TrialWaveFunction();

    vector<NewTimer*> myTimers;
    
  };
  /**@}*/
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


/////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
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
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Utilities/NewTimer.h"
/**@defgroup MBWfs Many-body wave function group
 * @brief Classes to handle many-body trial wave functions
 */

namespace qmcplusplus
  {

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
  class TrialWaveFunction: public MPIObjectBase
    {

    public:

      typedef OrbitalBase::RealType           RealType;
      typedef OrbitalBase::ValueType          ValueType;
      typedef OrbitalBase::PosType            PosType;
      typedef OrbitalBase::GradType           GradType;
      typedef OrbitalBase::BufferType         BufferType;

      ///differential gradients
      ParticleSet::ParticleGradient_t G;
      ///differential laplacians
      ParticleSet::ParticleLaplacian_t L;

      TrialWaveFunction(Communicate* c);

      ~TrialWaveFunction();

      inline int size() const
        {
          return Z.size();
        }
      inline RealType getPhase() const
        {
          return PhaseValue;
        }
      inline RealType getPhaseDiff() const
        {
          return PhaseDiff;
        }
      inline void resetPhaseDiff()
      {
        PhaseDiff=0.0;
      }
      inline RealType getLogPsi() const
        {
          return LogValue;
        }

      ///Add an OrbitalBase
      //void addOrbital(OrbitalBase* aterm);
      void addOrbital(OrbitalBase* aterm, const string& aname);

      ///read from xmlNode
      bool put(xmlNodePtr cur);
      ///implement the virtual function
      void reset();
      /** set OrbitalBase::IsOptimizing to true */
      void startOptimization();
      /** set OrbitalBase::IsOptimizing to flase */
      void stopOptimization();
      /** check in an optimizable parameter
       * * @param o a super set of optimizable variables
       *
       * Update myOptIndex if o is found among the "active" paramemters.
       */
      void checkInVariables(opt_variables_type& o);
      /** check out optimizable variables
       */
      void checkOutVariables(const opt_variables_type& o);
      ///reset member data
      void resetParameters(const opt_variables_type& active);
      /** print out state of the trial wavefunction
       */
      void reportStatus(ostream& os);

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

      /** evalaute the values of the wavefunction, gradient and laplacian  for a walkers */
      RealType evaluateLogOnly(ParticleSet& P);

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

      /** for releasednode calculations */
      RealType alternateRatio(ParticleSet& P);

      void update(ParticleSet& P, int iat);

      RealType ratio(ParticleSet& P, int iat,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL);

      GradType evalGrad(ParticleSet& P, int iat);
      /** Returns the logarithmic gradient of the trial wave function
       *  with respect to the iat^th atom of the source ParticleSet. */
      GradType evalGradSource(ParticleSet& P, ParticleSet &source, int iat);
      /** Returns the logarithmic gradient of the w.r.t. the iat^th atom
       * of the source ParticleSet of the sum of laplacians w.r.t. the
       * electrons (target ParticleSet) of the trial wave function. */
      GradType evalGradSource
      (ParticleSet& P, ParticleSet &source, int iat,
       TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
       TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

      RealType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);


      void rejectMove(int iat);
      void acceptMove(ParticleSet& P, int iat);

      RealType registerData(ParticleSet& P, BufferType& buf);
      RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);
      void copyFromBuffer(ParticleSet& P, BufferType& buf);
      RealType evaluateLog(ParticleSet& P, BufferType& buf);

      void dumpToBuffer(ParticleSet& P, BufferType& buf);
      void dumpFromBuffer(ParticleSet& P, BufferType& buf);

      RealType KECorrection() const;

      void evaluateDerivatives(ParticleSet& P,
                               const opt_variables_type& optvars,
                               vector<RealType>& dlogpsi,
                               vector<RealType>& dhpsioverpsi);

      /** evalaute the values of the wavefunction, gradient and laplacian  for all the walkers */
      //void evaluate(WalkerSetRef& W, OrbitalBase::ValueVectorType& psi);

      void reverse();

      void resizeTempP(ParticleSet& P)
      {
        tempP = new ParticleSet(P);
      };

      TrialWaveFunction* makeClone(ParticleSet& tqp) const;

      void setMassTerm(ParticleSet& P)
      {
        SpeciesSet tspecies(P.getSpeciesSet());
        int massind=tspecies.addAttribute("mass");
        RealType mass = tspecies(massind,0);
        OneOverM = 1.0/mass;
      }

      vector<OrbitalBase*>& getOrbitals()
      {
        return Z;
      }

      void get_ratios(ParticleSet& P, vector<ValueType>& ratios);
      void setTwist(vector<RealType> t)
      {
        myTwist=t;
      }
      const vector<RealType> twist()
      {
        return myTwist;
      }

    private:

      ///control how ratio is calculated
      bool Ordered;
      ///the size of ParticleSet
      int NumPtcls;

      ///the size of gradient component (QMCTraits::DIM)*the number of particles
      int TotalDim;

      ///index of the active particle
      int WorkingPtcl;

      ///starting index of the buffer
      size_t BufferCursor;

      ///sign of the trial wave function
      RealType PhaseValue;

      ///diff of the phase of the trial wave function during ratio calls
      RealType PhaseDiff;

      ///log of the trial wave function
      RealType LogValue;

      ///One over mass of target particleset, needed for Local Energy Derivatives
      RealType OneOverM;

      ///a list of OrbitalBases constituting many-body wave functions
      vector<OrbitalBase*> Z;

      ///a list of single-particle-orbital set
      //vector<OhmmsElementBase*> SPOSet;

      ///differential gradients
      ParticleSet::ParticleGradient_t delta_G;

      ///differential laplacians
      ParticleSet::ParticleLaplacian_t delta_L;

      ///fake particleset
      ParticleSet* tempP;

      TrialWaveFunction();

      vector<NewTimer*> myTimers;
      vector<RealType> myTwist;
    };
  /**@}*/
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

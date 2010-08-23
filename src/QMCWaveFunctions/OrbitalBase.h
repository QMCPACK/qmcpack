//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ORBITALBASE_H
#define QMCPLUSPLUS_ORBITALBASE_H
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/RecordProperty.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Particle/MCWalkerConfiguration.h"
#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

/**@file OrbitalBase.h
 *@brief Declaration of OrbitalBase
 */
namespace qmcplusplus
  {

#ifdef QMC_CUDA
    struct NLjob {
    int walker;
    int elec;
    int numQuadPoints;
    NLjob(int w, int e, int n) :
      walker(w), elec(e), numQuadPoints(n)
    { }
  };
#endif

  ///forward declaration of OrbitalBase
  class OrbitalBase;
  ///forward declaration of DiffOrbitalBase
  class DiffOrbitalBase;

#if defined(ENABLE_SMARTPOINTER)
  typedef boost::shared_ptr<OrbitalBase>     OrbitalBasePtr;
  typedef boost::shared_ptr<DiffOrbitalBase> DiffOrbitalBasePtr;
#else
  typedef OrbitalBase*                       OrbitalBasePtr;
  typedef DiffOrbitalBase*                   DiffOrbitalBasePtr;
#endif

  /**@defgroup OrbitalComponent Orbital group
   * @brief Classes which constitute a many-body trial wave function
   *
   * A many-body trial wave function is
   * \f[
   * \Psi(\{ {\bf R}\}) = \prod_i \psi_{i}(\{ {\bf R}\}),
   * \f]
   * where \f$\psi\f$s are represented by
   * the derived classes from OrbtialBase.
   */
  /** @ingroup OrbitalComponent
   * @brief An abstract class for a component of a many-body trial wave function
   */
  struct OrbitalBase: public QMCTraits
    {

      ///recasting enum of DistanceTableData to maintain consistency
      enum {SourceIndex  = DistanceTableData::SourceIndex,
            VisitorIndex = DistanceTableData::VisitorIndex,
            WalkerIndex  = DistanceTableData::WalkerIndex
           };

      /** enum for a update mode */
      enum
      {
        ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
        ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
        ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
        ORB_WALKER,    /*!< walker update */
        ORB_ALLWALKER  /*!< all walkers update */
      };

      typedef ParticleAttrib<ValueType> ValueVectorType;
      typedef ParticleAttrib<GradType>  GradVectorType;
      typedef PooledData<RealType>      BufferType;
      typedef ParticleSet::Walker_t     Walker_t;
      typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
      typedef OrbitalSetTraits<ValueType>::GradMatrix_t  GradMatrix_t;

      /** flag to set the optimization mode */
      bool IsOptimizing;
      /** boolean to set optimization
       *
       * If true, this object is actively modified during optimization
       */
      bool Optimizable;

      bool derivsDone;
      int parameterType;
      /** current update mode */
      int UpdateMode;
      /** current \f$\log\phi \f$
       */
      RealType LogValue;
      /** current phase
       */
      RealType PhaseValue;
      /** Pointer to the differential orbital of this object
       *
       * If dPsi=0, this orbital is constant with respect to the optimizable variables
       */
      DiffOrbitalBasePtr dPsi;
      /** A vector for \f$ \frac{\partial \nabla \log\phi}{\partial \alpha} \f$
       */
      GradVectorType dLogPsi;
      /** A vector for \f$ \frac{\partial \nabla^2 \log\phi}{\partial \alpha} \f$
       */
      ValueVectorType d2LogPsi;
      /** Name of this orbital
       */
      string OrbitalName;
      ///list of variables this orbital handles
      opt_variables_type myVars;

      /// default constructor
      OrbitalBase();
      //OrbitalBase(const OrbitalBase& old);

      ///default destructor
      virtual ~OrbitalBase() { }

      inline void setOptimizable(bool optimizeit)
      {
        Optimizable = optimizeit;
      }
      
      virtual RealType getAlternatePhaseDiff(){return 0.0;}
      virtual RealType getAlternatePhaseDiff(int iat){return 0.0;}

      ///assign a differential orbital
      virtual void setDiffOrbital(DiffOrbitalBasePtr d);

      /** check in optimizable parameters
       * @param active a super set of optimizable variables
       *
       * Add the paramemters this orbital manage to active.
       */
      virtual void checkInVariables(opt_variables_type& active)=0;

      /** check out optimizable variables
       *
       * Update myVars index map
       */
      virtual void checkOutVariables(const opt_variables_type& active)=0;

      /** reset the parameters during optimizations
       */
      virtual void resetParameters(const opt_variables_type& active)=0;

      /** print the state, e.g., optimizables */
      virtual void reportStatus(ostream& os)=0;

      /** reset properties, e.g., distance tables, for a new target ParticleSet
       * @param P ParticleSet
       */
      virtual void resetTargetParticleSet(ParticleSet& P)=0;

      /** evaluate the value of the orbital for a configuration P.R
       *@param P  active ParticleSet
       *@param G  Gradients
       *@param L  Laplacians
       *@return the value
       *
       *Mainly for walker-by-walker move. The initial stage of particle-by-particle
       *move also uses this.
       */
      virtual ValueType
      evaluate(ParticleSet& P,
               ParticleSet::ParticleGradient_t& G,
               ParticleSet::ParticleLaplacian_t& L) = 0;

      /** evaluate the value of the orbital
       * @param P active ParticleSet
       * @param G Gradients, \f$\nabla\ln\Psi\f$
       * @param L Laplacians, \f$\nabla^2\ln\Psi\f$
       *
       */
      virtual RealType
      evaluateLog(ParticleSet& P,
                  ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) = 0;

      /** evaluate the value of the orbital
       * @param P active ParticleSet
       * @param G Gradients, \f$\nabla\ln\Psi\f$
       * @param L Laplacians, \f$\nabla^2\ln\Psi\f$
       * @param buf Buffer, data for analytical derivative calculation
       *
       */
      virtual RealType
      evaluateLog(ParticleSet& P,
                  ParticleSet::ParticleGradient_t& G, 
                  ParticleSet::ParticleLaplacian_t& L,
                  PooledData<RealType>& buf,
                  bool fillBuffer )
      {
         return evaluateLog(P,G,L); 
      } 

      /** return the current gradient for the iat-th particle
       * @param Pquantum particle set
       * @param iat particle index
       * @return the gradient of the iat-th particle
       */
      virtual GradType evalGrad(ParticleSet& P, int iat)
      {
        APP_ABORT("OrbitalBase::evalGradient is not implemented");
        return GradType();
      }

      virtual GradType alternateEvalGrad(ParticleSet& P, int iat)
      {
        return GradType();
      }

      /** return the logarithmic gradient for the iat-th particle
       * of the source particleset
       * @param Pquantum particle set
       * @param iat particle index
       * @return the gradient of the iat-th particle
       */
      virtual GradType evalGradSource(ParticleSet& P,
                                      ParticleSet& source,
                                      int iat)
      {
        // APP_ABORT("OrbitalBase::evalGradSource is not implemented");
        return GradType();
      }

      /** Adds the gradient w.r.t. the iat-th particle of the
       *  source particleset (ions) of the logarithmic gradient
       *  and laplacian w.r.t. the target paritlceset (electrons).
       * @param P quantum particle set (electrons)
       * @param source classical particle set (ions)
       * @param iat particle index of source (ion)
       * @param the ion gradient of the elctron gradient
       * @param the ion gradient of the elctron laplacian.
       * @return the log gradient of psi w.r.t. the source particle iat
       */
      virtual GradType evalGradSource
      (ParticleSet& P, ParticleSet& source, int iat,
       TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
       TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
      {
        return GradType();
      }



      /** evaluate the ratio of the new to old orbital value
       * @param P the active ParticleSet
       * @param iat the index of a particle
       * @param grad_iat Gradient for the active particle
       */
      virtual ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
      {
        APP_ABORT("OrbitalBase::ratioGrad is not implemented");
        return ValueType();
      }

      virtual ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
      {
        return 1.0;
      }

      virtual void alternateGrad(ParticleSet::ParticleGradient_t& G){}

      /** evaluate the ratio of the new to old orbital value
       *@param P the active ParticleSet
       *@param iat the index of a particle
       *@param dG the differential gradient
       *@param dL the differential laplacian
       *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'} \}) \f$
       *
       *Paired with acceptMove(ParticleSet& P, int iat).
       */
      virtual ValueType ratio(ParticleSet& P, int iat,
                              ParticleSet::ParticleGradient_t& dG,
                              ParticleSet::ParticleLaplacian_t& dL) = 0;

      /** a move for iat-th particle is accepted. Update the content for the next moves
       * @param P target ParticleSet
       * @param iat index of the particle whose new position was proposed
       */
      virtual void acceptMove(ParticleSet& P, int iat) =0;

      /** a move for iat-th particle is reject. Restore to the content.
       * @param iat index of the particle whose new position was proposed
       */
      virtual void restore(int iat) = 0;

      /** evalaute the ratio of the new to old orbital value
       *@param P the active ParticleSet
       *@param iat the index of a particle
       *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'}\})\f$
       *
       *Specialized for particle-by-particle move.
       */
      virtual ValueType ratio(ParticleSet& P, int iat) =0;

      virtual ValueType alternateRatio(ParticleSet& P)
      {
        return 1.0;
      };
      
      /** update the gradient and laplacian values by accepting a move
       *@param P the active ParticleSet
       *@param dG the differential gradients
       *@param dL the differential laplacians
       *@param iat the index of a particle
       *
       *Specialized for particle-by-particle move. Each Hamiltonian
       *updates its data for next update and evaluates differential gradients
       *and laplacians.
       */
      virtual void update(ParticleSet& P,
                          ParticleSet::ParticleGradient_t& dG,
                          ParticleSet::ParticleLaplacian_t& dL,
                          int iat) =0;


      /** equivalent to evaluateLog(P,G,L) with write-back function */
      virtual RealType evaluateLog(ParticleSet& P,BufferType& buf)=0;

      /** add temporary data reserved for particle-by-particle move.
       *
       * Return the log|psi|  like evalaute evaluateLog
       */
      virtual RealType registerData(ParticleSet& P, BufferType& buf) =0;

      /** add temporary (constant) data used to calculate analytical
       *  derivatives during linear optimization of parameters 
       */
      virtual void registerDataForDerivatives(ParticleSet& P, BufferType& buf)
      {
      } 

      /** re-evaluate the content and buffer data
       * @param P particle set
       * @param buf Anonymous storage
       *
       * This function is introduced to update the data periodically for particle-by-particle move.
       */
      virtual RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) =0;

      /** copy the internal data saved for particle-by-particle move.*/
      virtual void copyFromBuffer(ParticleSet& P, BufferType& buf)=0;

      /** dump the internal data to buf for optimizations
       *
       * Implments the default function that does nothing
       */
      virtual void dumpToBuffer(ParticleSet& P, BufferType& buf) {}

      /** copy the internal data from buf for optimizations
       *
       * Implments the default function that does nothing
       */
      virtual void dumpFromBuffer(ParticleSet& P, BufferType& buf) {}

      /** return a proxy orbital of itself
       */
      OrbitalBasePtr makeProxy(ParticleSet& tqp);
      /** make clone
       * @param tqp target Quantum ParticleSet
       * @param deepcopy if true, make a decopy
       *
       * If not true, return a proxy class
       */
      virtual OrbitalBasePtr makeClone(ParticleSet& tqp) const;

      /** Return the Chiesa kinetic energy correction
       */
      virtual RealType KECorrection();

      virtual void evaluateDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       vector<RealType>& dlogpsi,
                                       vector<RealType>& dhpsioverpsi) ;

//      virtual void evaluateDerivatives(ParticleSet& P,
//                                       const opt_variables_type& optvars,
//                                       vector<RealType>& dlogpsi,
//                                       vector<RealType>& dhpsioverpsi, 
//                                       PooledData<RealType>& buf) 
//      {  
//         evaluateDerivatives(P,optvars,dlogpsi,dhpsioverpsi);
//      } 

      virtual void finalizeOptimization() { }

      /** evaluate ratios to evaluate the momentum distribution
       * @param P reference particleset
       * @param ratios \f$ ratios[i]=\{{\bf R}\}\rightarrow {r_0,\cdots,r_i^p=pos,\cdots,r_{N-1}}\f$
       */
      virtual void get_ratios(ParticleSet& P, vector<ValueType>& ratios);

      ///** copy data members from old
      // * @param old existing OrbitalBase from which all the data members are copied.
      // *
      // * It is up to the derived classes to determine to use deep, shallow and mixed copy methods.
      // */
      //virtual void copyFrom(const OrbitalBase& old);

      /////////////////////////////////////////////////////
      // Functions for vectorized evaluation and updates //
      /////////////////////////////////////////////////////
#ifdef QMC_CUDA
      virtual void freeGPUmem() 
      { }
      
      virtual void recompute(MCWalkerConfiguration &W, bool firstTime)
      { }
      
      virtual void reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool)
      { }
      
      /** Evaluate the log of the WF for all walkers
       *  @param walkers   vector of all walkers
       *  @param logPsi    output vector of log(psi)
       */
      virtual void 
      addLog (MCWalkerConfiguration &W,
	      vector<RealType> &logPsi)
      {
	app_error() << "Need specialization of OrbitalBase::addLog for "
		    << OrbitalName << ".\n";
	abort();
      }
      
      /** Evaluate the wave-function ratio w.r.t. moving particle iat
       *  for all walkers
       *  @param walkers     vector of all walkers
       *  @param iat         particle which is moving
       *  @param psi_ratios  output vector with psi_new/psi_old
       */
      virtual void 
      ratio (MCWalkerConfiguration &W, int iat,
	     vector<ValueType> &psi_ratios)
      {
	app_error() << "Need specialization of OrbitalBase::ratio.\n";
	abort();
      }
      
      // Returns the WF ratio and gradient w.r.t. iat for each walker
      // in the respective vectors
      virtual void 
      ratio (MCWalkerConfiguration &W, int iat,
	     vector<ValueType> &psi_ratios,	vector<GradType>  &grad)
      {
	app_error() << "Need specialization of OrbitalBase::ratio.\n";
	abort();
      }
      
      virtual void 
      ratio (MCWalkerConfiguration &W, int iat,
	     vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
	     vector<ValueType> &lapl)
      {
	app_error() << "Need specialization of OrbitalBase::ratio.\n";
	abort();
      }
      
      
      virtual void 
      ratio (vector<Walker_t*> &walkers, vector<int> &iatList,
	     vector<PosType> &rNew,  vector<ValueType> &psi_ratios,	
	     vector<GradType>  &grad,  vector<ValueType> &lapl)
      {
	app_error() << "Need specialization of OrbitalBase::ratio.\n";
	abort();
      }
      
      
      virtual void 
      addGradient(MCWalkerConfiguration &W, int iat,
		  vector<GradType> &grad) 
      {
	app_error() << "Need specialization of OrbitalBase::addGradient for "
		    << OrbitalName << ".\n";
	abort();
      }
      
      virtual void 
      gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
		ValueMatrix_t &lapl)
      {
	app_error() << "Need specialization of OrbitalBase::gradLapl for "
		    << OrbitalName << ".\n";
	abort();
      }
      
      
      virtual void 
      update (vector<Walker_t*> &walkers, int iat)
      {
	app_error() << "Need specialization of OrbitalBase::update.\n";
	abort();
      }
      
      virtual void 
      update (const vector<Walker_t*> &walkers, 
	      const vector<int> &iatList)
      {
	app_error() << "Need specialization of OrbitalBase::update.\n";
	abort();
      }
      
      
      virtual void 
      NLratios (MCWalkerConfiguration &W,  vector<NLjob> &jobList,
		vector<PosType> &quadPoints, vector<ValueType> &psi_ratios)
      {
	app_error() << "Need specialization of OrbitalBase::NLRatios.\n";
	abort();
      }
      
      virtual void 
      NLratios (MCWalkerConfiguration &W,  gpu::device_vector<CUDA_PRECISION*> &Rlist,
		gpu::device_vector<int*> &ElecList, gpu::device_vector<int>             &NumCoreElecs,
		gpu::device_vector<CUDA_PRECISION*> &QuadPosList,
		gpu::device_vector<CUDA_PRECISION*> &RatioList,
		int numQuadPoints)
      {
	app_error() << "Need specialization of OrbitalBase::NLRatios.\n";
	abort();
      }
      
      virtual void
      evaluateDerivatives (MCWalkerConfiguration &W, 
			   const opt_variables_type& optvars,
			   ValueMatrix_t &dgrad_logpsi,
			   ValueMatrix_t &dhpsi_over_psi)
      {
	app_error() << "Need specialization of OrbitalBase::evaluateDerivatives.\n";
	abort();
      }
#endif
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/


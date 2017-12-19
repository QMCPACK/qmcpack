//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ORBITALBASE_H
#define QMCPLUSPLUS_ORBITALBASE_H
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"
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
struct NLjob
{
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
 \Psi(\{ {\bf R}\}) = \prod_i \psi_{i}(\{ {\bf R}\}),
 * \f]
 * where \f$\Psi\f$s are represented by
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
  typedef ParticleSet::Walker_t     Walker_t;
  typedef Walker_t::WFBuffer_t      WFBufferType;
  typedef Walker_t::Buffer_t        BufferType;
  typedef OrbitalSetTraits<RealType>::ValueMatrix_t       RealMatrix_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t      ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t       GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType           HessType;
  typedef OrbitalSetTraits<ValueType>::HessVector_t       HessVector_t;

  /** flag to set the optimization mode */
  bool IsOptimizing;
  /** boolean to set optimization
   *
   * If true, this object is actively modified during optimization
   */
  bool Optimizable;
  /** true, if FermionWF */
  bool IsFermionWF;
  /** true, if it is done with derivatives */
  bool derivsDone;
  /** true, if compute for the ratio instead of buffering */
  bool Need2Compute4PbyP;
  /** define the level of storage in derivative buffer **/
  int DerivStorageType;

  /** flag to calculate and return ionic derivatives */
  bool ionDerivs;

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
  std::string OrbitalName;
  ///list of variables this orbital handles
  opt_variables_type myVars;
  ///Bytes in WFBuffer
  size_t Bytes_in_WFBuffer;

  /// default constructor
  OrbitalBase();
  //OrbitalBase(const OrbitalBase& old);

  ///default destructor
  virtual ~OrbitalBase() { }

  inline void setOptimizable(bool optimizeit)
  {
    Optimizable = optimizeit;
  }

  inline void setIonDerivs(bool calcionderiv)
  {
    ionDerivs=calcionderiv;
  }

  virtual RealType getAlternatePhaseDiff()
  {
    return 0.0;
  }
  virtual RealType getAlternatePhaseDiff(int iat)
  {
    return 0.0;
  }

  virtual void resetPhaseDiff() {}

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
  virtual void reportStatus(std::ostream& os)=0;

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

  /** recompute the value of the orbitals which require critical accuracy.
   * needed for Slater Determinants but not needed for most types of orbitals
   */
  virtual void recompute(ParticleSet& P) {};

 // virtual void evaluateHessian(ParticleSet& P, IndexType iat, HessType& grad_grad_psi)
 // {
 //   APP_ABORT("OrbitalBase::evaluateHessian is not implemented");  
 // }
  
  virtual void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi_all)
  {
    APP_ABORT("OrbitalBase::evaluateHessian is not implemented");  
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

  virtual void alternateGrad(ParticleSet::ParticleGradient_t& G) {}

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

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void registerData(ParticleSet& P, WFBufferType& buf) =0;

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param P particle set
   * @param buf Anonymous storage
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   * @return log value of the wavefunction.
   */
  virtual RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false) =0;

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   *  The log value, P.G and P.L contribution from the objects
   *  of this class are also added.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf)=0;

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
                                   std::vector<RealType>& dlogpsi,
                                   std::vector<RealType>& dhpsioverpsi) ;
  virtual void multiplyDerivsByOrbR(std::vector<RealType>& dlogpsi)
  {
    RealType myrat = std::exp(LogValue)*std::cos(PhaseValue);
    for(int j=0; j<myVars.size(); j++)
    {
      int loc=myVars.where(j);
      dlogpsi[loc] *= myrat;
    }
  };

  /** Calculates the derivatives of \grad(\textrm{log}(\psi)) with respect to
      the optimizable parameters, and the dot product of this is then
      performed with the passed-in G_in gradient vector. This object is then
      returned as dgradlogpsi.
   */
  virtual void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                                       std::vector<RealType>& dgradlogpsi) {
    app_error() << "Need specialization of OrbitalBase::evaluateGradDerivatives.\n";
    abort();
  }

  virtual void finalizeOptimization() { }

  /** evaluate the ratios of one virtual move with respect to all the particles
   * @param P reference particleset
   * @param ratios \f$ ratios[i]=\{{\bf R}\}\rightarrow {r_0,\cdots,r_i^p=pos,\cdots,r_{N-1}}\f$
   */
  virtual void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.activePtcl
   */
  virtual void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.activePtcl
   * @param dratios \f$\partial_{\alpha}(\ln \Psi ({\bf R}^{\prime}) - \ln \Psi ({\bf R})) \f$
   */
  virtual void evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
      std::vector<ValueType>& ratios, Matrix<ValueType>& dratios);

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

  virtual void reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool)
  { }

  /** Evaluate the log of the WF for all walkers
   *  @param walkers   vector of all walkers
   *  @param logPsi    output vector of log(psi)
   */
  virtual void
  addLog (MCWalkerConfiguration &W,
          std::vector<RealType> &logPsi)
  {
    app_error() << "Need specialization of OrbitalBase::addLog for "
                << OrbitalName << ".\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
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
         std::vector<ValueType> &psi_ratios)
  {
    app_error() << "Need specialization of OrbitalBase::ratio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  // Returns the WF ratio and gradient w.r.t. iat for each walker
  // in the respective vectors
  virtual void
  ratio (MCWalkerConfiguration &W, int iat,
         std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad)
  {
    app_error() << "Need specialization of OrbitalBase::ratio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  ratio (MCWalkerConfiguration &W, int iat,
         std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
         std::vector<ValueType> &lapl)
  {
    app_error() << "Need specialization of OrbitalBase::ratio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  calcRatio (MCWalkerConfiguration &W, int iat,
             std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
             std::vector<ValueType> &lapl)
  {
    app_error() << "Need specialization of OrbitalBase::calcRatio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  addRatio (MCWalkerConfiguration &W, int iat,
            std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
            std::vector<ValueType> &lapl)
  {
    app_error() << "Need specialization of OrbitalBase::addRatio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
         std::vector<PosType> &rNew,  std::vector<ValueType> &psi_ratios,
         std::vector<GradType>  &grad,  std::vector<ValueType> &lapl)
  {
    app_error() << "Need specialization of OrbitalBase::ratio.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }


  virtual void
  addGradient(MCWalkerConfiguration &W, int iat,
              std::vector<GradType> &grad)
  {
    app_error() << "Need specialization of OrbitalBase::addGradient for "
                << OrbitalName << ".\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  calcGradient(MCWalkerConfiguration &W, int iat,
               std::vector<GradType> &grad)
  {
    app_error() << "Need specialization of OrbitalBase::calcGradient for "
                << OrbitalName << ".\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
            ValueMatrix_t &lapl)
  {
    app_error() << "Need specialization of OrbitalBase::gradLapl for "
                << OrbitalName << ".\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }


  virtual void
  update (std::vector<Walker_t*> &walkers, int iat)
  {
    app_error() << "Need specialization of OrbitalBase::update.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  update (const std::vector<Walker_t*> &walkers,
          const std::vector<int> &iatList)
  {
    app_error() << "Need specialization of OrbitalBase::update.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }


  virtual void
  NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
            std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
  {
    app_error() << "Need specialization of OrbitalBase::NLRatios.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
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
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }

  virtual void
  evaluateDerivatives (MCWalkerConfiguration &W,
                       const opt_variables_type& optvars,
                       RealMatrix_t &dgrad_logpsi,
                       RealMatrix_t &dhpsi_over_psi)
  {
    app_error() << "Need specialization of OrbitalBase::evaluateDerivatives.\n";
    app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
    abort();
  }
#endif
};
}
#endif


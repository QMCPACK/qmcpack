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


#ifndef QMCPLUSPLUS_WAVEFUNCTIONCOMPONENT_H
#define QMCPLUSPLUS_WAVEFUNCTIONCOMPONENT_H
#include "Message/Communicate.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/RecordProperty.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Particle/MCWalkerConfiguration.h"
#include "type_traits/template_types.hpp"
#ifdef QMC_CUDA
#include "type_traits/CUDATypes.h"
#endif

/**@file WaveFunctionComponent.h
 *@brief Declaration of WaveFunctionComponent
 */
namespace qmcplusplus
{
#ifdef QMC_CUDA
struct NLjob
{
  int walker;
  int elec;
  int numQuadPoints;
  NLjob(int w, int e, int n) : walker(w), elec(e), numQuadPoints(n) {}
};
#endif

///forward declaration of WaveFunctionComponent
class WaveFunctionComponent;
///forward declaration of DiffWaveFunctionComponent
class DiffWaveFunctionComponent;

typedef WaveFunctionComponent* WaveFunctionComponentPtr;
typedef DiffWaveFunctionComponent* DiffWaveFunctionComponentPtr;

/**@defgroup WaveFunctionComponent group
 * @brief Classes which constitute a many-body trial wave function
 *
 * A many-body trial wave function is
 * \f[
 \Psi(\{ {\bf R}\}) = \prod_i \psi_{i}(\{ {\bf R}\}),
 * \f]
 * where \f$\Psi\f$s are represented by
 * the derived classes from WaveFunctionComponent.
 */
/** @ingroup WaveFunctionComponent
 * @brief An abstract class for a component of a many-body trial wave function
 *
 * mw_ prefix is a function name signature indicating it is for handling a batch of WaveFunctionComponent objects
 * which are required to be base class pointers of the same derived class type.
 * all the mw_ routines must be implemented in a way either stateless or maintains states of every walker.
 */
struct WaveFunctionComponent : public QMCTraits
{
  /** enum for a update mode */
  enum
  {
    ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
    ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
    ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
    ORB_WALKER,       /*!< walker update */
    ORB_ALLWALKER     /*!< all walkers update */
  };

  typedef ParticleAttrib<ValueType> ValueVectorType;
  typedef ParticleAttrib<GradType> GradVectorType;
  typedef ParticleSet::Walker_t Walker_t;
  typedef Walker_t::WFBuffer_t WFBufferType;
  typedef Walker_t::Buffer_t BufferType;
  typedef OrbitalSetTraits<RealType>::ValueMatrix_t RealMatrix_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType HessType;
  typedef OrbitalSetTraits<ValueType>::HessVector_t HessVector_t;

  // the value type for log(psi)
  using LogValueType = std::complex<QTFull::RealType>;
  // the value type for psi(r')/psi(r)
  using PsiValueType = QTFull::ValueType;

  /** flag to set the optimization mode */
  bool IsOptimizing;
  /** boolean to set optimization
   *
   * If true, this object is actively modified during optimization
   */
  bool Optimizable;
  /** true, if this component is fermionic */
  bool is_fermionic;

  /** current update mode */
  int UpdateMode;
  /** current \f$\log\phi \f$
   */
  LogValueType LogValue;
  /** Pointer to the differential WaveFunctionComponent of this object
   *
   * If dPsi=0, this WaveFunctionComponent is constant with respect to the optimizable variables
   */
  DiffWaveFunctionComponentPtr dPsi;
  /** A vector for \f$ \frac{\partial \nabla \log\phi}{\partial \alpha} \f$
   */
  GradVectorType dLogPsi;
  /** A vector for \f$ \frac{\partial \nabla^2 \log\phi}{\partial \alpha} \f$
   */
  ValueVectorType d2LogPsi;
  /** Name of the class derived from WaveFunctionComponent
   */
  std::string ClassName;
  ///list of variables this WaveFunctionComponent handles
  opt_variables_type myVars;
  ///Bytes in WFBuffer
  size_t Bytes_in_WFBuffer;

  /// default constructor
  WaveFunctionComponent();
  //WaveFunctionComponent(const WaveFunctionComponent& old);

  ///default destructor
  virtual ~WaveFunctionComponent() {}

  inline void setOptimizable(bool optimizeit) { Optimizable = optimizeit; }

  ///assign a differential WaveFunctionComponent
  virtual void setDiffOrbital(DiffWaveFunctionComponentPtr d);

  ///assembles the full value
  PsiValueType getValue() const
  {
    return LogToValue<PsiValueType>::convert(LogValue);
  }

  /** check in optimizable parameters
   * @param active a super set of optimizable variables
   *
   * Add the paramemters this WaveFunctionComponent manage to active.
   */
  virtual void checkInVariables(opt_variables_type& active) = 0;

  /** check out optimizable variables
   *
   * Update myVars index map
   */
  virtual void checkOutVariables(const opt_variables_type& active) = 0;

  /** reset the parameters during optimizations
   */
  virtual void resetParameters(const opt_variables_type& active) = 0;

  /** print the state, e.g., optimizables */
  virtual void reportStatus(std::ostream& os) = 0;

  /** reset properties, e.g., distance tables, for a new target ParticleSet
   * @param P ParticleSet
   */
  virtual void resetTargetParticleSet(ParticleSet& P) = 0;

  /** evaluate the value of the WaveFunctionComponent from scratch
   * @param P  active ParticleSet
   * @param G Gradients, \f$\nabla\ln\Psi\f$
   * @param L Laplacians, \f$\nabla^2\ln\Psi\f$
   * @return the log value
   *
   * Mainly for walker-by-walker move. The initial stage of particle-by-particle
   * move also uses this.
   */
  virtual LogValueType evaluateLog(ParticleSet& P,
                               ParticleSet::ParticleGradient_t& G,
                               ParticleSet::ParticleLaplacian_t& L) = 0;

  /** evaluate from scratch the same type WaveFunctionComponent of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param G_list the list of Gradients pointers in a walker batch, \f$\nabla\ln\Psi\f$
   * @param L_list the list of Laplacians pointers in a walker batch, \f$\nabla^2\ln\Psi\f$
   * @@param values the log WF values of walkers in a batch
   */
  virtual void mw_evaluateLog(const std::vector<WaveFunctionComponent*>& WFC_list,
                              const std::vector<ParticleSet*>& P_list,
                              const std::vector<ParticleSet::ParticleGradient_t*>& G_list,
                              const std::vector<ParticleSet::ParticleLaplacian_t*>& L_list)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->evaluateLog(*P_list[iw], *G_list[iw], *L_list[iw]);
  }

  /** recompute the value of the WaveFunctionComponents which require critical accuracy.
   * needed for Slater Determinants but not needed for most types of WaveFunctionComponents
   */
  virtual void recompute(ParticleSet& P) {}

  // virtual void evaluateHessian(ParticleSet& P, IndexType iat, HessType& grad_grad_psi)
  // {
  //   APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented");
  // }

  virtual void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi_all)
  {
    APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented in " + ClassName + " class.");
  }

  /** return the current gradient for the iat-th particle
   * @param P quantum particle set
   * @param iat particle index
   * @return the gradient of the iat-th particle
   */
  virtual GradType evalGrad(ParticleSet& P, int iat)
  {
    APP_ABORT("WaveFunctionComponent::evalGradient is not implemented in " + ClassName + " class.");
    return GradType();
  }

  /** compute the current gradients for the iat-th particle of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param grad_now the list of gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_evalGrad(const std::vector<WaveFunctionComponent*>& WFC_list,
                           const std::vector<ParticleSet*>& P_list,
                           int iat,
                           std::vector<GradType>& grad_now)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      grad_now[iw] = WFC_list[iw]->evalGrad(*P_list[iw], iat);
  }

  /** compute the current gradients for the iat-th particle of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param grad_now the list of gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_evalGrad(const std::vector<std::reference_wrapper<WaveFunctionComponent>>& WFC_list,
                           const std::vector<std::reference_wrapper<ParticleSet>>& P_list,
                           int iat,
                           std::vector<GradType>& grad_now)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      grad_now[iw] = WFC_list[iw].get().evalGrad(P_list[iw].get(), iat);
  }

  /** return the logarithmic gradient for the iat-th particle
   * of the source particleset
   * @param Pquantum particle set
   * @param iat particle index
   * @return the gradient of the iat-th particle
   */
  virtual GradType evalGradSource(ParticleSet& P, ParticleSet& source, int iat)
  {
    // unit_test_hamiltonian calls this function incorrectly; do not abort for now
    //    APP_ABORT("WaveFunctionComponent::evalGradSource is not implemented");
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
  virtual GradType evalGradSource(ParticleSet& P,
                                  ParticleSet& source,
                                  int iat,
                                  TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                                  TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad)
  {
    return GradType();
  }


  /** evaluate the ratio of the new to old WaveFunctionComponent value and the new gradient
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @param grad_iat Gradient for the active particle
   */
  virtual ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    APP_ABORT("WaveFunctionComponent::ratioGrad is not implemented in " + ClassName + " class.");
    return ValueType();
  }

  /** compute the ratio of the new to old WaveFunctionComponent value and the new gradient of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   * @param grad_now the list of new gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_ratioGrad(const std::vector<WaveFunctionComponent*>& WFC_list,
                            const std::vector<ParticleSet*>& P_list,
                            int iat,
                            std::vector<PsiValueType>& ratios,
                            std::vector<GradType>& grad_new)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      ratios[iw] = WFC_list[iw]->ratioGrad(*P_list[iw], iat, grad_new[iw]);
  }

  /** compute the ratio of the new to old WaveFunctionComponent value and the new gradient of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   * @param grad_now the list of new gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_ratioGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                            const RefVector<ParticleSet>& P_list,
                            int iat,
                            std::vector<PsiValueType>& ratios,
                            std::vector<GradType>& grad_new)
  {
    //#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      ratios[iw] = WFC_list[iw].get().ratioGrad(P_list[iw], iat, grad_new[iw]);
  }

  /** a move for iat-th particle is accepted. Update the current content.
   * @param P target ParticleSet
   * @param iat index of the particle whose new position was proposed
   */
  virtual void acceptMove(ParticleSet& P, int iat) = 0;

  /** moves of the iat-th particle on some walkers in a batch is accepted. Update the current content.
   *  Note that all the lists only include accepted walkers.
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   */
  virtual void mw_acceptMove(const std::vector<WaveFunctionComponent*>& WFC_list,
                             const std::vector<ParticleSet*>& P_list,
                             int iat)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->acceptMove(*P_list[iw], iat);
  }

  /** complete all the delayed updates, must be called after each substep or step during pbyp move
   */
  virtual void completeUpdates() {}

  /** complete all the delayed updates for all the walkers in a batch
   * must be called after each substep or step during pbyp move
   */
  virtual void mw_completeUpdates(const std::vector<WaveFunctionComponent*>& WFC_list)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->completeUpdates();
  }

  /** If a move for iat-th particle is rejected, restore to the content.
   * @param iat index of the particle whose new position was proposed
   *
   * Ye: hopefully we can gradually move away from restore
   */
  virtual void restore(int iat) = 0;

  /** If a move for iat-th particle on some walkers in a batch is rejected, restore their contents
   *  Note that all the lists only include rejected walkers.
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param iat index of the particle whose new position was proposed
   *
   * Ye: hopefully we can gradually move away from restore
   */
  virtual void mw_restore(const std::vector<WaveFunctionComponent*>& WFC_list, int iat)
  {
    //#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->restore(iat);
  }

  /** evaluate the ratio of the new to old WaveFunctionComponent value
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}\})\f$
   *
   * Specialized for particle-by-particle move
   */
  virtual ValueType ratio(ParticleSet& P, int iat) = 0;

  /** compute the ratio of the new to old WaveFunctionComponent value of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   */
  virtual void mw_calcRatio(const std::vector<WaveFunctionComponent*>& WFC_list,
                            const std::vector<ParticleSet*>& P_list,
                            int iat,
                            std::vector<PsiValueType>& ratios)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      ratios[iw] = WFC_list[iw]->ratio(*P_list[iw], iat);
  }

  /** compute the ratio of the new to old WaveFunctionComponent value of multiple walkers
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   */
  virtual void mw_calcRatio(const RefVector<WaveFunctionComponent>& WFC_list,
                            const RefVector<ParticleSet>& P_list,
                            int iat,
                            std::vector<PsiValueType>& ratios)
  {
    //#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      ratios[iw] = WFC_list[iw].get().ratio(P_list[iw], iat);
  }

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void registerData(ParticleSet& P, WFBufferType& buf) = 0;

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param buf_list Anonymous storage
   */
  virtual void mw_registerData(const std::vector<WaveFunctionComponent*>& WFC_list,
                               const std::vector<ParticleSet*>& P_list,
                               const std::vector<WFBufferType*>& buf_list)
  {
    // We can't make this static but we can use a lambda with no capture to
    // restrict access to *this scope
    auto registerComponentData = [](WaveFunctionComponent& wfc, ParticleSet& pset, WFBufferType& wfb) {
      wfc.registerData(pset, wfb);
    };
    for (int iw = 0; iw < WFC_list.size(); iw++)
      registerComponentData(*(WFC_list[iw]), *(P_list[iw]), *(buf_list[iw]));
  }

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param P particle set
   * @param buf Anonymous storage
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   * @return log value of the wavefunction.
   */
  virtual LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) = 0;

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param WFC_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param P_list the list of ParticleSet pointers in a walker batch
   * @param buf_list Anonymous storage
   * @@param values the log WF values of walkers in a batch
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   */
  virtual void mw_updateBuffer(const RefVector<WaveFunctionComponent>& WFC_list,
                               const RefVector<ParticleSet>& P_list,
                               const RefVector<WFBufferType>& buf_list,
                               bool fromscratch = false)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw].get().updateBuffer(P_list[iw], buf_list[iw], fromscratch);
  }

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   *  The log value, P.G and P.L contribution from the objects
   *  of this class are also added.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) = 0;

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void mw_copyFromBuffer(const std::vector<WaveFunctionComponent*>& WFC_list,
                                 const std::vector<ParticleSet*>& P_list,
                                 const std::vector<WFBufferType*>& buf_list)
  {
#pragma omp parallel for
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->copyFromBuffer(*P_list[iw], *buf_list[iw]);
  }

  /** make clone
   * @param tqp target Quantum ParticleSet
   * @param deepcopy if true, make a decopy
   *
   * If not true, return a proxy class
   */
  virtual WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;

  /** Intended as a handle to break 
   *
   *  
   */
  //virtual WaveFunctionComponentPtr makeThrScope(std::vector<std::pair<int,int>>& ptcl_group_indexes) const = 0;

  /** Return the Chiesa kinetic energy correction
   */
  virtual RealType KECorrection();

  /** Compute derivatives of the wavefunction with respect to the optimizable
   *  parameters.
   *  @param P particle set
   *  @param optvars optimizable parameters
   *  @param dlogpsi array of derivatives of the log of the wavefunction
   *  @param dhpsioverpsi array of derivatives of the Laplacian of the wavefunction divided by the wavefunction.
   *          Note that this does not use the Laplacian of the log of the wavefunction, as in evaluateLog.
   *          Also the factor of -1/2 from the kinetic energy must be included here.  The 1/m
   *          factor is applied in TrialWaveFunction.
   */
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   std::vector<ValueType>& dlogpsi,
                                   std::vector<ValueType>& dhpsioverpsi);

  /** Compute derivatives of rhe wavefunction with respect to the optimizable 
   *  parameters
   *  @param P particle set
   *  @param optvars optimizable parameters
   *  @param dlogpsi array of derivatives of the log of the wavefunction
   *  Note: this function differs from the evaluateDerivatives function in the way that it only computes
   *        the derivative of the log of the wavefunction. 
  */
  virtual void evaluateDerivativesWF(ParticleSet& P,
                                     const opt_variables_type& optvars,
                                     std::vector<ValueType>& dlogpsi);

  virtual void multiplyDerivsByOrbR(std::vector<ValueType>& dlogpsi)
  {
    RealType myrat = std::real(LogToValue<PsiValueType>::convert(LogValue));
    for (int j = 0; j < myVars.size(); j++)
    {
      int loc = myVars.where(j);
      dlogpsi[loc] *= myrat;
    }
  }

  /** Calculates the derivatives of \grad(\textrm{log}(\psi)) with respect to
      the optimizable parameters, and the dot product of this is then
      performed with the passed-in G_in gradient vector. This object is then
      returned as dgradlogpsi.
   */

  virtual void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in, std::vector<ValueType>& dgradlogpsi)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::evaluateGradDerivatives in " + ClassName + " class.\n");
  }

  virtual void finalizeOptimization() {}

  /** evaluate the ratios of one virtual move with respect to all the particles
   * @param P reference particleset
   * @param ratios \f$ ratios[i]=\{{\bf R}\}\rightarrow {r_0,\cdots,r_i^p=pos,\cdots,r_{N-1}}\f$
   */
  virtual void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.refPtcl
   */
  virtual void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.refPtcl
   * @param dratios \f$\partial_{\alpha}(\ln \Psi ({\bf R}^{\prime}) - \ln \Psi ({\bf R})) \f$
   */
  virtual void evaluateDerivRatios(VirtualParticleSet& VP,
                                   const opt_variables_type& optvars,
                                   std::vector<ValueType>& ratios,
                                   Matrix<ValueType>& dratios);

  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
#ifdef QMC_CUDA
  using CTS = CUDAGlobalTypes;

  virtual void freeGPUmem() {}

  virtual void recompute(MCWalkerConfiguration& W, bool firstTime) {}

  virtual void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool, int kblocksize) {}

  /** Evaluate the log of the WF for all walkers
   *  @param walkers   vector of all walkers
   *  @param logPsi    output vector of log(psi)
   */
  virtual void addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::addLog for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  /** Evaluate the wave-function ratio w.r.t. moving particle iat
   *  for all walkers
   *  @param walkers     vector of all walkers
   *  @param iat         particle which is moving
   *  @param psi_ratios  output vector with psi_new/psi_old
   */
  virtual void ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::ratio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  // Returns the WF ratio and gradient w.r.t. iat for each walker
  // in the respective vectors
  virtual void ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios, std::vector<GradType>& grad)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::ratio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void ratio(MCWalkerConfiguration& W,
                     int iat,
                     std::vector<ValueType>& psi_ratios,
                     std::vector<GradType>& grad,
                     std::vector<ValueType>& lapl)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::ratio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void calcRatio(MCWalkerConfiguration& W,
                         int iat,
                         std::vector<ValueType>& psi_ratios,
                         std::vector<GradType>& grad,
                         std::vector<ValueType>& lapl)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::calcRatio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void addRatio(MCWalkerConfiguration& W,
                        int iat,
                        int k,
                        std::vector<ValueType>& psi_ratios,
                        std::vector<GradType>& grad,
                        std::vector<ValueType>& lapl)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::addRatio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void ratio(std::vector<Walker_t*>& walkers,
                     std::vector<int>& iatList,
                     std::vector<PosType>& rNew,
                     std::vector<ValueType>& psi_ratios,
                     std::vector<GradType>& grad,
                     std::vector<ValueType>& lapl)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::ratio for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }


  virtual void addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::addGradient for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::calcGradient for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void gradLapl(MCWalkerConfiguration& W, GradMatrix_t& grads, ValueMatrix_t& lapl)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::gradLapl for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void det_lookahead(MCWalkerConfiguration& W,
                             std::vector<ValueType>& psi_ratios,
                             std::vector<GradType>& grad,
                             std::vector<ValueType>& lapl,
                             int iat,
                             int k,
                             int kd,
                             int nw)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::det_lookahead for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void update(MCWalkerConfiguration* W, std::vector<Walker_t*>& walkers, int iat, std::vector<bool>* acc, int k)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::update for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::update for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }


  virtual void NLratios(MCWalkerConfiguration& W,
                        std::vector<NLjob>& jobList,
                        std::vector<PosType>& quadPoints,
                        std::vector<ValueType>& psi_ratios)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::NLRatios for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void NLratios(MCWalkerConfiguration& W,
                        gpu::device_vector<CUDA_PRECISION*>& Rlist,
                        gpu::device_vector<int*>& ElecList,
                        gpu::device_vector<int>& NumCoreElecs,
                        gpu::device_vector<CUDA_PRECISION*>& QuadPosList,
                        gpu::device_vector<CUDA_PRECISION*>& RatioList,
                        int numQuadPoints)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::NLRatios for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }

  virtual void evaluateDerivatives(MCWalkerConfiguration& W,
                                   const opt_variables_type& optvars,
                                   RealMatrix_t& dgrad_logpsi,
                                   RealMatrix_t& dhpsi_over_psi)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::evaluateDerivatives for " + ClassName +
              ".\n Required CUDA functionality not implemented. Contact developers.\n");
  }
#endif
};
} // namespace qmcplusplus
#endif

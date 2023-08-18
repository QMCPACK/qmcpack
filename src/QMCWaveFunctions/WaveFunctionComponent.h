//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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

#include <memory>
#include "Message/Communicate.h"
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/VirtualParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "OptimizableObject.h"
#include "Particle/MCWalkerConfiguration.h"
#include "type_traits/template_types.hpp"
#include "TWFGrads.hpp"

/**@file WaveFunctionComponent.h
 *@brief Declaration of WaveFunctionComponent
 */
namespace qmcplusplus
{
///forward declaration
class WaveFunctionComponent;
class ResourceCollection;
class TWFFastDerivWrapper;
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
class WaveFunctionComponent : public QMCTraits
{
public:
  /** enum for a update mode */
  enum
  {
    ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
    ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
    ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
    ORB_WALKER,       /*!< walker update */
    ORB_ALLWALKER     /*!< all walkers update */
  };

  using Walker_t     = ParticleSet::Walker_t;
  using WFBufferType = Walker_t::WFBuffer_t;
  using BufferType   = Walker_t::Buffer_t;
  using RealMatrix_t = OrbitalSetTraits<RealType>::ValueMatrix;
  using ValueMatrix  = OrbitalSetTraits<ValueType>::ValueMatrix;
  using GradMatrix   = OrbitalSetTraits<ValueType>::GradMatrix;
  using HessType     = OrbitalSetTraits<ValueType>::HessType;
  using HessVector   = OrbitalSetTraits<ValueType>::HessVector;

  // the value type for log(psi)
  using LogValueType = std::complex<QTFull::RealType>;
  // the value type for psi(r')/psi(r)
  using PsiValueType = QTFull::ValueType;

  /** current update mode */
  int UpdateMode;
  ///list of variables this WaveFunctionComponent handles
  opt_variables_type myVars;
  ///Bytes in WFBuffer
  size_t Bytes_in_WFBuffer;

protected:
  /** Name of the object
   * It is required to be different for objects of the same derived type like multiple J1.
   * It can be left empty for object which is unique per many-body WF.
   */
  const std::string my_name_;
  /** Current \f$\log\phi \f$.
   *  Exception! Slater Determinant most of the time has inconsistent a log_value inconsistent with the determinants
   *  it contains dduring a move sequence. That case the log_value_ would be more safely calculated on the fly.
   *
   *  There could be others.
   */
  LogValueType log_value_;

public:
  const LogValueType& get_log_value() const { return log_value_; }

  /// default constructor
  WaveFunctionComponent(const std::string& obj_name = "");
  ///default destructor
  virtual ~WaveFunctionComponent();

  /// Validate the internal consistency of the object
  virtual void checkSanity() const {}

  /// return object name
  const std::string& getName() const { return my_name_; }

  /// return class name
  virtual std::string getClassName() const = 0;

  ///assembles the full value
  PsiValueType getValue() const { return LogToValue<PsiValueType>::convert(log_value_); }

  /** true, if this component is fermionic */
  virtual bool isFermionic() const { return false; }

  /** check out variational optimizable variables
   * @param active a super set of optimizable variables
   */
  virtual void checkOutVariables(const opt_variables_type& active);

  /** Register the component with the TWFFastDerivWrapper wrapper.  
   */
  virtual void registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const;

  /** evaluate the value of the WaveFunctionComponent from scratch
   * \param[in] P  active ParticleSet
   * \param[out] G Gradients, \f$\nabla\ln\Psi\f$
   * \param[out] L Laplacians, \f$\nabla^2\ln\Psi\f$
   * \return the log value
   *
   * Mainly for walker-by-walker move. The initial stage of particle-by-particle
   * move also uses this. causes complete state update in WFC's
   */
  virtual LogValueType evaluateLog(const ParticleSet& P,
                                   ParticleSet::ParticleGradient& G,
                                   ParticleSet::ParticleLaplacian& L) = 0;

  /** evaluate from scratch the same type WaveFunctionComponent of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param G_list the list of Gradients pointers in a walker batch, \f$\nabla\ln\Psi\f$
   * @param L_list the list of Laplacians pointers in a walker batch, \f$\nabla^2\ln\Psi\f$
   * @param values the log WF values of walkers in a batch
   */
  virtual void mw_evaluateLog(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const RefVector<ParticleSet::ParticleGradient>& G_list,
                              const RefVector<ParticleSet::ParticleLaplacian>& L_list) const;

  /** recompute the value of the WaveFunctionComponents which require critical accuracy.
   * needed for Slater Determinants but not needed for most types of WaveFunctionComponents
   */
  virtual void recompute(const ParticleSet& P);

  virtual void mw_recompute(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            const std::vector<bool>& recompute) const;

  // virtual void evaluateHessian(ParticleSet& P, IndexType iat, HessType& grad_grad_psi)
  // {
  //   APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented");
  // }

  virtual void evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi_all)
  {
    APP_ABORT("WaveFunctionComponent::evaluateHessian is not implemented in " + getClassName() + " class.");
  }

  /** Prepare internal data for updating WFC correspond to a particle group
   * It should be called before moving particles of a given group.
   * This call can be used to handle the precomputation of data used for moving this group of particle.
   * Such data should be static with respect to the moves of particles within this group.
   * Particle groups usually correspond to determinants of different spins.
   * @param P quantum particle set
   * @param ig particle group index
   */
  virtual void prepareGroup(ParticleSet& P, int ig) {}

  virtual void mw_prepareGroup(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                               const RefVectorWithLeader<ParticleSet>& p_list,
                               int ig) const;

  /** return the current gradient for the iat-th particle
   * @param P quantum particle set
   * @param iat particle index
   * @return the gradient of the iat-th particle
   */
  virtual GradType evalGrad(ParticleSet& P, int iat)
  {
    APP_ABORT("WaveFunctionComponent::evalGradient is not implemented in " + getClassName() + " class.");
    return GradType();
  }


  /** return the current spin gradient for the iat-th particle
   * Default implementation assumes that WaveFunctionComponent does not explicitly depend on Spin.
   * @param P quantum particle set
   * @param iat particle index
   * @return the spin gradient of the iat-th particle
   */
  virtual GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad) { return evalGrad(P, iat); }

  /** compute the current gradients for the iat-th particle of multiple walkers
   * @param[out] grad_now the list of gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  template<CoordsType CT>
  void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                   const RefVectorWithLeader<ParticleSet>& p_list,
                   const int iat,
                   TWFGrads<CT>& grads_now) const;

  /** compute the current gradients for the iat-th particle of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param grad_now the list of gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_evalGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<GradType>& grad_now) const;

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
                                  TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                                  TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad)
  {
    return GradType();
  }


  /** evaluate the ratio of the new to old WaveFunctionComponent value and the new gradient
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @param grad_iat Gradient for the active particle
   */
  virtual PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  /** evaluate the ratio of the new to old WaveFunctionComponent value and the new spin gradient
   * Default implementation assumes that WaveFunctionComponent does not explicitly depend on Spin.
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @param grad_iat realspace gradient for the active particle
   * @param spingrad_iat spin gradient for the active particle
   */
  virtual PsiValueType ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat)
  {
    return ratioGrad(P, iat, grad_iat);
  }

  template<CoordsType CT>
  void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                    const RefVectorWithLeader<ParticleSet>& p_list,
                    int iat,
                    std::vector<PsiValueType>& ratios,
                    TWFGrads<CT>& grad_new) const;

  /** compute the ratio of the new to old WaveFunctionComponent value and the new gradient of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   * @param grad_now the list of new gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValueType>& ratios,
                            std::vector<GradType>& grad_new) const;

  /** a move for iat-th particle is accepted. Update the current content.
   * @param P target ParticleSet
   * @param iat index of the particle whose new position was proposed
   * @param safe_to_delay if true, delayed accept is safe.
   */
  virtual void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) = 0;

  /** moves of the iat-th particle on some walkers in a batch is accepted. Update the current content.
   *  Note that all the lists only include accepted walkers.
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param safe_to_delay if true, delayed accept is safe.
   */
  virtual void mw_accept_rejectMove(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    int iat,
                                    const std::vector<bool>& isAccepted,
                                    bool safe_to_delay = false) const;

  /** complete all the delayed or asynchronous operations before leaving the p-by-p move region.
   * Must be called at the end of each substep if p-by-p move is used.
   * This function was initially introduced for determinant delayed updates to complete all the delayed operations.
   * It has been extended to handle asynchronous operations on accellerators before leaving the p-by-p move region.
   */
  virtual void completeUpdates() {}

  /** complete all the delayed or asynchronous operations for all the walkers in a batch before leaving the p-by-p move region.
   */
  virtual void mw_completeUpdates(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const;

  /** If a move for iat-th particle is rejected, restore to the content.
   * @param iat index of the particle whose new position was proposed
   *
   * Ye: hopefully we can gradually move away from restore
   */
  virtual void restore(int iat) = 0;

  /** evaluate the ratio of the new to old WaveFunctionComponent value
   * @param P the active ParticleSet
   * @param iat the index of a particle
   * @return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}\})\f$
   *
   * Specialized for particle-by-particle move
   */
  virtual PsiValueType ratio(ParticleSet& P, int iat) = 0;

  /** compute the ratio of the new to old WaveFunctionComponent value of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   */
  virtual void mw_calcRatio(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            int iat,
                            std::vector<PsiValueType>& ratios) const;

  /** compute gradients and laplacian of the TWF with respect to each particle.
   * @param P particle set
   * @param G Gradients, \f$\nabla\ln\Psi\f$
   * @param L Laplacians, \f$\nabla^2\ln\Psi\f$
   * @param fromscratch if true and this WFC is sensitive to numeical error accumulation,
   *        all the internal data are recomputed from scratch.
   * @return log(psi)
   */
  virtual LogValueType evaluateGL(const ParticleSet& P,
                                  ParticleSet::ParticleGradient& G,
                                  ParticleSet::ParticleLaplacian& L,
                                  bool fromscratch);

  /** evaluate gradients and laplacian of the same type WaveFunctionComponent of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param G_list the list of Gradients pointers in a walker batch, \f$\nabla\ln\Psi\f$
   * @param L_list the list of Laplacians pointers in a walker batch, \f$\nabla^2\ln\Psi\f$
   * @param fromscratch if true and this WFC is sensitive to numerical error accumulation,
   *        all the internal data are recomputed from scratch.
   */
  virtual void mw_evaluateGL(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             const RefVector<ParticleSet::ParticleGradient>& G_list,
                             const RefVector<ParticleSet::ParticleLaplacian>& L_list,
                             bool fromscratch) const;

  /** For particle-by-particle move. Requests space in the buffer
   *  based on the data type sizes of the objects in this class.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void registerData(ParticleSet& P, WFBufferType& buf) = 0;

  /** For particle-by-particle move. Put the objects of this class
   *  in the walker buffer or forward the memory cursor.
   * @param P particle set
   * @param buf Anonymous storage
   * @param fromscratch request recomputing the precision critical
   *        pieces of wavefunction from scratch
   * @return log value of the wavefunction.
   */
  virtual LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) = 0;

  /** For particle-by-particle move. Copy data or attach memory
   *  from a walker buffer to the objects of this class.
   *  The log value, P.G and P.L contribution from the objects
   *  of this class are also added.
   * @param P particle set
   * @param buf Anonymous storage
   */
  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf) = 0;

  /** initialize a shared resource and hand it to a collection
   */
  virtual void createResource(ResourceCollection& collection) const {}

  /** acquire a shared resource from a collection
   */
  virtual void acquireResource(ResourceCollection& collection,
                               const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
  {}

  /** return a shared resource to a collection
   */
  virtual void releaseResource(ResourceCollection& collection,
                               const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
  {}

  /** make clone
   * @param tqp target Quantum ParticleSet
   * @param deepcopy if true, make a decopy
   *
   * If not true, return a proxy class
   */
  virtual std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const;

  /** Return the Chiesa kinetic energy correction
   */
  virtual RealType KECorrection();

  /** if true, this contains optimizable components
   */
  virtual bool isOptimizable() const { return false; }

  /** extract underlying OptimizableObject references
   * @param opt_obj_refs aggregated list of optimizable object references
   */
  virtual void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs);

  /** Compute the derivatives of both the log of the wavefunction and kinetic energy
   * with respect to optimizable parameters.
   *  @param P particle set
   *  @param optvars optimizable parameters
   *  @param dlogpsi array of derivatives of the log of the wavefunction.
   *         Add the contribution from this component.
   *  @param dhpsioverpsi array of Hamiltonian derivatives.
   *         Add the kinetic energy derivatives contribution from this component.
   *         \f$ -\frac{1}{2}{\partial}_\alpha \tilde L - G \cdot {\partial}_\alpha \tilde G \f$.
   *         \f$ \tilde L \f$ and \f$ \tilde G \f$ are from this WaveFunctionComponent.
   *         \f$ G \f$ is from TrialWaveFunction. The 1/m factor is applied in TrialWaveFunction.
   *         This is a bug when the particle set doesn't hold equal mass particles.
   */
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   Vector<ValueType>& dhpsioverpsi) = 0;

  /** Compute the derivatives of the log of the wavefunction with respect to optimizable parameters.
   *  parameters
   *  @param P particle set
   *  @param optvars optimizable parameters
   *  @param dlogpsi array of derivatives of the log of the wavefunction.
   *  Note: this function differs from the evaluateDerivatives function in the way that it only computes
   *        the derivative of the log of the wavefunction.
  */
  virtual void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, Vector<ValueType>& dlogpsi);

  /** Calculates the derivatives of \f$ \nabla \textnormal{log} \psi_f \f$ with respect to
      the optimizable parameters, and the dot product of this is then
      performed with the passed-in G_in gradient vector. This object is then
      returned as dgradlogpsi.
   */

  virtual void evaluateGradDerivatives(const ParticleSet::ParticleGradient& G_in, std::vector<ValueType>& dgradlogpsi)
  {
    APP_ABORT("Need specialization of WaveFunctionComponent::evaluateGradDerivatives in " + getClassName() +
              " class.\n");
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
  virtual void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  /** evaluate ratios to evaluate the non-local PP multiple walkers
   * @param wfc_list the list of WaveFunctionComponent references of the same component in a walker batch
   * @param vp_list the list of VirtualParticleSet references in a walker batch
   * @param ratios of all the virtual moves of all the walkers
   */
  virtual void mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                 const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                 std::vector<std::vector<ValueType>>& ratios) const;

  /** evaluate ratios to evaluate the non-local PP
   * @param VP VirtualParticleSet
   * @param ratios ratios with new positions VP.R[k] the VP.refPtcl
   * @param dratios Nq x Num_param matrix. \f$\partial_{\alpha}(\ln \Psi ({\bf R}^{\prime}) - \ln \Psi ({\bf R})) \f$
   */
  virtual void evaluateDerivRatios(const VirtualParticleSet& VP,
                                   const opt_variables_type& optvars,
                                   std::vector<ValueType>& ratios,
                                   Matrix<ValueType>& dratios);

private:
  /** compute the current gradients and spin gradients for the iat-th particle of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param grad_now the list of gradients in a walker batch, \f$\nabla\ln\Psi\f$
   * @param spingrad_now the list of spin gradients in a walker batch, \f$\nabla_s\ln\Psi\f$
   *
   */
  virtual void mw_evalGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   int iat,
                                   std::vector<GradType>& grad_now,
                                   std::vector<ComplexType>& spingrad_now) const;

  /** compute the ratio of the new to old WaveFunctionComponent value and the new gradient/spingradient of multiple walkers
   * @param wfc_list the list of WaveFunctionComponent pointers of the same component in a walker batch
   * @param p_list the list of ParticleSet pointers in a walker batch
   * @param iat particle index
   * @param ratios the list of WF ratios of a walker batch, \f$ \Psi( \{ {\bf R}^{'} \} )/ \Psi( \{ {\bf R}\})\f$
   * @param grad_now the list of new gradients in a walker batch, \f$\nabla\ln\Psi\f$
   */
  virtual void mw_ratioGradWithSpin(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    int iat,
                                    std::vector<PsiValueType>& ratios,
                                    std::vector<GradType>& grad_new,
                                    std::vector<ComplexType>& spingrad_new) const;
};
} // namespace qmcplusplus
#endif

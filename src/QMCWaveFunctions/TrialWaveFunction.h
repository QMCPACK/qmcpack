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
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file TrialWaveFunction.h
 *@brief Declaration of a TrialWaveFunction
 */
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_H

#include "Message/MPIObjectBase.h"
#include "Particle/VirtualParticleSet.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "Utilities/TimerManager.h"
#include "type_traits/template_types.hpp"
#include "Containers/MinimalContainers/RecordArray.hpp"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#ifdef QMC_CUDA
#include "type_traits/CUDATypes.h"
#endif

/**@defgroup MBWfs Many-body wave function group
 * @brief Classes to handle many-body trial wave functions
 */

namespace qmcplusplus
{
/** @ingroup MBWfs
 * @brief Class to represent a many-body trial wave function
 *
 *A many-body trial wave function is represented by
 *\f$\Psi({\bf R}) = \prod_i \psi_i({\bf R})\f$,
 *where each function \f$\psi_i({\bf R})\f$ is an WaveFunctionComponent
 (see WaveFunctionComponent).
 *A Composite Pattern is used to handle \f$\prod\f$ operations.
 *Each WaveFunctionComponent should provide proper evaluate functions
 *for the value, gradient and laplacian values.
 *
 * mw_ prefix is a function name signature indicating it is for handling
 * a batch of TrialWaveFunction objects in a lock-step fashion. These functions
 * are defined statically because they should not have access to a
 * concrete TWF object except through the passed RefVectorWithLeader<TWF>&.
 * It dispatches to mw_ functions of WaveFunctionComponent
 */
class TrialWaveFunction
{
public:
  // derived types from WaveFunctionComponent
  using RealType         = WaveFunctionComponent::RealType;
  using ComplexType      = WaveFunctionComponent::ComplexType;
  using FullPrecRealType = WaveFunctionComponent::FullPrecRealType;
  using ValueType        = WaveFunctionComponent::ValueType;
  using PosType          = WaveFunctionComponent::PosType;
  using GradType         = WaveFunctionComponent::GradType;
  using BufferType       = WaveFunctionComponent::BufferType;
  using WFBufferType     = WaveFunctionComponent::WFBufferType;
  using HessType         = WaveFunctionComponent::HessType;
  using HessVector       = WaveFunctionComponent::HessVector;
  using LogValueType     = WaveFunctionComponent::LogValueType;
  using PsiValueType     = WaveFunctionComponent::PsiValueType;

#ifdef QMC_CUDA
  using CTS          = CUDAGlobalTypes;
  using RealMatrix_t = WaveFunctionComponent::RealMatrix_t;
  using ValueMatrix  = WaveFunctionComponent::ValueMatrix;
  using GradMatrix   = WaveFunctionComponent::GradMatrix;
  using Walker_t     = ParticleSet::Walker_t;
#endif

  /// enum type for computing partial WaveFunctionComponents
  enum class ComputeType
  {
    ALL,
    FERMIONIC,
    NONFERMIONIC
  };

  ///differential gradients
  ParticleSet::ParticleGradient G;
  ///differential laplacians
  ParticleSet::ParticleLaplacian L;

  TrialWaveFunction(const std::string& aname = "psi0", bool tasking = false, bool create_local_resource = true);

  // delete copy constructor
  TrialWaveFunction(const TrialWaveFunction&) = delete;
  // deleteFassign operator
  TrialWaveFunction& operator=(const TrialWaveFunction&) = delete;

  ~TrialWaveFunction();

  inline int size() const { return Z.size(); }
  inline RealType getPhase() const { return PhaseValue; }

  inline void setPhase(RealType PhaseValue_new) { PhaseValue = PhaseValue_new; }
  void getLogs(std::vector<RealType>& lvals);
  void getPhases(std::vector<RealType>& pvals);

  inline RealType getPhaseDiff() const { return PhaseDiff; }
  inline void resetPhaseDiff() { PhaseDiff = 0.0; }
  inline RealType getLogPsi() const { return log_real_; }
  inline void setLogPsi(RealType LogPsi_new) { log_real_ = LogPsi_new; }

  /** add a WaveFunctionComponent
   * @param aterm a WaveFunctionComponent pointer
   */
  void addComponent(std::unique_ptr<WaveFunctionComponent>&& aterm);

  ///read from xmlNode
  bool put(xmlNodePtr cur);
  ///implement the virtual function
  void reset();
  /** set WaveFunctionComponent::IsOptimizing to true */
  void startOptimization();
  /** set WaveFunctionComponent::IsOptimizing to flase */
  void stopOptimization();

  // Wavefunction Parameter Optimization
  //
  // The wavefunction consists of a set of components (derived from WaveFunctionComponent).
  // Each of these components may or may not have variational parameters.
  // In order to perform optimization:
  //  1. Parameters are collected into a single list.
  //  2. Optimization algorithm computes new values for those parameters.
  //  3. Changed parameters are propagated back to each of the components.
  //
  // The collection of variables is of type VariableSet (opt_variables_type is a typedef).
  // The variables local to each component are stored in WaveFunctionComponent::myVars, which
  // is set to the local parameters when the component is set up.
  // The call to checkInVariables collects all the local parameters into a global list (step 1).
  // The resetIndex function on VariableSet then computes indices for this global list.
  // The call to checkOutVariables sets up the mapping from global index to local index in each component's 'myVars'.
  // Finally, the call to resetParameters progates the new values (step 3).
  // The call to checkOutVariables is a prerequisite for resetParameters to set the local values successfully.

  /** Check in an optimizable parameter
   * * @param o aggregated list of optimizable variables
   *
   * Gather all the optimizable parameters from wavefunction components into a single list
   */
  void checkInVariables(opt_variables_type& o);

  /** Check out optimizable variables
   * Assign index mappings from global list (o) to local values in each component
   */
  void checkOutVariables(const opt_variables_type& o);

  /**  Set values of parameters in each component from the global list
   */
  void resetParameters(const opt_variables_type& active);


  /** print out state of the trial wavefunction
   */
  void reportStatus(std::ostream& os);

  /** Initialize a TWF wrapper for fast derivative evaluation
   */
  void initializeTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const;
  /** evalaute the log (internally gradients and laplacian) of the trial wavefunction. gold reference */
  RealType evaluateLog(ParticleSet& P);

  /** batched version of evaluateLog. gold reference */
  static void mw_evaluateLog(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                             const RefVectorWithLeader<ParticleSet>& p_list);

  /** recompute the value of the orbitals which require critical accuracy */
  void recompute(const ParticleSet& P);

  /** batched version of recompute*/
  static void mw_recompute(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           const std::vector<bool>& recompute);

  /** evaluate the log value of a many-body wave function
   * @param P input configuration containing N particles
   * @param recomputeall recompute all orbitals from scratch
   * @return the value of \f$ \log( \Pi_i \Psi_i) \f$  many-body wave function
   *
   * @if recompute == true
   *   all orbitals have "evaluateLog" called on them, including the non-optimized ones.
   * @else
   *   default value.  call evaluateLog only on optimizable orbitals.  OK if nonlocal pp's aren't used.
   *
   * To save time, logpsi, G, and L are only computed for orbitals that change over the course of the optimization.
   * It is assumed that the fixed components are stored elsewhere.  See evaluateDeltaLog(P,logpsi_fixed_r,logpsi_opt,fixedG,fixedL)
   * defined below.  Nonlocal pseudopotential evaluation requires temporary information like matrix inverses, so while
   * the logpsi, G, and L don't change, evaluateLog is called anyways to compute these auxiliary quantities from scratch.
   * logpsi, G, and L associated with these non-optimizable orbitals are discarded explicitly and with dummy variables.
   */

  RealType evaluateDeltaLog(ParticleSet& P, bool recompute = false);

  /** evaluate the sum of log value of optimizable many-body wavefunctions
   * @param P  input configuration containing N particles
   * @param logpsi_fixed log(std::abs(psi)) of the invariant orbitals
   * @param logpsi_opt log(std::abs(psi)) of the variable orbitals
   * @param fixedG gradients of log(psi) of the fixed wave functions
   * @param fixedL laplacians of log(psi) of the fixed wave functions
   *
   * This function is introduced for optimization only.
   * fixedG and fixedL save the terms coming from the wave functions
   * that are invariant during optimizations.
   * It is expected that evaluateDeltaLog(P,false) is called later
   * and the external object adds the varying G and L and the fixed terms.
   */
  void evaluateDeltaLog(ParticleSet& P,
                        RealType& logpsi_fixed,
                        RealType& logpsi_opt,
                        ParticleSet::ParticleGradient& fixedG,
                        ParticleSet::ParticleLaplacian& fixedL);

  /** evaluate the sum of log value of optimizable many-body wavefunctions
   * @param wf_list vector of wavefunctions
   * @param p_list vector of input particle configurations
   * @param logpsi_fixed_list vector of log(std::abs(psi)) of the invariant orbitals
   * @param logpsi_opt_list vector of log(std::abs(psi)) of the variable orbitals
   * @param fixedG_list vector of gradients of log(psi) of the fixed wave functions
   * @param fixedL_list vector of laplacians of log(psi) of the fixed wave functions
   *
   * For wavefunction optimization, it can speed evaluation to split the log value,
   * the gradient, and the laplacian computed from wavefunction components with optimizable
   * parameters from components that do not.  This function computes the log value of
   * both parts, and the gradient and laplacian of the fixed components.
   * During correlated sampling steps only the components with optimizable
   * parameters need to have the gradient and laplacian re-evaluated.
   *
   * Parameters fixedG_list and fixedL_list save the terms coming from the components
   * that do not have optimizable parameters.
   * It is expected that mw_evaluateDeltaLog(P,false) is called later
   * and the external object adds the varying G and L and the fixed terms.
   */
  static void mw_evaluateDeltaLogSetup(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                       std::vector<RealType>& logpsi_fixed_list,
                                       std::vector<RealType>& logpsi_opt_list,
                                       RefVector<ParticleSet::ParticleGradient>& fixedG_list,
                                       RefVector<ParticleSet::ParticleLaplacian>& fixedL_list);

  /** evaluate the log value for optimizable parts of a many-body wave function
   * @param wf_list vector of wavefunctions
   * @param p_list vector of input particle configurations
   * @param logpsi_list vector of log(std::abs(psi)) of the variable orbitals
   * @param dummyG_list vector of gradients of log(psi) of the fixed wave functions.
   * @param dummyL_list vector of laplacians of log(psi) of the fixed wave functions
   *
   * The dummyG_list and dummyL_list are only referenced if recompute is true.
   * If recompute is false, the storage of these lists are needed, but the values can be discarded.
   *
   * @if recompute == true
   *   all orbitals have "evaluateLog" called on them, including the non-optimized ones.
   * @else
   *   default value.  call evaluateLog only on optimizable orbitals.  OK if nonlocal pp's aren't used.
   *
   * To save time, logpsi, G, and L are only computed for orbitals that change over the course of the optimization.
   * It is assumed that the fixed components are stored elsewhere.  See mw_evaluateDeltaLogSetup defined above.
   * Nonlocal pseudopotential evaluation requires temporary information like matrix inverses, so while
   * the logpsi, G, and L don't change, evaluateLog is called anyways to compute these auxiliary quantities from scratch.
   * logpsi, G, and L associated with these non-optimizable orbitals are discarded explicitly and with dummy variables.
   */


  static void mw_evaluateDeltaLog(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                  std::vector<RealType>& logpsi_list,
                                  RefVector<ParticleSet::ParticleGradient>& dummyG_list,
                                  RefVector<ParticleSet::ParticleLaplacian>& dummyL_list,
                                  bool recompute = false);


  /** compute psi(R_new) / psi(R_current) ratio
   * It returns a complex value if the wavefunction is complex.
   * @param P the active ParticleSet
   * @param iat the index of a particle moved to the new position.
   * @param ct select ComputeType
   * @return ratio value
   */
  ValueType calcRatio(ParticleSet& P, int iat, ComputeType ct = ComputeType::ALL);

  /** batched version of calcRatio */
  static void mw_calcRatio(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                           const RefVectorWithLeader<ParticleSet>& p_list,
                           int iat,
                           std::vector<PsiValueType>& ratios,
                           ComputeType ct = ComputeType::ALL);

  /** compulte multiple ratios to handle non-local moves and other virtual moves
   */
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios, ComputeType ct = ComputeType::ALL);
  /** batched version of evaluateRatios
   * Note: unlike other mw_ static functions, *this is the batch leader instead of wf_list[0].
   */
  static void mw_evaluateRatios(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                const RefVectorWithLeader<const VirtualParticleSet>& Vp_list,
                                const RefVector<std::vector<ValueType>>& ratios_list,
                                ComputeType ct = ComputeType::ALL);

  /** compute both ratios and deriatives of ratio with respect to the optimizables*/
  void evaluateDerivRatios(const VirtualParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& ratios,
                           Matrix<ValueType>& dratio);

  void printGL(ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L, std::string tag = "GL");

  /** Returns the logarithmic gradient of the trial wave function
   *  with respect to the iat^th atom of the source ParticleSet. */
  GradType evalGradSource(ParticleSet& P, ParticleSet& source, int iat);
  /** Returns the logarithmic gradient of the w.r.t. the iat^th atom
   * of the source ParticleSet of the sum of laplacians w.r.t. the
   * electrons (target ParticleSet) of the trial wave function. */
  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& source,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad);

  /** compute psi(R_new) / psi(R_current) ratio and \nabla ln(psi(R_new)) gradients
   * It returns a complex value if the wavefunction is complex.
   * @param P the active ParticleSet
   * @param iat the index of a particle moved to the new position.
   * @param grad_iat gradients
   * @return ratio value
   */
  ValueType calcRatioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  /** compute psi(R_new) / psi(R_current) ratio and d/ds ln(psi(R_new)) spin gradient
   * It returns a complex value if the wavefunction is complex.
   * @param P the active ParticleSet
   * @param iat the index of a particle moved to the new position.
   * @param grad_iat real space gradient for iat
   * @param spingrad_iat spin gradient for iat
   * @return ratio value
   */
  ValueType calcRatioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat);

  /** batched version of ratioGrad
   *
   *  all vector sizes must match
   */
  static void mw_calcRatioGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                               const RefVectorWithLeader<ParticleSet>& p_list,
                               int iat,
                               std::vector<PsiValueType>& ratios,
                               std::vector<GradType>& grad_new);

  /** batched version of ratioGradWithSpin
   *
   *  all vector sizes must match
   */
  static void mw_calcRatioGradWithSpin(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                       int iat,
                                       std::vector<PsiValueType>& ratios,
                                       std::vector<GradType>& grad_new,
                                       std::vector<ComplexType>& spingrad_new);

  /** Prepare internal data for updating WFC correspond to a particle group
   *  Particle groups usually correspond to determinants of different spins.
   *  This call can be used to handle precomputation for PbyP moves.
   * @param P quantum particle set
   * @param ig particle group index
   */
  void prepareGroup(ParticleSet& P, int ig);
  /** batched version of prepareGroup
   *
   *  all vector sizes must match
   */
  static void mw_prepareGroup(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              int ig);

  GradType evalGrad(ParticleSet& P, int iat);

  /** compute d/ds ln(psi) spin gradient at current particle position for iat electron
   *
   * @param P active particle set.
   * @param iat index of the particle moved to the new position.
   * @param spingrad spingrad value.  Zeroed out first, then filled with d/ds ln(psi).
   * @return \nabla ln(psi) (complex)
   *
   */
  GradType evalGradWithSpin(ParticleSet& P, int iat, ComplexType& spingrad);

  /** batched version of evalGrad
    *
    * This is static because it should have no direct access
    * to any TWF.
    */
  static void mw_evalGrad(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                          const RefVectorWithLeader<ParticleSet>& p_list,
                          int iat,
                          std::vector<GradType>& grad_now);

  /** batched version of evalGradWithSpin
    *
    * This is static because it should have no direct access
    * to any TWF.
    */
  static void mw_evalGradWithSpin(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                  int iat,
                                  std::vector<GradType>& grad_now,
                                  std::vector<ComplexType>& spingrad_now);

  void rejectMove(int iat);

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false);
  /* batched version of acceptMove */
  static void mw_accept_rejectMove(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   int iat,
                                   const std::vector<bool>& isAccepted,
                                   bool safe_to_delay = false);

  /** complete all the delayed or asynchronous operations before leaving the p-by-p move region.
   *  See WaveFunctionComponent::completeUpdates for more detail */
  void completeUpdates();
  /* batched version of completeUpdates.  */
  static void mw_completeUpdates(const RefVectorWithLeader<TrialWaveFunction>& wf_list);

  /** compute gradients and laplacian of the TWF with respect to each particle.
   *  See WaveFunctionComponent::evaluateGL for more detail */
  LogValueType evaluateGL(ParticleSet& P, bool fromscratch);
  /* batched version of evaluateGL.
   */
  static void mw_evaluateGL(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                            const RefVectorWithLeader<ParticleSet>& p_list,
                            bool fromscratch);

  /** register all the wavefunction components in buffer.
   *  See WaveFunctionComponent::registerData for more detail */
  void registerData(ParticleSet& P, WFBufferType& buf);

  /** update all the wavefunction components in buffer.
   *  See WaveFunctionComponent::updateBuffer for more detail */
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);

  /** copy all the wavefunction components from buffer.
   *  See WaveFunctionComponent::updateBuffer for more detail */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /// initialize a shared resource and hand it to a collection
  void createResource(ResourceCollection& collection) const;
  /** acquire external resource
   * Note: use RAII ResourceCollectionLock whenever possible
   */
  static void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<TrialWaveFunction>& wf_list);
  /** release external resource
   * Note: use RAII ResourceCollectionLock whenever possible
   */
  static void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<TrialWaveFunction>& wf_list);

  RealType KECorrection() const;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi,
                           bool project = false);

  static void mw_evaluateParameterDerivatives(const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                              const RefVectorWithLeader<ParticleSet>& p_list,
                                              const opt_variables_type& optvars,
                                              RecordArray<ValueType>& dlogpsi,
                                              RecordArray<ValueType>& dhpsioverpsi);

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, std::vector<ValueType>& dlogpsi);

  void evaluateGradDerivatives(const ParticleSet::ParticleGradient& G_in, std::vector<ValueType>& dgradlogpsi);

  /** evaluate the hessian w.r.t. electronic coordinates of particle iat **/
  // void evaluateHessian(ParticleSet & P, int iat, HessType& grad_grad_psi);
  /** evaluate the hessian hessian w.r.t. electronic coordinates of particle iat **/
  void evaluateHessian(ParticleSet& P, HessVector& all_grad_grad_psi);

  std::unique_ptr<TrialWaveFunction> makeClone(ParticleSet& tqp) const;

  std::vector<std::unique_ptr<WaveFunctionComponent>> const& getOrbitals() { return Z; }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  void setTwist(std::vector<RealType> t) { myTwist = t; }
  const std::vector<RealType> twist() { return myTwist; }

  inline void setMassTerm(ParticleSet& P)
  {
    OneOverM = 1.0 / P.Mass[0];
    //SpeciesSet tspecies(P.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //OneOverM = 1.0/mass;
  }

  RealType getReciprocalMass() { return OneOverM; }

  const std::string& getName() const { return myName; }

  bool use_tasking() const { return use_tasking_; }

private:
  static void debugOnlyCheckBuffer(WFBufferType& buffer);

  ///getName is in the way
  const std::string myName;

  ///starting index of the buffer
  size_t BufferCursor;

  ///starting index of the scalar buffer
  size_t BufferCursor_scalar;

  ///sign of the trial wave function
  RealType PhaseValue;

  ///diff of the phase of the trial wave function during ratio calls
  RealType PhaseDiff;

  ///real part of trial wave function log
  RealType log_real_;

  ///One over mass of target particleset, needed for Local Energy Derivatives
  RealType OneOverM;

  /// if true, using internal tasking implementation
  const bool use_tasking_;

  ///a list of WaveFunctionComponents constituting many-body wave functions
  std::vector<std::unique_ptr<WaveFunctionComponent>> Z;

  /// For now, TrialWaveFunction will own the wrapper.
  TWFFastDerivWrapper twf_prototype;
  /// timers at TrialWaveFunction function call level
  TimerList_t TWF_timers_;
  /// timers at WaveFunctionComponent function call level
  TimerList_t WFC_timers_;
  std::vector<RealType> myTwist;

  /** @{
   *  @brief helper function for extracting a list of WaveFunctionComponent from a list of TrialWaveFunction
   */
  static std::vector<WaveFunctionComponent*> extractWFCPtrList(const UPtrVector<TrialWaveFunction>& wf_list, int id);

  static RefVectorWithLeader<WaveFunctionComponent> extractWFCRefList(
      const RefVectorWithLeader<TrialWaveFunction>& wf_list,
      int id);
  /** }@ */

  // helper function for extrating a list of gradients from a list of TrialWaveFunction
  static RefVector<ParticleSet::ParticleGradient> extractGRefList(
      const RefVectorWithLeader<TrialWaveFunction>& wf_list);

  // helper function for extracting a list of laplacian from a list of TrialWaveFunction
  static RefVector<ParticleSet::ParticleLaplacian> extractLRefList(
      const RefVectorWithLeader<TrialWaveFunction>& wf_list);

  ///////////////////////////////////////////
  // Vectorized version for GPU evaluation //
  ///////////////////////////////////////////
#ifdef QMC_CUDA
private:
  gpu::device_host_vector<CTS::ValueType> GPUratios;
  gpu::device_host_vector<CTS::GradType> GPUgrads;
  gpu::device_host_vector<CTS::ValueType> GPUlapls;
  int ndelay; // delay rank

public:
  void freeGPUmem();

  void recompute(MCWalkerConfiguration& W, bool firstTime = true);

  void reserve(PointerPool<gpu::device_vector<CTS::ValueType>>& pool, bool onlyOptimizable = false, int kblocksize = 1);
  void getGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad);
  void calcGradient(MCWalkerConfiguration& W, int iat, int k, std::vector<GradType>& grad);
  void calcGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad) { calcGradient(W, iat, 0, grad); }
  void addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad);
  void evaluateLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi);
  void ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios);
  void ratio(MCWalkerConfiguration& W, int iat, std::vector<ValueType>& psi_ratios, std::vector<GradType>& newG);
  void ratio(MCWalkerConfiguration& W,
             int iat,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& newG,
             std::vector<ValueType>& newL);
  void calcRatio(MCWalkerConfiguration& W,
                 int iat,
                 std::vector<ValueType>& psi_ratios,
                 std::vector<GradType>& newG,
                 std::vector<ValueType>& newL);
  void addRatio(MCWalkerConfiguration& W,
                int iat,
                int k,
                std::vector<ValueType>& psi_ratios,
                std::vector<GradType>& newG,
                std::vector<ValueType>& newL);
  void addRatio(MCWalkerConfiguration& W,
                int iat,
                std::vector<ValueType>& psi_ratios,
                std::vector<GradType>& newG,
                std::vector<ValueType>& newL)
  {
    addRatio(W, iat, 0, psi_ratios, newG, newL);
  }
  void det_lookahead(MCWalkerConfiguration& W,
                     std::vector<ValueType>& psi_ratios,
                     std::vector<GradType>& grad,
                     std::vector<ValueType>& lapl,
                     int iat,
                     int k,
                     int kd,
                     int nw);

#ifdef QMC_COMPLEX
  void convertRatiosFromComplexToReal(std::vector<ValueType>& psi_ratios, std::vector<RealType>& psi_ratios_real);
#endif
  void ratio(std::vector<Walker_t*>& walkers,
             std::vector<int>& iatList,
             std::vector<PosType>& rNew,
             std::vector<ValueType>& psi_ratios,
             std::vector<GradType>& newG,
             std::vector<ValueType>& newL);

  void NLratios(MCWalkerConfiguration& W,
                gpu::device_vector<CUDA_PRECISION*>& Rlist,
                gpu::device_vector<int*>& ElecList,
                gpu::device_vector<int>& NumCoreElecs,
                gpu::device_vector<CUDA_PRECISION*>& QuadPosList,
                gpu::device_vector<CUDA_PRECISION*>& RatioList,
                int numQuadPoints,
                ComputeType ct = ComputeType::ALL);

  void NLratios(MCWalkerConfiguration& W,
                std::vector<NLjob>& jobList,
                std::vector<PosType>& quadPoints,
                std::vector<ValueType>& psi_ratios,
                ComputeType ct = ComputeType::ALL);

  void update(MCWalkerConfiguration* W, std::vector<Walker_t*>& walkers, int iat, std::vector<bool>* acc, int k);
  void update(std::vector<Walker_t*>& walkers, int iat) { update(NULL, walkers, iat, NULL, 0); }
  void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList);

  void gradLapl(MCWalkerConfiguration& W, GradMatrix& grads, ValueMatrix& lapl);


  void evaluateDeltaLog(MCWalkerConfiguration& W, std::vector<RealType>& logpsi_opt);

  void evaluateDeltaLog(MCWalkerConfiguration& W,
                        std::vector<RealType>& logpsi_fixed,
                        std::vector<RealType>& logpsi_opt,
                        GradMatrix& fixedG,
                        ValueMatrix& fixedL);

  void evaluateOptimizableLog(MCWalkerConfiguration& W,
                              std::vector<RealType>& logpsi_opt,
                              GradMatrix& optG,
                              ValueMatrix& optL);

  void evaluateDerivatives(MCWalkerConfiguration& W,
                           const opt_variables_type& optvars,
                           RealMatrix_t& dlogpsi,
                           RealMatrix_t& dhpsioverpsi);


  void setndelay(int delay) { ndelay = delay; }

  int getndelay() { return ndelay; }
#endif
};
/**@}*/
} // namespace qmcplusplus
#endif

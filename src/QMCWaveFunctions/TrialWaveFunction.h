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
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"
#include "Utilities/NewTimer.h"
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
 * flex_ prefix is a function name signature indicating it is for handling
 * a batch of TrialWaveFunction objects in a lock-step fashion.
 * It dispatched mw_ functions of WaveFunctionComponent or single walker functions
 * based on the number of objects in the list to faciliate different parallelization
 * schemes to maximize performance.
 */
class TrialWaveFunction : public MPIObjectBase
{
public:
  // derived types from WaveFunctionComponent
  typedef WaveFunctionComponent::RealType RealType;
  typedef WaveFunctionComponent::ValueType ValueType;
  typedef WaveFunctionComponent::PosType PosType;
  typedef WaveFunctionComponent::GradType GradType;
  typedef WaveFunctionComponent::BufferType BufferType;
  typedef WaveFunctionComponent::WFBufferType WFBufferType;
  typedef WaveFunctionComponent::HessType HessType;
  typedef WaveFunctionComponent::HessVector_t HessVector_t;
  using LogValueType = WaveFunctionComponent::LogValueType;
  using PsiValueType = WaveFunctionComponent::PsiValueType;

#ifdef QMC_CUDA
  using CTS = CUDAGlobalTypes;
  typedef WaveFunctionComponent::RealMatrix_t RealMatrix_t;
  typedef WaveFunctionComponent::ValueMatrix_t ValueMatrix_t;
  typedef WaveFunctionComponent::GradMatrix_t GradMatrix_t;
  typedef ParticleSet::Walker_t Walker_t;
#endif

  /// enum type for computing partial WaveFunctionComponents
  enum class ComputeType
  {
    ALL,
    FERMIONIC,
    NONFERMIONIC
  };

  ///differential gradients
  ParticleSet::ParticleGradient_t G;
  ///differential laplacians
  ParticleSet::ParticleLaplacian_t L;

  TrialWaveFunction(Communicate* c);

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
  inline void resetPhaseDiff()
  {
    PhaseDiff = 0.0;
    for (int i = 0; i < Z.size(); i++)
      Z[i]->resetPhaseDiff();
  }
  inline RealType getLogPsi() const { return LogValue; }
  inline void setLogPsi(RealType LogPsi_new) { LogValue = LogPsi_new; }

  /** add a WaveFunctionComponent
   * @param aterm a WaveFunctionComponent pointer
   * @param aname a name to the added WaveFunctionComponent object for printing
   */
  void addComponent(WaveFunctionComponent* aterm, std::string aname);

  ///read from xmlNode
  bool put(xmlNodePtr cur);
  ///implement the virtual function
  void reset();
  /** set WaveFunctionComponent::IsOptimizing to true */
  void startOptimization();
  /** set WaveFunctionComponent::IsOptimizing to flase */
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
  void reportStatus(std::ostream& os);

  /** recursively change the ParticleSet whose G and L are evaluated */
  void resetTargetParticleSet(ParticleSet& P);

  /** evalaute the log (internally gradients and laplacian) of the trial wavefunction. gold reference */
  RealType evaluateLog(ParticleSet& P);

  /** batched version of evaluateLog. gold reference */
  void flex_evaluateLog(const std::vector<TrialWaveFunction*>& WF_list, const std::vector<ParticleSet*>& P_list) const;

  /** recompute the value of the orbitals which require critical accuracy */
  void recompute(ParticleSet& P);

  RealType evaluateDeltaLog(ParticleSet& P, bool recompute = false);

  void evaluateDeltaLog(ParticleSet& P,
                        RealType& logpsi_fixed,
                        RealType& logpsi_opt,
                        ParticleSet::ParticleGradient_t& fixedG,
                        ParticleSet::ParticleLaplacian_t& fixedL);

  /** functions to handle particle-by-particle update
   * both ratio and full_ratio will be replaced by calcRatio which will handle ValueType.
   */
  RealType ratio(ParticleSet& P, int iat);

  /** function that computes psi(R_new) / psi(R_current). It returns a complex value if the wavefunction 
  *   is complex. It differs from the ratio(ParticleSet& P, int iat) function in the way that the ratio
  *   function takes the absolute value of psi(R_new) / psi(R_current). */
  ValueType calcRatio(ParticleSet& P, int iat, ComputeType ct = ComputeType::ALL);
  /** batched version of calcRatio */
  void flex_calcRatio(const std::vector<TrialWaveFunction*>& WF_list,
                      const std::vector<ParticleSet*>& P_list,
                      int iat,
                      std::vector<PsiValueType>& ratios,
                      ComputeType ct = ComputeType::ALL) const;

  /** compulte multiple ratios to handle non-local moves and other virtual moves
   */
  void evaluateRatios(VirtualParticleSet& P, std::vector<ValueType>& ratios);
  /** compute both ratios and deriatives of ratio with respect to the optimizables*/
  void evaluateDerivRatios(VirtualParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& ratios,
                           Matrix<ValueType>& dratio);

  void printGL(ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L, std::string tag = "GL");

  /** Returns the logarithmic gradient of the trial wave function
   *  with respect to the iat^th atom of the source ParticleSet. */
  GradType evalGradSource(ParticleSet& P, ParticleSet& source, int iat);
  /** Returns the logarithmic gradient of the w.r.t. the iat^th atom
   * of the source ParticleSet of the sum of laplacians w.r.t. the
   * electrons (target ParticleSet) of the trial wave function. */
  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& source,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad);

  RealType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  /** batched version of ratioGrad */
  void flex_ratioGrad(const std::vector<TrialWaveFunction*>& WF_list,
                      const std::vector<ParticleSet*>& P_list,
                      int iat,
                      std::vector<PsiValueType>& ratios,
                      std::vector<GradType>& grad_new) const;

  GradType evalGrad(ParticleSet& P, int iat);
  /** batched version of evalGrad */
  void flex_evalGrad(const std::vector<TrialWaveFunction*>& WF_list,
                     const std::vector<ParticleSet*>& P_list,
                     int iat,
                     std::vector<GradType>& grad_now) const;

  void rejectMove(int iat);
  /* flexible batched version of rejectMove */
  void flex_rejectMove(const std::vector<TrialWaveFunction*>& WF_list, int iat) const;

  void acceptMove(ParticleSet& P, int iat);
  /* flexible batched version of acceptMove */
  void flex_acceptMove(const std::vector<TrialWaveFunction*>& WF_list,
                       const std::vector<ParticleSet*>& P_list,
                       int iat) const;
  void completeUpdates();
  /* flexible batched version of completeUpdates.  */
  void flex_completeUpdates(const std::vector<TrialWaveFunction*>& WF_list) const;

  /** register all the wavefunction components in buffer.
   *  See WaveFunctionComponent::registerData for more detail */
  void registerData(ParticleSet& P, WFBufferType& buf);
  /* flexible batched version of registerData. 
   * Ye: perhaps it doesn't need to be flexible but just operates on all the walkers
   */
  void flex_registerData(const std::vector<TrialWaveFunction*>& WF_list,
                         const std::vector<ParticleSet*>& P_list,
                         const std::vector<WFBufferType*>& buf_list) const;

  /** update all the wavefunction components in buffer.
   *  See WaveFunctionComponent::updateBuffer for more detail */
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false);
  /* flexible batched version of updateBuffer. 
   * Ye: perhaps it doesn't need to be flexible but just operates on all the walkers
   */
  void flex_updateBuffer(const std::vector<TrialWaveFunction*>& WF_list,
                         const std::vector<ParticleSet*>& P_list,
                         const std::vector<WFBufferType*>& buf_list,
                         bool fromscratch = false) const;

  /** copy all the wavefunction components from buffer.
   *  See WaveFunctionComponent::updateBuffer for more detail */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);
  /* flexible batched version of copyFromBuffer. 
   * Ye: perhaps it doesn't need to be flexible but just operates on all the walkers
   */
  void flex_copyFromBuffer(const std::vector<TrialWaveFunction*>& WF_list,
                           const std::vector<ParticleSet*>& P_list,
                           const std::vector<WFBufferType*>& buf_list) const;

  RealType KECorrection() const;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi,
                           bool project = false);

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& optvars, std::vector<ValueType>& dlogpsi);

  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in, std::vector<ValueType>& dgradlogpsi);

  /** evaluate the hessian w.r.t. electronic coordinates of particle iat **/
  // void evaluateHessian(ParticleSet & P, int iat, HessType& grad_grad_psi);
  /** evaluate the hessian hessian w.r.t. electronic coordinates of particle iat **/
  void evaluateHessian(ParticleSet& P, HessVector_t& all_grad_grad_psi);

  TrialWaveFunction* makeClone(ParticleSet& tqp) const;

  std::vector<WaveFunctionComponent*>& getOrbitals() { return Z; }

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

  /* flexible batched version of evaluateGL.
   * TODO: split the computation from updateBuffer to evaluateGL. Expected to be called by KE
   */
  void flex_evaluateGL(const std::vector<TrialWaveFunction*>& WF_list, const std::vector<ParticleSet*>& P_list) const;

private:
  ///starting index of the buffer
  size_t BufferCursor;

  ///starting index of the scalar buffer
  size_t BufferCursor_scalar;

  ///sign of the trial wave function
  RealType PhaseValue;

  ///diff of the phase of the trial wave function during ratio calls
  RealType PhaseDiff;

  ///log of the trial wave function
  RealType LogValue;

  ///One over mass of target particleset, needed for Local Energy Derivatives
  RealType OneOverM;

  ///a list of WaveFunctionComponents constituting many-body wave functions
  std::vector<WaveFunctionComponent*> Z;

  std::vector<NewTimer*> myTimers;
  std::vector<RealType> myTwist;

  // helper function for extracting a list of WaveFunctionComponent from a list of TrialWaveFunction
  std::vector<WaveFunctionComponent*> extract_WFC_list(const std::vector<TrialWaveFunction*>& WF_list, int id) const;

  // helper function for extracting a list of gradients from a list of TrialWaveFunction
  std::vector<ParticleSet::ParticleGradient_t*> extract_G_list(const std::vector<TrialWaveFunction*>& P_list) const;

  // helper function for extracting a list of laplacian from a list of TrialWaveFunction
  std::vector<ParticleSet::ParticleLaplacian_t*> extract_L_list(const std::vector<TrialWaveFunction*>& P_list) const;

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
                int numQuadPoints);

  void NLratios(MCWalkerConfiguration& W,
                std::vector<NLjob>& jobList,
                std::vector<PosType>& quadPoints,
                std::vector<ValueType>& psi_ratios);

  void update(MCWalkerConfiguration* W, std::vector<Walker_t*>& walkers, int iat, std::vector<bool>* acc, int k);
  void update(std::vector<Walker_t*>& walkers, int iat) { update(NULL, walkers, iat, NULL, 0); }
  void update(const std::vector<Walker_t*>& walkers, const std::vector<int>& iatList);

  void gradLapl(MCWalkerConfiguration& W, GradMatrix_t& grads, ValueMatrix_t& lapl);


  void evaluateDeltaLog(MCWalkerConfiguration& W, std::vector<RealType>& logpsi_opt);

  void evaluateDeltaLog(MCWalkerConfiguration& W,
                        std::vector<RealType>& logpsi_fixed,
                        std::vector<RealType>& logpsi_opt,
                        GradMatrix_t& fixedG,
                        ValueMatrix_t& fixedL);

  void evaluateOptimizableLog(MCWalkerConfiguration& W,
                              std::vector<RealType>& logpsi_opt,
                              GradMatrix_t& optG,
                              ValueMatrix_t& optL);

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

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file TrialWaveFunctionBatched.h
 *
 * TrialWaveFunction<Batching::BATCHED> is a trial wave function for batched walker evaulations.
 */

#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_BATCHED_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_BATCHED_H

#include "TrialWaveFunction.h"

namespace qmcplusplus
{

template<>
class TrialWaveFunction<Batching::BATCHED>: public TrialWaveFunction<Batching::SINGLE>
{
public:
  typedef WaveFunctionComponent::CudaValueType   CudaValueType;
  typedef WaveFunctionComponent::CudaGradType    CudaGradType;
  typedef WaveFunctionComponent::CudaRealType    CudaRealType;
  typedef WaveFunctionComponent::RealMatrix_t    RealMatrix_t;
  typedef WaveFunctionComponent::ValueMatrix_t   ValueMatrix_t;
  typedef WaveFunctionComponent::GradMatrix_t    GradMatrix_t;
  typedef ParticleSet::Walker_t        Walker_t;

private:
  gpu::device_host_vector<CudaValueType>   GPUratios;
  gpu::device_host_vector<CudaGradType>    GPUgrads;
  gpu::device_host_vector<CudaValueType>   GPUlapls;

public:
  TrialWaveFunction(Communicate* c) : TrialWaveFunction<Batching::SINGLE>(c) {}
  
  void freeGPUmem GPU_XRAY_TRACE ();

  void recompute GPU_XRAY_TRACE (MCWalkerConfiguration &W, bool firstTime=true);

  void reserve GPU_XRAY_TRACE (PointerPool<gpu::device_vector<CudaValueType> > &pool,
                bool onlyOptimizable=false);
  void getGradient GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);
  void calcGradient GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
                     std::vector<GradType> &grad);
  void addGradient GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);
  void evaluateLog GPU_XRAY_TRACE (MCWalkerConfiguration &W,
                    std::vector<RealType> &logPsi);
  void ratio GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios);
  void ratio GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG);
  void ratio(MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG,
              std::vector<ValueType> &newL);
  void calcRatio GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
                  std::vector<ValueType> &psi_ratios,
                  std::vector<GradType> &newG,
                  std::vector<ValueType> &newL);
  void addRatio GPU_XRAY_TRACE (MCWalkerConfiguration &W, int iat,
                 std::vector<ValueType> &psi_ratios,
                 std::vector<GradType> &newG,
                 std::vector<ValueType> &newL);
#ifdef QMC_COMPLEX
  void convertRatiosFromComplexToReal GPU_XRAY_TRACE (std::vector<ValueType> &psi_ratios,
                                       std::vector<RealType> &psi_ratios_real);
#endif
  void ratio GPU_XRAY_TRACE (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
              std::vector<PosType> &rNew,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG,
              std::vector<ValueType> &newL);

  void NLratios GPU_XRAY_TRACE (MCWalkerConfiguration &W,
                 gpu::device_vector<CUDA_PRECISION*> &Rlist,
                 gpu::device_vector<int*>            &ElecList,
                 gpu::device_vector<int>             &NumCoreElecs,
                 gpu::device_vector<CUDA_PRECISION*> &QuadPosList,
                 gpu::device_vector<CUDA_PRECISION*> &RatioList,
                 int numQuadPoints);

  void NLratios GPU_XRAY_TRACE (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);

  void update GPU_XRAY_TRACE GPU_XRAY_TRACE (std::vector<Walker_t*> &walkers, int iat);
  void update GPU_XRAY_TRACE (const std::vector<Walker_t*> &walkers,
               const std::vector<int> &iatList);

  void gradLapl GPU_XRAY_TRACE (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);


  void evaluateDeltaLog(MCWalkerConfiguration &W,
                        std::vector<RealType>& logpsi_opt);

  void evaluateDeltaLog GPU_XRAY_TRACE (MCWalkerConfiguration &W,
                        std::vector<RealType>& logpsi_fixed,
                        std::vector<RealType>& logpsi_opt,
                        GradMatrix_t&  fixedG,
                        ValueMatrix_t& fixedL);

  void evaluateOptimizableLog GPU_XRAY_TRACE (MCWalkerConfiguration &W,
                               std::vector<RealType>& logpsi_opt,
                               GradMatrix_t&  optG,
                               ValueMatrix_t& optL);

  void evaluateDerivatives GPU_XRAY_TRACE (MCWalkerConfiguration &W,
                            const opt_variables_type& optvars,
                            RealMatrix_t &dlogpsi,
                            RealMatrix_t &dhpsioverpsi);

};

}
#endif

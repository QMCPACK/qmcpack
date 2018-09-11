//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, Oak Ridge National Laboratory
//                      refactored from SPOSet.cpp
//
// File created by: Peter Doak, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "SPOSetEvaluation.h"
#include "SPOSet.h"

/** default implementation */

template<>
SPOSet::ValueType
SPOSetEvaluation<Batching::SINGLE>::RATIO(const ParticleSet& P, int iat, const ValueType* restrict arow)
{
  int ip=omp_get_thread_num();
  // YYYY to fix
  /*
  ValueVector_t psi(t_logpsi[ip],OrbitalSetSize);
  evaluate(P,iat,psi);
  return simd::dot(psi.data(),arow,OrbitalSetSize,ValueType());
  */
  return ValueType();
}

template<>
void SPOSetEvaluation<Batching::SINGLE>::evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
{
  APP_ABORT("SPOSet::evaluateVGL not implemented.");
}

template<>
void SPOSetEvaluation<Batching::SINGLE>::evaluateValues(const VirtualParticleSet& VP, ValueMatrix_t& psiM, ValueAlignedVector_t& SPOmem)
{
  for(int iat=0; iat<VP.getTotalNum(); ++iat)
  {
    ValueVector_t psi(psiM[iat],OrbitalSetSize);
    evaluate(VP,iat,psi);
  }
}

template<>
void SPOSetEvaluation<Batching::SINGLE>::evaluateThirdDeriv(const ParticleSet& P, int first, int last,
                                    GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluateThirdDeriv(). \n");
}

template<>
void SPOSetEvaluation<Batching::SINGLE>::evaluate_notranspose(const ParticleSet& P, int first, int last
                                      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
}

template<>
void SPOSetEvaluation<Batching::SINGLE>::evaluate_notranspose(const ParticleSet& P, int first, int last,
                                      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


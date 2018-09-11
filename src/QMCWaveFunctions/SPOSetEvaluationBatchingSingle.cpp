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
#include "Configuration.h"

namespace qmcplusplus
{
/** single walker implementation */

using QMCT = QMCTraits;
  
QMCT::ValueType
SPOSetEvaluation<Batching::SINGLE>::RATIO(const ParticleSet& P, int iat, const QMCT::ValueType* restrict arow)
{
  int ip=omp_get_thread_num();
  // YYYY to fix
  /*
  ValueVector_t psi(t_logpsi[ip],OrbitalSetSize);
  evaluate(P,iat,psi);
  return simd::dot(psi.data(),arow,OrbitalSetSize,ValueType());
  */
  return QMCT::ValueType();
}

void SPOSetEvaluation<Batching::SINGLE>::evaluateVGL(const ParticleSet& P, int iat, SSTA::VGLVector_t& vgl)
{
  APP_ABORT("SPOSet::evaluateVGL not implemented.");
}

void SPOSetEvaluation<Batching::SINGLE>::evaluateValues(const VirtualParticleSet& VP,
							SSTA::ValueMatrix_t& psiM,
							SSTA::ValueAlignedVector_t& SPOmem,
							const QMCT::IndexType orbital_set_size)
{
  for(int iat=0; iat<VP.getTotalNum(); ++iat)
  {
    SSTA::ValueVector_t psi(psiM[iat],orbital_set_size);
    evaluate(VP,iat,psi);
  }
}

void SPOSetEvaluation<Batching::SINGLE>::evaluateThirdDeriv(const ParticleSet& P, int first, int last,
							    SSTA::GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluateThirdDeriv(). \n");
}

void SPOSetEvaluation<Batching::SINGLE>::evaluate_notranspose(const ParticleSet& P, int first, int last
							      , SSTA::ValueMatrix_t& logdet,
							      SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSetEvaluation<Batching::SINGLE>::evaluate_notranspose(const ParticleSet& P, int first, int last,
							      SSTA::ValueMatrix_t& logdet,
							      SSTA::GradMatrix_t& dlogdet,
							      SSTA::HessMatrix_t& grad_grad_logdet,
							      SSTA::GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}

void SPOSetEvaluation<Batching::SINGLE>::evaluateGradSource (const ParticleSet &P
                                     , int first, int last, const ParticleSet &source
							     , int iat_src, SSTA::GradMatrix_t &gradphi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

void SPOSetEvaluation<Batching::SINGLE>::evaluateGradSource (const ParticleSet &P, int first, int last,
                                     const ParticleSet &source, int iat_src,
							     SSTA::GradMatrix_t &grad_phi,
							     SSTA::HessMatrix_t &grad_grad_phi,
							     SSTA::GradMatrix_t &grad_lapl_phi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

}

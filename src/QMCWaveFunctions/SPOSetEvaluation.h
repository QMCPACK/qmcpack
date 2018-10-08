//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   refactored from SPOSet.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_EVALUATION_H
#define QMCPLUSPLUS_SPOSET_EVALUATION_H

#include "Configuration.h"
#include "Batching.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
//! SPOSet evaluation interface depends on walker batching strategy
/*!
  SPOSet inherits this so the correct evaluation function signatures and translations
  for the walker batching strategy are specified.  The actual evaluations
  are implement by SplineAdoptor descendents.
*/
namespace qmcplusplus
{

class SPOSetEvaluationDefault
{
public:
  using QMCT = QMCTraits;
  using SSTA = SPOSetTypeAliases;

  void appAbort()
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }
  
  virtual void
  evaluate (const ParticleSet& P, QMCT::PosType &r, SSTA::ValueVector_t &psi)
  {
    APP_ABORT("Need specialization for single walker SPOSet::evaluate "
              "(const ParticleSet& P, PosType &r).\n");
  }

  virtual void
  evaluate(const ParticleSet& P, int iat, SSTA::ValueVector_t& psi)
  {
    APP_ABORT("Need specialization for single walker SPOSet::evaluate"
	      "(const ParticleSet& P, int iat);");
  }

  virtual QMCT::ValueType RATIO(const ParticleSet& P, int iat, const QMCT::ValueType*
      restrict arow)
  {
    APP_ABORT("Need specialization for single walker SPOSet::RATIO"
	      "(const ParticleSet& P, int iat);")
    return QMCT::ValueType();
  }

  virtual void
  evaluateVGL(const ParticleSet& P, int iat, SSTA::VGLVector_t& vgl)
  {
    APP_ABORT("Need specialization for single walker SPOSet::evaluateVGL"
	      "(const ParticleSet& P, int iat");
  }

  virtual void
  evaluateValues(const VirtualParticleSet& VP, SSTA::ValueMatrix_t& psiM, SSTA::ValueAlignedVector_t& SPOMem, const QMCT::IndexType orbital_set_size)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::ValueVector_t& d2psi)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::HessVector_t& grad_grad_psi)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , SSTA::GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet,
				    SSTA::GradMatrix_t& dlogdet, SSTA::ValueMatrix_t& d2logdet)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet, SSTA::GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src, SSTA::GradMatrix_t &gradphi)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src
                                   , SSTA::GradMatrix_t &grad_phi, SSTA::HessMatrix_t &grad_grad_phi, SSTA::GradMatrix_t &grad_lapl_phi)
  {
    APP_ABORT("SPOSetEvaluationDefault methods should not be called");
  }

  // Interface for Batched Evaluation
  virtual void initGPU() { appAbort(); }

};
  
template<Batching batching = Batching::SINGLE>
class SPOSetEvaluation;

  
template<>
class SPOSetEvaluation<Batching::SINGLE> : public SPOSetEvaluationDefault
{
public:

}  
#endif

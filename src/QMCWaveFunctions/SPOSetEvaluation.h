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
#include "QMCWaveFunctions/BsplineFactory/temp_batch_type.h"
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

};
  
template<Batching batching>
class SPOSetEvaluation : public SPOSetEvaluationDefault
{
};

  
template<>
class SPOSetEvaluation<Batching::SINGLE> : public SPOSetEvaluationDefault
{
public:
  using SSTA = SPOSetTypeAliases;
  using QMCT = QMCTraits;
  virtual void
  evaluate (const ParticleSet& P, QMCT::PosType &r, SSTA::ValueVector_t &psi)
  {
    app_error() << "Need specialization for SPOSet::evaluate "
                << "(const ParticleSet& P, PosType &r).\n";
    abort();
  }

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat, SSTA::ValueVector_t& psi)=0;

  /** compute dot_product of new row and old row */
  virtual QMCT::ValueType RATIO(const ParticleSet& P, int iat, const QMCT::ValueType*
      restrict arow);

  /** evaluate VGL of SPOs using SoA container for gl
   */
  virtual void
  evaluateVGL(const ParticleSet& P, int iat, SSTA::VGLVector_t& vgl);

  /** evaluate values for the virtual moves, e.g., sphere move for nonlocalPP
   * @param VP virtual particle set
   * @param psiM single-particle orbitals psiM(i,j) for the i-th particle and the j-th orbital
   * @param SPOMem scratch space for SPO value evaluation, alignment is required.
   */
  virtual void
  evaluateValues(const VirtualParticleSet& VP, SSTA::ValueMatrix_t& psiM, SSTA::ValueAlignedVector_t& SPOMem, const QMCT::IndexType orbital_set_size);


  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::ValueVector_t& d2psi)=0;

  /** evaluate the values, gradients and hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           SSTA::ValueVector_t& psi, SSTA::GradVector_t& dpsi, SSTA::HessVector_t& grad_grad_psi)=0;

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , SSTA::GGGMatrix_t& grad_grad_grad_logdet);

  
  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last) particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   */
  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::ValueMatrix_t& d2logdet)=0;

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet);

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , SSTA::ValueMatrix_t& logdet, SSTA::GradMatrix_t& dlogdet, SSTA::HessMatrix_t& grad_grad_logdet, SSTA::GGGMatrix_t& grad_grad_grad_logdet);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src, SSTA::GradMatrix_t &gradphi);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src
                                   , SSTA::GradMatrix_t &grad_phi, SSTA::HessMatrix_t &grad_grad_phi, SSTA::GradMatrix_t &grad_lapl_phi);
};

  
template<>
class SPOSetEvaluation<Batching::BATCHED> : public SPOSetEvaluationDefault
{
public:
  using SSTA = SPOSetTypeAliases;
  using QMCT = QMCTraits;
  virtual void initGPU() {  }

  //////////////////////////////////////////
  // Walker-parallel vectorized functions //
  //////////////////////////////////////////
  virtual void
  reserve (PointerPool<gpu::device_vector<QMCT::CudaValueType> > &pool) { }

  virtual void
  evaluate (std::vector<SSTA::Walker_t*> &walkers, int iat, gpu::device_vector<QMCT::CudaValueType*> &phi);

  virtual void evaluate (std::vector<SSTA::Walker_t*> &walkers, std::vector<QMCT::PosType> &new_pos
                         , gpu::device_vector<QMCT::CudaValueType*> &phi);

  virtual void
  evaluate (std::vector<SSTA::Walker_t*> &walkers,
            std::vector<QMCT::PosType> &new_pos,
            gpu::device_vector<QMCT::CudaValueType*> &phi,
            gpu::device_vector<QMCT::CudaValueType*> &grad_lapl_list,
            int row_stride);

  virtual void
  evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaRealType*> &phi);
  virtual void
  evaluate (std::vector<QMCT::PosType> &pos, gpu::device_vector<QMCT::CudaComplexType*> &phi);

};




}
#endif

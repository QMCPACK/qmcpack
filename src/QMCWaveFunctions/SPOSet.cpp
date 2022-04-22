//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SPOSet.h"
#include "Message/Communicate.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "CPU/SIMD/simd.hpp"
#include "Utilities/ProgressReportEngine.h"
#include "hdf/hdf_archive.h"
#include <limits>

namespace qmcplusplus
{
SPOSet::SPOSet(bool use_OMP_offload, bool ion_deriv, bool optimizable)
    : useOMPoffload(use_OMP_offload), ionDerivs(ion_deriv), Optimizable(optimizable), OrbitalSetSize(0)
{
  className = "invalid";
}

void SPOSet::evaluateDetRatios(const VirtualParticleSet& VP,
                               ValueVector& psi,
                               const ValueVector& psiinv,
                               std::vector<ValueType>& ratios)
{
  assert(psi.size() == psiinv.size());
  for (int iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    evaluateValue(VP, iat, psi);
    ratios[iat] = simd::dot(psi.data(), psiinv.data(), psi.size());
  }
}

void SPOSet::mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                  const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                  const RefVector<ValueVector>& psi_list,
                                  const std::vector<const ValueType*>& invRow_ptr_list,
                                  std::vector<std::vector<ValueType>>& ratios_list) const
{
  assert(this == &spo_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
  {
    Vector<ValueType> invRow(const_cast<ValueType*>(invRow_ptr_list[iw]), psi_list[iw].get().size());
    spo_list[iw].evaluateDetRatios(vp_list[iw], psi_list[iw], invRow, ratios_list[iw]);
  }
}

void SPOSet::evaluateVGL_spin(const ParticleSet& P,
                              int iat,
                              ValueVector& psi,
                              GradVector& dpsi,
                              ValueVector& d2psi,
                              ValueVector& dspin)
{
  throw std::runtime_error("Need specialization of SPOSet::evaluateVGL_spin");
}

void SPOSet::mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& spo_list,
                            const RefVectorWithLeader<ParticleSet>& P_list,
                            int iat,
                            const RefVector<ValueVector>& psi_v_list,
                            const RefVector<GradVector>& dpsi_v_list,
                            const RefVector<ValueVector>& d2psi_v_list) const
{
  assert(this == &spo_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluateVGL(P_list[iw], iat, psi_v_list[iw], dpsi_v_list[iw], d2psi_v_list[iw]);
}

void SPOSet::mw_evaluateValue(const RefVectorWithLeader<SPOSet>& spo_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector>& psi_v_list) const
{
  assert(this == &spo_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluateValue(P_list[iw], iat, psi_v_list[iw]);
}

void SPOSet::mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                    const RefVectorWithLeader<ParticleSet>& P_list,
                                    int iat,
                                    const RefVector<ValueVector>& psi_v_list,
                                    const RefVector<GradVector>& dpsi_v_list,
                                    const RefVector<ValueVector>& d2psi_v_list,
                                    const RefVector<ValueVector>& dspin_v_list) const
{
  assert(this == &spo_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluateVGL_spin(P_list[iw], iat, psi_v_list[iw], dpsi_v_list[iw], d2psi_v_list[iw], dspin_v_list[iw]);
}

void SPOSet::mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                            const RefVectorWithLeader<ParticleSet>& P_list,
                                            int iat,
                                            const std::vector<const ValueType*>& invRow_ptr_list,
                                            VGLVector& phi_vgl_v,
                                            std::vector<ValueType>& ratios,
                                            std::vector<GradType>& grads) const
{
  assert(this == &spo_list.getLeader());
  const size_t nw             = spo_list.size();
  const size_t norb_requested = phi_vgl_v.size() / nw;
#pragma omp parallel for
  for (int iw = 0; iw < nw; iw++)
  {
    ValueVector phi_v(phi_vgl_v.data() + norb_requested * iw, norb_requested);
    GradVector dphi_v(reinterpret_cast<GradType*>(phi_vgl_v.data(1)) + norb_requested * iw, norb_requested);
    ValueVector d2phi_v(phi_vgl_v.data(4) + norb_requested * iw, norb_requested);
    spo_list[iw].evaluateVGL(P_list[iw], iat, phi_v, dphi_v, d2phi_v);

    ratios[iw] = simd::dot(invRow_ptr_list[iw], phi_v.data(), norb_requested);
    grads[iw]  = simd::dot(invRow_ptr_list[iw], dphi_v.data(), norb_requested) / ratios[iw];
  }
}

void SPOSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix& grad_grad_grad_logdet)
{
  throw std::runtime_error("Need specialization of SPOSet::evaluateThirdDeriv(). \n");
}

void SPOSet::evaluate_notranspose_spin(const ParticleSet& P,
                                       int first,
                                       int last,
                                       ValueMatrix& logdet,
                                       GradMatrix& dlogdet,
                                       ValueMatrix& d2logdet,
                                       ValueMatrix& dspinlogdet)
{
  throw std::runtime_error("Need specialization of " + className +
                           "::evaluate_notranspose_spin(P,iat,psi,dpsi,d2logdet, dspin_logdet) (vector quantities)\n");
}

void SPOSet::mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                                     const RefVectorWithLeader<ParticleSet>& P_list,
                                     int first,
                                     int last,
                                     const RefVector<ValueMatrix>& logdet_list,
                                     const RefVector<GradMatrix>& dlogdet_list,
                                     const RefVector<ValueMatrix>& d2logdet_list) const
{
  assert(this == &spo_list.getLeader());
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluate_notranspose(P_list[iw], first, last, logdet_list[iw], dlogdet_list[iw], d2logdet_list[iw]);
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix& logdet,
                                  GradMatrix& dlogdet,
                                  HessMatrix& grad_grad_logdet)
{
  throw std::runtime_error("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix& logdet,
                                  GradMatrix& dlogdet,
                                  HessMatrix& grad_grad_logdet,
                                  GGGMatrix& grad_grad_grad_logdet)
{
  throw std::runtime_error("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


std::unique_ptr<SPOSet> SPOSet::makeClone() const
{
  throw std::runtime_error("Missing  SPOSet::makeClone for " + className);
  return std::unique_ptr<SPOSet>();
}

void SPOSet::basic_report(const std::string& pad) const
{
  app_log() << pad << "size = " << size() << std::endl;
  app_log() << pad << "state info:" << std::endl;
  //states.report(pad+"  ");
  app_log().flush();
}

void SPOSet::evaluateVGH(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, HessVector& grad_grad_psi)
{
  throw std::runtime_error("Need specialization of " + className +
                           "::evaluate(P,iat,psi,dpsi,dhpsi) (vector quantities)\n");
}

void SPOSet::evaluateVGHGH(const ParticleSet& P,
                           int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& grad_grad_psi,
                           GGGVector& grad_grad_grad_psi)
{
  throw std::runtime_error("Need specialization of " + className +
                           "::evaluate(P,iat,psi,dpsi,dhpsi,dghpsi) (vector quantities)\n");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix& gradphi)
{
  throw std::runtime_error("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix& grad_phi,
                                HessMatrix& grad_grad_phi,
                                GradMatrix& grad_lapl_phi)
{
  throw std::runtime_error("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluateGradSourceRow(const ParticleSet& P,
                                   int iel,
                                   const ParticleSet& source,
                                   int iat_src,
                                   GradVector& gradphi)
{
  throw std::runtime_error("SPOSetBase::evalGradSourceRow is not implemented");
}

void SPOSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector& psi, ValueVector& dpsi)
{
  throw std::runtime_error("Need specialization of " + className +
                           "::evaluate_spin(P,iat,psi,dpsi) (vector quantities)\n");
}

#ifdef QMC_CUDA

void SPOSet::evaluate(const ParticleSet& P, PosType& r, ValueVector& psi)
{
  throw std::runtime_error("Need specialization for SPOSet::evaluate(const ParticleSet& P, PosType &r)\n");
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers, int iat, gpu::device_vector<CTS::ValueType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi,
                      gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                      int row_stride)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<Walker_t*>& walkers,
                      std::vector<PosType>& new_pos,
                      gpu::device_vector<CTS::ValueType*>& phi,
                      gpu::device_vector<CTS::ValueType*>& grad_lapl_list,
                      int row_stride,
                      int k,
                      bool klinear)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::RealType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

void SPOSet::evaluate(std::vector<PosType>& pos, gpu::device_vector<CTS::ComplexType*>& phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSet.\n";
  app_error() << "Required CUDA functionality not implemented. Contact developers.\n";
  abort();
}

#endif
} // namespace qmcplusplus

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


#include "QMCWaveFunctions/SPOSet.h"
#include "Message/Communicate.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"
#include <simd/simd.hpp>
#include "Utilities/ProgressReportEngine.h"
#include <io/hdf_archive.h>
#include <limits>

namespace qmcplusplus
{
SPOSet::SPOSet(bool use_OMP_offload, bool ion_deriv, bool optimizable)
    :
      useOMPoffload(use_OMP_offload),
      ionDerivs(ion_deriv),
      Optimizable(optimizable),
      OrbitalSetSize(0)
{
  className = "invalid";
}

void SPOSet::evaluate(const ParticleSet& P, PosType& r, ValueVector_t& psi)
{
  APP_ABORT("Need specialization for SPOSet::evaluate(const ParticleSet& P, PosType &r)\n");
}

void SPOSet::mw_evaluateValue(const RefVector<SPOSet>& spo_list,
                              const RefVector<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector_t>& psi_v_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateValue(P_list[iw], iat, psi_v_list[iw]);
}

void SPOSet::evaluateDetRatios(const VirtualParticleSet& VP,
                               ValueVector_t& psi,
                               const ValueVector_t& psiinv,
                               std::vector<ValueType>& ratios)
{
  assert(psi.size() == psiinv.size());
  for (int iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    evaluateValue(VP, iat, psi);
    ratios[iat] = simd::dot(psi.data(), psiinv.data(), psi.size());
  }
}

void SPOSet::mw_evaluateDetRatios(const RefVector<SPOSet>& spo_list,
                                  const RefVector<const VirtualParticleSet>& vp_list,
                                  const RefVector<ValueVector_t>& psi_list,
                                  const RefVector<const ValueVector_t>& psiinv_list,
                                  std::vector<std::vector<ValueType>>& ratios_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateDetRatios(vp_list[iw], psi_list[iw], psiinv_list[iw], ratios_list[iw]);
}

void SPOSet::mw_evaluateVGL(const RefVector<SPOSet>& spo_list,
                            const RefVector<ParticleSet>& P_list,
                            int iat,
                            const RefVector<ValueVector_t>& psi_v_list,
                            const RefVector<GradVector_t>& dpsi_v_list,
                            const RefVector<ValueVector_t>& d2psi_v_list)
{
#pragma omp parallel for
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].get().evaluateVGL(P_list[iw], iat, psi_v_list[iw], dpsi_v_list[iw], d2psi_v_list[iw]);
}

void SPOSet::mw_evaluateVGLandDetRatioGrads(const RefVector<SPOSet>& spo_list,
                                            const RefVector<ParticleSet>& P_list,
                                            int iat,
                                            const std::vector<const ValueType*>& invRow_ptr_list,
                                            VGLVector_t& phi_vgl_v,
                                            std::vector<ValueType>& ratios,
                                            std::vector<GradType>& grads)
{
  const size_t nw             = spo_list.size();
  const size_t norb_requested = phi_vgl_v.size() / nw;
#pragma omp parallel for
  for (int iw = 0; iw < nw; iw++)
  {
    ValueVector_t phi_v(phi_vgl_v.data() + norb_requested * iw, norb_requested);
    GradVector_t dphi_v(reinterpret_cast<GradType*>(phi_vgl_v.data(1)) + norb_requested * iw, norb_requested);
    ValueVector_t d2phi_v(phi_vgl_v.data(4) + norb_requested * iw, norb_requested);
    spo_list[iw].get().evaluateVGL(P_list[iw], iat, phi_v, dphi_v, d2phi_v);

    ratios[iw] = simd::dot(invRow_ptr_list[iw], phi_v.data(), norb_requested);
    grads[iw]  = simd::dot(invRow_ptr_list[iw], dphi_v.data(), norb_requested) / ratios[iw];
  }
}

void SPOSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluateThirdDeriv(). \n");
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix_t& logdet,
                                  GradMatrix_t& dlogdet,
                                  HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix_t& logdet,
                                  GradMatrix_t& dlogdet,
                                  HessMatrix_t& grad_grad_logdet,
                                  GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


SPOSet* SPOSet::makeClone() const
{
  APP_ABORT("Missing  SPOSet::makeClone for " + className);
  return 0;
}

void SPOSet::basic_report(const std::string& pad)
{
  app_log() << pad << "size = " << size() << std::endl;
  app_log() << pad << "state info:" << std::endl;
  //states.report(pad+"  ");
  app_log().flush();
}

void SPOSet::evaluateVGH(const ParticleSet& P,
                         int iat,
                         ValueVector_t& psi,
                         GradVector_t& dpsi,
                         HessVector_t& grad_grad_psi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate(P,iat,psi,dpsi,dhpsi) (vector quantities)\n");
}

void SPOSet::evaluateVGHGH(const ParticleSet& P,
                           int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           HessVector_t& grad_grad_psi,
                           GGGVector_t& grad_grad_grad_psi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate(P,iat,psi,dpsi,dhpsi,dghpsi) (vector quantities)\n");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix_t& gradphi)
{
  APP_ABORT("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix_t& grad_phi,
                                HessMatrix_t& grad_grad_phi,
                                GradMatrix_t& grad_lapl_phi)
{
  APP_ABORT("SPOSetBase::evalGradSource is not implemented");
}

void SPOSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector_t& psi, ValueVector_t& dpsi)
{
  APP_ABORT("Need specialization of " + className + "::evaluate_spin(P,iat,psi,dpsi) (vector quantities)\n");
}

#ifdef QMC_CUDA

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

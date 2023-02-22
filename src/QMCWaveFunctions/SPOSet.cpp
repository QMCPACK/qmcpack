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
SPOSet::SPOSet(const std::string& my_name) : my_name_(my_name), OrbitalSetSize(0) {}

void SPOSet::extractOptimizableObjectRefs(UniqueOptObjRefs&)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::extractOptimizableObjectRefs "
                           "must be overloaded when the SPOSet is optimizable.");
}

void SPOSet::checkOutVariables(const opt_variables_type& active)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::checkOutVariables "
                           "must be overloaded when the SPOSet is optimizable.");
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
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluateVGL(P_list[iw], iat, psi_v_list[iw], dpsi_v_list[iw], d2psi_v_list[iw]);
}

void SPOSet::mw_evaluateValue(const RefVectorWithLeader<SPOSet>& spo_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector>& psi_v_list) const
{
  assert(this == &spo_list.getLeader());
  for (int iw = 0; iw < spo_list.size(); iw++)
    spo_list[iw].evaluateValue(P_list[iw], iat, psi_v_list[iw]);
}

void SPOSet::mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                    const RefVectorWithLeader<ParticleSet>& P_list,
                                    int iat,
                                    const RefVector<ValueVector>& psi_v_list,
                                    const RefVector<GradVector>& dpsi_v_list,
                                    const RefVector<ValueVector>& d2psi_v_list,
                                    OffloadMatrix<ComplexType>& mw_dspin) const
{
  throw std::runtime_error(getClassName() + "::mw_evaluateVGLWithSpin() is not supported. \n");
}

void SPOSet::mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                            const RefVectorWithLeader<ParticleSet>& P_list,
                                            int iat,
                                            const std::vector<const ValueType*>& invRow_ptr_list,
                                            OffloadMWVGLArray& phi_vgl_v,
                                            std::vector<ValueType>& ratios,
                                            std::vector<GradType>& grads) const
{
  assert(this == &spo_list.getLeader());
  assert(phi_vgl_v.size(0) == DIM_VGL);
  assert(phi_vgl_v.size(1) == spo_list.size());
  const size_t nw             = spo_list.size();
  const size_t norb_requested = phi_vgl_v.size(2);
  GradVector dphi_v(norb_requested);
  for (int iw = 0; iw < nw; iw++)
  {
    ValueVector phi_v(phi_vgl_v.data_at(0, iw, 0), norb_requested);
    ValueVector d2phi_v(phi_vgl_v.data_at(4, iw, 0), norb_requested);
    spo_list[iw].evaluateVGL(P_list[iw], iat, phi_v, dphi_v, d2phi_v);

    ratios[iw] = simd::dot(invRow_ptr_list[iw], phi_v.data(), norb_requested);
    grads[iw]  = simd::dot(invRow_ptr_list[iw], dphi_v.data(), norb_requested) / ratios[iw];

    // transpose the array of gradients to SoA in phi_vgl_v
    for (size_t idim = 0; idim < DIM; idim++)
    {
      ValueType* phi_g = phi_vgl_v.data_at(idim + 1, iw, 0);
      for (size_t iorb = 0; iorb < norb_requested; iorb++)
        phi_g[iorb] = dphi_v[iorb][idim];
    }
  }
  phi_vgl_v.updateTo();
}

void SPOSet::mw_evaluateVGLandDetRatioGradsWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                                    const RefVectorWithLeader<ParticleSet>& P_list,
                                                    int iat,
                                                    const std::vector<const ValueType*>& invRow_ptr_list,
                                                    OffloadMWVGLArray& phi_vgl_v,
                                                    std::vector<ValueType>& ratios,
                                                    std::vector<GradType>& grads,
                                                    std::vector<ValueType>& spingrads) const
{
  throw std::runtime_error("Need specialization of " + getClassName() +
                           "::mw_evaluateVGLandDetRatioGradsWithSpin(). \n");
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
  throw std::runtime_error("Need specialization of " + getClassName() +
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
  throw std::runtime_error("Missing  SPOSet::makeClone for " + getClassName());
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
  throw std::runtime_error("Need specialization of " + getClassName() +
                           "::evaluate(P,iat,psi,dpsi,dhpsi) (vector quantities)\n");
}

void SPOSet::evaluateVGHGH(const ParticleSet& P,
                           int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& grad_grad_psi,
                           GGGVector& grad_grad_grad_psi)
{
  throw std::runtime_error("Need specialization of " + getClassName() +
                           "::evaluate(P,iat,psi,dpsi,dhpsi,dghpsi) (vector quantities)\n");
}

void SPOSet::applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy)
{
  if (isRotationSupported())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::applyRotation "
                           "must be overloaded when the SPOSet supports rotation.");
}

void SPOSet::evaluateDerivatives(ParticleSet& P,
                                 const opt_variables_type& optvars,
                                 Vector<ValueType>& dlogpsi,
                                 Vector<ValueType>& dhpsioverpsi,
                                 const int& FirstIndex,
                                 const int& LastIndex)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateDerivatives "
                           "must be overloaded when the SPOSet is optimizable.");
}

void SPOSet::evaluateDerivativesWF(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   int FirstIndex,
                                   int LastIndex)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateDerivativesWF "
                           "must be overloaded when the SPOSet is optimizable.");
}

void SPOSet::evaluateDerivRatios(const VirtualParticleSet& VP,
                                 const opt_variables_type& optvars,
                                 ValueVector& psi,
                                 const ValueVector& psiinv,
                                 std::vector<ValueType>& ratios,
                                 Matrix<ValueType>& dratios,
                                 int FirstIndex,
                                 int LastIndex)
{
  // Match the fallback in WaveFunctionComponent that evaluates just the ratios
  evaluateDetRatios(VP, psi, psiinv, ratios);

  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateDerivRatios "
                           "must be overloaded when the SPOSet is optimizable.");
}


/** Evaluate the derivative of the optimized orbitals with respect to the parameters
   *  this is used only for MSD, to be refined for better serving both single and multi SD
   */
void SPOSet::evaluateDerivatives(ParticleSet& P,
                                 const opt_variables_type& optvars,
                                 Vector<ValueType>& dlogpsi,
                                 Vector<ValueType>& dhpsioverpsi,
                                 const ValueType& psiCurrent,
                                 const std::vector<ValueType>& Coeff,
                                 const std::vector<size_t>& C2node_up,
                                 const std::vector<size_t>& C2node_dn,
                                 const ValueVector& detValues_up,
                                 const ValueVector& detValues_dn,
                                 const GradMatrix& grads_up,
                                 const GradMatrix& grads_dn,
                                 const ValueMatrix& lapls_up,
                                 const ValueMatrix& lapls_dn,
                                 const ValueMatrix& M_up,
                                 const ValueMatrix& M_dn,
                                 const ValueMatrix& Minv_up,
                                 const ValueMatrix& Minv_dn,
                                 const GradMatrix& B_grad,
                                 const ValueMatrix& B_lapl,
                                 const std::vector<int>& detData_up,
                                 const size_t N1,
                                 const size_t N2,
                                 const size_t NP1,
                                 const size_t NP2,
                                 const std::vector<std::vector<int>>& lookup_tbl)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateDerivatives "
                           "must be overloaded when the SPOSet is optimizable.");
}

/** Evaluate the derivative of the optimized orbitals with respect to the parameters
   *  this is used only for MSD, to be refined for better serving both single and multi SD
   */
void SPOSet::evaluateDerivativesWF(ParticleSet& P,
                                   const opt_variables_type& optvars,
                                   Vector<ValueType>& dlogpsi,
                                   const QTFull::ValueType& psiCurrent,
                                   const std::vector<ValueType>& Coeff,
                                   const std::vector<size_t>& C2node_up,
                                   const std::vector<size_t>& C2node_dn,
                                   const ValueVector& detValues_up,
                                   const ValueVector& detValues_dn,
                                   const ValueMatrix& M_up,
                                   const ValueMatrix& M_dn,
                                   const ValueMatrix& Minv_up,
                                   const ValueMatrix& Minv_dn,
                                   const std::vector<int>& detData_up,
                                   const std::vector<std::vector<int>>& lookup_tbl)
{
  if (isOptimizable())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateDerivativesWF "
                           "must be overloaded when the SPOSet is optimizable.");
}


void SPOSet::evaluateGradSource(const ParticleSet& P,
                                int first,
                                int last,
                                const ParticleSet& source,
                                int iat_src,
                                GradMatrix& gradphi)
{
  if (hasIonDerivs())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateGradSource "
                           "must be overloaded when the SPOSet has ion derivatives.");
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
  if (hasIonDerivs())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateGradSource "
                           "must be overloaded when the SPOSet has ion derivatives.");
}

void SPOSet::evaluateGradSourceRow(const ParticleSet& P,
                                   int iel,
                                   const ParticleSet& source,
                                   int iat_src,
                                   GradVector& gradphi)
{
  if (hasIonDerivs())
    throw std::logic_error("Bug!! " + getClassName() +
                           "::evaluateGradSourceRow "
                           "must be overloaded when the SPOSet has ion derivatives.");
}

void SPOSet::evaluate_spin(const ParticleSet& P, int iat, ValueVector& psi, ValueVector& dpsi)
{
  throw std::runtime_error("Need specialization of " + getClassName() +
                           "::evaluate_spin(P,iat,psi,dpsi) (vector quantities)\n");
}

} // namespace qmcplusplus

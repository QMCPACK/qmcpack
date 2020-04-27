//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"
#include "simd/simd.hpp"

namespace qmcplusplus
{
/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
template<typename DET_ENGINE_TYPE>
DiracDeterminantBatched<DET_ENGINE_TYPE>::DiracDeterminantBatched(SPOSetPtr const spos, int first)
    : DiracDeterminantBase(spos, first)
{
  ClassName = "DiracDeterminantBatched";
}

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::set(int first, int nel, int delay)
{
  FirstIndex = first;

  resize(nel, nel);

  if (Optimizable)
    Phi->buildOptVariables(nel);
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::invertPsiM(const ValueMatrix_t& logdetT,
                                                          OffloadPinnedValueMatrix_t& invMat)
{
  InverseTimer.start();
  det_engine_.invert_transpose(logdetT, invMat, LogValue);
  InverseTimer.stop();
}


///reset the size: with the number of particles and number of orbtials
template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::resize(int nel, int morb)
{
  int norb = morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiMinv.resize(nel, norb);
  psiMinv_dev_ptr = getOffloadDevicePtr(psiMinv.data());
  psiM_vgl.resize(nel * norb);
  psiM_vgl_dev_ptr = getOffloadDevicePtr(psiM_vgl.data());
  psiM_temp.attachReference(psiM_vgl.data(0), nel, norb);
  dpsiM.attachReference(reinterpret_cast<GradType*>(psiM_vgl.data(1)), nel, norb);
  d2psiM.attachReference(psiM_vgl.data(4), nel, norb);

  LastIndex   = FirstIndex + nel;
  NumPtcls    = nel;
  NumOrbitals = norb;

  psiV.resize(NumOrbitals);
  psiV_host_view.attachReference(psiV.data(), NumOrbitals);
  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
}

template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::GradType DiracDeterminantBatched<DET_ENGINE_TYPE>::evalGrad(
    ParticleSet& P,
    int iat)
{
  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  GradType g = simd::dot(psiMinv[WorkingIndex], dpsiM[WorkingIndex], psiMinv.rows());
  RatioTimer.stop();
  return g;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_evalGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                           const RefVector<ParticleSet>& P_list,
                           int iat,
                           std::vector<GradType>& grad_now)
{
  RatioTimer.start();

  const int nw = WFC_list.size();
  std::vector<const ValueType*> invRow_list(nw, nullptr);
  std::vector<const ValueType*> dpsiM_row_list(nw, nullptr);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    invRow_list[iw] = det.psiMinv_dev_ptr + NumOrbitals * WorkingIndex;
    dpsiM_row_list[iw] = det.psiM_vgl_dev_ptr + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
  }

  det_engine_.mw_evalGrad(invRow_list, dpsiM_row_list, NumOrbitals, grad_now);

  RatioTimer.stop();
}

template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::PsiValueType DiracDeterminantBatched<DET_ENGINE_TYPE>::ratioGrad(
    ParticleSet& P,
    int iat,
    GradType& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  SPOVGLTimer.start();
  Phi->evaluateVGL(P, iat, psiV_host_view, dpsiV, d2psiV);
  SPOVGLTimer.stop();

  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  curRatio               = simd::dot(psiMinv[WorkingIndex], psiV.data(), psiV.size());
  grad_iat += static_cast<ValueType>(static_cast<PsiValueType>(1.0) / curRatio) *
      simd::dot(psiMinv[WorkingIndex], dpsiV.data(), psiMinv.rows());
  RatioTimer.stop();
  return curRatio;
}


template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_ratioGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                                                            const RefVector<ParticleSet>& P_list,
                                                            int iat,
                                                            std::vector<PsiValueType>& ratios,
                                                            std::vector<GradType>& grad_new)
{
  SPOVGLTimer.start();
  RefVector<SPOSet> phi_list;
  phi_list.reserve(WFC_list.size());

  const int WorkingIndex = iat - FirstIndex;
  std::vector<const ValueType*> psiMinv_row_dev_ptr_list(WFC_list.size(), nullptr);
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    phi_list.push_back(*det.Phi);
    if (Phi->isOMPoffload())
      psiMinv_row_dev_ptr_list[iw] = det.psiMinv_dev_ptr + NumOrbitals * WorkingIndex;
    else
    {
      psiMinv_row_dev_ptr_list[iw] = det.psiMinv.data() + NumOrbitals * WorkingIndex;
      auto* Ainv_ptr = det.psiMinv.data();
      PRAGMA_OFFLOAD("omp target update from(Ainv_ptr[NumOrbitals*WorkingIndex:NumOrbitals])")
    }
  }

  resizeMultiWalkerScratch(psiMinv.cols(), WFC_list.size());
  ratios_local.resize(WFC_list.size());
  grad_new_local.resize(WFC_list.size());

  VectorSoaContainer<ValueType, DIM + 2> phi_vgl_v_view(phi_vgl_v.data(), phi_vgl_v.size(), phi_vgl_v.capacity());
  Phi->mw_evaluateVGLandDetRatioGrads(phi_list, P_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v_view, ratios_local,
                                      grad_new_local);
  SPOVGLTimer.stop();

  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det      = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    det.UpdateMode = ORB_PBYP_PARTIAL;
    ratios[iw] = det.curRatio = ratios_local[iw];
    grad_new[iw] += grad_new_local[iw];
  }
}


/** move was accepted, update the real container
*/
template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  const int WorkingIndex = iat - FirstIndex;
  LogValue += convertValueToLog(curRatio);
  UpdateTimer.start();
  det_engine_.updateRow(psiMinv, WorkingIndex, psiV, curRatio);
  if (UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsiM[WorkingIndex], dpsiV.data(), NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  }
  UpdateTimer.stop();
  curRatio = 1.0;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_accept_rejectMove(const RefVector<WaveFunctionComponent>& WFC_list,
                                                                    const RefVector<ParticleSet>& P_list,
                                                                    int iat,
                                                                    const std::vector<bool>& isAccepted,
                                                                    bool safe_to_delay)
{
  UpdateTimer.start();

  const int nw = WFC_list.size();
  int count    = 0;
  for (int iw = 0; iw < nw; iw++)
    if (isAccepted[iw])
      count++;
  const int n_accepted = count;
  std::vector<ValueType*> psiMinv_dev_ptr_list(n_accepted, nullptr);
  std::vector<ValueType*> psiM_g_dev_ptr_list(n_accepted, nullptr);
  std::vector<ValueType*> psiM_l_dev_ptr_list(n_accepted, nullptr);

  const int WorkingIndex = iat - FirstIndex;
  for (int iw = 0, count = 0; iw < nw; iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    if (isAccepted[iw])
    {
      psiMinv_dev_ptr_list[count] = det.psiMinv_dev_ptr;
      psiM_g_dev_ptr_list[count] = det.psiM_vgl_dev_ptr + psiM_vgl.capacity() + NumOrbitals * WorkingIndex * DIM;
      psiM_l_dev_ptr_list[count] = det.psiM_vgl_dev_ptr + psiM_vgl.capacity() * 4 + NumOrbitals * WorkingIndex;
      det.LogValue += convertValueToLog(det.curRatio);
      count++;
    }
    det.curRatio = 1.0;
  }

  if (!Phi->isOMPoffload() && n_accepted > 0)
  {
    auto* phi_vgl_v_ptr = phi_vgl_v.data();
    PRAGMA_OFFLOAD("omp target update to(phi_vgl_v_ptr[:phi_vgl_v.capacity()*5])")
  }

  det_engine_.mw_updateRow(psiMinv_dev_ptr_list, psiM_g_dev_ptr_list, psiM_l_dev_ptr_list, psiMinv.rows(), WorkingIndex, isAccepted, phi_vgl_v_dev_ptr, phi_vgl_v.capacity(), ratios_local);

  UpdateTimer.stop();
}

/** move was rejected. copy the real container to the temporary to move on
*/
template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::restore(int iat)
{
  curRatio = 1.0;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::completeUpdates()
{
  UpdateTimer.start();
  /// no action here because single walker code path keep Ainv, dpsiM, d2psiM up to date on the host.
  UpdateTimer.stop();
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_completeUpdates(const RefVector<WaveFunctionComponent>& WFC_list)
{
  UpdateTimer.start();
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    auto* Ainv_ptr = det.psiMinv.data();
    PRAGMA_OFFLOAD("omp target update from(Ainv_ptr[:psiMinv.size()])")

    auto& my_psiM_vgl = det.psiM_vgl;
    auto* psiM_vgl_ptr = my_psiM_vgl.data();
    PRAGMA_OFFLOAD("omp target update from(psiM_vgl_ptr[my_psiM_vgl.capacity():my_psiM_vgl.capacity()*4])")
  }
  UpdateTimer.stop();
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::updateAfterSweep(ParticleSet& P,
                                                                ParticleSet::ParticleGradient_t& G,
                                                                ParticleSet::ParticleLaplacian_t& L)
{
  if (UpdateMode == ORB_PBYP_RATIO)
  { //need to compute dpsiM and d2psiM. Do not touch psiM!
    SPOVGLTimer.start();
    Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
    //FIXME maybe need the same transfer as recompute
    SPOVGLTimer.stop();
  }

  if (NumPtcls == 1)
  {
    ValueType y = psiMinv(0, 0);
    GradType rv = y * dpsiM(0, 0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y * d2psiM(0, 0) - dot(rv, rv);
  }
  else
  {
    for (size_t i = 0, iat = FirstIndex; i < NumPtcls; ++i, ++iat)
    {
      mValueType dot_temp = simd::dot(psiMinv[i], d2psiM[i], NumOrbitals);
      mGradType rv        = simd::dot(psiMinv[i], dpsiM[i], NumOrbitals);
      G[iat] += rv;
      L[iat] += dot_temp - dot(rv, rv);
    }
  }
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::registerData(ParticleSet& P, WFBufferType& buf)
{
  buf.add(LogValue);
}

template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::LogValueType DiracDeterminantBatched<DET_ENGINE_TYPE>::updateBuffer(
    ParticleSet& P,
    WFBufferType& buf,
    bool fromscratch)
{
  if (fromscratch)
  {
    LogValue = evaluateLog(P, P.G, P.L);
  }
  else
  {
    updateAfterSweep(P, P.G, P.L);
  }
  BufferTimer.start();
  buf.put(LogValue);
  BufferTimer.stop();
  return LogValue;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BufferTimer.start();
  recompute(P);
  buf.get(LogValue);
  BufferTimer.stop();
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::PsiValueType DiracDeterminantBatched<DET_ENGINE_TYPE>::ratio(
    ParticleSet& P,
    int iat)
{
  UpdateMode             = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - FirstIndex;
  SPOVTimer.start();
  Phi->evaluateValue(P, iat, psiV_host_view);
  SPOVTimer.stop();
  RatioTimer.start();
  curRatio = simd::dot(psiMinv[WorkingIndex], psiV.data(), psiV.size());
  RatioTimer.stop();
  return curRatio;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_calcRatio(const RefVector<WaveFunctionComponent>& WFC_list,
                                                            const RefVector<ParticleSet>& P_list,
                                                            int iat,
                                                            std::vector<PsiValueType>& ratios)
{
  SPOVTimer.start();
  RefVector<SPOSet> phi_list;
  phi_list.reserve(WFC_list.size());

  const int WorkingIndex = iat - FirstIndex;
  std::vector<const ValueType*> psiMinv_row_dev_ptr_list(WFC_list.size(), nullptr);
  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    phi_list.push_back(*det.Phi);
    if (Phi->isOMPoffload())
      psiMinv_row_dev_ptr_list[iw] = det.psiMinv_dev_ptr + NumOrbitals * WorkingIndex;
    else
    {
      psiMinv_row_dev_ptr_list[iw] = det.psiMinv.data() + NumOrbitals * WorkingIndex;
      auto* Ainv_ptr = det.psiMinv.data();
      PRAGMA_OFFLOAD("omp target update from(Ainv_ptr[NumOrbitals*WorkingIndex:NumOrbitals])")
    }
  }

  resizeMultiWalkerScratch(psiMinv.cols(), WFC_list.size());
  ratios_local.resize(WFC_list.size());
  grad_new_local.resize(WFC_list.size());

  VectorSoaContainer<ValueType, DIM + 2> phi_vgl_v_view(phi_vgl_v.data(), phi_vgl_v.size(), phi_vgl_v.capacity());

  // calling Phi->mw_evaluateVGLandDetRatioGrads is a temporary workaround.
  // We may implement mw_evaluateVandDetRatio in the future.
  Phi->mw_evaluateVGLandDetRatioGrads(phi_list, P_list, iat, psiMinv_row_dev_ptr_list, phi_vgl_v_view, ratios_local,
                                      grad_new_local);
  SPOVTimer.stop();

  for (int iw = 0; iw < WFC_list.size(); iw++)
  {
    auto& det      = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(WFC_list[iw].get());
    det.UpdateMode = ORB_PBYP_PARTIAL;
    ratios[iw] = det.curRatio = ratios_local[iw];
  }
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::evaluateRatios(const VirtualParticleSet& VP,
                                                              std::vector<ValueType>& ratios)
{
  RatioTimer.start();
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  std::copy_n(psiMinv[WorkingIndex], d2psiV.size(), d2psiV.data());
  RatioTimer.stop();
  SPOVTimer.start();
  Phi->evaluateDetRatios(VP, psiV_host_view, d2psiV, ratios);
  SPOVTimer.stop();
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::mw_evaluateRatios(const RefVector<WaveFunctionComponent>& wfc_list,
                                                                 const RefVector<const VirtualParticleSet>& vp_list,
                                                                 std::vector<std::vector<ValueType>>& ratios)
{
  RatioTimer.start();
  const size_t nw = wfc_list.size();

  RefVector<SPOSet> phi_list;
  RefVector<ValueVector_t> psiV_list;
  RefVector<const ValueVector_t> invRow_list;
  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  invRow_list.reserve(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det = static_cast<DiracDeterminantBatched<DET_ENGINE_TYPE>&>(wfc_list[iw].get());
    const VirtualParticleSet& vp(vp_list[iw]);
    const int WorkingIndex = vp.refPtcl - FirstIndex;
    std::copy_n(det.psiMinv[WorkingIndex], det.d2psiV.size(), det.d2psiV.data());
    // build lists
    phi_list.push_back(*det.Phi);
    psiV_list.push_back(det.psiV_host_view);
    invRow_list.push_back(det.d2psiV);
  }
  RatioTimer.stop();

  SPOVTimer.start();
  Phi->mw_evaluateDetRatios(phi_list, vp_list, psiV_list, invRow_list, ratios);
  SPOVTimer.stop();
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluateValue(P, -1, psiV_host_view);
  SPOVTimer.stop();
  ValueMatrix_t psiMinv_host_view(psiMinv.data(), psiMinv.rows(), psiMinv.cols());
  MatrixOperators::product(psiMinv_host_view, psiV.data(), &ratios[FirstIndex]);
}


template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::resizeScratchObjectsForIonDerivs()
{
  grad_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_lapl_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_grad_source_psiM.resize(NumPtcls, NumOrbitals);
  phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_Minv.resize(NumPtcls, NumOrbitals);
  lapl_phi_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::resizeMultiWalkerScratch(int norb, int nw)
{
  const size_t total_size = norb * nw;
  if (phi_vgl_v.size() < total_size)
  {
    phi_vgl_v.resize(total_size);
    phi_vgl_v_dev_ptr = getOffloadDevicePtr(phi_vgl_v.data());
  }
}

template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::GradType DiracDeterminantBatched<DET_ENGINE_TYPE>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat)
{
  GradType g(0.0);
  if (Phi->hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    Phi->evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
    g = simd::dot(psiMinv.data(), grad_source_psiM.data(), psiMinv.size());
  }

  return g;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiMinv.rows(), psiMinv.cols());
  //IM A HACK.  Assumes evaluateLog has already been executed.
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, grad_grad_source_psiM);
  invertPsiM(psiM_temp, psiMinv);

  phi_alpha_Minv      = 0.0;
  grad_phi_Minv       = 0.0;
  lapl_phi_Minv       = 0.0;
  grad_phi_alpha_Minv = 0.0;
  //grad_grad_psi.resize(NumPtcls);

  for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
  {
    GradType rv = simd::dot(psiMinv[i], dpsiM[i], NumOrbitals);
    //  HessType hess_tmp=simd::dot(psiMinv[i],grad_grad_source_psiM[i],NumOrbitals);
    HessType hess_tmp;
    hess_tmp           = 0.0;
    hess_tmp           = simd::dot(psiMinv[i], grad_grad_source_psiM[i], NumOrbitals);
    grad_grad_psi[iat] = hess_tmp - outerProduct(rv, rv);
  }
}

template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::GradType DiracDeterminantBatched<DET_ENGINE_TYPE>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat,
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad)
{
  GradType gradPsi(0.0);
  if (Phi->hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    Phi->evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM, grad_grad_source_psiM,
                            grad_lapl_source_psiM);
    // HACK HACK HACK
    // Phi->evaluateVGL(P, FirstIndex, LastIndex, psiMinv, dpsiM, d2psiM);
    // psiM_temp = psiMinv;
    // LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,
    // 			   WorkSpace.data(),Pivot.data(),PhaseValue);
    // for (int i=0; i<NumPtcls; i++)
    //   for (int j=0; j<NumPtcls; j++) {
    // 	double val = 0.0;
    // 	for (int k=0; k<NumPtcls; k++)
    // 	  val += psiMinv(i,k) * psiM_temp(k,j);
    // 	val -= (i == j) ? 1.0 : 0.0;
    // 	if (std::abs(val) > 1.0e-12)
    // 	  std::cerr << "Error in inverse.\n";
    //   }
    // for (int i=0; i<NumPtcls; i++) {
    //   P.G[FirstIndex+i] = GradType();
    //   for (int j=0; j<NumOrbitals; j++)
    // 	P.G[FirstIndex+i] += psiMinv(i,j)*dpsiM(i,j);
    // }
    // Compute matrices
    phi_alpha_Minv      = 0.0;
    grad_phi_Minv       = 0.0;
    lapl_phi_Minv       = 0.0;
    grad_phi_alpha_Minv = 0.0;
    for (int i = 0; i < NumPtcls; i++)
      for (int j = 0; j < NumOrbitals; j++)
      {
        lapl_phi_Minv(i, j) = 0.0;
        for (int k = 0; k < NumOrbitals; k++)
          lapl_phi_Minv(i, j) += d2psiM(i, k) * psiMinv(j, k);
      }
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      for (int i = 0; i < NumPtcls; i++)
        for (int j = 0; j < NumOrbitals; j++)
        {
          for (int k = 0; k < NumOrbitals; k++)
          {
            phi_alpha_Minv(i, j)[dim] += grad_source_psiM(i, k)[dim] * psiMinv(j, k);
            grad_phi_Minv(i, j)[dim] += dpsiM(i, k)[dim] * psiMinv(j, k);
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              grad_phi_alpha_Minv(i, j)(dim, dim_el) += grad_grad_source_psiM(i, k)(dim, dim_el) * psiMinv(j, k);
          }
        }
    }
    for (int i = 0, iel = FirstIndex; i < NumPtcls; i++, iel++)
    {
      HessType dval(0.0);
      GradType d2val(0.0);
      for (int dim = 0; dim < OHMMS_DIM; dim++)
        for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
          dval(dim, dim_el) = grad_phi_alpha_Minv(i, i)(dim, dim_el);
      for (int j = 0; j < NumOrbitals; j++)
      {
        gradPsi += grad_source_psiM(i, j) * psiMinv(i, j);
        for (int dim = 0; dim < OHMMS_DIM; dim++)
          for (int k = 0; k < OHMMS_DIM; k++)
            dval(dim, k) -= phi_alpha_Minv(j, i)[dim] * grad_phi_Minv(i, j)[k];
      }
      for (int dim = 0; dim < OHMMS_DIM; dim++)
      {
        for (int k = 0; k < OHMMS_DIM; k++)
          grad_grad[dim][iel][k] += dval(dim, k);
        for (int j = 0; j < NumOrbitals; j++)
        {
          // First term, eq 9
          lapl_grad[dim][iel] += grad_lapl_source_psiM(i, j)[dim] * psiMinv(i, j);
          // Second term, eq 9
          if (j == i)
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              lapl_grad[dim][iel] -=
                  (RealType)2.0 * grad_phi_alpha_Minv(j, i)(dim, dim_el) * grad_phi_Minv(i, j)[dim_el];
          // Third term, eq 9
          // First term, eq 10
          lapl_grad[dim][iel] -= phi_alpha_Minv(j, i)[dim] * lapl_phi_Minv(i, j);
          // Second term, eq 11
          for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
            lapl_grad[dim][iel] +=
                (RealType)2.0 * phi_alpha_Minv(j, i)[dim] * grad_phi_Minv(i, i)[dim_el] * grad_phi_Minv(i, j)[dim_el];
        }
      }
    }
  }
  return gradPsi;
}


/** Calculate the log value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
template<typename DET_ENGINE_TYPE>
typename DiracDeterminantBatched<DET_ENGINE_TYPE>::LogValueType DiracDeterminantBatched<DET_ENGINE_TYPE>::evaluateLog(
    ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L)
{
  recompute(P);

  if (NumPtcls == 1)
  {
    ValueType y = psiMinv(0, 0);
    GradType rv = y * dpsiM(0, 0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y * d2psiM(0, 0) - dot(rv, rv);
  }
  else
  {
    for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
    {
      mGradType rv   = simd::dot(psiMinv[i], dpsiM[i], NumOrbitals);
      mValueType lap = simd::dot(psiMinv[i], d2psiM[i], NumOrbitals);
      G[iat] += rv;
      L[iat] += lap - dot(rv, rv);
    }
  }
  return LogValue;
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::recompute(ParticleSet& P)
{
  SPOVGLTimer.start();
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
  // mw_evaluate_notranspose will be needed. if Phi supports offload, it only guarantees device ready in results.
  // if Phi is not offload. A transfer of dpsiM, d2psiM to device is needed.
  // now evaluate_notranspose only guarantees host. So always transfer.
  //if(!Phi->isOMPoffload())
  {
    auto* psiM_vgl_ptr = psiM_vgl.data();
    PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[psiM_vgl.capacity():psiM_vgl.capacity()*4])")
  }

  SPOVGLTimer.stop();
  if (NumPtcls == 1)
  {
    ValueType det = psiM_temp(0, 0);
    psiMinv(0, 0) = RealType(1) / det;
    LogValue      = convertValueToLog(det);
  }
  else
  {
    invertPsiM(psiM_temp, psiMinv);
  }
}

template<typename DET_ENGINE_TYPE>
void DiracDeterminantBatched<DET_ENGINE_TYPE>::evaluateDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& active,
                                                                   std::vector<ValueType>& dlogpsi,
                                                                   std::vector<ValueType>& dhpsioverpsi)
{
  Phi->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi, FirstIndex, LastIndex);
}

template<typename DET_ENGINE_TYPE>
DiracDeterminantBatched<DET_ENGINE_TYPE>* DiracDeterminantBatched<DET_ENGINE_TYPE>::makeCopy(SPOSetPtr spo) const
{
  DiracDeterminantBatched<DET_ENGINE_TYPE>* dclone = new DiracDeterminantBatched<DET_ENGINE_TYPE>(spo);
  dclone->set(FirstIndex, LastIndex - FirstIndex);
  return dclone;
}

template class DiracDeterminantBatched<>;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
template class DiracDeterminantBatched<MatrixUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus

// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "DiracDeterminant.h"
#include <stdexcept>
#include "CPU/BLAS.hpp"
#include "CPU/SIMD/simd.hpp"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
template<typename DU_TYPE>
DiracDeterminant<DU_TYPE>::DiracDeterminant(SPOSetPtr const spos, int first)
    : DiracDeterminantBase("DiracDeterminant", spos, first), ndelay(1), invRow_id(-1)
{}

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::set(int first, int nel, int delay)
{
  FirstIndex = first;
  ndelay     = delay;

  resize(nel, nel);

  if (Optimizable)
    Phi->buildOptVariables(nel);

  if (Phi->getOrbitalSetSize() < nel)
  {
    std::ostringstream err_msg;
    err_msg << "The SPOSet " << Phi->getName() << " only has " << Phi->getOrbitalSetSize() << " orbitals "
            << "but this determinant needs at least " << nel << std::endl;
    throw std::runtime_error(err_msg.str());
  }
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat)
{
  InverseTimer.start();
  updateEng.invert_transpose(logdetT, invMat, LogValue);
  InverseTimer.stop();
}


///reset the size: with the number of particles and number of orbtials
template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::resize(int nel, int morb)
{
  if (Bytes_in_WFBuffer > 0)
    throw std::runtime_error("DiracDeterimnant just went out of sync with buffer");
  int norb = morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)
  updateEng.resize(norb, ndelay);
  psiM.resize(nel, norb);
  dpsiM.resize(nel, norb);
  d2psiM.resize(nel, norb);
  psiV.resize(norb);
  invRow.resize(norb);
  psiM_temp.resize(nel, norb);
  LastIndex   = FirstIndex + nel;
  NumPtcls    = nel;
  NumOrbitals = norb;

  dpsiV.resize(NumOrbitals);
  dspin_psiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  FirstAddressOfdV = &(dpsiM(0, 0)[0]); //(*dpsiM.begin())[0]);
  LastAddressOfdV  = FirstAddressOfdV + NumPtcls * NumOrbitals * DIM;
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::GradType DiracDeterminant<DU_TYPE>::evalGrad(ParticleSet& P, int iat)
{
  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  invRow_id = WorkingIndex;
  updateEng.getInvRow(psiM, WorkingIndex, invRow);
  GradType g = simd::dot(invRow.data(), dpsiM[WorkingIndex], invRow.size());
  RatioTimer.stop();
#ifndef NDEBUG
  ValueType g_norm = simd::dot(g.data(), g.data(), g.Size);
  if (std::abs(g_norm) < std::abs(std::numeric_limits<QMCTraits::ValueType>::epsilon()))
  {
    std::cerr << "evalGrad gradient is " << g[0] << ' ' << g[1] << ' ' << g[2] << '\n';
    throw std::runtime_error("gradient of zero");
  }
#endif
  return g;
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::GradType DiracDeterminant<DU_TYPE>::evalGradWithSpin(ParticleSet& P,
                                                                                         int iat,
                                                                                         ComplexType& spingrad)
{
  Phi->evaluate_spin(P, iat, psiV, dspin_psiV);
  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  invRow_id = WorkingIndex;
  updateEng.getInvRow(psiM, WorkingIndex, invRow);
  GradType g         = simd::dot(invRow.data(), dpsiM[WorkingIndex], invRow.size());
  ComplexType spin_g = simd::dot(invRow.data(), dspin_psiV.data(), invRow.size());
  RatioTimer.stop();

  spingrad += spin_g;
  return g;
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::PsiValueType DiracDeterminant<DU_TYPE>::ratioGrad(ParticleSet& P,
                                                                                      int iat,
                                                                                      GradType& grad_iat)
{
  SPOVGLTimer.start();
  Phi->evaluateVGL(P, iat, psiV, dpsiV, d2psiV);
  SPOVGLTimer.stop();
  return ratioGrad_compute(iat, grad_iat);
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::PsiValueType DiracDeterminant<DU_TYPE>::ratioGrad_compute(int iat,
                                                                                              GradType& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;
  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  // This is an optimization.
  // check invRow_id against WorkingIndex to see if getInvRow() has been called already
  // Some code paths call evalGrad before calling ratioGrad.
  if (invRow_id != WorkingIndex)
  {
    invRow_id = WorkingIndex;
    updateEng.getInvRow(psiM, WorkingIndex, invRow);
  }
  curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
  grad_iat += static_cast<ValueType>(static_cast<PsiValueType>(1.0) / curRatio) *
      simd::dot(invRow.data(), dpsiV.data(), invRow.size());
  RatioTimer.stop();
  return curRatio;
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::PsiValueType DiracDeterminant<DU_TYPE>::ratioGradWithSpin(ParticleSet& P,
                                                                                              int iat,
                                                                                              GradType& grad_iat,
                                                                                              ComplexType& spingrad_iat)
{
  SPOVGLTimer.start();
  Phi->evaluateVGL(P, iat, psiV, dpsiV, d2psiV);
  Phi->evaluate_spin(P, iat, psiV, dspin_psiV);
  SPOVGLTimer.stop();

  UpdateMode = ORB_PBYP_PARTIAL;
  RatioTimer.start();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  // This is an optimization.
  // check invRow_id against WorkingIndex to see if getInvRow() has been called already
  // Some code paths call evalGrad before calling ratioGrad.
  if (invRow_id != WorkingIndex)
  {
    invRow_id = WorkingIndex;
    updateEng.getInvRow(psiM, WorkingIndex, invRow);
  }
  curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
  grad_iat += static_cast<ValueType>(static_cast<PsiValueType>(1.0) / curRatio) *
      simd::dot(invRow.data(), dpsiV.data(), invRow.size());

  spingrad_iat += static_cast<ValueType>(static_cast<PsiValueType>(1.0) / curRatio) *
      simd::dot(invRow.data(), dspin_psiV.data(), invRow.size());
  RatioTimer.stop();

  return curRatio;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::mw_ratioGrad(const RefVector<WaveFunctionComponent>& WFC_list,
                                             const RefVector<ParticleSet>& P_list,
                                             int iat,
                                             std::vector<PsiValueType>& ratios,
                                             std::vector<GradType>& grad_new)
{
  SPOVGLTimer.start();
  RefVector<SPOSet> phi_list;
  phi_list.reserve(WFC_list.size());
  RefVector<ValueVector_t> psi_v_list;
  psi_v_list.reserve(WFC_list.size());
  RefVector<GradVector_t> dpsi_v_list;
  dpsi_v_list.reserve(WFC_list.size());
  RefVector<ValueVector_t> d2psi_v_list;
  d2psi_v_list.reserve(WFC_list.size());

  for (WaveFunctionComponent& wfc : WFC_list)
  {
    auto& det = static_cast<DiracDeterminant<DU_TYPE>&>(wfc);
    phi_list.push_back(*det.Phi);
    psi_v_list.push_back(det.psiV);
    dpsi_v_list.push_back(det.dpsiV);
    d2psi_v_list.push_back(det.d2psiV);
  }

  Phi->mw_evaluateVGL(phi_list, P_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list);
  SPOVGLTimer.stop();

  for (int iw = 0; iw < WFC_list.size(); iw++)
    ratios[iw] = static_cast<DiracDeterminant<DU_TYPE>&>(WFC_list[iw].get()).ratioGrad_compute(iat, grad_new[iw]);
}


/** move was accepted, update the real container
*/
template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  LogValue += convertValueToLog(curRatio);
  UpdateTimer.start();
  updateEng.acceptRow(psiM, WorkingIndex, psiV);
  if (!safe_to_delay)
    updateEng.updateInvMat(psiM);
  // invRow becomes invalid after accepting a move
  invRow_id = -1;
  if (UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsiM[WorkingIndex], dpsiV.data(), NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  }
  UpdateTimer.stop();
  curRatio = 1.0;
}

/** move was rejected. copy the real container to the temporary to move on
*/
template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::restore(int iat)
{
  curRatio = 1.0;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::completeUpdates()
{
  UpdateTimer.start();
  // invRow becomes invalid after updating the inverse matrix
  invRow_id = -1;
  updateEng.updateInvMat(psiM);
  UpdateTimer.stop();
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::updateAfterSweep(ParticleSet& P,
                                                 ParticleSet::ParticleGradient_t& G,
                                                 ParticleSet::ParticleLaplacian_t& L)
{
  if (UpdateMode == ORB_PBYP_RATIO)
  { //need to compute dpsiM and d2psiM. Do not touch psiM!
    SPOVGLTimer.start();
    Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
    SPOVGLTimer.stop();
  }

  if (NumPtcls == 1)
  {
    ValueType y = psiM(0, 0);
    GradType rv = y * dpsiM(0, 0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y * d2psiM(0, 0) - dot(rv, rv);
  }
  else
  {
    for (size_t i = 0, iat = FirstIndex; i < NumPtcls; ++i, ++iat)
    {
      mValueType dot_temp = simd::dot(psiM[i], d2psiM[i], NumOrbitals);
      mGradType rv        = simd::dot(psiM[i], dpsiM[i], NumOrbitals);
      G[iat] += rv;
      L[iat] += dot_temp - dot(rv, rv);
    }
  }
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::registerData(ParticleSet& P, WFBufferType& buf)
{
  if (Bytes_in_WFBuffer == 0)
  {
    //add the data: inverse, gradient and laplacian
    Bytes_in_WFBuffer = buf.current();
    buf.add(psiM.first_address(), psiM.last_address());
    buf.add(FirstAddressOfdV, LastAddressOfdV);
    buf.add(d2psiM.first_address(), d2psiM.last_address());
    Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
    // free local space
    psiM.free();
    dpsiM.free();
    d2psiM.free();
  }
  else
  {
    buf.forward(Bytes_in_WFBuffer);
#ifndef NDEBUG
    // this causes too much output in the legacy code.
    // \todo turn this back on after legacy is dropped,
    // I don't think it should print at all in the new design
    // std::cerr << ("You really should know whether you have registered this objects data previously!, consider this an error in the unified code");
#endif
  }
  buf.add(LogValue);
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::LogValueType DiracDeterminant<DU_TYPE>::updateBuffer(ParticleSet& P,
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
  buf.forward(Bytes_in_WFBuffer);
  buf.put(LogValue);
  BufferTimer.stop();
  return LogValue;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BufferTimer.start();
  psiM.attachReference(buf.lendReference<ValueType>(psiM.size()));
  dpsiM.attachReference(buf.lendReference<GradType>(dpsiM.size()));
  d2psiM.attachReference(buf.lendReference<ValueType>(d2psiM.size()));
  buf.get(LogValue);
  // start with invRow labelled invalid
  invRow_id = -1;
  updateEng.initializeInv(psiM);
  BufferTimer.stop();
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::PsiValueType DiracDeterminant<DU_TYPE>::ratio(ParticleSet& P, int iat)
{
  UpdateMode             = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  SPOVTimer.start();
  Phi->evaluateValue(P, iat, psiV);
  SPOVTimer.stop();
  RatioTimer.start();
  // This is an optimization.
  // check invRow_id against WorkingIndex to see if getInvRow() has been called
  // This is intended to save redundant compuation in TM1 and TM3
  if (invRow_id != WorkingIndex)
  {
    invRow_id = WorkingIndex;
    updateEng.getInvRow(psiM, WorkingIndex, invRow);
  }
  curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
  RatioTimer.stop();
  return curRatio;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  RatioTimer.start();
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  assert(WorkingIndex >= 0);
  std::copy_n(psiM[WorkingIndex], invRow.size(), invRow.data());
  RatioTimer.stop();
  SPOVTimer.start();
  Phi->evaluateDetRatios(VP, psiV, invRow, ratios);
  SPOVTimer.stop();
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::mw_evaluateRatios(const RefVector<WaveFunctionComponent>& wfc_list,
                                                  const RefVector<const VirtualParticleSet>& vp_list,
                                                  std::vector<std::vector<ValueType>>& ratios)
{
  RatioTimer.start();
  const size_t nw = wfc_list.size();

  RefVector<SPOSet> phi_list;
  RefVector<ValueVector_t> psiV_list;
  std::vector<const ValueType*> invRow_ptr_list;
  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  invRow_ptr_list.reserve(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det = static_cast<DiracDeterminant<DU_TYPE>&>(wfc_list[iw].get());
    const VirtualParticleSet& vp(vp_list[iw]);
    const int WorkingIndex = vp.refPtcl - FirstIndex;
    assert(WorkingIndex >= 0);
    // If DiracDeterminant is in a valid state this copy_n is not necessary.
    // That is at minimum a call to evaluateLog and ...
    // std::copy_n(det.psiM[WorkingIndex], det.invRow.s.ize(), det.invRow.data());
    // build lists
    phi_list.push_back(*det.Phi);
    psiV_list.push_back(det.psiV);
    invRow_ptr_list.push_back(det.psiM[WorkingIndex]);
  }
  RatioTimer.stop();

  SPOVTimer.start();
  // Phi->isOMPoffload() requires device invRow pointers for mw_evaluateDetRatios.
  // evaluateDetRatios only requires host invRow pointers.
  if (Phi->isOMPoffload())
    for (int iw = 0; iw < phi_list.size(); iw++)
    {
      Vector<ValueType> invRow(const_cast<ValueType*>(invRow_ptr_list[iw]), psiV_list[iw].get().size());
      phi_list[iw].get().evaluateDetRatios(vp_list[iw], psiV_list[iw], invRow, ratios[iw]);
    }
  else
    Phi->mw_evaluateDetRatios(phi_list, vp_list, psiV_list, invRow_ptr_list, ratios);
  SPOVTimer.stop();
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluateValue(P, -1, psiV);
  SPOVTimer.stop();
  MatrixOperators::product(psiM, psiV.data(), &ratios[FirstIndex]);
}


template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::resizeScratchObjectsForIonDerivs()
{
  grad_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_lapl_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_grad_source_psiM.resize(NumPtcls, NumOrbitals);
  phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_Minv.resize(NumPtcls, NumOrbitals);
  lapl_phi_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::GradType DiracDeterminant<DU_TYPE>::evalGradSource(ParticleSet& P,
                                                                                       ParticleSet& source,
                                                                                       int iat)
{
  GradType g(0.0);
  if (Phi->hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    Phi->evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
    g = simd::dot(psiM.data(), grad_source_psiM.data(), psiM.size());
  }

  return g;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiM.rows(), psiM.cols());
  //IM A HACK.  Assumes evaluateLog has already been executed.
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, grad_grad_source_psiM);
  invertPsiM(psiM_temp, psiM);

  phi_alpha_Minv      = 0.0;
  grad_phi_Minv       = 0.0;
  lapl_phi_Minv       = 0.0;
  grad_phi_alpha_Minv = 0.0;
  //grad_grad_psi.resize(NumPtcls);

  for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
  {
    GradType rv = simd::dot(psiM[i], dpsiM[i], NumOrbitals);
    //  HessType hess_tmp=simd::dot(psiM[i],grad_grad_source_psiM[i],NumOrbitals);
    HessType hess_tmp;
    hess_tmp           = 0.0;
    hess_tmp           = simd::dot(psiM[i], grad_grad_source_psiM[i], NumOrbitals);
    grad_grad_psi[iat] = hess_tmp - outerProduct(rv, rv);
  }
}

template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::GradType DiracDeterminant<DU_TYPE>::evalGradSource(
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
    // Phi->evaluateVGL(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
    // psiM_temp = psiM;
    // LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
    // 			   WorkSpace.data(),Pivot.data(),PhaseValue);
    // for (int i=0; i<NumPtcls; i++)
    //   for (int j=0; j<NumPtcls; j++) {
    // 	double val = 0.0;
    // 	for (int k=0; k<NumPtcls; k++)
    // 	  val += psiM(i,k) * psiM_temp(k,j);
    // 	val -= (i == j) ? 1.0 : 0.0;
    // 	if (std::abs(val) > 1.0e-12)
    // 	  std::cerr << "Error in inverse.\n";
    //   }
    // for (int i=0; i<NumPtcls; i++) {
    //   P.G[FirstIndex+i] = GradType();
    //   for (int j=0; j<NumOrbitals; j++)
    // 	P.G[FirstIndex+i] += psiM(i,j)*dpsiM(i,j);
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
          lapl_phi_Minv(i, j) += d2psiM(i, k) * psiM(j, k);
      }
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      for (int i = 0; i < NumPtcls; i++)
        for (int j = 0; j < NumOrbitals; j++)
        {
          for (int k = 0; k < NumOrbitals; k++)
          {
            phi_alpha_Minv(i, j)[dim] += grad_source_psiM(i, k)[dim] * psiM(j, k);
            grad_phi_Minv(i, j)[dim] += dpsiM(i, k)[dim] * psiM(j, k);
            for (int dim_el = 0; dim_el < OHMMS_DIM; dim_el++)
              grad_phi_alpha_Minv(i, j)(dim, dim_el) += grad_grad_source_psiM(i, k)(dim, dim_el) * psiM(j, k);
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
        gradPsi += grad_source_psiM(i, j) * psiM(i, j);
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
          lapl_grad[dim][iel] += grad_lapl_source_psiM(i, j)[dim] * psiM(i, j);
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
template<typename DU_TYPE>
typename DiracDeterminant<DU_TYPE>::LogValueType DiracDeterminant<DU_TYPE>::evaluateLog(
    ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L)
{
  recompute(P);

  if (NumPtcls == 1)
  {
    ValueType y = psiM(0, 0);
    GradType rv = y * dpsiM(0, 0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y * d2psiM(0, 0) - dot(rv, rv);
  }
  else
  {
    for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
    {
      mGradType rv   = simd::dot(psiM[i], dpsiM[i], NumOrbitals);
      mValueType lap = simd::dot(psiM[i], d2psiM[i], NumOrbitals);
      G[iat] += rv;
      L[iat] += lap - dot(rv, rv);
    }
  }
  return LogValue;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::recompute(ParticleSet& P)
{
  SPOVGLTimer.start();
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
  SPOVGLTimer.stop();
  if (NumPtcls == 1)
  {
    ValueType det = psiM_temp(0, 0);
    psiM(0, 0)    = RealType(1) / det;
    LogValue      = convertValueToLog(det);
  }
  else
  {
    invertPsiM(psiM_temp, psiM);
  }

  // invRow becomes invalid after updating the inverse matrix
  invRow_id = -1;
}

template<typename DU_TYPE>
void DiracDeterminant<DU_TYPE>::evaluateDerivatives(ParticleSet& P,
                                                    const opt_variables_type& active,
                                                    std::vector<ValueType>& dlogpsi,
                                                    std::vector<ValueType>& dhpsioverpsi)
{
  Phi->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi, FirstIndex, LastIndex);
}

template<typename DU_TYPE>
DiracDeterminant<DU_TYPE>* DiracDeterminant<DU_TYPE>::makeCopy(SPOSetPtr spo) const
{
  DiracDeterminant<DU_TYPE>* dclone = new DiracDeterminant<DU_TYPE>(spo);
  dclone->set(FirstIndex, LastIndex - FirstIndex, ndelay);
  return dclone;
}

template class DiracDeterminant<>;
#if defined(ENABLE_CUDA)
template class DiracDeterminant<DelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>;
#endif

} // namespace qmcplusplus

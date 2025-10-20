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
#include "CPU/SIMD/inner_product.hpp"
#include "Numerics/MatrixOperators.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "QMCWaveFunctions/RotatedSPOs.h"

namespace qmcplusplus
{
/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
template<PlatformKind PL, typename VT, typename FPVT>
DiracDeterminant<PL, VT, FPVT>::DiracDeterminant(SPOSet& phi,
                                                 int first,
                                                 int last,
                                                 int ndelay,
                                                 DetMatInvertor matrix_inverter_kind)
    : DiracDeterminantBase(getClassName(), phi, first, last),
      ndelay_(ndelay),
      invRow_id(-1),
      matrix_inverter_kind_(matrix_inverter_kind)
{
  resize(NumPtcls, NumPtcls);

  RotatedSPOs* rot_spo = dynamic_cast<RotatedSPOs*>(&phi_);
  if (rot_spo)
    rot_spo->buildOptVariables(NumPtcls);
}

template<PlatformKind PL, typename VT, typename FPVT>
DiracDeterminant<PL, VT, FPVT>::~DiracDeterminant()
{
  if (!psiM.isAttached())
    accel_engine_.update_eng_.releaseFromDeviceCopy(psiM);
  accel_engine_.update_eng_.releaseFromDeviceCopy(psiM_temp);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::invertPsiM(const ValueMatrix& logdetT, ValueMatrix& invMat)
{
  ScopedTimer local_timer(InverseTimer);
  if (matrix_inverter_kind_ == DetMatInvertor::ACCEL && AccelEngine<PL, VT, FPVT>::inverter_supported)
  {
    bool success        = false;
    int failure_counter = 0;
    do
    {
      try
      {
        if constexpr (AccelEngine<PL, VT, FPVT>::inverter_supported)
          accel_engine_.update_eng_.invert_transpose(accel_engine_.inverter_, logdetT, invMat, log_value_);
        else
          throw std::runtime_error("DiracDeterminan::invertPsiM report bug invalid code path!");
        if (failure_counter > 0)
        {
          std::ostringstream success_msg;
          success_msg << "Successful rerun matrix inversion on Rank " << OHMMS::Controller->rank() << " Thread "
                      << omp_get_thread_num() << std::endl;
          std::cerr << success_msg.str();
        }
        success = true;
      }
      catch (const std::exception& e)
      {
        failure_counter++;
        std::ostringstream err_msg;
        err_msg << failure_counter << "th matrix inversion on Rank " << OHMMS::Controller->rank() << " Thread "
                << omp_get_thread_num() << " which failed earlier with an error:\n  " << e.what() << std::endl;
        std::cerr << err_msg.str();
        if (failure_counter == 1)
        {
          //record the bad matrix to a file at the first failure
          std::ostringstream matfname;
          matfname << "badmatrix.r" << OHMMS::Controller->rank() << "t" << omp_get_thread_num() << ".txt";
          std::ofstream matfile(matfname.str().c_str(), std::ios::app);
          matfile << std::setprecision(14) << std::scientific;
          for (size_t i = 0; i < logdetT.rows(); i++)
          {
            for (size_t j = 0; j < logdetT.cols(); j++)
              matfile << "  " << logdetT[i][j];
            matfile << std::endl;
          }
        }
      }
    } while (!success && failure_counter < 5); // try 5 times at maximum
    if (!success)
      throw std::runtime_error("Matrix inversion failed after " + std::to_string(failure_counter) + " attempts.\n");
  }
  else
  {
    host_inverter_.invert_transpose(logdetT, invMat, log_value_);
    accel_engine_.update_eng_.initializeInv(invMat);
  }
}


///reset the size: with the number of particles and number of orbtials
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::resize(int nel, int morb)
{
  if (Bytes_in_WFBuffer > 0)
    throw std::runtime_error("DiracDeterimnant just went out of sync with buffer");
  int norb = morb;
  if (norb <= 0)
    norb = nel; // for morb == -1 (default)

  accel_engine_.update_eng_.releaseFromDeviceCopy(psiM);
  accel_engine_.update_eng_.releaseFromDeviceCopy(psiM_temp);

  accel_engine_.update_eng_.resize(norb, ndelay_);
  psiM.resize(nel, norb);
  dpsiM.resize(nel, norb);
  d2psiM.resize(nel, norb);
  psiV.resize(norb);
  invRow.resize(norb);
  psiM_temp.resize(nel, norb);

  accel_engine_.update_eng_.prepareForDeviceCopy(psiM);
  accel_engine_.update_eng_.prepareForDeviceCopy(psiM_temp);

  dpsiV.resize(NumOrbitals);
  dspin_psiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  FirstAddressOfdV = &(dpsiM(0, 0)[0]); //(*dpsiM.begin())[0]);
  LastAddressOfdV  = FirstAddressOfdV + NumPtcls * NumOrbitals * DIM;
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::GradType DiracDeterminant<PL, VT, FPVT>::evalGrad(ParticleSet& P, int iat)
{
  ScopedTimer local_timer(RatioTimer);
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  invRow_id = WorkingIndex;
  accel_engine_.update_eng_.getInvRow(psiM, WorkingIndex, invRow);
  GradType g = simd::dot(invRow.data(), dpsiM[WorkingIndex], invRow.size());
  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::GradType DiracDeterminant<PL, VT, FPVT>::evalGradWithSpin(
    ParticleSet& P,
    int iat,
    ComplexType& spingrad)
{
  phi_.evaluate_spin(P, iat, psiV, dspin_psiV);
  ScopedTimer local_timer(RatioTimer);
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  invRow_id = WorkingIndex;
  accel_engine_.update_eng_.getInvRow(psiM, WorkingIndex, invRow);
  GradType g         = simd::dot(invRow.data(), dpsiM[WorkingIndex], invRow.size());
  ComplexType spin_g = simd::dot(invRow.data(), dspin_psiV.data(), invRow.size());
  spingrad += spin_g;

  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::PsiValue DiracDeterminant<PL, VT, FPVT>::ratioGrad(ParticleSet& P,
                                                                                            int iat,
                                                                                            GradType& grad_iat)
{
  {
    ScopedTimer local_timer(SPOVGLTimer);
    phi_.evaluateVGL(P, iat, psiV, dpsiV, d2psiV);
  }
  return ratioGrad_compute(iat, grad_iat);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::PsiValue DiracDeterminant<PL, VT, FPVT>::ratioGrad_compute(int iat,
                                                                                                    GradType& grad_iat)
{
  ScopedTimer local_timer(RatioTimer);

  UpdateMode             = ORB_PBYP_PARTIAL;
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  // This is an satefy mechanism.
  // check invRow_id against WorkingIndex to see if getInvRow() has been called already
  // when evalGrad has not been called already or the particle id is not consistent,
  // invRow is recomputed.
  if (invRow_id != WorkingIndex)
  {
    invRow_id = WorkingIndex;
    accel_engine_.update_eng_.getInvRow(psiM, WorkingIndex, invRow);
  }
  curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
  grad_iat += static_cast<ValueType>(static_cast<PsiValue>(1.0) / curRatio) *
      simd::dot(invRow.data(), dpsiV.data(), invRow.size());
  return curRatio;
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::PsiValue DiracDeterminant<PL, VT, FPVT>::ratioGradWithSpin(
    ParticleSet& P,
    int iat,
    GradType& grad_iat,
    ComplexType& spingrad_iat)
{
  {
    ScopedTimer local_timer(SPOVGLTimer);
    phi_.evaluateVGL_spin(P, iat, psiV, dpsiV, d2psiV, dspin_psiV);
  }

  {
    ScopedTimer local_timer(RatioTimer);
    UpdateMode             = ORB_PBYP_PARTIAL;
    const int WorkingIndex = iat - FirstIndex;
    assert(WorkingIndex >= 0);
    // This is an optimization.
    // check invRow_id against WorkingIndex to see if getInvRow() has been called already
    // Some code paths call evalGrad before calling ratioGrad.
    if (invRow_id != WorkingIndex)
    {
      invRow_id = WorkingIndex;
      accel_engine_.update_eng_.getInvRow(psiM, WorkingIndex, invRow);
    }
    curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
    grad_iat += static_cast<ValueType>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(invRow.data(), dpsiV.data(), invRow.size());

    spingrad_iat += static_cast<ValueType>(static_cast<PsiValue>(1.0) / curRatio) *
        simd::dot(invRow.data(), dspin_psiV.data(), invRow.size());
  }

  return curRatio;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::mw_ratioGrad(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                  const RefVectorWithLeader<ParticleSet>& p_list,
                                                  int iat,
                                                  std::vector<PsiValue>& ratios,
                                                  std::vector<GradType>& grad_new) const
{
  {
    ScopedTimer local_timer(SPOVGLTimer);
    RefVectorWithLeader<SPOSet> phi_list(phi_);
    phi_list.reserve(wfc_list.size());
    RefVector<ValueVector> psi_v_list;
    psi_v_list.reserve(wfc_list.size());
    RefVector<GradVector> dpsi_v_list;
    dpsi_v_list.reserve(wfc_list.size());
    RefVector<ValueVector> d2psi_v_list;
    d2psi_v_list.reserve(wfc_list.size());

    for (WaveFunctionComponent& wfc : wfc_list)
    {
      auto& det = static_cast<DiracDeterminant<PL, VT, FPVT>&>(wfc);
      phi_list.push_back(det.phi_);
      psi_v_list.push_back(det.psiV);
      dpsi_v_list.push_back(det.dpsiV);
      d2psi_v_list.push_back(det.d2psiV);
    }

    phi_.mw_evaluateVGL(phi_list, p_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list);
  }

  for (int iw = 0; iw < wfc_list.size(); iw++)
    ratios[iw] = wfc_list.getCastedElement<DiracDeterminant<PL, VT, FPVT>>(iw).ratioGrad_compute(iat, grad_new[iw]);
}

/** move was accepted, update the real container
*/
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  if (curRatio == PsiValue(0))
  {
    std::ostringstream msg;
    msg << "DiracDeterminant::acceptMove curRatio is " << curRatio << "! Report a bug." << std::endl;
    throw std::runtime_error(msg.str());
  }
  ScopedTimer local_timer(UpdateTimer);
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  log_value_ += convertValueToLog(curRatio);
  accel_engine_.update_eng_.acceptRow(psiM, WorkingIndex, psiV, curRatio);
  if (!safe_to_delay)
    accel_engine_.update_eng_.updateInvMat(psiM);
  // invRow becomes invalid after accepting a move
  invRow_id = -1;
  if (UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsiM[WorkingIndex], dpsiV.data(), NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  }
  curRatio = 1.0;
}

/** move was rejected. copy the real container to the temporary to move on
*/
template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::restore(int iat)
{
  curRatio = 1.0;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::completeUpdates()
{
  ScopedTimer local_timer(UpdateTimer);
  // invRow becomes invalid after updating the inverse matrix
  invRow_id = -1;
  accel_engine_.update_eng_.updateInvMat(psiM);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::updateAfterSweep(const ParticleSet& P,
                                                      ParticleSet::ParticleGradient& G,
                                                      ParticleSet::ParticleLaplacian& L)
{
  if (UpdateMode == ORB_PBYP_RATIO)
  { //need to compute dpsiM and d2psiM. Do not touch psiM!
    ScopedTimer local_timer(SPOVGLTimer);
    phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
    UpdateMode = ORB_WALKER;
  }

  for (size_t i = 0, iat = FirstIndex; i < NumPtcls; ++i, ++iat)
  {
    mValueType dot_temp = simd::dot(psiM[i], d2psiM[i], NumOrbitals);
    mGradType rv        = simd::dot(psiM[i], dpsiM[i], NumOrbitals);
    G[iat] += rv;
    L[iat] += dot_temp - dot(rv, rv);
  }
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::registerData(ParticleSet& P, WFBufferType& buf)
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
    accel_engine_.update_eng_.releaseFromDeviceCopy(psiM);
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
  buf.add(log_value_);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::LogValue DiracDeterminant<PL, VT, FPVT>::evaluateGL(
    const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L,
    bool fromscratch)
{
  if (fromscratch)
    evaluateLog(P, G, L);
  else
    updateAfterSweep(P, G, L);
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::LogValue DiracDeterminant<PL, VT, FPVT>::updateBuffer(ParticleSet& P,
                                                                                               WFBufferType& buf,
                                                                                               bool fromscratch)
{
  evaluateGL(P, P.G, P.L, fromscratch);
  {
    ScopedTimer local_timer(BufferTimer);
    buf.forward(Bytes_in_WFBuffer);
    buf.put(log_value_);
  }
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  ScopedTimer local_timer(BufferTimer);
  psiM.attachReference(buf.lendReference<ValueType>(NumPtcls * NumPtcls), NumPtcls, NumPtcls);
  dpsiM.attachReference(buf.lendReference<GradType>(NumPtcls * NumPtcls), NumPtcls, NumPtcls);
  d2psiM.attachReference(buf.lendReference<ValueType>(NumPtcls * NumPtcls), NumPtcls, NumPtcls);
  buf.get(log_value_);
  // start with invRow labelled invalid
  invRow_id = -1;
  accel_engine_.update_eng_.initializeInv(psiM);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::registerTWFFastDerivWrapper(const ParticleSet& P, TWFFastDerivWrapper& twf) const
{
  twf.addGroup(P, P.getGroupID(FirstIndex), &phi_);
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::PsiValue DiracDeterminant<PL, VT, FPVT>::ratio(ParticleSet& P, int iat)
{
  UpdateMode             = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0);
  {
    ScopedTimer local_timer(SPOVTimer);
    phi_.evaluateValue(P, iat, psiV);
  }
  {
    ScopedTimer local_timer(RatioTimer);
    // This is an optimization.
    // check invRow_id against WorkingIndex to see if getInvRow() has been called
    // This is intended to save redundant compuation in TM1 and TM3
    if (invRow_id != WorkingIndex)
    {
      invRow_id = WorkingIndex;
      accel_engine_.update_eng_.getInvRow(psiM, WorkingIndex, invRow);
    }
    curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
  }
  return curRatio;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  ScopedTimer local_timer(SPOVTimer);
  Vector<ValueType> inv_row(psiM[VP.refPtcl - FirstIndex], psiV.size());
  phi_.evaluateDetRatios(VP, psiV, inv_row, ratios);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateSpinorRatios(const VirtualParticleSet& VP,
                                                          const std::pair<ValueVector, ValueVector>& spinor_multiplier,
                                                          std::vector<ValueType>& ratios)
{
  ScopedTimer local_timer(SPOVTimer);
  Vector<ValueType> inv_row(psiM[VP.refPtcl - FirstIndex], psiV.size());
  phi_.evaluateDetSpinorRatios(VP, psiV, spinor_multiplier, inv_row, ratios);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateDerivRatios(const VirtualParticleSet& VP,
                                                         const opt_variables_type& optvars,
                                                         std::vector<ValueType>& ratios,
                                                         Matrix<ValueType>& dratios)
{
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  assert(WorkingIndex >= 0);
  std::copy_n(psiM[WorkingIndex], invRow.size(), invRow.data());
  phi_.evaluateDerivRatios(VP, optvars, psiV, invRow, ratios, dratios, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateSpinorDerivRatios(
    const VirtualParticleSet& VP,
    const std::pair<ValueVector, ValueVector>& spinor_multiplier,
    const opt_variables_type& optvars,
    std::vector<ValueType>& ratios,
    Matrix<ValueType>& dratios)
{
  const int WorkingIndex = VP.refPtcl - FirstIndex;
  assert(WorkingIndex >= 0);
  std::copy_n(psiM[WorkingIndex], invRow.size(), invRow.data());
  phi_.evaluateSpinorDerivRatios(VP, spinor_multiplier, optvars, psiV, invRow, ratios, dratios, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  ScopedTimer local_timer(SPOVTimer);
  phi_.evaluateValue(P, -1, psiV);
  Vector<ValueType> ratios_this_det(ratios.data() + FirstIndex, NumPtcls);
  MatrixOperators::product(psiM, psiV, ratios_this_det);
}


template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::resizeScratchObjectsForIonDerivs()
{
  grad_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_lapl_source_psiM.resize(NumPtcls, NumOrbitals);
  grad_grad_source_psiM.resize(NumPtcls, NumOrbitals);
  phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_Minv.resize(NumPtcls, NumOrbitals);
  lapl_phi_Minv.resize(NumPtcls, NumOrbitals);
  grad_phi_alpha_Minv.resize(NumPtcls, NumOrbitals);
}

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::GradType DiracDeterminant<PL, VT, FPVT>::evalGradSource(ParticleSet& P,
                                                                                                 ParticleSet& source,
                                                                                                 int iat)
{
  GradType g(0.0);
  if (phi_.hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    phi_.evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
    g = simd::dot(psiM.data(), grad_source_psiM.data(), psiM.size());
  }

  return g;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateHessian(ParticleSet& P, HessVector& grad_grad_psi)
{
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiM.rows(), psiM.cols());
  //IM A HACK.  Assumes evaluateLog has already been executed.
  phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, grad_grad_source_psiM);
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

template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::GradType DiracDeterminant<PL, VT, FPVT>::evalGradSource(
    ParticleSet& P,
    ParticleSet& source,
    int iat,
    TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad)
{
  GradType gradPsi(0.0);
  if (phi_.hasIonDerivs())
  {
    resizeScratchObjectsForIonDerivs();
    phi_.evaluateGradSource(P, FirstIndex, LastIndex, source, iat, grad_source_psiM, grad_grad_source_psiM,
                            grad_lapl_source_psiM);
    // HACK HACK HACK
    // phi_.evaluateVGL(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
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
template<PlatformKind PL, typename VT, typename FPVT>
typename DiracDeterminant<PL, VT, FPVT>::LogValue DiracDeterminant<PL, VT, FPVT>::evaluateLog(
    const ParticleSet& P,
    ParticleSet::ParticleGradient& G,
    ParticleSet::ParticleLaplacian& L)
{
  recompute(P);

  for (int i = 0, iat = FirstIndex; i < NumPtcls; i++, iat++)
  {
    mGradType rv   = simd::dot(psiM[i], dpsiM[i], NumOrbitals);
    mValueType lap = simd::dot(psiM[i], d2psiM[i], NumOrbitals);
    G[iat] += rv;
    L[iat] += lap - dot(rv, rv);
  }
  return log_value_;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::recompute(const ParticleSet& P)
{
  {
    ScopedTimer local_timer(SPOVGLTimer);
    UpdateMode = ORB_WALKER;
    phi_.evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
  }

  invertPsiM(psiM_temp, psiM);

  // invRow becomes invalid after updating the inverse matrix
  invRow_id = -1;
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateDerivatives(ParticleSet& P,
                                                         const opt_variables_type& active,
                                                         Vector<ValueType>& dlogpsi,
                                                         Vector<ValueType>& dhpsioverpsi)
{
  phi_.evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
void DiracDeterminant<PL, VT, FPVT>::evaluateDerivativesWF(ParticleSet& P,
                                                           const opt_variables_type& active,
                                                           Vector<ValueType>& dlogpsi)
{
  phi_.evaluateDerivativesWF(P, active, dlogpsi, FirstIndex, LastIndex);
}

template<PlatformKind PL, typename VT, typename FPVT>
std::unique_ptr<DiracDeterminantBase> DiracDeterminant<PL, VT, FPVT>::makeCopy(SPOSet& phi) const
{
  return std::make_unique<DiracDeterminant<PL, VT, FPVT>>(phi, FirstIndex, LastIndex, ndelay_, matrix_inverter_kind_);
}

template class DiracDeterminant<>;
#if defined(ENABLE_CUDA)
template class DiracDeterminant<PlatformKind::CUDA, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif
#if defined(ENABLE_SYCL)
template class DiracDeterminant<PlatformKind::SYCL, QMCTraits::ValueType, QMCTraits::QTFull::ValueType>;
#endif

} // namespace qmcplusplus

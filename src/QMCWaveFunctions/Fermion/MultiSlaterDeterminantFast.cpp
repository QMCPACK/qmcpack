//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{
MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl,
                                                       std::vector<std::unique_ptr<MultiDiracDeterminant>>&& dets,
                                                       bool use_pre_computing)
    : WaveFunctionComponent("MultiSlaterDeterminantFast"),
      RatioTimer(*timer_manager.createTimer(ClassName + "::ratio")),
      EvalGradTimer(*timer_manager.createTimer(ClassName + "::evalGrad")),
      RatioGradTimer(*timer_manager.createTimer(ClassName + "::ratioGrad")),
      PrepareGroupTimer(*timer_manager.createTimer(ClassName + "::prepareGroup")),
      UpdateTimer(*timer_manager.createTimer(ClassName + "::updateBuffer")),
      AccRejTimer(*timer_manager.createTimer(ClassName + "::Accept_Reject")),
      EvaluateTimer(*timer_manager.createTimer(ClassName + "::evaluate")),
      CI_Optimizable(false),
      use_pre_computing_(use_pre_computing)
{
  registerTimers();
  //Optimizable=true;
  Optimizable  = false;
  is_fermionic = true;
  usingCSF     = false;
  Dets         = std::move(dets);
  C_otherDs.resize(Dets.size());
  int NP = targetPtcl.getTotalNum();
  myG.resize(NP);
  myL.resize(NP);
  myG_temp.resize(NP);
  myL_temp.resize(NP);

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i) - 1;

  usingBF = false;
  BFTrans = 0;
}

void MultiSlaterDeterminantFast::initialize()
{
  C2node       = std::make_shared<std::vector<std::vector<size_t>>>(Dets.size());
  C            = std::make_shared<std::vector<ValueType>>();
  CSFcoeff     = std::make_shared<std::vector<ValueType>>();
  DetsPerCSF   = std::make_shared<std::vector<size_t>>();
  CSFexpansion = std::make_shared<std::vector<RealType>>();
  myVars       = std::make_shared<opt_variables_type>();
}

MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() = default;

WaveFunctionComponentPtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
{
  std::vector<std::unique_ptr<MultiDiracDeterminant>> dets_clone;
  for (auto& det : Dets)
    dets_clone.emplace_back(std::make_unique<MultiDiracDeterminant>(*det));

  MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp, std::move(dets_clone), use_pre_computing_);
  if (usingBF)
  {
    BackflowTransformation* tr = BFTrans->makeClone(tqp);
    clone->setBF(tr);
  }

  clone->C2node = C2node;
  clone->C      = C;
  clone->myVars = myVars;

  clone->Optimizable = Optimizable;
  clone->usingCSF    = usingCSF;
  clone->usingBF     = usingBF;

  clone->CI_Optimizable = CI_Optimizable;

  if (usingCSF)
  {
    clone->CSFcoeff     = CSFcoeff;
    clone->CSFexpansion = CSFexpansion;
    clone->DetsPerCSF   = DetsPerCSF;
  }
  return clone;
}

void MultiSlaterDeterminantFast::testMSD(ParticleSet& P, int iat)
{
  //     APP_ABORT("Testing disabled for safety");
  app_log() << "Testing MSDFast. \n";
  int n = P.getTotalNum();
  ParticleSet::ParticleGradient_t G(n), G0(n);
  ParticleSet::ParticleLaplacian_t L(n), L0(n);
  ValueType log0;
  GradType G1;
  //     log = msd->evaluate(P,G,L);
  log0 = evaluate(P, G0, L0);
  /*
       app_log() <<"Testing evaluate(P,G,L). \n";
       std::cout << std::endl << std::endl;
       std::cout <<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 << std::endl;

       for(int i=0; i<n; i++) {
         std::cout <<i  <<"\n"
             <<"  x: " <<G[i][0]-G0[i][0] <<"\n"
             <<"  y: " <<G[i][1]-G0[i][1] <<"\n"
             <<"  z: " <<G[i][2]-G0[i][2] <<"\n"
             <<"  d2: " <<L(i)-L0(i) <<"\n"
             << std::endl;
       }
       std::cout << std::endl << std::endl;
       APP_ABORT("end of test 1");
  */
  Walker_t::WFBuffer_t wbuffer;
  wbuffer.clear();
  registerData(P, wbuffer);
  //     log = msd->evaluate(P,G,L);
  log0 = evaluate(P, G0, L0);
  PosType dr;
  dr[0] = 0.1;
  dr[1] = 0.05;
  dr[2] = -0.01;
  P.makeMove(iat, dr);
  app_log() << "Testing ratio(P,dG,dL). \n";
  G       = 0;
  G0      = 0;
  L       = 0;
  log0    = ratioGrad(P, iat, G1);
  G0[iat] = G1;
  std::cout << "Psi: " << log0 << std::endl;
  for (int i = 0; i < n; i++)
  {
    std::cout << i << "\n"
              << "  x: " << G[i][0] - G0[i][0] << "  " << G[i][0] << "\n"
              << "  y: " << G[i][1] - G0[i][1] << "  " << G[i][1] << "\n"
              << "  z: " << G[i][2] - G0[i][2] << "  " << G[i][2] << "\n"
              << std::endl;
  }
  std::cout << std::endl << std::endl;
  APP_ABORT("After MultiSlaterDeterminantFast::testMSD()");
}

/** Compute VGL of this MultiSlaterDeterminantFast
 *
 * THis is introduced to remove redundant code in 
 * - evaluate(P,G,L)
 * - evaluateLog(P,G,L,buf,fillbuffer)
 * Miguel's note: can this change over time??? I don't know yet
 */
WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::evaluate_vgl_impl(
    ParticleSet& P,
    ParticleSet::ParticleGradient_t& g_tmp,
    ParticleSet::ParticleLaplacian_t& l_tmp)
{
  const ValueVector_t& detValues_up = Dets[0]->detValues;
  const ValueVector_t& detValues_dn = Dets[1]->detValues;
  const GradMatrix_t& grads_up      = Dets[0]->grads;
  const GradMatrix_t& grads_dn      = Dets[1]->grads;
  const ValueMatrix_t& lapls_up     = Dets[0]->lapls;
  const ValueMatrix_t& lapls_dn     = Dets[1]->lapls;
  const size_t N1                   = Dets[0]->FirstIndex;
  const size_t N2                   = Dets[1]->FirstIndex;
  const size_t NP1                  = Dets[0]->NumPtcls;
  const size_t NP2                  = Dets[1]->NumPtcls;
  const ValueType czero(0);
  PsiValueType psi = czero;
  g_tmp            = czero;
  l_tmp            = czero;

  const ValueType* restrict cptr = C->data();
  const size_t nc                = C->size();
  const auto& upC                = (*C2node)[0];
  const auto& dnC                = (*C2node)[1];
  for (size_t i = 0; i < nc; ++i)
  {
    const ValueType c = cptr[i];
    const size_t up   = upC[i];
    const size_t down = dnC[i];
    psi += c * detValues_up[up] * detValues_dn[down];
    const ValueType c_up = c * detValues_dn[down];
    const ValueType c_dn = c * detValues_up[up];
    for (int k = 0, n = N1; k < NP1; k++, n++)
    {
      g_tmp[n] += c_up * grads_up(up, k); //myG[n] += c*grads_up(up,k)*detValues_dn[down];
      l_tmp[n] += c_up * lapls_up(up, k); //myL[n] += c*lapls_up(up,k)*detValues_dn[down];
    }
    for (int k = 0, n = N2; k < NP2; k++, n++)
    {
      g_tmp[n] += c_dn * grads_dn(down, k); //myG[n] += c*grads_dn(down,k)*detValues_up[up];
      l_tmp[n] += c_dn * lapls_dn(down, k); //myL[n] += c*lapls_dn(down,k)*detValues_up[up];
    }
  }
  ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psi);
  g_tmp *= psiinv;
  l_tmp *= psiinv;

  return psi;
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::evaluate(ParticleSet& P,
                                                                         ParticleSet::ParticleGradient_t& G,
                                                                         ParticleSet::ParticleLaplacian_t& L)
{
  ScopedTimer local_timer(EvaluateTimer);
  Dets[0]->evaluateForWalkerMove(P);
  Dets[1]->evaluateForWalkerMove(P);

  psiCurrent = evaluate_vgl_impl(P, myG, myL);

  G += myG;
  for (size_t i = 0; i < L.size(); i++)
    L[i] += myL[i] - dot(myG[i], myG[i]);

  return psiCurrent;
}

WaveFunctionComponent::LogValueType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P,
                                                                            ParticleSet::ParticleGradient_t& G,
                                                                            ParticleSet::ParticleLaplacian_t& L)
{
  return LogValue = convertValueToLog(evaluate(P, G, L));
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::evalGrad_impl(ParticleSet& P,
                                                                              int iat,
                                                                              bool newpos,
                                                                              GradType& g_at)
{
  const int det_id = getDetID(iat);

  if (newpos)
    Dets[det_id]->evaluateDetsAndGradsForPtclMove(P, iat);
  else
    Dets[det_id]->evaluateGrads(P, iat);

  const GradMatrix_t& grads            = (newpos) ? Dets[det_id]->new_grads : Dets[det_id]->grads;
  const ValueType* restrict detValues0 = (newpos) ? Dets[det_id]->new_detValues.data() : Dets[det_id]->detValues.data();
  const size_t noffset                 = Dets[det_id]->FirstIndex;

  PsiValueType psi(0);
  for (size_t i = 0; i < Dets[det_id]->NumDets; i++)
  {
    psi += detValues0[i] * C_otherDs[det_id][i];
    g_at += C_otherDs[det_id][i] * grads(i, iat - noffset);
  }
  return psi;
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::evalGrad_impl_no_precompute(ParticleSet& P,
                                                                                            int iat,
                                                                                            bool newpos,
                                                                                            GradType& g_at)
{
  const bool upspin = getDetID(iat) == 0;
  const int spin0   = (upspin) ? 0 : 1;
  const int spin1   = (upspin) ? 1 : 0;

  if (newpos)
    Dets[spin0]->evaluateDetsAndGradsForPtclMove(P, iat);
  else
    Dets[spin0]->evaluateGrads(P, iat);

  const GradMatrix_t& grads            = (newpos) ? Dets[spin0]->new_grads : Dets[spin0]->grads;
  const ValueType* restrict detValues0 = (newpos) ? Dets[spin0]->new_detValues.data() : Dets[spin0]->detValues.data();
  const ValueType* restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t* restrict det0          = (*C2node)[spin0].data();
  const size_t* restrict det1          = (*C2node)[spin1].data();
  const ValueType* restrict cptr       = C->data();
  const size_t nc                      = C->size();
  const size_t noffset                 = Dets[spin0]->FirstIndex;
  PsiValueType psi(0);
  for (size_t i = 0; i < nc; ++i)
  {
    const size_t d0 = det0[i];
    //const size_t d1=det1[i];
    //psi +=  cptr[i]*detValues0[d0]        * detValues1[d1];
    //g_at += cptr[i]*grads(d0,iat-noffset) * detValues1[d1];
    const ValueType t = cptr[i] * detValues1[det1[i]];
    psi += t * detValues0[d0];
    g_at += t * grads(d0, iat - noffset);
  }
  return psi;
}

WaveFunctionComponent::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: evalGrad not implemented. \n");
  }

  ScopedTimer local_timer(EvalGradTimer);

  GradType grad_iat;
  PsiValueType psi;
  if (use_pre_computing_)
    psi = evalGrad_impl(P, iat, false, grad_iat);
  else
    psi = evalGrad_impl_no_precompute(P, iat, false, grad_iat);

  grad_iat *= (PsiValueType(1.0) / psi);
  return grad_iat;
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratioGrad not implemented. \n");
  }

  ScopedTimer local_timer(RatioGradTimer);
  UpdateMode = ORB_PBYP_PARTIAL;

  GradType dummy;
  PsiValueType psiNew;
  if (use_pre_computing_)
    psiNew = evalGrad_impl(P, iat, true, dummy);
  else
    psiNew = evalGrad_impl_no_precompute(P, iat, true, dummy);

  grad_iat += static_cast<ValueType>(PsiValueType(1.0) / psiNew) * dummy;
  curRatio = psiNew / psiCurrent;
  return curRatio;
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::ratio_impl(ParticleSet& P, int iat)
{
  const int det_id = getDetID(iat);

  Dets[det_id]->evaluateDetsForPtclMove(P, iat);

  const ValueType* restrict detValues0 = Dets[det_id]->new_detValues.data(); //always new

  PsiValueType psi = 0;
  // This function computes
  // psi=Det_Coeff[i]*Det_Value[unique_det_up]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // Since only one electron group is moved at the time, identified by det_id, We precompute:
  // C_otherDs[det_id][i]=Det_Coeff[i]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  for (size_t i = 0; i < Dets[det_id]->NumDets; i++)
    psi += detValues0[i] * C_otherDs[det_id][i];

  return psi;
}

WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::ratio_impl_no_precompute(ParticleSet& P, int iat)
{
  const bool upspin = getDetID(iat) == 0;
  const int spin0   = (upspin) ? 0 : 1;
  const int spin1   = (upspin) ? 1 : 0;

  Dets[spin0]->evaluateDetsForPtclMove(P, iat);

  const ValueType* restrict detValues0 = Dets[spin0]->new_detValues.data(); //always new
  const ValueType* restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t* restrict det0          = (*C2node)[spin0].data();
  const size_t* restrict det1          = (*C2node)[spin1].data();
  const ValueType* restrict cptr       = C->data();
  const size_t nc                      = C->size();

  PsiValueType psi = 0;
  for (size_t i = 0; i < nc; ++i)
    psi += cptr[i] * detValues0[det0[i]] * detValues1[det1[i]];
  return psi;
}

// use ci_node for this routine only
WaveFunctionComponent::PsiValueType MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio not implemented. \n");
  }

  ScopedTimer local_timer(RatioTimer);
  UpdateMode = ORB_PBYP_RATIO;

  PsiValueType psiNew;
  if (use_pre_computing_)
    psiNew = ratio_impl(P, iat);
  else
    psiNew = ratio_impl_no_precompute(P, iat);

  curRatio = psiNew / psiCurrent;
  return curRatio;
}

void MultiSlaterDeterminantFast::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  // this should depend on the type of update, ratio / ratioGrad
  // for now is incorrect fot ratio(P,iat,dG,dL) updates
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: acceptMove not implemented. \n");
  }

  ScopedTimer local_timer(AccRejTimer);
  // update psiCurrent,myG_temp,myL_temp
  psiCurrent *= curRatio;
  curRatio = 1.0;

  Dets[getDetID(iat)]->acceptMove(P, iat, safe_to_delay);
}

void MultiSlaterDeterminantFast::restore(int iat)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }

  ScopedTimer local_timer(AccRejTimer);
  Dets[getDetID(iat)]->restore(iat);
  curRatio = 1.0;
}

void MultiSlaterDeterminantFast::registerData(ParticleSet& P, WFBufferType& buf)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }

  Dets[0]->registerData(P, buf);
  Dets[1]->registerData(P, buf);

  buf.add(psiCurrent);
}

// FIX FIX FIX
WaveFunctionComponent::LogValueType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P,
                                                                             WFBufferType& buf,
                                                                             bool fromscratch)
{
  ScopedTimer local_timer(UpdateTimer);

  Dets[0]->updateBuffer(P, buf, fromscratch);
  Dets[1]->updateBuffer(P, buf, fromscratch);

  psiCurrent = evaluate_vgl_impl(P, myG, myL);

  P.G += myG;
  for (int i = 0; i < P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i], myG[i]);

  buf.put(psiCurrent);

  return LogValue = convertValueToLog(psiCurrent);
}

void MultiSlaterDeterminantFast::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  if (usingBF)
  {
    APP_ABORT("Fast MSD+BF: copyFromBuffer not implemented. \n");
  }
  Dets[0]->copyFromBuffer(P, buf);
  Dets[1]->copyFromBuffer(P, buf);

  buf.get(psiCurrent);
}


void MultiSlaterDeterminantFast::checkInVariables(opt_variables_type& active)
{
  if (CI_Optimizable)
  {
    if (myVars->size())
      active.insertFrom(*myVars);
    else
      Optimizable = false;
  }
  if (Dets[0]->Optimizable && Dets[1]->Optimizable)
  {
    Dets[0]->checkInVariables(active);
    Dets[1]->checkInVariables(active);
  }
}

void MultiSlaterDeterminantFast::checkOutVariables(const opt_variables_type& active)
{
  if (CI_Optimizable)
    myVars->getIndex(active);

  if (Dets[0]->Optimizable && Dets[1]->Optimizable)
  {
    Dets[0]->checkOutVariables(active);
    Dets[1]->checkOutVariables(active);
  }
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
{
  if (CI_Optimizable)
  {
    if (usingCSF)
    {
      ValueType* restrict CSFcoeff_p = CSFcoeff->data();
      for (int i = 0; i < CSFcoeff->size() - 1; i++)
      {
        int loc = myVars->where(i);
        if (loc >= 0)
        {
          CSFcoeff_p[i + 1] = (*myVars)[i] = active[loc];
        }
      }
      int cnt                                 = 0;
      ValueType* restrict C_p                 = C->data();
      const RealType* restrict CSFexpansion_p = CSFexpansion->data();
      for (int i = 0; i < DetsPerCSF->size(); i++)
      {
        for (int k = 0; k < (*DetsPerCSF)[i]; k++)
        {
          C_p[cnt] = CSFcoeff_p[i] * CSFexpansion_p[cnt];
          cnt++;
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
    else
    {
      ValueType* restrict C_p = C->data();
      for (int i = 0; i < C->size() - 1; i++)
      {
        int loc = myVars->where(i);
        if (loc >= 0)
        {
          C_p[i + 1] = (*myVars)[i] = active[loc];
        }
      }
      //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
    }
  }
  if (Dets[0]->Optimizable && Dets[1]->Optimizable)
  {
    Dets[0]->resetParameters(active);
    Dets[1]->resetParameters(active);
  }
}
void MultiSlaterDeterminantFast::reportStatus(std::ostream& os) {}


void MultiSlaterDeterminantFast::evaluateDerivatives(ParticleSet& P,
                                                     const opt_variables_type& optvars,
                                                     std::vector<ValueType>& dlogpsi,
                                                     std::vector<ValueType>& dhpsioverpsi)
{
  evaluateDerivativesWF(P, optvars, dlogpsi);
  if (CI_Optimizable)
  {
    bool recalculate(false);
    for (int k = 0; k < myVars->size(); ++k)
    {
      int kk = myVars->where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
    }
    // need to modify for CSF later on, right now assume Slater Det basis
    if (recalculate)
    {
      if (usingCSF)
      {
        if (laplSum_up.size() == 0)
          laplSum_up.resize(Dets[0]->detValues.size());
        if (laplSum_dn.size() == 0)
          laplSum_dn.resize(Dets[1]->detValues.size());
        // assume that evaluateLog has been called in opt routine before
        //   Dets[0]->evaluateForWalkerMove(P);
        //   Dets[1]->evaluateForWalkerMove(P);
        ValueVector_t& detValues_up = Dets[0]->detValues;
        ValueVector_t& detValues_dn = Dets[1]->detValues;
        GradMatrix_t& grads_up      = Dets[0]->grads;
        GradMatrix_t& grads_dn      = Dets[1]->grads;
        ValueMatrix_t& lapls_up     = Dets[0]->lapls;
        ValueMatrix_t& lapls_dn     = Dets[1]->lapls;
        const size_t N1             = Dets[0]->FirstIndex;
        const size_t N2             = Dets[1]->FirstIndex;
        const size_t NP1            = Dets[0]->NumPtcls;
        const size_t NP2            = Dets[1]->NumPtcls;
        // myG,myL should already be calculated
        const size_t n = P.getTotalNum();

        ValueType psiinv   = static_cast<ValueType>(PsiValueType(1.0) / psiCurrent);
        ValueType lapl_sum = 0.0;
        ValueType gg = 0.0, ggP = 0.0;
        myG_temp = 0.0;
        int num  = laplSum_up.size();
        ValueVector_t::iterator it(laplSum_up.begin());
        ValueVector_t::iterator last(laplSum_up.end());
        ValueType* ptr0 = lapls_up[0];
        int nels_up     = P.last(0) - P.first(0);
        while (it != last)
        {
          (*it) = 0.0;
          for (int k = 0; k < nels_up; k++, ptr0++)
            (*it) += *ptr0;
          it++;
        }
        it          = laplSum_dn.begin();
        last        = laplSum_dn.end();
        ptr0        = lapls_dn[0];
        int nels_dn = P.last(1) - P.first(1);
        while (it != last)
        {
          (*it) = 0.0;
          for (int k = 0; k < nels_dn; k++, ptr0++)
            (*it) += *ptr0;
          it++;
        }

        const ValueType* restrict C_p = C->data();
        for (size_t i = 0; i < C->size(); i++)
        {
          size_t upC     = (*C2node)[0][i];
          size_t dnC     = (*C2node)[1][i];
          ValueType tmp1 = C_p[i] * detValues_dn[dnC] * psiinv;
          ValueType tmp2 = C_p[i] * detValues_up[upC] * psiinv;
          lapl_sum += tmp1 * laplSum_up[upC] + tmp2 * laplSum_dn[dnC];
          for (size_t k = 0, j = N1; k < NP1; k++, j++)
            myG_temp[j] += tmp1 * grads_up(upC, k);
          for (size_t k = 0, j = N2; k < NP2; k++, j++)
            myG_temp[j] += tmp2 * grads_dn(dnC, k);
        }
        gg = ggP = 0.0;
        for (size_t i = 0; i < n; i++)
        {
          gg += dot(myG_temp[i], myG_temp[i]) - dot(P.G[i], myG_temp[i]);
        }
        //       for(int i=0; i<C.size(); i++)
        num     = CSFcoeff->size() - 1;
        int cnt = 0;
        //        this one is not optable
        cnt += (*DetsPerCSF)[0];
        int ip(1);
        for (int i = 0; i < num; i++, ip++)
        {
          int kk = myVars->where(i);
          if (kk < 0)
          {
            cnt += (*DetsPerCSF)[ip];
            continue;
          }
          ValueType q0 = 0.0, v1 = 0.0, v2 = 0.0;
          const RealType* restrict CSFexpansion_p = CSFexpansion->data();
          for (int k = 0; k < (*DetsPerCSF)[ip]; k++)
          {
            size_t upC     = (*C2node)[0][cnt];
            size_t dnC     = (*C2node)[1][cnt];
            ValueType tmp1 = CSFexpansion_p[cnt] * detValues_dn[dnC] * psiinv;
            ValueType tmp2 = CSFexpansion_p[cnt] * detValues_up[upC] * psiinv;
            q0 += (tmp1 * laplSum_up[upC] + tmp2 * laplSum_dn[dnC]);
            for (size_t l = 0, j = N1; l < NP1; l++, j++)
              v1 += tmp1 * static_cast<ValueType>(dot(P.G[j], grads_up(upC, l)) - dot(myG_temp[j], grads_up(upC, l)));
            for (size_t l = 0, j = N2; l < NP2; l++, j++)
              v2 += tmp2 * static_cast<ValueType>(dot(P.G[j], grads_dn(dnC, l)) - dot(myG_temp[j], grads_dn(dnC, l)));
            cnt++;
          }
          ValueType dhpsi  = (RealType)-0.5 * (q0 - dlogpsi[kk] * lapl_sum) - dlogpsi[kk] * gg - v1 - v2;
          dhpsioverpsi[kk] = dhpsi;
        }
      }
      else
      //usingDETS
      {
        if (laplSum_up.size() == 0)
          laplSum_up.resize(Dets[0]->detValues.size());
        if (laplSum_dn.size() == 0)
          laplSum_dn.resize(Dets[1]->detValues.size());
        // assume that evaluateLog has been called in opt routine before
        //   Dets[0]->evaluateForWalkerMove(P);
        //   Dets[1]->evaluateForWalkerMove(P);
        ValueVector_t& detValues_up = Dets[0]->detValues;
        ValueVector_t& detValues_dn = Dets[1]->detValues;
        GradMatrix_t& grads_up      = Dets[0]->grads;
        GradMatrix_t& grads_dn      = Dets[1]->grads;
        ValueMatrix_t& lapls_up     = Dets[0]->lapls;
        ValueMatrix_t& lapls_dn     = Dets[1]->lapls;
        int N1                      = Dets[0]->FirstIndex;
        int N2                      = Dets[1]->FirstIndex;
        int NP1                     = Dets[0]->NumPtcls;
        int NP2                     = Dets[1]->NumPtcls;
        int n                       = P.getTotalNum();
        ValueType psiinv            = static_cast<ValueType>(PsiValueType(1.0) / psiCurrent);
        ValueType lapl_sum          = 0.0;
        ValueType gg = 0.0, ggP = 0.0;
        myG_temp = 0.0;
        ValueVector_t::iterator it(laplSum_up.begin());
        ValueVector_t::iterator last(laplSum_up.end());
        ValueType* ptr0 = lapls_up[0];
        int nels_up     = P.last(0) - P.first(0);
        while (it != last)
        {
          (*it) = 0.0;
          for (int k = 0; k < nels_up; k++, ptr0++)
            (*it) += *ptr0;
          it++;
        }
        it          = laplSum_dn.begin();
        last        = laplSum_dn.end();
        ptr0        = lapls_dn[0];
        int nels_dn = P.last(1) - P.first(1);
        while (it != last)
        {
          (*it) = 0.0;
          for (size_t k = 0; k < nels_dn; k++, ptr0++)
            (*it) += *ptr0;
          it++;
        }
        const ValueType* restrict C_p = C->data();
        for (size_t i = 0; i < C->size(); i++)
        {
          size_t upC     = (*C2node)[0][i];
          size_t dnC     = (*C2node)[1][i];
          ValueType tmp1 = C_p[i] * detValues_dn[dnC] * psiinv;
          ValueType tmp2 = C_p[i] * detValues_up[upC] * psiinv;
          lapl_sum += tmp1 * laplSum_up[upC] + tmp2 * laplSum_dn[dnC];
          for (size_t k = 0, j = N1; k < NP1; k++, j++)
            myG_temp[j] += tmp1 * grads_up(upC, k);
          for (size_t k = 0, j = N2; k < NP2; k++, j++)
            myG_temp[j] += tmp2 * grads_dn(dnC, k);
        }
        gg = ggP = 0.0;
        for (size_t i = 0; i < n; i++)
        {
          gg += dot(myG_temp[i], myG_temp[i]) - dot(P.G[i], myG_temp[i]);
        }
        for (size_t i = 1; i < C->size(); i++)
        {
          int kk = myVars->where(i - 1);
          if (kk < 0)
            continue;
          const size_t upC = (*C2node)[0][i];
          const size_t dnC = (*C2node)[1][i];
          ValueType tmp1   = detValues_dn[dnC] * psiinv;
          ValueType tmp2   = detValues_up[upC] * psiinv;
          ValueType v1 = 0.0, v2 = 0.0;
          for (size_t k = 0, j = N1; k < NP1; k++, j++)
            v1 += (dot(P.G[j], grads_up(upC, k)) - dot(myG_temp[j], grads_up(upC, k)));
          for (size_t k = 0, j = N2; k < NP2; k++, j++)
            v2 += (dot(P.G[j], grads_dn(dnC, k)) - dot(myG_temp[j], grads_dn(dnC, k)));
          ValueType dhpsi =
              (RealType)-0.5 * (tmp1 * laplSum_up[upC] + tmp2 * laplSum_dn[dnC] - dlogpsi[kk] * lapl_sum) -
              dlogpsi[kk] * gg - (tmp1 * v1 + tmp2 * v2);
          dhpsioverpsi[kk] = dhpsi;
        }
      }
    }
  }

  Dets[0]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, *Dets[1], static_cast<ValueType>(psiCurrent), *C,
                               (*C2node)[0], (*C2node)[1]);
  Dets[1]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, *Dets[0], static_cast<ValueType>(psiCurrent), *C,
                               (*C2node)[1], (*C2node)[0]);
}

void MultiSlaterDeterminantFast::evaluateDerivativesWF(ParticleSet& P,
                                                       const opt_variables_type& optvars,
                                                       std::vector<ValueType>& dlogpsi)
{
  if (CI_Optimizable)
  {
    bool recalculate(false);
    for (int k = 0; k < myVars->size(); ++k)
    {
      int kk = myVars->where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
    }
    // need to modify for CSF later on, right now assume Slater Det basis
    if (recalculate)
    {
      if (usingCSF)
      {
        // assume that evaluateLog has been called in opt routine before
        ValueVector_t& detValues_up = Dets[0]->detValues;
        ValueVector_t& detValues_dn = Dets[1]->detValues;

        ValueType psiinv = static_cast<ValueType>(PsiValueType(1.0) / psiCurrent);

        int num = CSFcoeff->size() - 1;
        int cnt = 0;
        //        this one is not optable
        cnt += (*DetsPerCSF)[0];
        int ip(1);
        for (int i = 0; i < num; i++, ip++)
        {
          int kk = myVars->where(i);
          if (kk < 0)
          {
            cnt += (*DetsPerCSF)[ip];
            continue;
          }
          ValueType cdet                          = 0.0;
          const RealType* restrict CSFexpansion_p = CSFexpansion->data();
          for (int k = 0; k < (*DetsPerCSF)[ip]; k++)
          {
            size_t upC = (*C2node)[0][cnt];
            size_t dnC = (*C2node)[1][cnt];
            cdet += CSFexpansion_p[cnt] * detValues_up[upC] * detValues_dn[dnC] * psiinv;
            cnt++;
          }
          dlogpsi[kk] = cdet;
        }
      }
      else
      //usingDETS
      {
        ValueVector_t& detValues_up = Dets[0]->detValues;
        ValueVector_t& detValues_dn = Dets[1]->detValues;
        ValueType psiinv            = static_cast<ValueType>(PsiValueType(1.0) / psiCurrent);
        for (size_t i = 1; i < C->size(); i++)
        {
          int kk = myVars->where(i - 1);
          if (kk < 0)
            continue;
          const size_t upC = (*C2node)[0][i];
          const size_t dnC = (*C2node)[1][i];
          ValueType cdet   = detValues_up[upC] * detValues_dn[dnC] * psiinv;
          dlogpsi[kk]      = cdet;
        }
      }
    }
  }

  // FIXME this needs to be fixed by SPF to separate evaluateDerivatives and evaluateDerivativesWF for otbital rotation matrix
  Dets[0]->evaluateDerivativesWF(P, optvars, dlogpsi, *Dets[1], psiCurrent, *C, (*C2node)[0], (*C2node)[1]);
  Dets[1]->evaluateDerivativesWF(P, optvars, dlogpsi, *Dets[0], psiCurrent, *C, (*C2node)[1], (*C2node)[0]);
}

void MultiSlaterDeterminantFast::buildOptVariables()
{
  Dets[0]->buildOptVariables((*C2node)[0]);
  Dets[1]->buildOptVariables((*C2node)[1]);
}

void MultiSlaterDeterminantFast::registerTimers()
{
  RatioTimer.reset();
  EvalGradTimer.reset();
  RatioGradTimer.reset();
  PrepareGroupTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
}

void MultiSlaterDeterminantFast::prepareGroup(ParticleSet& P, int ig)
{
  if (!use_pre_computing_)
    return;
  // This function computes
  // C_otherDs[det_id][i]=Det_Coeff[i]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // Since only one electron group is moved at the time, identified by det_id, We precompute C_otherDs[det_id][i]:
  // psi=Det_Coeff[i]*Det_Value[unique_det_up]*Det_Value[unique_det_dn]*Det_Value[unique_det_AnyOtherType]
  // becomes:
  // psi=Det_Value[unique_det_up]*C_otherDs[det_id][i]
  // ig is the id of the group electron being moved. In this function, we compute the other groups
  // of electrons.
  // We loop over the number of type of determinants (up, diwn, positrons, etc), but we only multiply for ll types BUT ig
  // C_otherDs(0, :) stores C x D_dn x D_pos
  // C_otherDs(1, :) stores C x D_up x D_pos
  // C_otherDs(2, :) stores C x D_up x D_dn

  ScopedTimer local_timer(PrepareGroupTimer);
  C_otherDs[ig].resize(Dets[ig]->NumDets);
  std::fill(C_otherDs[ig].begin(), C_otherDs[ig].end(), ValueType(0));
  for (size_t i = 0; i < C->size(); i++)
  {
    auto product = (*C)[i];
    for (size_t id = 0; id < Dets.size(); id++)
      if (id != ig)
        product *= Dets[id]->detValues[(*C2node)[id][i]];
    C_otherDs[ig][(*C2node)[ig][i]] += product;
  }
}


} // namespace qmcplusplus

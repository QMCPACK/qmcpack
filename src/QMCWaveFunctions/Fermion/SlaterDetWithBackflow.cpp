//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SlaterDetWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
SlaterDetWithBackflow::SlaterDetWithBackflow(ParticleSet& targetPtcl, BackflowTransformation* BF)
    : SlaterDet(targetPtcl, "SlaterDetWithBackflow"), BFTrans(BF)
{
  Optimizable = false;
}

///destructor
SlaterDetWithBackflow::~SlaterDetWithBackflow()
{
  ///clean up SPOSet
}

void SlaterDetWithBackflow::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->evaluateRatiosAlltoOne(P, ratios);
}

SlaterDetWithBackflow::LogValueType SlaterDetWithBackflow::evaluateLog(ParticleSet& P,
                                                                       ParticleSet::ParticleGradient_t& G,
                                                                       ParticleSet::ParticleLaplacian_t& L)
{
  BFTrans->evaluate(P);
  LogValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->evaluateLog(P, G, L);
  return LogValue;
}

void SlaterDetWithBackflow::registerData(ParticleSet& P, WFBufferType& buf)
{
  BFTrans->registerData(P, buf);
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
}

SlaterDetWithBackflow::LogValueType SlaterDetWithBackflow::updateBuffer(ParticleSet& P,
                                                                        WFBufferType& buf,
                                                                        bool fromscratch)
{
  //BFTrans->updateBuffer(P,buf,fromscratch);
  BFTrans->updateBuffer(P, buf, fromscratch);
  //BFTrans->evaluate(P);
  LogValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->updateBuffer(P, buf, fromscratch);
  return LogValue;
}

void SlaterDetWithBackflow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BFTrans->copyFromBuffer(P, buf);
  //BFTrans->evaluate(P);
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
}

WaveFunctionComponentPtr SlaterDetWithBackflow::makeClone(ParticleSet& tqp) const
{
  BackflowTransformation* tr     = BFTrans->makeClone(tqp);
  SlaterDetWithBackflow* myclone = new SlaterDetWithBackflow(tqp, tr);
  myclone->Optimizable           = Optimizable;
  for (int i = 0; i < Dets.size(); ++i)
  {
    DiracDeterminantBase* dclne = Dets[i]->makeCopy(std::unique_ptr<SPOSet>(Dets[i]->getPhi()->makeClone()));
    myclone->add(dclne, i);
  }
  myclone->setBF(tr);
  return myclone;
}

void SlaterDetWithBackflow::testDerivGL(ParticleSet& P)
{
  // testing derivatives of G and L
  app_log() << "testing derivatives of G and L \n";
  opt_variables_type wfVars, wfvar_prime;
  checkInVariables(wfVars);
  checkOutVariables(wfVars);
  int Nvars   = wfVars.size();
  wfvar_prime = wfVars;
  wfVars.print(std::cout);
  std::vector<RealType> dlogpsi;
  std::vector<RealType> dhpsi;
  dlogpsi.resize(Nvars);
  dhpsi.resize(Nvars);
  ParticleSet::ParticleGradient_t G0, G1, G2;
  ParticleSet::ParticleLaplacian_t L0, L1, L2;
  G0.resize(P.getTotalNum());
  G1.resize(P.getTotalNum());
  G2.resize(P.getTotalNum());
  L0.resize(P.getTotalNum());
  L1.resize(P.getTotalNum());
  L2.resize(P.getTotalNum());
  LogValueType psi1 = 1.0;
  LogValueType psi2 = 1.0;
  RealType dh       = 0.00001;
  for (int k = 0; k < Dets.size(); k++)
  {
    DiracDeterminantWithBackflow* Dets_ = dynamic_cast<DiracDeterminantWithBackflow*>(Dets[k].get());
    Dets_->testGGG(P);
    for (int i = 0; i < Nvars; i++)
    {
      Dets_->testDerivFjj(P, i);
      Dets_->testDerivLi(P, i);
    }
  }
  app_log() << "Nvars: " << Nvars << std::endl;
  for (int i = 0; i < Nvars; i++)
  {
    for (int j = 0; j < Nvars; j++)
      wfvar_prime[j] = wfVars[j];
    resetParameters(wfvar_prime);
    BFTrans->evaluateDerivatives(P);
    G0 = 0.0;
    G1 = 0.0;
    G2 = 0.0;
    L0 = 0.0;
    L1 = 0.0;
    L2 = 0.0;
    for (int k = 0; k < Dets.size(); k++)
    {
      DiracDeterminantWithBackflow* Dets_ = dynamic_cast<DiracDeterminantWithBackflow*>(Dets[k].get());
      Dets_->evaluateDerivatives(P, wfVars, dlogpsi, dhpsi, &G0, &L0, i);
    }
    for (int j = 0; j < Nvars; j++)
      wfvar_prime[j] = wfVars[j];
    wfvar_prime[i] = wfVars[i] + dh;
    resetParameters(wfvar_prime);
    BFTrans->evaluate(P);
    for (int k = 0; k < Dets.size(); k++)
      psi1 += Dets[k]->evaluateLog(P, G1, L1);
    for (int j = 0; j < Nvars; j++)
      wfvar_prime[j] = wfVars[j];
    wfvar_prime[i] = wfVars[i] - dh;
    resetParameters(wfvar_prime);
    BFTrans->evaluate(P);
    for (int k = 0; k < Dets.size(); k++)
      psi2 += Dets[k]->evaluateLog(P, G2, L2);
    ParticleSet::SingleParticleValue_t tmp = 0.0;
    for (int q = 0; q < P.getTotalNum(); q++)
      tmp += (L1[q] - L2[q]) / (2.0 * dh);
    app_log() << i << "\n"
              << "Ldiff : " << L0[0] << "  " << tmp << "  " << L0[0] - tmp << std::endl;
    for (int k = 0; k < P.getTotalNum(); k++)
    {
      app_log() << G0[k] << std::endl
                << (G1[k] - G2[k]) / (2.0 * dh) << std::endl
                << "Gdiff: " << G0[k] - (G1[k] - G2[k]) / (2.0 * dh) << std::endl
                << std::endl;
    }
  }
  resetParameters(wfVars);
  APP_ABORT("Testing bF derivs \n");
}


void SlaterDetWithBackflow::evaluateDerivatives(ParticleSet& P,
                                                const opt_variables_type& optvars,
                                                std::vector<ValueType>& dlogpsi,
                                                std::vector<ValueType>& dhpsioverpsi)
{
  //testDerivGL(P);
  if (BFTrans->isOptimizable())
  {
    // build QP,Amat,Bmat_full,Xmat,Cmat,Ymat
    BFTrans->evaluateDerivatives(P);
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  }
}


} // namespace qmcplusplus

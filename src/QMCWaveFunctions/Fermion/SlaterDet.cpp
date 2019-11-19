//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
SlaterDet::SlaterDet(ParticleSet& targetPtcl)
{
  Optimizable = false;
  is_fermionic = true;
  ClassName   = "SlaterDet";

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i) - 1;

  Dets.resize(targetPtcl.groups(), nullptr);
}

///destructor
SlaterDet::~SlaterDet()
{
  ///clean up SPOSet
}

///add a new SPOSet to the list of determinants
void SlaterDet::add(SPOSet* sposet, const std::string& aname)
{
  if (mySPOSet.find(aname) == mySPOSet.end())
  {
    mySPOSet[aname]    = sposet;
    sposet->objectName = aname;
  }
  else
  {
    APP_ABORT(" SlaterDet::add(SPOSet*, const std::string&) Cannot reuse the " + aname);
  }
}

///add a new DiracDeterminant to the list of determinants
void SlaterDet::add(Determinant_t* det, int ispin)
{
  if (Dets[ispin] != nullptr)
  {
    APP_ABORT("SlaterDet::add(Determinant_t* det, int ispin) is alreaded instantiated.");
  }
  else
    Dets[ispin] = det;
  Optimizable = Optimizable || det->Optimizable;
}

void SlaterDet::checkInVariables(opt_variables_type& active)
{
  myVars.clear();
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
    {
      Dets[i]->checkInVariables(active);
      Dets[i]->checkInVariables(myVars);
    }
}

void SlaterDet::checkOutVariables(const opt_variables_type& active)
{
  myVars.clear();
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
    {
      Dets[i]->checkOutVariables(active);
      myVars.insertFrom(Dets[i]->myVars);
    }
  myVars.getIndex(active);
}

///reset all the Dirac determinants, Optimizable is true
void SlaterDet::resetParameters(const opt_variables_type& active)
{
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->resetParameters(active);
}

void SlaterDet::reportStatus(std::ostream& os) {}

void SlaterDet::resetTargetParticleSet(ParticleSet& P)
{
  std::map<std::string, SPOSetPtr>::iterator sit(mySPOSet.begin());
  while (sit != mySPOSet.end())
  {
    (*sit).second->resetTargetParticleSet(P);
    ++sit;
  }
  //BasisSet->resetTargetParticleSet(P);
  //LOGMSG("\nSlaterDet::resetTargetParticleSet")
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->resetTargetParticleSet(P);
}

void SlaterDet::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->evaluateRatiosAlltoOne(P, ratios);
}

SlaterDet::LogValueType SlaterDet::evaluateLog(ParticleSet& P,
                                           ParticleSet::ParticleGradient_t& G,
                                           ParticleSet::ParticleLaplacian_t& L)
{
  LogValue   = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->evaluateLog(P, G, L);
  return LogValue;
}

void SlaterDet::mw_evaluateLog(const std::vector<WaveFunctionComponent*>& WFC_list,
                               const std::vector<ParticleSet*>& P_list,
                               const std::vector<ParticleSet::ParticleGradient_t*>& G_list,
                               const std::vector<ParticleSet::ParticleLaplacian_t*>& L_list)
{
  constexpr RealType czero(0);

  for (int iw = 0; iw < WFC_list.size(); iw++)
    WFC_list[iw]->LogValue   = czero;

  for (int i = 0; i < Dets.size(); ++i)
  {
    const std::vector<WaveFunctionComponent*> Det_list(extract_Det_list(WFC_list, i));
    Dets[i]->mw_evaluateLog(Det_list, P_list, G_list, L_list);
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw]->LogValue += Det_list[iw]->LogValue;
  }
}

void SlaterDet::recompute(ParticleSet& P)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->recompute(P);
}

void SlaterDet::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  grad_grad_psi.resize(P.getTotalNum());
  HessVector_t tmp;
  tmp.resize(P.getTotalNum());
  for (int i = 0; i < Dets.size(); ++i)
  {
    tmp = 0;
    Dets[i]->evaluateHessian(P, tmp);
    //  app_log()<<"squee ----- "<<i<< std::endl;
    //  app_log()<<"grad_grad_psi = "<<grad_grad_psi<< std::endl;
    //  app_log()<<"tmp = "<<tmp<< std::endl<< std::endl;
    grad_grad_psi += tmp;
  }
}

void SlaterDet::registerData(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
}

SlaterDet::LogValueType SlaterDet::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  LogValue   = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->updateBuffer(P, buf, fromscratch);
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  return LogValue;
}

void SlaterDet::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
}

WaveFunctionComponentPtr SlaterDet::makeClone(ParticleSet& tqp) const
{
  SlaterDet* myclone   = new SlaterDet(tqp);
  myclone->Optimizable = Optimizable;
  if (mySPOSet.size() > 1)
  {
    std::map<std::string, SPOSetPtr>::const_iterator Mit, Lit;
    Mit = mySPOSet.begin();
    Lit = mySPOSet.end();
    while (Mit != Lit)
    {
      SPOSetPtr spo = (*Mit).second;
      SPOSetPtr spo_clone;
      spo_clone = spo->makeClone();
      spo_clone->resetTargetParticleSet(tqp);
      myclone->add(spo_clone, spo->objectName);
      for (int i = 0; i < Dets.size(); ++i)
      {
        if (spo == Dets[i]->getPhi())
        {
          Determinant_t* newD = Dets[i]->makeCopy(spo_clone);
          newD->resetTargetParticleSet(tqp);
          myclone->add(newD, i);
        }
      }
      Mit++;
    }
  }
  else
  {
    SPOSetPtr spo       = Dets[0]->getPhi();
    SPOSetPtr spo_clone = spo->makeClone();
    spo_clone->resetTargetParticleSet(tqp);
    myclone->add(spo_clone, spo->objectName);
    for (int i = 0; i < Dets.size(); ++i)
    {
      Determinant_t* newD = Dets[i]->makeCopy(spo_clone);
      newD->resetTargetParticleSet(tqp);
      myclone->add(newD, i);
    }
  }
  return myclone;
}

} // namespace qmcplusplus

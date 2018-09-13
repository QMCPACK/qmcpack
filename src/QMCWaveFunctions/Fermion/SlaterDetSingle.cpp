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
    
    
#include "QMCWaveFunctions/Fermion/SlaterDetSingle.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

SlaterDetSingle::SlaterDetSingle(ParticleSet& targetPtcl)
{
  Optimizable = false;
  OrbitalName = "SlaterDetSingle";

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i)-1;

  Dets.resize(targetPtcl.groups(), nullptr);
}

///destructor
SlaterDetSingle::~SlaterDetSingle()
{
  ///clean up SPOSet
}

///add a new SPOSet to the list of determinants
void SlaterDetSingle::add(SPOSet* sposet, const std::string& aname)
{
  if (mySPOSet.find(aname) == mySPOSet.end())
  {
    mySPOSet[aname] = sposet;
    sposet->objectName = aname;
  }
  else
  {
    APP_ABORT(" SlaterDetSingle::add(SPOSet*, const std::string&) Cannot reuse the " + aname )
    ;
  }
}

///add a new DiracDeterminant to the list of determinants
void SlaterDetSingle::add(Determinant_t* det, int ispin)
{
  if (Dets[ispin]!=nullptr)
  {
    APP_ABORT("SlaterDetSingle::add(Determinant_t* det, int ispin) is alreaded instantiated.");
  }
  else
    Dets[ispin] = det;
  Optimizable = Optimizable || det->Optimizable;
}

void SlaterDetSingle::checkInVariables(opt_variables_type& active)
{
  myVars.clear();
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
    {
      Dets[i]->checkInVariables(active);
      Dets[i]->checkInVariables(myVars);
    }
}

void SlaterDetSingle::checkOutVariables(const opt_variables_type& active)
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
void SlaterDetSingle::resetParameters(const opt_variables_type& active)
{
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->resetParameters(active);
}

void SlaterDetSingle::reportStatus(std::ostream& os)
{
}

void SlaterDetSingle::resetTargetParticleSet(ParticleSet& P)
{
  std::map<std::string, SPOSet*>::iterator sit(mySPOSet.begin());
  while (sit != mySPOSet.end())
  {
    (*sit).second->resetTargetParticleSet(P);
    ++sit;
  }
  //BasisSet->resetTargetParticleSet(P);
  //LOGMSG("\nSlaterDetSingle::resetTargetParticleSet")
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->resetTargetParticleSet(P);
}

void SlaterDetSingle::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->evaluateRatiosAlltoOne(P, ratios);
}

SlaterDetSingle::RealType SlaterDetSingle::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  //ValueType psi = 1.0;
  //for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
  //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  LogValue = 0.0;
  PhaseValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
  {
    LogValue += Dets[i]->evaluateLog(P, G, L);
    PhaseValue += Dets[i]->PhaseValue;
  }
  return LogValue;
}

void SlaterDetSingle::recompute(ParticleSet& P)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->recompute(P);
}

void SlaterDetSingle::evaluateHessian(ParticleSet & P, HessVector_t& grad_grad_psi)
{
	grad_grad_psi.resize(P.getTotalNum());
	HessVector_t tmp;
	tmp.resize(P.getTotalNum());
	for (int i = 0; i < Dets.size(); ++i)
    {
	  tmp=0;
      Dets[i]->evaluateHessian(P, tmp);
    //  app_log()<<"squee ----- "<<i<< std::endl;
    //  app_log()<<"grad_grad_psi = "<<grad_grad_psi<< std::endl;
    //  app_log()<<"tmp = "<<tmp<< std::endl<< std::endl;
      grad_grad_psi += tmp;

    }
	
}

void SlaterDetSingle::registerData(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDetSingle::registerData ",buf.current());
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
  DEBUG_PSIBUFFER(" SlaterDetSingle::registerData ",buf.current());
}

void SlaterDetSingle::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
{
  for (size_t i = 0, n=Dets.size(); i < n; ++i)
  {
    Dets[i]->updateAfterSweep(P,G,L);
  }
}

SlaterDetSingle::RealType SlaterDetSingle::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  DEBUG_PSIBUFFER(" SlaterDetSingle::updateBuffer ",buf.current());
  //ValueType psi = 1.0;
  //for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->updateBuffer(P,buf,fromscratch);
  //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  LogValue = 0.0;
  PhaseValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
  {
    LogValue += Dets[i]->updateBuffer(P, buf, fromscratch);
    PhaseValue += Dets[i]->PhaseValue;
  }
  DEBUG_PSIBUFFER(" SlaterDetSingle::updateBuffer ",buf.current());
  return LogValue;
}

void SlaterDetSingle::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDetSingle::copyFromBuffer ",buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDetSingle::copyFromBuffer ",buf.current());
}

WaveFunctionComponentPtr SlaterDetSingle::makeClone(ParticleSet& tqp) const
{
  SlaterDetSingle* myclone = new SlaterDetSingle(tqp);
  myclone->Optimizable=Optimizable;
  if (mySPOSet.size() > 1)
  {
    std::map<std::string,SPOSet*>::const_iterator Mit,Lit;
    Mit= mySPOSet.begin();
    Lit= mySPOSet.end();
    while (Mit!=Lit)
    {
      SPOSet* spo = (*Mit).second;
      SPOSet* spo_clone;
      spo_clone = spo->makeClone();
      spo_clone->resetTargetParticleSet(tqp);
      myclone->add(spo_clone,spo->objectName);
      for (int i = 0; i < Dets.size(); ++i)
      {
        if (spo == Dets[i]->getPhi())
        {
          Determinant_t* newD=Dets[i]->makeCopy(dynamic_cast<SPOSetSingle*>(spo_clone));
          newD->resetTargetParticleSet(tqp);
          myclone->add(newD, i);
        }
      }
      Mit++;
    }
  }
  else
  {
    SPOSet* spo = Dets[0]->getPhi();
    SPOSet* spo_clone = spo->makeClone();
    spo_clone->resetTargetParticleSet(tqp);
    myclone->add(spo_clone, spo->objectName);
    for (int i = 0; i < Dets.size(); ++i)
    {
      Determinant_t* newD=Dets[i]->makeCopy(dynamic_cast<SPOSetSingle*>(spo_clone));
      newD->resetTargetParticleSet(tqp);
      myclone->add(newD, i);
    }
  }
  return myclone;
}


}

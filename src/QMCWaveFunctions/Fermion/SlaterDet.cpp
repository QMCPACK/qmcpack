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
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBaseAlternate.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

SlaterDet::SlaterDet(ParticleSet& targetPtcl)
{
  Optimizable = false;
  OrbitalName = "SlaterDet";

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i)-1;

  Dets.resize(targetPtcl.groups(), nullptr);
}

///destructor
SlaterDet::~SlaterDet()
{
  ///clean up SPOSet
}

///add a new SPOSet to the list of determinants
void SlaterDet::add(SPOSetBase* sposet, const std::string& aname)
{
  if (mySPOSet.find(aname) == mySPOSet.end())
  {
    mySPOSet[aname] = sposet;
    sposet->objectName = aname;
  }
  else
  {
    APP_ABORT(" SlaterDet::add(SPOSetBase*, const std::string&) Cannot reuse the " + aname )
    ;
  }
}

///add a new DiracDeterminant to the list of determinants
void SlaterDet::add(Determinant_t* det, int ispin)
{
  if (Dets[ispin]!=nullptr)
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

void SlaterDet::reportStatus(std::ostream& os)
{
}

void SlaterDet::resetTargetParticleSet(ParticleSet& P)
{
  std::map<std::string, SPOSetBasePtr>::iterator sit(mySPOSet.begin());
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

SlaterDet::ValueType SlaterDet::evaluate(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = 1.0;
  for (int i = 0; i < Dets.size(); i++)
    psi *= Dets[i]->evaluate(P, G, L);
  return psi;
}

SlaterDet::RealType SlaterDet::evaluateLog(ParticleSet& P,
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

void SlaterDet::recompute(ParticleSet& P)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->recompute(P);
}

void SlaterDet::evaluateHessian(ParticleSet & P, HessVector_t& grad_grad_psi)
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

void SlaterDet::registerData(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::registerData ",buf.current());
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::registerData ",buf.current());
}

void SlaterDet::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
{
  for (size_t i = 0, n=Dets.size(); i < n; ++i)
  {
    Dets[i]->updateAfterSweep(P,G,L);
  }
}

SlaterDet::RealType SlaterDet::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ",buf.current());
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
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ",buf.current());
  return LogValue;
}

void SlaterDet::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ",buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ",buf.current());
}

OrbitalBasePtr SlaterDet::makeClone(ParticleSet& tqp) const
{
  SlaterDet* myclone = new SlaterDet(tqp);
  myclone->Optimizable=Optimizable;
  if (mySPOSet.size() > 1)
  {
    std::map<std::string,SPOSetBasePtr>::const_iterator Mit,Lit;
    Mit= mySPOSet.begin();
    Lit= mySPOSet.end();
    while (Mit!=Lit)
    {
      SPOSetBasePtr spo = (*Mit).second;
      SPOSetBasePtr spo_clone;
      spo_clone = spo->makeClone();
      spo_clone->resetTargetParticleSet(tqp);
      myclone->add(spo_clone,spo->objectName);
      for (int i = 0; i < Dets.size(); ++i)
      {
        if (spo == Dets[i]->getPhi())
        {
          Determinant_t* newD=Dets[i]->makeCopy(spo_clone);
          newD->resetTargetParticleSet(tqp);
          myclone->add(newD, i);
        }
      }
      Mit++;
    }
  }
//         {
//       for (int i = 0; i < Dets.size(); ++i)
//       {
//         SPOSetBasePtr spo = Dets[i]->getPhi();
//         // Check to see if this determinants SPOSet has already been
//         // cloned
//         bool found = false;
//         SPOSetBasePtr spo_clone;
//         for (int j = 0; j < i; j++)
//           if (spo == Dets[j]->getPhi())
//           {
//             found = true;
//             spo_clone = myclone->Dets[j]->getPhi();
//             spo_clone->resetTargetParticleSet(tqp);
//           }
//         // If it hasn't, clone it now
//         if (!found)
//         {
//           spo_clone = spo->makeClone();
//           spo_clone->resetTargetParticleSet(tqp);
//           myclone->add(spo_clone, spo->objectName);
//         }
//         // Make a copy of the determinant.
//         Determinant_t* newD=Dets[i]->makeCopy(spo_clone);
//         newD->resetTargetParticleSet(tqp);
//         myclone->add(newD, i);
//       }
//     }
  else
  {
    SPOSetBasePtr spo = Dets[0]->getPhi();
    SPOSetBasePtr spo_clone = spo->makeClone();
    spo_clone->resetTargetParticleSet(tqp);
    myclone->add(spo_clone, spo->objectName);
    for (int i = 0; i < Dets.size(); ++i)
    {
      Determinant_t* newD=Dets[i]->makeCopy(spo_clone);
      newD->resetTargetParticleSet(tqp);
      myclone->add(newD, i);
    }
  }
  //map<SPOSetBase*,SPOSetBase*> spomap;
  //SlaterDet* myclone= new SlaterDet(*this);
  //myclone->releasedNode=releasedNode;
  //for(int i=0; i<Dets.size(); i++)
  //{
  //  std::map<SPOSetBase*,SPOSetBase*>::iterator it=spomap.find(Dets[i]->Phi);
  //  if (releasedNode==1)
  //  {
  //    RNDiracDeterminantBase* adet=new RNDiracDeterminantBase(static_cast<RNDiracDeterminantBase&>(*Dets[i]));
  //    adet->NP=0;
  //    if(it == spomap.end())
  //    {
  //      SPOSetBase* newspo=Dets[i]->clonePhi();
  //      spomap[Dets[i]->Phi]=newspo;
  //      adet->Phi=newspo;//assign a new SPOSet
  //    }
  //    else
  //    {
  //      adet->Phi=(*it).second;//safe to transfer
  //    }
  //    adet->resetTargetParticleSet(tqp);
  //    myclone->Dets[i]=adet;
  //  }
  //  else if (releasedNode==2)
  //  {
  //    RNDiracDeterminantBaseAlternate* adet=new RNDiracDeterminantBaseAlternate(static_cast<RNDiracDeterminantBaseAlternate&>(*Dets[i]));
  //    adet->NP=0;
  //    if(it == spomap.end())
  //    {
  //      SPOSetBase* newspo=Dets[i]->clonePhi();
  //      spomap[Dets[i]->Phi]=newspo;
  //      adet->Phi=newspo;//assign a new SPOSet
  //    }
  //    else
  //    {
  //      adet->Phi=(*it).second;//safe to transfer
  //    }
  //    adet->resetTargetParticleSet(tqp);
  //    myclone->Dets[i]=adet;
  //  }
  //  else
  //  {
  //  Determinant_t* adet=new Determinant_t(*Dets[i]);
  //  adet->NP=0;
  //  if(it == spomap.end())
  //  {
  //    SPOSetBase* newspo=Dets[i]->clonePhi();
  //    spomap[Dets[i]->Phi]=newspo;
  //    adet->Phi=newspo;//assign a new SPOSet
  //  }
  //  else
  //  {
  //    adet->Phi=(*it).second;//safe to transfer
  //  }
  //  adet->resetTargetParticleSet(tqp);
  //  myclone->Dets[i]=adet;
  // }
  //}
  return myclone;
}

}

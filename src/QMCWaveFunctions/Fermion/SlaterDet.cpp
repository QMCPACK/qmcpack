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
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

SlaterDet::SlaterDet(ParticleSet& targetPtcl)
{
  Optimizable = false;
  OrbitalName = "SlaterDet";
  M.resize(targetPtcl.groups() + 1, 0);
  for (int i = 0; i < M.size(); ++i)
    M[i] = targetPtcl.first(i);
  DetID.resize(targetPtcl.getTotalNum());
  for (int i = 0; i < targetPtcl.groups(); ++i)
    for (int j = targetPtcl.first(i); j < targetPtcl.last(i); ++j)
      DetID[j] = i;
  Dets.resize(targetPtcl.groups(), 0);
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
  if (Dets[ispin])
  {
    APP_ABORT("SlaterDet::add(Determinant_t* det, int ispin) is alreaded instantiated.");
  }
  else
    Dets[ispin] = det;
  Optimizable = Optimizable || det->Optimizable;
  //int last=Dets.size();
  //Dets.push_back(det);
  //M[last+1]=M[last]+Dets[last]->rows();
  //DetID.insert(DetID.end(),det->rows(),last);
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

void SlaterDet::get_ratios(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->get_ratios(P, ratios);
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

void SlaterDet::registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerDataForDerivatives(P,buf,storageType);
}

SlaterDet::RealType SlaterDet::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L,
    PooledData<RealType>& buf,
    bool fillBuffer )
{
  LogValue = 0.0;
  PhaseValue = 0.0;
  if(fillBuffer)
  {
    for (int i = 0; i < Dets.size(); ++i)
    {
      LogValue +=Dets[i]->evaluateLogForDerivativeBuffer(P, buf);
      Dets[i]->copyToDerivativeBuffer(P, buf);
      PhaseValue += Dets[i]->PhaseValue;
    }
  }
  else
  {
    for (int i = 0; i < Dets.size(); ++i)
    {
      Dets[i]->copyFromDerivativeBuffer(P,buf);
      LogValue += Dets[i]->evaluateLogFromDerivativeBuffer(P, buf);
      PhaseValue += Dets[i]->PhaseValue;
    }
  }
  return LogValue;
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

SlaterDet::RealType SlaterDet::registerData(ParticleSet& P,
    PooledData<RealType>& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::registerData ",buf.current());
  //ValueType psi = 1.0;
  //for(int i=0; i<Dets.size(); i++)
  //  psi *= Dets[i]->registerData(P,buf);
  //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
  LogValue = 0.0;
  PhaseValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
  {
    LogValue += Dets[i]->registerData(P, buf);
    PhaseValue += Dets[i]->PhaseValue;
  }
  DEBUG_PSIBUFFER(" SlaterDet::registerData ",buf.current());
  return LogValue;
}

SlaterDet::RealType SlaterDet::updateBuffer(ParticleSet& P,
    PooledData<RealType>& buf, bool fromscratch)
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

void SlaterDet::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ",buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ",buf.current());
}

/** reimplements the virtual function
 *
 * The DiractDeterminants of SlaterDet need to save the inverse
 * of the determinant matrix to evaluate ratio
 */
void SlaterDet::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->dumpToBuffer(P, buf);
}

/** reimplements the virtual function
 *
 * Matching function to dumpToBuffer.
 */
void SlaterDet::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->dumpFromBuffer(P, buf);
}

SlaterDet::RealType SlaterDet::evaluateLog(ParticleSet& P,
    PooledData<RealType>& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::evaluateLog ",buf.current());
  LogValue = 0.0;
  PhaseValue = 0.0;
  for (int i = 0; i < Dets.size(); i++)
  {
    LogValue += Dets[i]->evaluateLog(P, buf);
    PhaseValue += Dets[i]->PhaseValue;
  }
  DEBUG_PSIBUFFER(" SlaterDet::evaluateLog ",buf.current());
  return LogValue;
}
//SlaterDet::ValueType
//  SlaterDet::evaluate(ParticleSet& P, PooledData<RealType>& buf)
//  {
//    ValueType r=1.0;
//    for(int i=0; i<Dets.size(); i++) 	r *= Dets[i]->evaluate(P,buf);
//    return r;
//  }

OrbitalBasePtr SlaterDet::makeClone(ParticleSet& tqp) const
{
  SlaterDet* myclone = new SlaterDet(tqp);
  myclone->Optimizable=Optimizable;
  myclone->RecomputeNeedsDistanceTable=RecomputeNeedsDistanceTable;
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

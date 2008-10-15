//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  SlaterDet::SlaterDet() 
  {
    Optimizable=false;
    OrbitalName="SlaterDet";
    M.resize(3,0);
  }
  ///destructor
  SlaterDet::~SlaterDet() { }

  ///add a new DiracDeterminant to the list of determinants
  void SlaterDet::add(Determinant_t* det) 
  { 
    int last=Dets.size();
    Dets.push_back(det);
    M[last+1]=M[last]+Dets[last]->rows();
    DetID.insert(DetID.end(),det->rows(),last);
  }

  void SlaterDet::checkInVariables(opt_variables_type& active)
  {
  }

  void SlaterDet::checkOutVariables(const opt_variables_type& active)
  {
  }

  ///reset all the Dirac determinants, Optimizable is true
  void SlaterDet::resetParameters(const opt_variables_type& active) 
  {  
    //if(Optimizable) 
    //  for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
  }

  void SlaterDet::reportStatus(ostream& os)
  {
  }

  void SlaterDet::resetTargetParticleSet(ParticleSet& P) 
  {
    //BasisSet->resetTargetParticleSet(P);
    //LOGMSG("\nSlaterDet::resetTargetParticleSet")
    for(int i=0; i<Dets.size(); i++) Dets[i]->resetTargetParticleSet(P);
  }

  SlaterDet::ValueType 
    SlaterDet::evaluate(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L) 
    {
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      return psi;
    }

  SlaterDet::RealType 
    SlaterDet::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L) 
    {
      //ValueType psi = 1.0;
      //for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
      LogValue=0.0;
      PhaseValue=0.0;
      for(int i=0; i<Dets.size(); ++i)
      {
        LogValue+=Dets[i]->evaluateLog(P,G,L);
        PhaseValue += Dets[i]->PhaseValue;
      }
      return LogValue;
    }

  SlaterDet::RealType SlaterDet::registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    //ValueType psi = 1.0;
    //for(int i=0; i<Dets.size(); i++) 
    //  psi *= Dets[i]->registerData(P,buf);
    //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
    LogValue=0.0;
    PhaseValue=0.0;
    for(int i=0; i<Dets.size(); ++i)
    {
      LogValue+=Dets[i]->registerData(P,buf);
      PhaseValue += Dets[i]->PhaseValue;
    }
    return LogValue;
  }

  SlaterDet::RealType SlaterDet::updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
      bool fromscratch)
  {
    //ValueType psi = 1.0;
    //for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->updateBuffer(P,buf,fromscratch);
    //return LogValue = evaluateLogAndPhase(psi,PhaseValue);
    LogValue=0.0;
    PhaseValue=0.0;
    for(int i=0; i<Dets.size(); ++i)
    {
      LogValue+=Dets[i]->updateBuffer(P,buf,fromscratch);
      PhaseValue+=Dets[i]->PhaseValue;
    }
    return LogValue;
  }

  void SlaterDet::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->copyFromBuffer(P,buf);
  }

  /** reimplements the virtual function
   *
   * The DiractDeterminants of SlaterDet need to save the inverse
   * of the determinant matrix to evaluate ratio
   */
  void SlaterDet::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpToBuffer(P,buf);
  }

  /** reimplements the virtual function
   *
   * Matching function to dumpToBuffer.
   */
  void SlaterDet::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpFromBuffer(P,buf);
  }

  SlaterDet::RealType 
    SlaterDet::evaluateLog(ParticleSet& P, PooledData<RealType>& buf) 
    {
      LogValue=0.0;
      PhaseValue=0.0;
      for(int i=0; i<Dets.size(); i++) 	
      { 
        LogValue += Dets[i]->evaluateLog(P,buf);
        PhaseValue +=Dets[i]->PhaseValue;
      }
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
    map<SPOSetBase*,SPOSetBase*> spomap;
    SlaterDet* myclone= new SlaterDet(*this);
    for(int i=0; i<Dets.size(); i++) 
    {
      map<SPOSetBase*,SPOSetBase*>::iterator it=spomap.find(Dets[i]->Phi);
      Determinant_t* adet=new Determinant_t(*Dets[i]);
      adet->NP=0;
      if(it == spomap.end())
      {
        SPOSetBase* newspo=Dets[i]->clonePhi();
        spomap[Dets[i]->Phi]=newspo;
        adet->Phi=newspo;//assign a new SPOSet
      }
      else
      {
        adet->Phi=(*it).second;//safe to transfer
      }
      adet->resetTargetParticleSet(tqp);
      myclone->Dets[i]=adet;
    }
    return myclone;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

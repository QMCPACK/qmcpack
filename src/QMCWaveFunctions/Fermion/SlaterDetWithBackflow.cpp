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
#include "QMCWaveFunctions/Fermion/SlaterDetWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/RNDiracDeterminantBaseAlternate.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {

  SlaterDetWithBackflow::SlaterDetWithBackflow(ParticleSet& targetPtcl, BackflowTransformation *BF):SlaterDet(targetPtcl),BFTrans(BF)
  {
    Optimizable=true;
    OrbitalName="SlaterDetWithBackflow";
  }

  ///destructor
  SlaterDetWithBackflow::~SlaterDetWithBackflow() 
  { 
    ///clean up SPOSet
  }

  void SlaterDetWithBackflow::get_ratios(ParticleSet& P, vector<ValueType>& ratios)
  {
    for(int i=0; i<Dets.size(); ++i) Dets[i]->get_ratios(P,ratios);
  }


  SlaterDetWithBackflow::ValueType 
    SlaterDetWithBackflow::evaluate(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L) 
    {
      BFTrans->evaluate(P);

      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      return psi;
    }

  SlaterDetWithBackflow::RealType 
    SlaterDetWithBackflow::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L) 
    {
      BFTrans->evaluate(P);

      LogValue=0.0;
      PhaseValue=0.0;
      for(int i=0; i<Dets.size(); ++i)
      {
        LogValue+=Dets[i]->evaluateLog(P,G,L);
        PhaseValue += Dets[i]->PhaseValue;
      }
      return LogValue;
    }

  SlaterDetWithBackflow::RealType SlaterDetWithBackflow::registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    BFTrans->evaluate(P);

    LogValue=0.0;
    PhaseValue=0.0;
    for(int i=0; i<Dets.size(); ++i)
    {
      LogValue+=Dets[i]->registerData(P,buf);
      PhaseValue += Dets[i]->PhaseValue;
    }
    return LogValue;
  }

  SlaterDetWithBackflow::RealType SlaterDetWithBackflow::updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
      bool fromscratch)
  {
    BFTrans->evaluate(P);

    LogValue=0.0;
    PhaseValue=0.0;
    for(int i=0; i<Dets.size(); ++i)
    {
      LogValue+=Dets[i]->updateBuffer(P,buf,fromscratch);
      PhaseValue+=Dets[i]->PhaseValue;
    }
    return LogValue;
  }

  void SlaterDetWithBackflow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->copyFromBuffer(P,buf);
  }

  void SlaterDetWithBackflow::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpToBuffer(P,buf);
  }

  void SlaterDetWithBackflow::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) 
  {
    for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpFromBuffer(P,buf);
  }

  SlaterDetWithBackflow::RealType 
    SlaterDetWithBackflow::evaluateLog(ParticleSet& P, PooledData<RealType>& buf) 
    {
//      BFTrans->evaluate(P);

      LogValue=0.0;
      PhaseValue=0.0;
      for(int i=0; i<Dets.size(); i++) 	
      { 
        LogValue += Dets[i]->evaluateLog(P,buf);
        PhaseValue +=Dets[i]->PhaseValue;
      }
      return LogValue;
    }
  //SlaterDetWithBackflow::ValueType 
  //  SlaterDetWithBackflow::evaluate(ParticleSet& P, PooledData<RealType>& buf) 
  //  {
  //    ValueType r=1.0;
  //    for(int i=0; i<Dets.size(); i++) 	r *= Dets[i]->evaluate(P,buf);
  //    return r;
  //  }

  OrbitalBasePtr SlaterDetWithBackflow::makeClone(ParticleSet& tqp) const
  {
    BackflowTransformation *tr = BFTrans->makeClone();
    SlaterDetWithBackflow* myclone=new SlaterDetWithBackflow(tqp,tr);
    if(mySPOSet.size()>1)//each determinant owns its own set
    {
      for(int i=0; i<Dets.size(); ++i)
      {
        SPOSetBasePtr spo=Dets[i]->getPhi();
	// Check to see if this determinants SPOSet has already been
	// cloned
	bool found = false;
        SPOSetBasePtr spo_clone;
	for (int j=0; j<i; j++)
	  if (spo == Dets[j]->getPhi()) {
	    found = true;
	    spo_clone = myclone->Dets[j]->getPhi();
	  }
	// If it hasn't, clone it now
	if (!found) {
	  spo_clone=spo->makeClone();
	  myclone->add(spo_clone,spo->objectName);
	}
	// Make a copy of the determinant.
        myclone->add(Dets[i]->makeCopy(spo_clone),i);
      }
    }
    else
    {
      SPOSetBasePtr spo=Dets[0]->getPhi();
      SPOSetBasePtr spo_clone=spo->makeClone();
      myclone->add(spo_clone,spo->objectName);
      for(int i=0; i<Dets.size(); ++i)
        myclone->add(Dets[i]->makeCopy(spo_clone),i);
    }

    return myclone;
  }

}
/***************************************************************************
 * $RCSfile$   $Author: kpesler $
 * $Revision: 4721 $   $Date: 2010-03-12 17:11:47 -0600 (Fri, 12 Mar 2010) $
 * $Id: SlaterDetWithBackflow.cpp 4721 2010-03-12 23:11:47Z kpesler $ 
 ***************************************************************************/

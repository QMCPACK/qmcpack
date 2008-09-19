//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim 
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
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  LocalECPotential::LocalECPotential(const ParticleSet& ions, ParticleSet& els):
    IonConfig(ions)
  { 
    NumIons=ions.getTotalNum();
    myTableIndex=els.addTable(ions);
    //allocate null
    PPset.resize(ions.getSpeciesSet().getTotalNum(),0);
    PP.resize(NumIons,0);
    Zeff.resize(NumIons,0.0);
    gZeff.resize(ions.getSpeciesSet().getTotalNum(),0);
  } 

  ///destructor
  LocalECPotential::~LocalECPotential() { 
    delete_iter(PPset.begin(),PPset.end());
    //map<int,RadialPotentialType*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    //while(pit != pit_end) {
    //  delete (*pit).second; ++pit;
    //}
  }

  void LocalECPotential::resetTargetParticleSet(ParticleSet& P) {
    int tid=P.addTable(IonConfig);
    if(tid != myTableIndex)
    {
      APP_ABORT("  LocalECPotential::resetTargetParticleSet found a different distance table index.");
    }
  }

  void LocalECPotential::add(int groupID, RadialPotentialType* ppot, RealType z) {
    PPset[groupID]=ppot;
    gZeff[groupID]=z;
    for(int iat=0; iat<PP.size(); iat++) {
      if(IonConfig.GroupID[iat]==groupID) {
        PP[iat]=ppot;
        Zeff[iat]=z;
      }
    }
  }

  LocalECPotential::Return_t 
    LocalECPotential::evaluate(ParticleSet& P) 
    {
      const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
      Value=0.0;
      //loop over all the ions
      for(int iat=0; iat<NumIons; iat++) {
        RadialPotentialType* ppot(PP[iat]);
        if(ppot==0) continue;
        Return_t esum(0.0);
        for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
          esum += ppot->splint(d_table.r(nn))*d_table.rinv(nn);
        //count the sign and effective charge
        Value -= esum*Zeff[iat];
      }
      return Value;
    }

  LocalECPotential::Return_t 
    LocalECPotential::registerData(ParticleSet& P, BufferType& buffer) 
    {
      PPart.resize(P.getTotalNum());
      NewValue=Value=evaluateForPbyP(P);
      buffer.add(PPart.begin(),PPart.end());
      buffer.add(Value);
      return Value;
    }

  LocalECPotential::Return_t 
    LocalECPotential::updateBuffer(ParticleSet& P, BufferType& buffer) 
    {
      NewValue=Value=evaluateForPbyP(P);
      buffer.put(PPart.begin(),PPart.end());
      buffer.put(Value);
      return Value;
    }

  void LocalECPotential::copyFromBuffer(ParticleSet& P, BufferType& buffer) 
  {
    buffer.get(PPart.begin(),PPart.end());
    buffer.get(Value);
  }

  void LocalECPotential::copyToBuffer(ParticleSet& P, BufferType& buffer) 
  {
    buffer.put(PPart.begin(),PPart.end());
    buffer.put(Value);
  }

  LocalECPotential::Return_t
    LocalECPotential::evaluateForPbyP(ParticleSet& P)
    {
      const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
      PPart=0.0;
      Return_t res=0.0;
      for(int iat=0; iat<NumIons; ++iat) {
        RadialPotentialType* ppot(PP[iat]);
        if(ppot) {
          Return_t z=-Zeff[iat];
          for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
          {
            Return_t e= z*ppot->splint(d_table.r(nn))*d_table.rinv(nn);
            PPart[d_table.J[nn]]+=e;
            res+=e;
          }
        }
      }
      return res;
    }

  LocalECPotential::Return_t
    LocalECPotential::evaluatePbyP(ParticleSet& P, int active)
    {
      const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[myTableIndex]->Temp);
      PPtmp=0.0;
      for(int iat=0; iat<NumIons; ++iat) 
      {
        if(PP[iat]) 
          PPtmp -= Zeff[iat]*PP[iat]->splint(temp[iat].r1)*temp[iat].rinv1;
      }
      return NewValue=Value+PPtmp-PPart[active];
    }

  void LocalECPotential::acceptMove(int active)
  {
    PPart[active]=PPtmp;
    Value=NewValue;
  }

  QMCHamiltonianBase* LocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    LocalECPotential* myclone=new LocalECPotential(IonConfig,qp);
    for(int ig=0; ig<PPset.size(); ++ig)
    {
      if(PPset[ig]) 
      {
        RadialPotentialType* ppot=PPset[ig]->makeClone();
        myclone->add(ig,ppot,gZeff[ig]);
      }
    }
    return myclone;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


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
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  LocalECPotential::LocalECPotential(const ParticleSet& ions, ParticleSet& els):
    IonConfig(ions), d_table(0)
  { 
    NumIons=ions.getTotalNum();
    d_table = DistanceTable::add(ions,els);
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
    d_table = DistanceTable::add(IonConfig,P);
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
    LocalECPotential::evaluate(ParticleSet& P) {
      Value=0.0;
      //loop over all the ions
      for(int iat=0; iat<NumIons; iat++) {
        RadialPotentialType* ppot(PP[iat]);
        if(ppot) {
          Return_t esum(0.0);
          for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++){
            esum += ppot->splint(d_table->r(nn))*d_table->rinv(nn);
          }
          //count the sign and effective charge
          Value -= esum*Zeff[iat];
        }
      }
      return Value;
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


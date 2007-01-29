//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/NonLocalECPotential.h"

namespace qmcplusplus {

  void NonLocalECPotential::resetTargetParticleSet(ParticleSet& P) {
    d_table = DistanceTable::add(IonConfig,P);
  }

  /** constructor
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\param psi trial wavefunction
   */
  NonLocalECPotential::NonLocalECPotential(ParticleSet& ions, ParticleSet& els,
      TrialWaveFunction& psi): IonConfig(ions), d_table(0), Psi(psi)
  { 
    d_table = DistanceTable::add(ions,els);
    NumIons=ions.getTotalNum();
    //els.resizeSphere(NumIons);
    PP.resize(NumIons,0);
  }

  ///destructor
  NonLocalECPotential::~NonLocalECPotential() { 
    map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    while(pit != pit_end) {
       delete (*pit).second; ++pit;
    }
  }

  NonLocalECPotential::Return_t
  NonLocalECPotential::evaluate(ParticleSet& P) { 
    Value=0.0;
    //loop over all the ions
    for(int iat=0; iat<NumIons; iat++) {
      if(PP[iat]) {
        PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
        Value += PP[iat]->evaluate(P,iat,Psi);
      }
    }
    return Value;
  }

  NonLocalECPotential::Return_t
  NonLocalECPotential::evaluate(ParticleSet& P, vector<NonLocalData>& Txy) { 
    Value=0.0;
    //loop over all the ions
    for(int iat=0; iat<NumIons; iat++) {
      if(PP[iat]) {
        PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
        Value += PP[iat]->evaluate(P,Psi,iat,Txy);
      }
    }
    return Value;
  }

  void 
  NonLocalECPotential::add(int groupID, NonLocalECPComponent* ppot) {
    map<int,NonLocalECPComponent*>::iterator pit(PPset.find(groupID));
    ppot->myTable=d_table;
    if(pit  == PPset.end()) {
      for(int iat=0; iat<PP.size(); iat++) {
        if(IonConfig.GroupID[iat]==groupID) PP[iat]=ppot;
      }
      PPset[groupID]=ppot;
    }
  }

  QMCHamiltonianBase* NonLocalECPotential::clone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new NonLocalECPotential(IonConfig,qp,psi);
  }


  void NonLocalECPotential::setRandomGenerator(RandomGenerator_t* rng)
  {
    map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    while(pit != pit_end) {
      (*pit).second->setRandomGenerator(rng);
      ++pit;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

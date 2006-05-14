//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim 
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
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/LocalECPotential.h"

namespace qmcplusplus {

  LocalECPotential::LocalECPotential(ParticleSet& ions, ParticleSet& els):
    IonConfig(ions), d_table(0)
  { 
    NumIons=ions.getTotalNum();
    d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
    //allocate null
    PP.resize(NumIons,0);
  } 

  ///destructor
  LocalECPotential::~LocalECPotential() { 
    map<int,RadialPotentialType*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    while(pit != pit_end) {
      delete (*pit).second; ++pit;
    }
  }

  void LocalECPotential::resetTargetParticleSet(ParticleSet& P) {
    d_table = DistanceTable::getTable(DistanceTable::add(IonConfig,P));
  }

  void LocalECPotential::add(int groupID, RadialPotentialType* ppot) {
    map<int,RadialPotentialType*>::iterator pit(PPset.find(groupID));
    if(pit  == PPset.end()) {
      for(int iat=0; iat<PP.size(); iat++) {
        if(IonConfig.GroupID[iat]==groupID) PP[iat]=ppot;
      }
      PPset[groupID]=ppot;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


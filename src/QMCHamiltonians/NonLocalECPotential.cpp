//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus {

  void NonLocalECPotential::resetTargetParticleSet(ParticleSet& P) 
  {
    d_table = DistanceTable::add(IonConfig,P);
  }

  /** constructor
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\param psi trial wavefunction
   */
  NonLocalECPotential::NonLocalECPotential(ParticleSet& ions, ParticleSet& els,
					   TrialWaveFunction& psi, 
					   bool computeForces): 
    IonConfig(ions), d_table(0), Psi(psi), 
    ComputeForces(computeForces), ForceBase(ions,els)
  { 
    d_table = DistanceTable::add(ions,els);
    NumIons=ions.getTotalNum();
    //els.resizeSphere(NumIons);
    PP.resize(NumIons,0);
    prefix="FNL";
    PPset.resize(IonConfig.getSpeciesSet().getTotalNum(),0);
    PulayTerm.resize(NumIons);
  }

  ///destructor
  NonLocalECPotential::~NonLocalECPotential() 
  { 
    delete_iter(PPset.begin(),PPset.end());
    //map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    //while(pit != pit_end) {
    //   delete (*pit).second; ++pit;
    //}
  }

  NonLocalECPotential::Return_t
  NonLocalECPotential::evaluate(ParticleSet& P) { 
    Value=0.0;
    //loop over all the ions
    if (ComputeForces) {
      for(int iat=0; iat<NumIons; iat++) 
	if(PP[iat]) {
	  PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
	  //Value += PP[iat]->evaluate(P,iat,Psi, forces[iat]);
	  Value += PP[iat]->evaluate(P,IonConfig,iat,Psi, forces[iat], 
	   			     PulayTerm[iat]);
	}
    }
    else {
      for(int iat=0; iat<NumIons; iat++) 
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
    if (ComputeForces) {
      for(int iat=0; iat<NumIons; iat++) 
	if(PP[iat]) {
	  PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
	  Value += PP[iat]->evaluate(P,Psi,iat,Txy, forces[iat]);
	}
    }
    else {
      for(int iat=0; iat<NumIons; iat++) {
	if(PP[iat]) {
	  PP[iat]->randomize_grid(*(P.Sphere[iat]),UpdateMode[PRIMARY]);
	  Value += PP[iat]->evaluate(P,Psi,iat,Txy);
	}
      }
    }
    return Value;
  }

  void 
  NonLocalECPotential::add(int groupID, NonLocalECPComponent* ppot) {
    //map<int,NonLocalECPComponent*>::iterator pit(PPset.find(groupID));
    //ppot->myTable=d_table;
    //if(pit  == PPset.end()) {
    //  for(int iat=0; iat<PP.size(); iat++) {
    //    if(IonConfig.GroupID[iat]==groupID) PP[iat]=ppot;
    //  }
    //  PPset[groupID]=ppot;
    //}
    ppot->myTable=d_table;
    for(int iat=0; iat<PP.size(); iat++) 
      if(IonConfig.GroupID[iat]==groupID) PP[iat]=ppot;
    PPset[groupID]=ppot;
  }

  QMCHamiltonianBase* NonLocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    NonLocalECPotential* myclone=new NonLocalECPotential(IonConfig,qp,psi,
							 ComputeForces);

    for(int ig=0; ig<PPset.size(); ++ig)
    {
      if(PPset[ig]) myclone->add(ig,PPset[ig]->makeClone());
    }

    //resize sphere
    qp.resizeSphere(IonConfig.getTotalNum());
    for(int ic=0; ic<IonConfig.getTotalNum(); ic++) {
      if(PP[ic] && PP[ic]->nknot) qp.Sphere[ic]->resize(PP[ic]->nknot);
    }
    return myclone;
  }


  void NonLocalECPotential::setRandomGenerator(RandomGenerator_t* rng)
  {
    for(int ig=0; ig<PPset.size(); ++ig)
      if(PPset[ig]) PPset[ig]->setRandomGenerator(rng);
    //map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
    //while(pit != pit_end) {
    //  (*pit).second->setRandomGenerator(rng);
    //  ++pit;
    //}
  }


  void NonLocalECPotential::addObservables(PropertySetType& plist
      , BufferType& collectables)
  {
    QMCHamiltonianBase::addValue(plist);
    if (ComputeForces) {
      if(FirstForceIndex<0) 
        FirstForceIndex=plist.size();
      for(int iat=0; iat<Nnuc; iat++) {
        for(int x=0; x<OHMMS_DIM; x++) {
          ostringstream obsName1, obsName2;
          obsName1 << "FNL" << "_" << iat << "_" << x;
          plist.add(obsName1.str());
          obsName2 << "FNL_Pulay" << "_" << iat << "_" << x;
          plist.add(obsName2.str());
        }
      }
    }
  }

  void 
  NonLocalECPotential::registerObservables(vector<observable_helper*>& h5list,
				  hid_t gid) const
  {
    QMCHamiltonianBase::registerObservables(h5list, gid);
    if (ComputeForces) 
    {
      vector<int> ndim(2);
      ndim[0]=Nnuc;
      ndim[1]=OHMMS_DIM;
      observable_helper* h5o1 = new observable_helper("FNL");
      h5o1->set_dimensions(ndim,FirstForceIndex);
      h5o1->open(gid);
      h5list.push_back(h5o1);

      observable_helper* h5o2 = new observable_helper("FNL_Pulay");
      h5o2->set_dimensions(ndim,FirstForceIndex+Nnuc*OHMMS_DIM);
      h5o2->open(gid);
      h5list.push_back(h5o2);
    }
  }

  void 
  NonLocalECPotential::setObservables(QMCTraits::PropertySetType& plist) 
  {
    QMCHamiltonianBase::setObservables(plist);
    if (ComputeForces) {
      int index = FirstForceIndex;
      for(int iat=0; iat<Nnuc; iat++) {
	for(int x=0; x<OHMMS_DIM; x++) {
	  plist[index++] = forces[iat][x];
	  plist[index++] = PulayTerm[iat][x];
	}
      }
    }
  }
   

  void 
  NonLocalECPotential::setParticlePropertyList(QMCTraits::PropertySetType& plist, 
				      int offset)
  {
    QMCHamiltonianBase::setParticlePropertyList (plist, offset);
    if (ComputeForces) {
      int index = FirstForceIndex + offset;
      for(int iat=0; iat<Nnuc; iat++) {
	for(int x=0; x<OHMMS_DIM; x++) {
	  plist[index++] = forces[iat][x];
	  plist[index++] = PulayTerm[iat][x];
	}
      }
    }
  }


}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QuantumSystem/STOBuilder.h"
#include "QuantumSystem/CompositeSlaterDeterminant.h"
#include "QuantumSystem/JastrowFunction.h"
#include "QuantumSystem/TrialWaveFunction.h"
using namespace qmcplusplus;

void 
initialize(TrialWaveFunction& Psi, MCWalkerConfiguration& el, 
	   ParticleBase& ion) {

  enum {Up = 0, Down};

  ///add a distance table for ee IndexTypeeractions
  IndexType iee = DistanceTable::add(el,"ee");
  DistanceTableData* d_ee = DistanceTable::getTable(iee);
  ///add a distance table for ie interactions
  int iei = DistanceTable::add(ion,el,"ie");
  DistanceTableData* d_ei = DistanceTable::getTable(iei);

  ///create molecular orbital
  MolecularOrbitals<SlaterTypeOrbitals_t> *MO =
    new MolecularOrbitals<SlaterTypeOrbitals_t>;

  AtomicOrbitalBuilder<SlaterTypeOrbitals_t> orbitalbuilder;
  for(int iat=0; iat<ion.getTotalNum(); iat++) {
    orbitalbuilder.add("He",*MO);
  }
  MO->setTable(d_ei);

  typedef SlaterDeterminant<MolecularOrbitals<SlaterTypeOrbitals_t> > Det_t;
  Det_t *DetU = new Det_t(*MO,el.first(Up));
  DetU->set(el.first(Up),el.last(Up)-el.first(Up));

  Det_t* DetD = new Det_t(*MO,el.first(Down));
  DetD->set(el.first(Down),el.last(Down)-el.first(Down));

  LinearSlaterDeterminant<MolecularOrbitals<SlaterTypeOrbitals_t> >
    *asymmpsi 
    = new LinearSlaterDeterminant<MolecularOrbitals<SlaterTypeOrbitals_t> >;

  asymmpsi->add(DetU);
  asymmpsi->add(DetD);

  Psi.add(asymmpsi,d_ei);

  /*
  OneBodyJastrow<PadeJastrow<ValueType>,FastWalkerIndex> *Jie 
    = new OneBodyJastrow<PadeJastrow<ValueType>,FastWalkerIndex>;
  Jie->F.push_back(PadeJastrow<ValueType>(1.0,1.0));
  Psi.add(Jie,d_ei);
  */

  TwoBodyJastrow<PadeJastrow<ValueType>,FastWalkerIndex> *Jee
    = new TwoBodyJastrow<PadeJastrow<ValueType>,FastWalkerIndex>;

  Jee->F.push_back(PadeJastrow<double>(0.5,1.0));
  Jee->F.push_back(PadeJastrow<double>(1.0/6.0,1.0));
  Jee->F.push_back(PadeJastrow<double>(1.0/6.0,1.0));
  Jee->F.push_back(PadeJastrow<double>(0.5,1.0));
  Psi.add(Jee,d_ee);

}

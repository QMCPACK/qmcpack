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
//#include "QuantumSystem/STOBuilder.h"
//#include "QuantumSystem/CompositeSlaterDeterminant.h"
//#include "QuantumSystem/JastrowFunction.h"
#include "QuantumSystem/TrialWaveFunctionBuilder.h"
#include "QMCBase/RandomSeqGenerator.h"
using namespace ohmmsqmc;

//extern void 
//initialize(TrialWaveFunction&, MCWalkerConfiguration&, ParticleSet&);

int main(int argc, char **argv) {

  IndexType nblocks = 1;
  IndexType nsteps = 100;
  IndexType nup = 1;
  IndexType ndown = 1;
  IndexType nw = 10;
  RealType Tau = 10;

  //Random.init(0,1,-1);
  ///create electrons: using MCWalkerConfiguration
  MCWalkerConfiguration el;
  TinyVector<IndexType,2> N(nup,ndown);
  el.create(N);
  makeGaussRandom(el.R);

  ///create walkers and sample new configurations
  el.createWalkers(nw);
  for (MCWalkerConfiguration::iterator it = el.begin(); it != el.end(); ++it) {
    makeGaussRandom((*it)->R);
  }

  //el.R[0] = PosType(-0.0692796, 0.0541703, 0.0816451);
  //el.R[1] = PosType(0.103747, -0.330282, 0.237803);
  ///create ions
  ParticleSet ion;
  ion.create(1);
  ion.R[0] = 0.0;

  IndexType iee = DistanceTable::add(el,"ee");
  IndexType iei = DistanceTable::add(ion,el,"ie");

  ///create a trial wave function
  TrialWaveFunction Psi;
  QMCDOMNode<TrialWaveFunction> wfs_reader(Psi);
  //initialize(Psi,el,ion);

  const ValueType delta = 0.0001;
  ValueType deltainv = 1.0/delta/2.0;
  ValueType dh2 = 1.0/delta/delta;

  DistanceTable::update(el);
  ValueType psi = Psi.evaluate(el);    

  IndexType nat = el.getTotalNum();
  ParticleSet::ParticlePos_t R(nat);
  ParticleSet::ParticleGradient_t G0(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L0(nat),L1(nat);
  G0 = el.G;
  L0 = el.L;
  R = el.R;

  for(IndexType iat=0; iat<el.getTotalNum(); iat++) {
    ValueType lap = 0.0;
    for(int idim=0; idim<DIM; idim++) {
      el.R[iat][idim] += delta;
      DistanceTable::update(el);
      ValueType psip =  0.5*log(pow(Psi.evaluate(el),2));
      el.R[iat][idim] -= 2.0*delta;
      DistanceTable::update(el);
      ValueType psim = 0.5*log(pow(Psi.evaluate(el),2));

      G1[iat][idim] = (psip-psim);
      lap += psip+psim;    
      el.R[iat][idim] += delta;
    }    
    G1[iat] *= deltainv;
    L1[iat]  = dh2*(lap-3.0*log(pow(psi,2)));
  }

  cout << endl;
  cout << "Difference between the analytic and finite-difference method " << endl;
  cout << "gradient(ln(|det|)) " << endl;
  for(IndexType iat=0; iat<nat; iat++)  
    cout << G0[iat] << G0[iat]-G1[iat] << endl;

  cout << "laplacian(ln(|det|)) " << endl;
  for(IndexType iat=0; iat<nat; iat++)  
    cout << L0[iat] << " " << L1[iat] << " " << L0[iat]-L1[iat] << endl;

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

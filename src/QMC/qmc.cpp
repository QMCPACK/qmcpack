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
#include "QMCWaveFunctions/TrialWaveFunctionBuilder.h"
#include "QMCHamiltonians/HartreePotential.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"
#include "QMCBase/RandomSeqGenerator.h"
#include "QMC/VMC.h"
#include "QMC/DMC.h"
#include "QMC/EstimatorManager.h"
#include "OhmmsData/DOMProcessor.h"

namespace ohmmsqmc {
  extern bool initQuantumParticle(MCWalkerConfiguration& el, DOMProcessor& reader);
  extern bool initTrialWaveFunction(TrialWaveFunction& Psi, DOMProcessor& reader);
}

using namespace ohmmsqmc;

int main(int argc, char **argv) {

  int nblocks = 10;
  int nsteps = 10000;
  int nup = 1;
  int ndown = 1;
  int nw = 1;
  RealType Tau = 0.01;

  int iargc = 0;
  while(iargc<argc) {
    if(!strcmp(argv[iargc],"--blocks")) {
      nblocks = atoi(argv[++iargc]);
    } else if(!strcmp(argv[iargc],"--steps")) {
      nsteps = atoi(argv[++iargc]);
    } else if(!strcmp(argv[iargc],"--walkers")) {
      nw = atoi(argv[++iargc]);
    } else if(!strcmp(argv[iargc],"--up")) {
      nup = atoi(argv[++iargc]);
    } else if(!strcmp(argv[iargc],"--down")) {
      ndown = atoi(argv[++iargc]);
    } else if(!strcmp(argv[iargc],"--tau")) {
      Tau = atof(argv[++iargc]);
    }
    iargc++;
  }

  Random.init(0,1,0);
  DOMProcessor reader(argv[1]);

  MCWalkerConfiguration el;
  initQuantumParticle(el,reader);

  ///create ions
  ParticleBase ion;
  ion.create(1);
  ion.R[0] = 0.0;

  IndexType iee = DistanceTable::add(el,"ee");
  IndexType iei = DistanceTable::add(ion,el,"ie");

  ///create a trial wave function
  TrialWaveFunction Psi;
  initTrialWaveFunction(Psi,reader);

  nup = el.last(0);
  QMCHamiltonian H;

  DistanceTableData* d_ee = DistanceTable::getTable(iee);
  DistanceTableData* d_ei = DistanceTable::getTable(iei);
  H.add(new CoulombPotential(el.getTotalNum()), d_ei);
  H.add(new HartreePotential, d_ee);
  H.add(new BareKineticEnergy, NULL);

  EstimatorManager Estimators;

  Estimators.reset("vmc");
  VMC vmc(el,Psi,H,Estimators);
  vmc.initialize(nw);
  vmc.run(nblocks,nsteps,Tau);

/*
  Estimators.reset("dmc");
  DMC dmc(el,Psi,H,Estimators);
  dmc.run(nblocks,nsteps,Tau);
*/
}



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
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Message/CommCreate.h"
#include "Utilities/Clock.h"
#include "QMC/WaveFunctionTester.h"
#include "Utilities/OhmmsInform.h"
#include "math.h"
using namespace ohmmsqmc;


WaveFunctionTester::WaveFunctionTester(MCWalkerConfiguration& w, 
				       TrialWaveFunction& psi, 
				       QMCHamiltonian& h, 
				       xmlNodePtr q):
  QMCDriver(w,psi,h,q) { }


/*!
 * \brief Test the evaluation of the wavefunction, gradient and laplacian
 by comparing to the numerical evaluation.
 *
 Use the finite difference formulas formulas
 \f[ 
 \nabla_i f({\bf R}) = \frac{f({\bf R+\Delta r_i}) - f({\bf R})}{2\Delta r_i}
 \f]
 and
 \f[
 \nabla_i^2 f({\bf R}) = \sum_{x,y,z} \frac{f({\bf R}+\Delta x_i) 
 - 2 f({\bf R}) + f({\bf R}-\Delta x_i)}{2\Delta x_i^2},
 \f]
 where \f$ f = \ln \Psi \f$ and \f$ \Delta r_i \f$ is a 
 small displacement for the ith particle.
*/

bool 
WaveFunctionTester::run() {

  LogOut->getStream() << "Starting a Wavefucntion tester" << endl;

  DistanceTable::create(1);

  put(qmc_node);

  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 0.001;
  RealType delta2 = 2*delta;
  ValueType c1 = 1.0/delta/2.0;
  ValueType c2 = 1.0/delta/delta;

  int nat = W.getTotalNum();

  ParticleSet::ParticlePos_t deltaR(nat);
  MCWalkerConfiguration::PropertyContainer_t Properties;
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

  //copy the properties of the working walker
  Properties = awalker->Properties;
  
  //sample a new walker configuration and copy to ParticleSet::R
  makeGaussRandom(deltaR);
 
  W.R = awalker->R;

  W.R += deltaR;

  DistanceTable::update(W);
  ValueType psi =log(fabs(Psi.evaluate(W)));

  ParticleSet::ParticlePos_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;

  for(int iat=0; iat<nat; iat++) {
    PosType r0 = W.R[iat];
    PosType g0;  
    ValueType lap = 0.0;
    for(int idim=0; idim<3; idim++) {
   
      W.R[iat][idim] = r0[idim]+delta;         
      DistanceTable::update(W);
      ValueType psi_p = log(fabs(Psi.evaluate(W)));

      W.R[iat][idim] = r0[idim]-delta;         
      DistanceTable::update(W);
      ValueType psi_m = log(fabs(Psi.evaluate(W)));

      lap += psi_m + psi_p;
      g0[idim] = psi_p - psi_m;

      W.R[iat] = r0;
    }

    G1[iat] = c1*g0;
    L1[iat] = c2*(lap-6.0*psi);
  }
  
  cout << "Gradients:" << endl;
  cout.precision(15);
  for(int iat=0; iat<nat; iat++) {
    cout.precision(15);
    cout << "For particle #" << iat << " at coordinate " << W.R[iat] << endl;
    cout << "Code result =         " << setw(12) << G[iat] << endl 
	 << "Finite diff. result = " << setw(12) << G1[iat] << endl << endl;
  }

  cout << "Laplacians:" << endl;
  cout.precision(15);
  for(int iat=0; iat<nat; iat++) {
    cout << "For particle #" << iat << " at coordinate " << W.R[iat] << endl;
    cout << "Code result =         " << setw(12) << L[iat] << endl 
	 << "Finite diff. result = " << setw(12) << L1[iat] << endl << endl;
  }

  return true;
} 


bool 
WaveFunctionTester::put(xmlNodePtr q){
  return putQMCInfo(q);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

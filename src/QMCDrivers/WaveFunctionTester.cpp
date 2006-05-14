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
#include "Message/Communicate.h"
#include "Utilities/Clock.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "Utilities/OhmmsInform.h"
#include "math.h"
using namespace qmcplusplus;


WaveFunctionTester::WaveFunctionTester(MCWalkerConfiguration& w, 
				       TrialWaveFunction& psi, 
				       QMCHamiltonian& h):
  QMCDriver(w,psi,h),checkRatio("no") { 
    m_param.add(checkRatio,"ratio","string");
  }


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

  app_log() << "Starting a Wavefunction tester" << endl;

  DistanceTable::create(1);

  put(qmcNode);

  if(checkRatio == "yes") {
    runRatioTest();
  } else {
    runBasicTest();
  }
  return true;
}

void WaveFunctionTester::runBasicTest() {
  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 0.0001;
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
  //makeGaussRandom(deltaR);
 
  W.R = awalker->R;

  //W.R += deltaR;

  W.update();
  //DistanceTable::update(W);
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
      W.update();
      ValueType psi_p = log(fabs(Psi.evaluate(W)));

      W.R[iat][idim] = r0[idim]-delta;         
      W.update();
      ValueType psi_m = log(fabs(Psi.evaluate(W)));
      lap += psi_m + psi_p;
      g0[idim] = psi_p - psi_m;

      W.R[iat] = r0;
    }

    G1[iat] = c1*g0;
    L1[iat] = c2*(lap-6.0*psi);
  }
  
  cout.precision(15);
  for(int iat=0; iat<nat; iat++) {
    cout.precision(15);
    cout << "For particle #" << iat << " at " << W.R[iat] << endl;
    cout << "Gradient      = " << setw(12) << G[iat] << endl 
	 << "  Finite diff = " << setw(12) << G1[iat] << endl 
	 << "  Error       = " << setw(12) << G[iat]-G1[iat] << endl << endl;
    cout << "Laplacian     = " << setw(12) << L[iat] << endl 
	 << "  Finite diff = " << setw(12) << L1[iat] << endl 
	 << "  Error       = " << setw(12) << L[iat]-L1[iat] << endl << endl;
  }

  makeGaussRandom(deltaR);
  //testing ratio alone
  for(int iat=0; iat<nat; iat++) {
    W.update();
    ValueType psi_p = log(fabs(Psi.evaluate(W)));

    W.makeMove(iat,deltaR[iat]);
    RealType aratio=Psi.ratio(W,iat);
    W.rejectMove(iat);

    W.R[iat] += deltaR[iat];         
    W.update();
    ValueType psi_m = log(fabs(Psi.evaluate(W)));

    cout << iat << " ratio " << aratio << " " << exp(psi_m-psi_p) << endl;
  }
} 

void WaveFunctionTester::runRatioTest() {

  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t G(nat), dG(nat);
  ParticleSet::ParticleLaplacian_t L(nat), dL(nat);

  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while(it != it_end) {
    (*it)->DataSet.rewind();
    W.registerData(**it,(*it)->DataSet);
    ValueType logpsi=Psi.registerData(W,(*it)->DataSet);
    RealType vsq = Dot(W.G,W.G);
    (*it)->Drift = Tau*W.G;

    RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi.getSign(),ene);
    H.saveProperty((*it)->getPropertyBase());
    ++it;
  } 

  it=W.begin();
  while(it != it_end) {

    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

    w_buffer.rewind();
    W.copyFromBuffer(w_buffer);
    Psi.copyFromBuffer(W,w_buffer);

    ValueType eold(thisWalker.Properties(LOCALENERGY));
    ValueType logpsi(thisWalker.Properties(LOGPSI));
    ValueType emixed(eold), enew(eold);

    makeGaussRandom(deltaR);

    //mave a move
    RealType ratio_accum=1.0;
    for(int iat=0; iat<nat; iat++) {
      PosType dr(Tau*deltaR[iat]);

      PosType newpos(W.makeMove(iat,dr));

      RealType ratio=Psi.ratio(W,iat,dG,dL);

      //if(ratio > 0 && iat%2 == 1) {//accept only the even index
      if(ratio > 0) {
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);

        W.G += dG;
        W.L += dL;
        //G=W.G+dG;
        //L=W.L+dL;
        ratio_accum *= ratio;
      } else {
        cout << " Rejecting a move for " << iat << endl;
        W.rejectMove(iat); Psi.rejectMove(iat);
      }
    }

    G = W.G;
    L = W.L;

    w_buffer.rewind();
    W.copyToBuffer(w_buffer);
    ValueType psi = Psi.evaluate(W,w_buffer);

    W.update();
    ValueType newlogpsi=Psi.evaluateLog(W);

    cout << " Ratio " << ratio_accum*ratio_accum 
      << " | " << exp(2.0*(newlogpsi-logpsi)) << " " << ratio_accum*ratio_accum/exp(2.0*(newlogpsi-logpsi)) << endl
      << " new log(psi) " << newlogpsi 
      << " old log(psi) " << logpsi << endl;

    cout << " Gradients " << endl;
    for(int iat=0; iat<nat; iat++) {
      cout << W.G[iat]-G[iat] << W.G[iat] << endl; //W.G[iat] << G[iat] << endl;
    }
    cout << " Laplacians " << endl;
    for(int iat=0; iat<nat; iat++) {
      cout << W.L[iat]-L[iat] << endl;
    }
    ++it;
  }
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

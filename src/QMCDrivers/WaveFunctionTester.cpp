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
using namespace std;
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/OhmmsInform.h"
#include "LongRange/StructFact.h"
namespace qmcplusplus {


WaveFunctionTester::WaveFunctionTester(MCWalkerConfiguration& w, 
				       TrialWaveFunction& psi, 
				       QMCHamiltonian& h,
				       ParticleSetPool &ptclPool):
  QMCDriver(w,psi,h),checkRatio("no"),checkClone("no"), checkHamPbyP("no"),
  PtclPool(ptclPool)
  { 
    m_param.add(checkRatio,"ratio","string");
    m_param.add(checkClone,"clone","string");
    m_param.add(checkHamPbyP,"hamiltonianpbyp","string");
    m_param.add(sourceName,"source","string");
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

  //DistanceTable::create(1);

  put(qmcNode);

  if(checkRatio == "yes") 
  {
    runRatioTest();
    runRatioTest2();
  }
  else if (checkClone == "yes") 
    runCloneTest();
  else if (sourceName.size() != 0)
    runGradSourceTest();
  else
    runBasicTest();

  RealType ene = H.evaluate(W);

  cout << " Energy " << ene << endl;

  return true;
}

  void WaveFunctionTester::runCloneTest() {

    for(int iter=0; iter<4; ++iter)
    {
    ParticleSet* w_clone = new MCWalkerConfiguration(W);
    TrialWaveFunction *psi_clone = Psi.makeClone(*w_clone);
    QMCHamiltonian *h_clone = H.makeClone(*w_clone,*psi_clone);
    h_clone->setPrimary(false);

    IndexType nskipped = 0;
    RealType sig2Enloc=0, sig2Drift=0;
    RealType delta = 0.00001;
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
    
    W.R = awalker->R;
    W.update();
    ValueType logpsi1 = Psi.evaluateLog(W);
    RealType eloc1  = H.evaluate(W);

    w_clone->R=awalker->R;
    w_clone->update();
    ValueType logpsi2 = psi_clone->evaluateLog(*w_clone); 
    RealType eloc2  = h_clone->evaluate(*w_clone);

    app_log() << "Testing walker-by-walker functions " << endl;
    app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << endl;
    app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << endl;

        
    app_log() << "Testing pbyp functions " << endl;
    Walker_t::Buffer_t &wbuffer(awalker->DataSet);
    wbuffer.clear();
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;
    W.registerData(wbuffer);
    logpsi1=Psi.registerData(W,wbuffer);
    eloc1= H.evaluate(W);
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;

    wbuffer.clear();
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;
    w_clone->registerData(wbuffer);
    logpsi2=psi_clone->registerData(W,wbuffer);
    eloc2= H.evaluate(*w_clone);
    app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;

    app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << endl;
    app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << endl;

    delete h_clone;
    delete psi_clone;
    delete w_clone;
    }
  }

void WaveFunctionTester::runBasicTest() {
  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 0.00001;
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
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);

  app_log() << "  HamTest " << "  Total " <<  eloc << endl;
  for(int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;

//#if defined(QMC_COMPLEX)
//  ValueType logpsi(std::log(psi));
//#else
//  ValueType logpsi(std::log(std::abs(psi)));
//#endif

  for(int iat=0; iat<nat; iat++) {
    PosType r0 = W.R[iat];
    GradType g0;  
    ValueType lap = 0.0;
    for(int idim=0; idim<3; idim++) {
   
      W.R[iat][idim] = r0[idim]+delta;         
      W.update();
      ValueType logpsi_p = Psi.evaluateLog(W);
      
      W.R[iat][idim] = r0[idim]-delta;         
      W.update();
      ValueType logpsi_m = Psi.evaluateLog(W);
      lap += logpsi_m+logpsi_p;
      g0[idim] = logpsi_p-logpsi_m;
//#if defined(QMC_COMPLEX)
//      lap += std::log(psi_m) + std::log(psi_p);
//      g0[idim] = std::log(psi_p)-std::log(psi_m);
//#else
//      lap += std::log(std::abs(psi_m)) + std::log(abs(psi_p));
//      g0[idim] = std::log(std::abs(psi_p/psi_m));
//#endif
      W.R[iat] = r0;
    }
    G1[iat] = c1*g0;
    L1[iat] = c2*(lap-6.0*logpsi);
    cerr << "G1 = " << G1[iat] << endl;
    cerr << "L1 = " << L1[iat] << endl;
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
    //ValueType psi_p = log(fabs(Psi.evaluate(W)));
    RealType psi_p = Psi.evaluateLog(W);
    RealType phase_p=Psi.getPhase();

    W.makeMove(iat,deltaR[iat]);
    RealType aratio = Psi.ratio(W,iat);
    W.rejectMove(iat);
    Psi.rejectMove(iat);

    W.R[iat] += deltaR[iat];         
    W.update();
    //ValueType psi_m = log(fabs(Psi.evaluate(W)));
    RealType psi_m = Psi.evaluateLog(W);
    RealType phase_m=Psi.getPhase();

    RealType ratDiff=std::exp(psi_m-psi_p)*std::cos(phase_m-phase_p) ;
    cout << iat << " ratio " << aratio/ratDiff << " " << ratDiff << endl;
  }
} 

void WaveFunctionTester::runRatioTest() {

  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
  ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);

  bool checkHam=(checkHamPbyP == "yes");
  Tau=0.025;
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  while(it != it_end) {
    makeGaussRandom(deltaR);
    Walker_t::Buffer_t tbuffer;
    W.R = (**it).R+Tau*deltaR;
    (**it).R=W.R;
    //W.registerData(**it,tbuffer);
    W.registerData(tbuffer);
    RealType logpsi=Psi.registerData(W,tbuffer);
    RealType ene;
    if(checkHam) 
      ene = H.registerData(W,tbuffer);
    else
      ene = H.evaluate(W);
    (*it)->DataSet=tbuffer;

    //RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
    H.saveProperty((*it)->getPropertyBase());
    ++it;

    app_log() << "  HamTest " << "  Total " <<  ene << endl;
    for(int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;
  } 

  cout << "  Update using drift " << endl;
  for(int iter=0; iter<4;++iter)
  {
    int iw=0;
    it=W.begin();
    while(it != it_end) {

      cout << "\nStart Walker " << iw++ << endl;
      Walker_t& thisWalker(**it);
      W.R = thisWalker.R;
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);
      H.copyFromBuffer(W,w_buffer);

      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);

      makeGaussRandom(deltaR);

      //mave a move
      RealType ratio_accum(1.0);
      for(int iat=0; iat<nat; iat++) {
        PosType dr(Tau*deltaR[iat]);

        PosType newpos(W.makeMove(iat,dr));

        //RealType ratio=Psi.ratio(W,iat,dGp,dLp);
        RealType ratio=Psi.ratio(W,iat,W.dG,W.dL);
        Gp = W.G + W.dG;
        //RealType enew = H.evaluatePbyP(W,iat);
        if(checkHam) enew = H.evaluatePbyP(W,iat);

        if(ratio > Random()) {
          cout << " Accepting a move for " << iat << endl;
          cout << " Energy after a move " << enew << endl;
          W.G += W.dG;
          W.L += W.dL;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          if(checkHam) H.acceptMove(iat);
          ratio_accum *= ratio;
        } else {
          cout << " Rejecting a move for " << iat << endl;
          W.rejectMove(iat); 
          Psi.rejectMove(iat);
          //H.rejectMove(iat);
        }
      }

      cout << " Energy after pbyp = " << H.getLocalEnergy() << endl;
      thisWalker.R=W.R;
      w_buffer.rewind();
      W.copyToBuffer(w_buffer);
      RealType newlogpsi_up = Psi.evaluateLog(W,w_buffer);
      RealType ene_up;
      if(checkHam)
        ene_up= H.evaluate(W,w_buffer);
      else
        ene_up = H.evaluate(W);

      Gp=W.G;
      Lp=W.L;
      W.R=thisWalker.R;
      W.update();
      RealType newlogpsi=Psi.evaluateLog(W);
      RealType ene = H.evaluate(W);
      thisWalker.resetProperty(newlogpsi,Psi.getPhase(),ene);
      //thisWalker.resetProperty(std::log(psi),Psi.getPhase(),ene);

      cout << iter << "  Energy by update = "<< ene_up << " " << ene << " "  << ene_up-ene << endl;
      cout << iter << " Ratio " << ratio_accum*ratio_accum 
        << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " " 
        << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << endl
        << " new log(psi) updated " << newlogpsi_up
        << " new log(psi) calculated " << newlogpsi  
        << " old log(psi) " << logpsi << endl;

      cout << " Gradients " << endl;
      for(int iat=0; iat<nat; iat++) 
        cout << W.G[iat]-Gp[iat] << W.G[iat] << endl; //W.G[iat] << G[iat] << endl;
      cout << " Laplacians " << endl;
      for(int iat=0; iat<nat; iat++) 
        cout << W.L[iat]-Lp[iat] << " " << W.L[iat] << endl;

      ++it;
    }
  }

  cout << "  Update without drift : for VMC useDrift=\"no\"" << endl;
  for(int iter=0; iter<4;++iter)
  {
    it=W.begin();
    int iw=0;
    while(it != it_end) {

      cout << "\nStart Walker " << iw++ << endl;
      Walker_t& thisWalker(**it);
      W.R = thisWalker.R;
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);

      //mave a move
      RealType ratio_accum(1.0);

      for(int substep=0; substep<3; ++substep)
      {
        makeGaussRandom(deltaR);

        for(int iat=0; iat<nat; iat++) {
          PosType dr(Tau*deltaR[iat]);

          PosType newpos(W.makeMove(iat,dr));

          RealType ratio=Psi.ratio(W,iat);
          RealType prob = ratio*ratio;
          if(prob > Random()) {
            cout << " Accepting a move for " << iat << endl;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            ratio_accum *= ratio;
          } else {
            cout << " Rejecting a move for " << iat << endl;
            W.rejectMove(iat); 
            Psi.rejectMove(iat);
          }
        }

        thisWalker.R=W.R;
        w_buffer.rewind();
        W.updateBuffer(w_buffer);
        RealType logpsi_up = Psi.updateBuffer(W,w_buffer,false);
        RealType ene = H.evaluate(W);
        thisWalker.resetProperty(logpsi_up,Psi.getPhase(),ene);
      }

      Gp=W.G;
      Lp=W.L;

      W.update();
      RealType newlogpsi=Psi.evaluateLog(W);
      cout << iter << " Ratio " << ratio_accum*ratio_accum 
        << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " " 
        << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << endl
        << " new log(psi) " << newlogpsi 
        << " old log(psi) " << logpsi << endl;

      cout << " Gradients " << endl;
      for(int iat=0; iat<nat; iat++) {
        cout << W.G[iat]-Gp[iat] << W.G[iat] << endl; //W.G[iat] << G[iat] << endl;
      }
      cout << " Laplacians " << endl;
      for(int iat=0; iat<nat; iat++) {
        cout << W.L[iat]-Lp[iat] << " " << W.L[iat] << endl;
      }
      ++it;
    }
  }

  //for(it=W.begin();it != it_end; ++it)
  //{
  //  Walker_t& thisWalker(**it);
  //  Walker_t::Buffer_t& w_buffer((*it)->DataSet);
  //  w_buffer.rewind();
  //  W.updateBuffer(**it,w_buffer);
  //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
  //}


}

void WaveFunctionTester::runRatioTest2() 
{

  int nat = W.getTotalNum();
  ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
  ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);

  Tau=0.025;
  MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
  for(; it != it_end; ++it)
  {
    makeGaussRandom(deltaR);
    Walker_t::Buffer_t tbuffer;
    W.R = (**it).R+Tau*deltaR;
    (**it).R=W.R;
    //W.registerData(**it,tbuffer);
    W.registerData(tbuffer);
    RealType logpsi=Psi.registerData(W,tbuffer);
    RealType ene = H.evaluate(W);
    (*it)->DataSet=tbuffer;

    //RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
    H.saveProperty((*it)->getPropertyBase());

    app_log() << "  HamTest " << "  Total " <<  ene << endl;
    for(int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;
  } 

  for(int iter=0; iter<20;++iter)
  {
    int iw=0;
    it=W.begin();
    //while(it != it_end) 
    for(; it != it_end; ++it)
    {
      cout << "\nStart Walker " << iw++ << endl;
      Walker_t& thisWalker(**it);
      W.R = thisWalker.R;
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);

      RealType eold(thisWalker.Properties(LOCALENERGY));
      RealType logpsi(thisWalker.Properties(LOGPSI));
      RealType emixed(eold), enew(eold);

      makeGaussRandom(deltaR);

      //mave a move
      RealType ratio_accum(1.0);
      for(int iat=0; iat<nat; iat++) 
      {
        GradType grad_now=Psi.evalGrad(W,iat), grad_new;
        PosType dr(Tau*deltaR[iat]);
        PosType newpos(W.makeMove(iat,dr));

        RealType ratio2 = Psi.ratioGrad(W,iat,grad_new);
        W.rejectMove(iat); 
        Psi.rejectMove(iat);

        newpos=W.makeMove(iat,dr);
        RealType ratio1 = Psi.ratio(W,iat);
        W.rejectMove(iat); 
        cout << " ratio1 = " << ratio1 << " ration2 = " << ratio2 << endl;
      }
    }
  }

  //for(it=W.begin();it != it_end; ++it)
  //{
  //  Walker_t& thisWalker(**it);
  //  Walker_t::Buffer_t& w_buffer((*it)->DataSet);
  //  w_buffer.rewind();
  //  W.updateBuffer(**it,w_buffer);
  //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
  //}


}


void WaveFunctionTester::runGradSourceTest() {
  ParticleSetPool::PoolType::iterator p;
  for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
    app_log() << "ParticelSet = " << p->first << endl;

  // Find source ParticleSet
  ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
  app_log() << pit->first << endl;
  // if(pit == PtclPool.getPool().end()) 
  //   APP_ABORT("Unknown source \"" + sourceName + "\" WaveFunctionTester.");
  
  ParticleSet& source = *((*pit).second);

  IndexType nskipped = 0;
  RealType sig2Enloc=0, sig2Drift=0;
  RealType delta = 0.00001;
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
  //ValueType psi = Psi.evaluate(W);
  ValueType logpsi = Psi.evaluateLog(W);
  RealType eloc=H.evaluate(W);

  app_log() << "  HamTest " << "  Total " <<  eloc << endl;
  for(int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

  //RealType psi = Psi.evaluateLog(W);
  ParticleSet::ParticleGradient_t G(nat), G1(nat);
  ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
  G = W.G;
  L = W.L;

  for (int isrc=0; isrc < 1/*source.getTotalNum()*/; isrc++) {
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad;
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad_FD;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_FD;
    for (int dim=0; dim<OHMMS_DIM; dim++) {
      grad_grad[dim].resize(nat);    lapl_grad[dim].resize(nat);
      grad_grad_FD[dim].resize(nat); lapl_grad_FD[dim].resize(nat);
    }
    Psi.evaluateLog(W);
    GradType grad_log = Psi.evalGradSource (W, source, isrc, grad_grad, lapl_grad);
    ValueType log = Psi.evaluateLog(W);
    //grad_log = Psi.evalGradSource (W, source, isrc);

    for(int iat=0; iat<nat; iat++) {
      PosType r0 = W.R[iat];
      GradType gFD[OHMMS_DIM];  
      GradType lapFD = 0.0;
      for(int eldim=0; eldim<3; eldim++) {
	W.R[iat][eldim] = r0[eldim]+delta;         
	W.update();
	//ValueType log_p = Psi.evaluateLog(W);
	GradType gradlogpsi_p =  Psi.evalGradSource(W, source, isrc);
	W.R[iat][eldim] = r0[eldim]-delta;         
	W.update();
	//ValueType log_m = Psi.evaluateLog(W);
	GradType gradlogpsi_m = Psi.evalGradSource(W, source, isrc);
	lapFD    += gradlogpsi_m + gradlogpsi_p;
	gFD[eldim] = gradlogpsi_p - gradlogpsi_m;
	W.R[iat] = r0;
	W.update();
	//Psi.evaluateLog(W);
      }
      for (int iondim=0; iondim<OHMMS_DIM; iondim++) {
	for (int eldim=0; eldim<OHMMS_DIM; eldim++)
	  grad_grad_FD[iondim][iat][eldim] = c1*gFD[eldim][iondim];
	lapl_grad_FD[iondim][iat] = c2*(lapFD[iondim]-6.0*grad_log[iondim]);
      }
    }
    cout.precision(15);
    for (int dimsrc=0; dimsrc<OHMMS_DIM; dimsrc++) {
      for(int iat=0; iat<nat; iat++) {
	cout.precision(15);
	cout << "For particle #" << iat << " at " << W.R[iat] << endl;
	cout << "Gradient      = " << setw(12) << grad_grad[dimsrc][iat] << endl 
	     << "  Finite diff = " << setw(12) << grad_grad_FD[dimsrc][iat] << endl 
	     << "  Error       = " << setw(12) 
	     <<  grad_grad_FD[dimsrc][iat] - grad_grad[dimsrc][iat] << endl << endl;
	// cout << "Laplacian     = " << setw(12) << lapl_grad[dimsrc][iat] << endl 
	//      << "  Finite diff = " << setw(12) << lapl_grad_FD[dimsrc][iat] << endl 
	//      << "  Error       = " << setw(12) 
	//      << lapl_grad_FD[dimsrc][iat] - lapl_grad[dimsrc][iat] << endl << endl;
      }
    }
  }
} 



bool 
WaveFunctionTester::put(xmlNodePtr q){
  return putQMCInfo(q);
}
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/QMCUpdateBase.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/// Constructor.
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, RandomGenerator_t& rg)
  : W(w),Psi(psi),Guide(guide),H(h), RandomGen(rg), branchEngine(0), Estimators(0), Traces(0), csoffset(0)
{
  setDefaults();
}

/// Constructor.
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
  : W(w),Psi(psi),H(h),Guide(psi), RandomGen(rg), branchEngine(0), Estimators(0), Traces(0), csoffset(0)
{
  setDefaults();
}

///copy constructor
QMCUpdateBase::QMCUpdateBase(const QMCUpdateBase& a)
  : W(a.W), Psi(a.Psi), Guide(a.Guide), H(a.H), RandomGen(a.RandomGen)
  , branchEngine(0), Estimators(0), Traces(0)
{
  APP_ABORT("QMCUpdateBase::QMCUpdateBase(const QMCUpdateBase& a) Not Allowed");
}

/// destructor
QMCUpdateBase::~QMCUpdateBase()
{
}

void QMCUpdateBase::setDefaults()
{
  UpdatePbyP=true;
  UseTMove=false;
  NumPtcl=0;
  nSubSteps=1;
  MaxAge=10;
  m_r2max=-1;
  myParams.add(m_r2max,"maxDisplSq","double"); //maximum displacement
  //store 1/mass per species
  SpeciesSet tspecies(W.getSpeciesSet());
  int massind=tspecies.addAttribute("mass");
  MassInvS.resize(tspecies.getTotalNum());
  for(int ig=0; ig<tspecies.getTotalNum(); ++ig)
    MassInvS[ig]=1.0/tspecies(massind,ig);
  MassInvP.resize(W.getTotalNum());
  for(int ig=0; ig<W.groups(); ++ig)
  {
    for(int iat=W.first(ig); iat<W.last(ig); ++iat)
      MassInvP[iat]=MassInvS[ig];
  }
}

bool QMCUpdateBase::put(xmlNodePtr cur)
{
  //nonlocal operator is very light
  UseTMove = nonLocalOps.put(cur);
  bool s=myParams.put(cur);
  if (branchEngine)
    branchEngine->put(cur);
  return s;
}

void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManager* est)
{
  Estimators=est;
  branchEngine=brancher;
  branchEngine->setEstimatorManager(est);
  NumPtcl = W.getTotalNum();
  deltaR.resize(NumPtcl);
  drift.resize(NumPtcl);
  G.resize(NumPtcl);
  dG.resize(NumPtcl);
  L.resize(NumPtcl);
  dL.resize(NumPtcl);
  //set the default tau-mass related values with electrons
  Tau=brancher->getTau();
  m_tauovermass = Tau*MassInvS[0];
  m_oneover2tau = 0.5/(m_tauovermass);
  m_sqrttau = std::sqrt(m_tauovermass);
  if(!UpdatePbyP)
  {
    // store sqrt(tau/mass)
    SqrtTauOverMass.resize(W.getTotalNum());
    for(int iat=0; iat<W.getTotalNum(); ++iat)
      SqrtTauOverMass[iat]=std::sqrt(Tau*MassInvP[iat]);
  }
  //app_log() << "  QMCUpdateBase::resetRun m/tau=" << m_tauovermass << endl;
  if (m_r2max<0)
    m_r2max =  W.Lattice.LR_rc* W.Lattice.LR_rc;
  //app_log() << "  Setting the bound for the displacement max(r^2) = " <<  m_r2max << endl;
}


void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManager* est, TraceManager* traces)
{
  resetRun(brancher,est);
  Traces = traces;
}


void QMCUpdateBase::resetEtrial(RealType et)
{
  //branchEngine->E_T=et;
  branchEngine->setTrialEnergy(et,1.0);
  branchEngine->flush(0);
}

void QMCUpdateBase::startRun(int blocks, bool record)
{
  Estimators->start(blocks,record);
  H.initialize_traces(*Traces,W);
  Estimators->initialize_traces(*Traces);
  Traces->initialize_traces();
}

void QMCUpdateBase::stopRun()
{
  Estimators->stop();
}

//ugly, but will use until general usage of stopRun is clear
//  DMCOMP and VMCSingleOMP do not use stopRun anymore
void QMCUpdateBase::stopRun2()
{
  H.finalize_traces();
  Estimators->finalize_traces();
  Traces->finalize_traces();
}

void QMCUpdateBase::startBlock(int steps)
{
  Estimators->startBlock(steps);
  Traces->startBlock(steps);
  nAccept = 0;
  nReject=0;
  nAllRejected=0;
  nNodeCrossing=0;
  NonLocalMoveAccepted=0;
}

void QMCUpdateBase::stopBlock(bool collectall)
{
  Estimators->stopBlock(acceptRatio(),collectall);
  Traces->stopBlock();
}

void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=false;
  //ignore different mass
  //RealType tauovermass = Tau*MassInv[0];
  for (; it != it_end; ++it)
  {
    W.R = (*it)->R;
    W.update();
    RealType logpsi(Psi.evaluateLog(W));
    (*it)->G=W.G;
    (*it)->L=W.L;
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    RealType ene = H.evaluate(W);
    H.auxHevaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0, nodecorr);
    (*it)->Weight=1;
    H.saveProperty((*it)->getPropertyBase());
  }
}


void QMCUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=true;
  for (; it != it_end; ++it)
  {
    Walker_t& awalker(**it);
    W.R=awalker.R;
    W.update();
    //W.loadWalker(awalker,UpdatePbyP);
    if (awalker.DataSet.size())
      awalker.DataSet.clear();
    awalker.DataSet.rewind();
    RealType logpsi=Psi.registerData(W,awalker.DataSet);
    RealType logpsi2=Psi.updateBuffer(W,awalker.DataSet,false);
    awalker.G=W.G;
    awalker.L=W.L;
    randomize(awalker);
  }
}

/** randomize a walker with a diffusion MC using gradients */
void QMCUpdateBase::randomize(Walker_t& awalker)
{
  BadState=false;
  //Walker_t::Buffer_t& w_buffer(awalker.DataSet);
  //W.loadWalker(awalker,true);
  //Psi.copyFromBuffer(W,w_buffer);
  RealType eloc_tot=0.0;
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR,RandomGen);
  for(int ig=0; ig<W.groups(); ++ig) //loop over species
  {
    RealType tauovermass = Tau*MassInvS[ig];
    RealType oneover2tau = 0.5/(tauovermass);
    RealType sqrttau = std::sqrt(tauovermass);
    for (int iat=W.first(ig); iat<W.last(ig); ++iat)
    {
      GradType grad_now=Psi.evalGrad(W,iat), grad_new;
      PosType dr;
      getScaledDrift(tauovermass,grad_now,dr);
      dr += sqrttau*deltaR[iat];
      if (!W.makeMoveAndCheck(iat,dr))
        continue;
      //PosType newpos = W.makeMove(iat,dr);
      RealType ratio = Psi.ratioGrad(W,iat,grad_new);
      RealType prob = ratio*ratio;
      //zero is always rejected
      if (prob<numeric_limits<RealType>::epsilon())
      {
        ++nReject;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
        continue;
      }
      RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
      getScaledDrift(tauovermass,grad_new,dr);
      dr = awalker.R[iat]-W.R[iat]-dr;
      RealType logGb = -oneover2tau*dot(dr,dr);
      if (RandomGen() < prob*std::exp(logGb-logGf))
      {
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
      }
      else
      {
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
    }
  }
  //for subSteps must update thiswalker
  awalker.R=W.R;
  awalker.G=W.G;
  awalker.L=W.L;
  RealType logpsi = Psi.updateBuffer(W,awalker.DataSet,false);
  W.saveWalker(awalker);
  RealType eloc=H.evaluate(W);
  BadState |= isnan(eloc);
  //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
  awalker.resetProperty(logpsi,Psi.getPhase(), eloc);
  H.auxHevaluate(W,awalker);
  H.saveProperty(awalker.getPropertyBase());
  awalker.ReleasedNodeAge=0;
  awalker.ReleasedNodeWeight=0;
  awalker.Weight=1;
}

QMCUpdateBase::RealType
QMCUpdateBase::getNodeCorrection(const ParticleSet::ParticleGradient_t& g, ParticleSet::ParticlePos_t& gscaled)
{
  //setScaledDrift(m_tauovermass,g,gscaled);
  //RealType vsq=Dot(g,g);
  //RealType x=m_tauovermass*vsq;
  //return (vsq<numeric_limits<RealType>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
  return  setScaledDriftPbyPandNodeCorr(Tau,MassInvP,g,gscaled);
}

void QMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,UpdatePbyP);
    Walker_t::Buffer_t& w_buffer((*it)->DataSet);
    RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
    W.saveWalker(thisWalker);
  }
}

void QMCUpdateBase::setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    RealType M=std::abs((*it)->Weight);
    (*it)->Multiplicity = std::floor(M + RandomGen());
  }
}

void QMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    RealType M=(*it)->Weight;
    if ((*it)->Age>MaxAge)
      M = std::min(0.5,M);
    else
      if ((*it)->Age > 0)
        M = std::min(1.0,M);
    (*it)->Multiplicity = M + RandomGen();
  }
}

void QMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip)
{
  char fname[16];
  sprintf(fname,"test.%i",ip);
  ofstream fout(fname,ios::app);
  int i=0;
  fout << "benchMark started." << endl;
  for (; it != it_end; ++it,++i)
  {
    Walker_t& thisWalker(**it);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    W.R = m_sqrttau*deltaR+ thisWalker.R;
    W.update();
    ValueType logpsi(Psi.evaluateLog(W));
    RealType e = H.evaluate(W);
    fout << W.R[0] << W.G[0] << endl;
    fout <<  i << " " << logpsi << " " << e << endl;
  }
  fout << "benchMark completed." << endl;
}

/** advance a walker: walker move, use drift and vmc
 */
void QMCUpdateBase::advanceWalker(Walker_t& thisWalker)
{
  W.loadWalker(thisWalker,false);
  RealType logpsi0=thisWalker.Properties(LOGPSI);
  RealType phase0=thisWalker.Properties(SIGN);
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    if (!W.makeMoveWithDrift(thisWalker,drift ,deltaR, m_sqrttau))
    {
      ++nReject;
      continue;
    }
    RealType logpsi(Psi.evaluateLog(W));
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    // setScaledDrift(m_tauovermass,W.G,drift);
    nodecorr=setScaledDriftPbyPandNodeCorr(Tau,MassInvP,W.G,drift);
    //backward GreenFunction needs \f$d{\bf R} = {\bf R}_{old} - {\bf R}_{new} - {\bf V}_d\f$
    deltaR = thisWalker.R - W.R - drift;
    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
    RealType g= std::exp(logGb-logGf+2.0*(logpsi-thisWalker.Properties(LOGPSI)));
    if (RandomGen() > g)
    {
      thisWalker.Age++;
      ++nReject;
    }
    else
    {
      thisWalker.R=W.R;
      thisWalker.G=W.G;
      thisWalker.L=W.L;
      //skip energy
      thisWalker.resetProperty(logpsi,Psi.getPhase(),0);
      //update logpsi0, phase0
      logpsi0=logpsi;
      phase0=Psi.getPhase();
      ++nAccept;
    }
  }
  //measure energy
  W.loadWalker(thisWalker,true);
  RealType eloc=H.evaluate(W);
  thisWalker.resetProperty(logpsi0,phase0,eloc);
  H.auxHevaluate(W,thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());
}

/** advance of a walker using VMC+drift */
void QMCUpdateBase::advancePbyP(Walker_t& thisWalker)
{
  Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
  W.loadWalker(thisWalker,true);
  Psi.copyFromBuffer(W,w_buffer);
  bool moved = false;
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    for(int ig=0; ig<W.groups(); ++ig) //loop over species
    {
      RealType tauovermass = Tau*MassInvS[ig];
      RealType oneover2tau = 0.5/(tauovermass);
      RealType sqrttau = std::sqrt(tauovermass);
      for (int iat=W.first(ig); iat<W.last(ig); ++iat)
      {
        GradType grad_now=Psi.evalGrad(W,iat), grad_new;
        PosType dr;
        getScaledDrift(tauovermass,grad_now,dr);
        dr += sqrttau*deltaR[iat];
        if (!W.makeMoveAndCheck(iat,dr))
        {
          ++nReject;
          continue;
        }
        //PosType newpos = W.makeMove(iat,dr);
        RealType ratio = Psi.ratioGrad(W,iat,grad_new);
        RealType prob = ratio*ratio;
        //zero is always rejected
        if (prob<numeric_limits<RealType>::epsilon())
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
          continue;
        }
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
        getScaledDrift(tauovermass,grad_new,dr);
        dr = thisWalker.R[iat]-W.R[iat]-dr;
        RealType logGb = -oneover2tau*dot(dr,dr);
        //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
        if (RandomGen() < prob*std::exp(logGb-logGf))
        {
          moved = true;
          ++nAccept;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
        }
        else
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
    //for subSteps must update thiswalker
    thisWalker.R=W.R;
    thisWalker.G=W.G;
    thisWalker.L=W.L;
  }
  //Always compute the energy
  {
    RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
    W.saveWalker(thisWalker);
    RealType eloc=H.evaluate(W);
    //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
    thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
    H.auxHevaluate(W,thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());
  }
  if(!moved)
    ++nAllRejected;
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1618 $   $Date: 2007-01-14 18:10:10 -0600 (Sun, 14 Jan 2007) $
 * $Id: QMCUpdateBase.cpp 1618 2007-01-15 00:10:10Z jnkim $
 ***************************************************************************/

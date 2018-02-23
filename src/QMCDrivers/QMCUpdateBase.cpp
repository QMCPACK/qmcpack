//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Platforms/sysutil.h"
#include "QMCDrivers/QMCUpdateBase.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

/// Constructor.
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, RandomGenerator_t& rg)
  : W(w), Psi(psi), Guide(guide), H(h), nonLocalOps(w.getTotalNum()), RandomGen(rg), branchEngine(0), Estimators(0), Traces(0), csoffset(0)
{
  setDefaults();
}

/// Constructor.
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
  : W(w), Psi(psi), H(h), nonLocalOps(w.getTotalNum()), Guide(psi), RandomGen(rg), branchEngine(0), Estimators(0), Traces(0), csoffset(0)
{
  setDefaults();
}

///copy constructor
QMCUpdateBase::QMCUpdateBase(const QMCUpdateBase& a)
  : W(a.W), Psi(a.Psi), Guide(a.Guide), H(a.H), nonLocalOps(a.W.getTotalNum()), RandomGen(a.RandomGen)
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
  UseDrift=true;
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

  InitWalkersTimer = TimerManager.createTimer("QMCUpdateBase::WalkerInit", timer_level_medium);
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

void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManagerBase* est)
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
  //app_log() << "  QMCUpdateBase::resetRun m/tau=" << m_tauovermass << std::endl;
  if (m_r2max<0)
    m_r2max =  W.Lattice.LR_rc* W.Lattice.LR_rc;
  //app_log() << "  Setting the bound for the displacement std::max(r^2) = " <<  m_r2max << std::endl;
}


void QMCUpdateBase::resetRun(BranchEngineType* brancher, EstimatorManagerBase* est, TraceManager* traces)
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
#if !defined(REMOVE_TRACEMANAGER)
  if(!Traces)
  {
    APP_ABORT("QMCUpdateBase::startRun\n  derived QMCDriver class has not setup trace clones properly\n  null TraceManager pointer encountered in derived QMCUpdateBase class\n  see VMCLinearOptOMP.cpp for a correct minimal interface (search on 'trace')\n  refer to changes made in SVN revision 6597 for further guidance");
  }
  H.initialize_traces(*Traces,W);
  Traces->initialize_traces();
#endif
}

void QMCUpdateBase::stopRun()
{
  Estimators->stop();
}

//ugly, but will use until general usage of stopRun is clear
//  DMCOMP and VMCSingleOMP do not use stopRun anymore
void QMCUpdateBase::stopRun2()
{
#if !defined(REMOVE_TRACEMANAGER)
  H.finalize_traces();
  Traces->finalize_traces();
#endif
}

void QMCUpdateBase::startBlock(int steps)
{
  Estimators->startBlock(steps);
#if !defined(REMOVE_TRACEMANAGER)
  Traces->startBlock(steps);
#endif
  nAccept = 0;
  nReject=0;
  nAllRejected=0;
  nNodeCrossing=0;
  NonLocalMoveAccepted=0;
}

void QMCUpdateBase::stopBlock(bool collectall)
{
  Estimators->stopBlock(acceptRatio(),collectall);
#if !defined(REMOVE_TRACEMANAGER)
  Traces->stopBlock();
#endif
}

void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=false;
  InitWalkersTimer->start();
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
    // cannot call auxHevalate() here because walkers are not initialized
    // for example, DensityEstimator needs the weights of the walkers
    //H.auxHevaluate(W);
    (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0, nodecorr);
    (*it)->Weight=1;
    H.saveProperty((*it)->getPropertyBase());
  }
  InitWalkersTimer->stop();
}

void QMCUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=true;
  InitWalkersTimer->start();
  for (; it != it_end; ++it)
  {
    Walker_t& awalker(**it);
    W.R=awalker.R;
    W.update(true);
    //W.loadWalker(awalker,UpdatePbyP);
    if (awalker.DataSet.size())
      awalker.DataSet.clear();
    awalker.DataSet.rewind();
    awalker.registerData();
    Psi.registerData(W,awalker.DataSet);
    awalker.DataSet.allocate();
    Psi.copyFromBuffer(W,awalker.DataSet);
    Psi.evaluateLog(W);
    RealType logpsi=Psi.updateBuffer(W,awalker.DataSet,false);
    awalker.G=W.G;
    awalker.L=W.L;
    randomize(awalker);
  }
  InitWalkersTimer->stop();
  #pragma omp master
  print_mem("Memory Usage after the buffer registration", app_log());
}

/** randomize a walker with a diffusion MC using gradients */
void QMCUpdateBase::randomize(Walker_t& awalker)
{
  BadState=false;
  //Walker_t::WFBuffer_t& w_buffer(awalker.DataSet);
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
      W.setActive(iat);
      GradType grad_now=Psi.evalGrad(W,iat), grad_new;
      mPosType dr;
      getScaledDrift(tauovermass,grad_now,dr);
      dr += sqrttau*deltaR[iat];
      if (!W.makeMoveAndCheck(iat,dr))
        continue;
      //PosType newpos = W.makeMove(iat,dr);
      RealType ratio = Psi.ratioGrad(W,iat,grad_new);
      RealType prob = ratio*ratio;
      //zero is always rejected
      if (prob<std::numeric_limits<RealType>::epsilon())
      {
        ++nReject;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
        continue;
      }
      RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
      getScaledDrift(tauovermass,grad_new,dr);
      dr = W.R[iat]-W.activePos-dr;
      RealType logGb = -oneover2tau*dot(dr,dr);
      if (RandomGen() < prob*std::exp(logGb-logGf))
      {
        Psi.acceptMove(W,iat);
        W.acceptMove(iat);
      }
      else
      {
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
    }
  }

  W.donePbyP();
  //for subSteps must update thiswalker
  //awalker.R=W.R;
  //awalker.G=W.G;
  //awalker.L=W.L;
  RealType logpsi = Psi.updateBuffer(W,awalker.DataSet,false);
  W.saveWalker(awalker);

  RealType eloc=H.evaluate(W);
#if (__cplusplus >= 201103L)
  BadState |= std::isnan(eloc);
#else
  BadState |= isnan(eloc);
#endif
  //thisWalker.resetProperty(std::log(std::abs(psi)), psi,eloc);
  awalker.resetProperty(logpsi,Psi.getPhase(), eloc);
  H.auxHevaluate(W,awalker);
  H.saveProperty(awalker.getPropertyBase());
  awalker.ReleasedNodeAge=0;
  awalker.ReleasedNodeWeight=0;
  awalker.Weight=1;

//  printf("energy  %.3f\n",eloc);
//  for(size_t i=0, n=W.getTotalNum(); i<n; ++i)
//    printf("evalGrad  %.3f %.3f %.3f %.3f\n",W.G[i][0],W.G[i][1],W.G[i][2],W.L[i]);
//
}

QMCUpdateBase::RealType
QMCUpdateBase::getNodeCorrection(const ParticleSet::ParticleGradient_t& g, ParticleSet::ParticlePos_t& gscaled)
{
  //setScaledDrift(m_tauovermass,g,gscaled);
  //RealType vsq=Dot(g,g);
  //RealType x=m_tauovermass*vsq;
  //return (vsq<std::numeric_limits<RealType>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
  return  setScaledDriftPbyPandNodeCorr(Tau,MassInvP,g,gscaled);
}

void QMCUpdateBase::updateWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,UpdatePbyP);
    //recompute distance tables
    W.update();
    Walker_t::WFBuffer_t& w_buffer((*it)->DataSet);
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
      M = std::min((RealType)0.5,M);
    else
      if ((*it)->Age > 0)
        M = std::min((RealType)1.0,M);
    (*it)->Multiplicity = M + RandomGen();
  }
}

void QMCUpdateBase::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool recompute)
{
  for (; it != it_end; ++it)
  {
    advanceWalker(**it,recompute);
  }
}

void QMCUpdateBase::benchMark(WalkerIter_t it, WalkerIter_t it_end, int ip)
{
  char fname[16];
  sprintf(fname,"test.%i",ip);
  std::ofstream fout(fname,std::ios::app);
  int i=0;
  fout << "benchMark started." << std::endl;
  for (; it != it_end; ++it,++i)
  {
    Walker_t& thisWalker(**it);
    makeGaussRandomWithEngine(deltaR,RandomGen);
    W.R = m_sqrttau*deltaR+ thisWalker.R;
    W.update();
    ValueType logpsi(Psi.evaluateLog(W));
    RealType e = H.evaluate(W);
    fout << W.R[0] << W.G[0] << std::endl;
    fout <<  i << " " << logpsi << " " << e << std::endl;
  }
  fout << "benchMark completed." << std::endl;
}

}


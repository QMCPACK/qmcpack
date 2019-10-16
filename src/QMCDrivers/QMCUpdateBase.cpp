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
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w,
                             TrialWaveFunction& psi,
                             TrialWaveFunction& guide,
                             QMCHamiltonian& h,
                             RandomGenerator_t& rg)
    : csoffset(0),
      Traces(0),
      W(w),
      Psi(psi),
      Guide(guide),
      H(h),
      RandomGen(rg),
      branchEngine(0),
      DriftModifier(0),
      Estimators(0)
{
  setDefaults();
}

/// Constructor.
QMCUpdateBase::QMCUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
    : csoffset(0),
      Traces(0),
      W(w),
      Psi(psi),
      Guide(psi),
      H(h),
      RandomGen(rg),
      branchEngine(0),
      DriftModifier(0),
      Estimators(0)
{
  setDefaults();
}

/// destructor
QMCUpdateBase::~QMCUpdateBase() {}

void QMCUpdateBase::setDefaults()
{
  UpdatePbyP = true;
  UseDrift   = true;
  NumPtcl    = 0;
  nSubSteps  = 1;
  MaxAge     = 10;
  m_r2max    = -1;
  myParams.add(m_r2max, "maxDisplSq", "double"); //maximum displacement
  //store 1/mass per species
  SpeciesSet tspecies(W.getSpeciesSet());
  int massind = tspecies.addAttribute("mass");
  MassInvS.resize(tspecies.getTotalNum());
  for (int ig = 0; ig < tspecies.getTotalNum(); ++ig)
    MassInvS[ig] = 1.0 / tspecies(massind, ig);
  MassInvP.resize(W.getTotalNum());
  for (int ig = 0; ig < W.groups(); ++ig)
  {
    for (int iat = W.first(ig); iat < W.last(ig); ++iat)
      MassInvP[iat] = MassInvS[ig];
  }

  InitWalkersTimer = TimerManager.createTimer("QMCUpdateBase::WalkerInit", timer_level_medium);
}

bool QMCUpdateBase::put(xmlNodePtr cur)
{
  H.setNonLocalMoves(cur);
  bool s = myParams.put(cur);
  return s;
}

void QMCUpdateBase::resetRun(BranchEngineType* brancher,
                             EstimatorManagerBase* est,
                             TraceManager* traces,
                             const DriftModifierBase* driftmodifer)
{
  Estimators    = est;
  branchEngine  = brancher;
  DriftModifier = driftmodifer;
  Traces        = traces;

  NumPtcl = W.getTotalNum();
  deltaR.resize(NumPtcl);
  drift.resize(NumPtcl);
  G.resize(NumPtcl);
  dG.resize(NumPtcl);
  L.resize(NumPtcl);
  dL.resize(NumPtcl);
  //set the default tau-mass related values with electrons
  Tau           = brancher->getTau();
  m_tauovermass = Tau * MassInvS[0];
  m_oneover2tau = 0.5 / (m_tauovermass);
  m_sqrttau     = std::sqrt(m_tauovermass);
  if (!UpdatePbyP)
  {
    // store sqrt(tau/mass)
    SqrtTauOverMass.resize(W.getTotalNum());
    for (int iat = 0; iat < W.getTotalNum(); ++iat)
      SqrtTauOverMass[iat] = std::sqrt(Tau * MassInvP[iat]);
  }
  //app_log() << "  QMCUpdateBase::resetRun m/tau=" << m_tauovermass << std::endl;
  if (m_r2max < 0)
    m_r2max = W.Lattice.LR_rc * W.Lattice.LR_rc;
  //app_log() << "  Setting the bound for the displacement std::max(r^2) = " <<  m_r2max << std::endl;
}

void QMCUpdateBase::startRun(int blocks, bool record)
{
  Estimators->start(blocks, record);
#if !defined(REMOVE_TRACEMANAGER)
  if (!Traces)
  {
    APP_ABORT(
        "QMCUpdateBase::startRun\n  derived QMCDriver class has not setup trace clones properly\n  null TraceManager "
        "pointer encountered in derived QMCUpdateBase class\n  see VMCLinearOptOMP.cpp for a correct minimal interface "
        "(search on 'trace')\n  refer to changes made in SVN revision 6597 for further guidance");
  }
  H.initialize_traces(*Traces, W);
  Traces->initialize_traces();
#endif
}

void QMCUpdateBase::stopRun() { Estimators->stop(); }

//ugly, but will use until general usage of stopRun is clear
//  DMC and VMC do not use stopRun anymore
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
  nAccept              = 0;
  nReject              = 0;
  nAllRejected         = 0;
  nNodeCrossing        = 0;
  NonLocalMoveAccepted = 0;
}

void QMCUpdateBase::stopBlock(bool collectall)
{
  Estimators->stopBlock(acceptRatio(), collectall);
#if !defined(REMOVE_TRACEMANAGER)
  Traces->stopBlock();
#endif
}

void QMCUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP = false;
  InitWalkersTimer->start();
  //ignore different mass
  //RealType tauovermass = Tau*MassInv[0];
  for (; it != it_end; ++it)
  {
    W.R = (*it)->R;
    W.update();
    RealType logpsi(Psi.evaluateLog(W));
    (*it)->G          = W.G;
    (*it)->L          = W.L;
    RealType nodecorr = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
    RealType ene      = H.evaluate(W);
    // cannot call auxHevalate() here because walkers are not initialized
    // for example, DensityEstimator needs the weights of the walkers
    //H.auxHevaluate(W);
    (*it)->resetProperty(logpsi, Psi.getPhase(), ene, 0.0, 0.0, nodecorr);
    (*it)->Weight = 1;
    H.saveProperty((*it)->getPropertyBase());
  }
  InitWalkersTimer->stop();
}

void QMCUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP = true;
  BadState   = false;
  InitWalkersTimer->start();
  if (it == it_end)
  {
    // a particular case, no walker enters in this call.
    // but need to free the memory of Psi.
    Walker_t dummy_walker(W.getTotalNum());
    dummy_walker.Properties = W.Properties;
    dummy_walker.registerData();
    Psi.registerData(W, dummy_walker.DataSet);
  }
  for (; it != it_end; ++it)
  {
    Walker_t& awalker(**it);
    W.R = awalker.R;
    W.update();
    if (awalker.DataSet.size())
      awalker.DataSet.clear();
    awalker.DataSet.rewind();
    awalker.registerData();
    Psi.registerData(W, awalker.DataSet);
    awalker.DataSet.allocate();
    // This from here on should happen in the scope of the block
    Psi.copyFromBuffer(W, awalker.DataSet);
    Psi.evaluateLog(W);
    RealType logpsi = Psi.updateBuffer(W, awalker.DataSet, false);
    W.saveWalker(awalker);
    RealType eloc = H.evaluate(W);
    BadState |= std::isnan(eloc);
    awalker.resetProperty(logpsi, Psi.getPhase(), eloc);
    H.auxHevaluate(W, awalker);
    H.saveProperty(awalker.getPropertyBase());
    awalker.ReleasedNodeAge    = 0;
    awalker.ReleasedNodeWeight = 0;
    awalker.Weight             = 1;
  }
  InitWalkersTimer->stop();
#pragma omp master
  print_mem("Memory Usage after the buffer registration", app_log());
}

QMCUpdateBase::RealType QMCUpdateBase::getNodeCorrection(const ParticleSet::ParticleGradient_t& g,
                                                         ParticleSet::ParticlePos_t& gscaled)
{
  //setScaledDrift(m_tauovermass,g,gscaled);
  //RealType vsq=Dot(g,g);
  //RealType x=m_tauovermass*vsq;
  //return (vsq<std::numeric_limits<RealType>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
  return setScaledDriftPbyPandNodeCorr(Tau, MassInvP, g, gscaled);
}

void QMCUpdateBase::setReleasedNodeMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    RealType M          = std::abs((*it)->Weight);
    (*it)->Multiplicity = std::floor(M + RandomGen());
  }
}

void QMCUpdateBase::setMultiplicity(WalkerIter_t it, WalkerIter_t it_end)
{
  for (; it != it_end; ++it)
  {
    RealType M = (*it)->Weight;
    if ((*it)->Age > MaxAge)
      M = std::min((RealType)0.5, M);
    else if ((*it)->Age > 0)
      M = std::min((RealType)1.0, M);
    (*it)->Multiplicity = M + RandomGen();
  }
}

void QMCUpdateBase::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool recompute)
{
  for (; it != it_end; ++it)
  {
    advanceWalker(**it, recompute);
  }
}

} // namespace qmcplusplus

#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#include "Configuration.h"
#include "Particle/Reptile.h"
#include <cmath>

//////////////////////////////////////////////////////////////////////////
//
//  This driver proposes all-electron moves like in DMCUpdateAll.cpp.
//
//  For H4 and H2 systems, it appears as though the Umrigar pbyp scaled drift has similar instabilities
//  as in the RQMCMultiple.cpp driver.  Hence, the bare propagator with filtered local energies is used.
//
//  The symmetric link action (Ceperley and Pierleoni) is hardcoded into this driver.
//
//////////////////////////////////////////////////////////////////////////////

namespace qmcplusplus
{

/// Constructor.
RMCUpdateAllWithDrift::RMCUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg,std::vector<int> act, std::vector<int> tp):  QMCUpdateBase(w,psi,h,rg), Action(act), TransProb(tp)
{
  scaleDrift=false;
  actionType=SYM_ACTION;
}

RMCUpdateAllWithDrift::~RMCUpdateAllWithDrift()
{
}

void RMCUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  IndexType direction = W.reptile->direction;
  IndexType forward =(1-direction)/2;
  IndexType backward=(1+direction)/2;
  Walker_t& curhead=W.reptile->getHead();
  W.loadWalker(curhead, false);
  RealType nodecorr=1;
  if (scaleDrift==true)
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,curhead.G,drift);
  else
    assignDrift(m_tauovermass,curhead.G,drift);
  //app_log()<<"Sign head = "<<curhead.Properties(SIGN)<<endl;
  //app_log()<<"Old phase = "<<Psi.getPhase()<<endl;
  makeGaussRandomWithEngine(deltaR,RandomGen);
  RealType r2proposed=Dot(deltaR,deltaR);
  RealType r2accept=0.0;
  W.reptile->r2prop += r2proposed;
  W.reptile->r2samp++;
  if (!W.makeMoveWithDrift(curhead,drift ,deltaR, m_sqrttau))
  {
    ++nReject;
    H.rejectedMove(W,curhead);
    curhead.Age+=1;
    W.reptile->flip();
    return;
  }
  RealType logpsi(Psi.evaluateLog(W));
  // app_log()<<"Sign newhead = "<<W.Properties(SIGN)<<endl;
  //RealType* restrict old_headProp ((*it)->getPropertyBase());
  //old_headProp[TransProb[forward]]= 0.5*Dot(deltaR,deltaR);
  curhead.Properties(W.reptile->TransProb[forward])=-0.5*Dot(deltaR,deltaR);
  curhead.Properties(W.reptile->Action[forward])= 0.5*0.5*Dot(deltaR,deltaR);
  Walker_t::ParticlePos_t fromdeltaR(deltaR);
  RealType nodecorrnew = 1.0;
  if (scaleDrift==true)
    nodecorrnew = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  else
    assignDrift(m_tauovermass,W.G,drift);
  fromdeltaR = curhead.R - W.R - drift;
  RealType* restrict new_headProp (W.getPropertyBase());
  W.Properties(W.reptile->TransProb[backward])= -m_oneover2tau*Dot(fromdeltaR,fromdeltaR);
  W.Properties(W.reptile->Action[backward])= 0.5*m_oneover2tau*Dot(fromdeltaR,fromdeltaR);
  Walker_t& lastbead(W.reptile->getTail()), nextlastbead(W.reptile->getNext());
  //Implementing the fixed-node approximation.  If phase difference is not a multiple of 2pi, bounce away from node.
  RealType newphase = Psi.getPhase();
  RealType phasediff = newphase - curhead.Properties(SIGN);
  //Reject & bounce if node crossed.
  if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
  {
    ++nReject;
    H.rejectedMove(W,curhead);
    curhead.Age+=1;
    W.reptile->flip();
    //app_log()<<"hit a node.  Bouncing...\n";
    return;
  }
  RealType eloc=H.evaluate(W);
  //new_headProp[Action[2]]= 0.5*Tau*eloc;
  ////////////////////////////////////////////////////////////////////////
  ///  Like DMC, this filters the local energy to ignore divergences near pathological points in phase space.
  ////////////////////////////////////////////////////////////////////////
  RealType eest = W.reptile->eest;
  RealType fbet = std::max(eest - curhead.Properties(LOCALENERGY), eest - eloc);
  //  app_log()<<"eval = "<<eest<<" estdev="<<stddev<<endl;
  RealType rawcutoff=100*std::sqrt(W.reptile->evar);
  RealType cutoffmax = 1.5*rawcutoff;
  RealType cutoff=1;
  if (fbet > rawcutoff)
    cutoff = 1-(fbet - rawcutoff)/(rawcutoff*0.5);
  if( fbet > cutoffmax )
    cutoff=0;
  //////////////////////////////////////////////////////////////////////////
  RealType tauscale = W.reptile->tauscale;
  W.Properties(W.reptile->Action[2])= 0.5*Tau*eloc*cutoff*tauscale;
  RealType dS = 0;
  RealType acceptProb=1;
  if (actionType==SYM_ACTION)
  {
    dS = +(curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI))
         + curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
         + curhead.Properties(W.reptile->Action[forward]) + W.Properties(W.reptile->Action[backward])
         - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]))
         - (lastbead.Properties(W.reptile->Action[forward]) + nextlastbead.Properties(W.reptile->Action[backward]));
    acceptProb=std::exp(-dS + (nextlastbead.Properties(W.reptile->TransProb[backward]) - curhead.Properties(W.reptile->TransProb[forward])));
  }
  else
  {
    dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
         - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
    acceptProb=std::min(1.0,std::exp(-dS ));
  }
  /*	RealType dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
  					- (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
  		RealType acceptProb=std::min(1.0,std::exp(-dS ));		*/
  if ((RandomGen() < acceptProb ) || (curhead.Age>=MaxAge))
  {
    r2accept=r2proposed;
    W.reptile->r2accept+=r2accept;
    MCWalkerConfiguration::Walker_t& overwriteWalker(W.reptile->getNewHead());
    if (curhead.Age>=MaxAge)
      app_log()<<"\tForce Acceptance...\n";
    W.saveWalker(overwriteWalker);
    //overwriteWalker.Properties(LOCALENERGY)=eloc;
    overwriteWalker.Properties(W.reptile->Action[forward])=0;
    overwriteWalker.Properties(W.reptile->Action[backward])=W.Properties(W.reptile->Action[backward]);
    overwriteWalker.Properties(W.reptile->Action[2])=W.Properties(W.reptile->Action[2]);
    overwriteWalker.Properties(W.reptile->TransProb[forward])=W.Properties(W.reptile->TransProb[forward]);
    overwriteWalker.Properties(W.reptile->TransProb[backward])=W.Properties(W.reptile->TransProb[backward]);
    overwriteWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
    H.auxHevaluate(W,overwriteWalker);
    H.saveProperty(overwriteWalker.getPropertyBase());
    overwriteWalker.Age=0;
    ++nAccept;
  }
  else
  {
    ++nReject;
    H.rejectedMove(W,curhead);
    curhead.Age+=1;
    W.reptile->flip();
    return;
  }
}


}

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
//#include "Utilities/OhmmsInfo.h"
//#include "Particle/MCWalkerConfiguration.h"
//#include "Particle/HDFWalkerIO.h"
//#include "ParticleBase/ParticleUtility.h"
//#include "ParticleBase/RandomSeqGenerator.h"
//#include "ParticleBase/ParticleAttribOps.h"
//#include "Message/Communicate.h"
//#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
CSRMCUpdateAll::CSRMCUpdateAll(MCWalkerConfiguration& w,
                               TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{ }

/**  Advance all the walkers one timstep.
 */
void CSRMCUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  int iwlk(0);
  int nPsi_minus_one(nPsi-1);
  while(it != it_end)
  {
    MCWalkerConfiguration::Walker_t &thisWalker(**it);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    if(useDrift)
      W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
    else
      W.R = m_sqrttau*deltaR + thisWalker.R;
    //update the distance table associated with W
    //DistanceTable::update(W);
    W.update();
    //Evaluate Psi and graidients and laplacians
    //\f$\sum_i \ln(|psi_i|)\f$ and catching the sign separately
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(W);
      Psi1[ipsi]->L=W.L;
      Psi1[ipsi]->G=W.G;
      sumratio[ipsi]=1.0;
    }
    // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    for(int ipsi=0; ipsi< nPsi_minus_one; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=avgNorm[ipsi]/avgNorm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
    }
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      invsumratio[ipsi]=1.0/sumratio[ipsi];
    RealType g = sumratio[0]/thisWalker.Multiplicity*
                 std::exp(2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
    if(useDrift)
    {
      //forward green function
      RealType logGf = -0.5*Dot(deltaR,deltaR);
      PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
      for(int ipsi=1; ipsi< nPsi ; ipsi++)
      {
        PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
      }
      setScaledDrift(Tau,drift);
      //backward green function
      deltaR = thisWalker.R - W.R - drift;
      RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      g *= std::exp(logGb-logGf);
    }
    //Original
    //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*
    //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));
    //Reuse Multiplicity to store the sumratio[0]
    //This is broken up into two pieces
    //RealType g = sumratio[0]/thisWalker.Multiplicity*
    // 	std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
    if(Random() > g)
    {
      thisWalker.Age++;
      ++nReject;
    }
    else
    {
      thisWalker.Age=0;
      thisWalker.Multiplicity=sumratio[0];
      thisWalker.R = W.R;
      thisWalker.Drift = drift;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        W.L=Psi1[ipsi]->L;
        W.G=Psi1[ipsi]->G;
        RealType et = H1[ipsi]->evaluate(W);
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
        //multiEstimator->updateSample(iwlk,ipsi,et,invsumratio[ipsi]);
        H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
      ++nAccept;
    }
    ++it;
    ++iwlk;
  }
}

void CSRMCUpdateAll:advanceWalkersRMC()
{
void RMCUpdateAllWithDrift::advanceWalkersVMC()
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
	//app_log()<<"Sign head = "<<curhead.Properties(SIGN)<< std::endl;
	//app_log()<<"Old phase = "<<Psi.getPhase()<< std::endl;
	makeGaussRandomWithEngine(deltaR,RandomGen);
	RealType r2proposed=Dot(deltaR,deltaR);
	RealType r2accept=0.0;
	W.reptile->r2prop += r2proposed;
	W.reptile->r2samp++;
	if (!W.makeMoveWithDrift(curhead,drift ,deltaR, m_sqrttau))
	{
	  ++nReject;
	  H.rejectedMove(W,curhead);
	  // curhead.Age+=1;
	  //W.reptile->flip();
	  return;
	}
	
	RealType logpsi(Psi.evaluateLog(W));
	RealType logGf=-0.5*Dot(deltaR,deltaR);
	RealType Action_forward=-0.5*logGf;
	curhead.Properties(W.reptile->TransProb[forward])=-0.5*Dot(deltaR,deltaR);
	curhead.Properties(W.reptile->Action[forward])= 0.5*0.5*Dot(deltaR,deltaR);
	
	W.reptile->saveTransProb(curhead,+1,logGf);
	W.reptile->saveAction(curhead,+1,Action_forward);
	
	Walker_t::ParticlePos_t fromdeltaR(deltaR);
	
	
	if (scaleDrift==true)
	  setScaledDrift(m_tauovermass,W.G,drift);
	else
	  assignDrift(m_tauovermass,W.G,drift);
	fromdeltaR = curhead.R - W.R - drift;
	RealType* restrict new_headProp (W.getPropertyBase());
	
	RealType logGb = -m_oneover2tau*Dot(fromdeltaR,fromdeltaR);
	RealType Action_backward=-0.5*logGb;
	
	W.reptile->saveTransProb(W,-1,logGb);
	//W.reptile->saveAction(W,-1,Action_backward);
	
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
	  // curhead.Age+=1;
	  // W.reptile->flip();
	  //app_log()<<"hit a node.  Bouncing...\n";
	  return;
	}
	RealType eloc=H.evaluate(W);
	new_headProp[Action[2]]= 0.5*Tau*eloc;

	////////////////////////////////////////////////////////////////////////
	///  Like DMC, this filters the local energy to ignore divergences near pathological points in phase space.
	////////////////////////////////////////////////////////////////////////
	RealType eest = W.reptile->eest;
	RealType fbet = std::max(eest - curhead.Properties(LOCALENERGY), eest - eloc);
	//   app_log()<<"eval = "<<eest<<" estdev="<<stddev<< std::endl;
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
	  RealType oldhead_logpsi=curhead.Properties(LOGPSI);
	  RealType oldtail_logpsi=lastbead.Properties(LOGPSI);
	  RealType newhead_logpsi=logpsi;
	  RealType newtail_logpsi=nextlastbead.Properties(LOGPSI);
	
  	  RealType oldhead_e=curhead.Properties(LOCALENERGY);
	  RealType oldtail_e=lastbead.Properties(LOCALENERGY);
	  RealType newhead_e=W.Properties(LOCALENERGY);
	  RealType newtail_e=nextlastbead.Properties(LOCALENERGY);
	
	  RealType head_forward=W.reptile->getTransProb(curhead,+1);
	  RealType head_backward=W.reptile->getTransProb(W,-1);
	  RealType tail_forward=W.reptile->getTransProb(lastbead,+1);
	  RealType tail_backward=W.reptile->getTransProb(nextlastbead,-1);
	
	  RealType dS_head=branchEngine->symLinkAction(head_forward, head_backward, newhead_e, oldhead_e);
	  RealType dS_tail=branchEngine->symLinkAction(tail_forward, tail_backward, newtail_e, oldtail_e);
	
	}
	else
	{
	  RealType dS_head=branchEngine->DMCLinkAction(eloc,curhead.Properties(LOCALENERGY));
	  RealType dS_tail=branchEngine->DMCLinkAction(lastbead.Properties(LOCALENERGY),nextlastbead.Properties(LOCALENERGY));
	}
	acceptProb=std::exp(logGb-logGf+2.0*(logpsi-curhead.Properties(LOGPSI)));

	if (RandomGen() < acceptProb )
	{
	  r2accept=r2proposed;
	  W.reptile->r2accept+=r2accept;
	  MCWalkerConfiguration::Walker_t& overwriteWalker(W.reptile->getNewHead());
	
	  W.saveWalker(overwriteWalker);
	  overwriteWalker.Properties(LOCALENERGY)=eloc;
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
	  return;
	}
}
}
}


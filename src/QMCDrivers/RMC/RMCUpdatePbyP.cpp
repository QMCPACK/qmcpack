#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"
#include "Configuration.h"
#include "Particle/Reptile.h"
#include <cmath>
//////////////////////////////////////////////////////////////////////////
//  This driver is strongly based off the method used by Lucas Wagner in QWalk.
//
//  This driver proposes an "all-electron" configuration by performing N single particle moves,
//  accepting or rejecting at each step.  Umrigar's scaled drift is used to generate the VMC configurations.
//
//  The configuration is then accepted/rejected based only on the symmetrized DMC action, which
//  amounts to Maroni/Baroni's method.  Note that energy filtering (like in DMC) is used here as well.
//
//  Until serious discrepencies are detected, it is strongly advised that this driver be used over the
//  RMCUpdateAll, as the time step errors seem to be really reduced using this method.
//////////////////////////////////////////////////////////////////////////////


namespace qmcplusplus
{


/*  void add_rmc_timers(vector<NewTimer*>& timers)
  {
    timers.push_back(new NewTimer("RMCUpdatePbyP::advance")); //timer for the walker loop
    timers.push_back(new NewTimer("RMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
   // timers.push_back(new NewTimer("RMCUpdatePbyP::updateMBO")); //timer for measurements
   // timers.push_back(new NewTimer("RMCUpdatePbyP::energy")); //timer for measurements
    for (int i=0; i<timers.size(); ++i) TimerManager.addTimer(timers[i]);
  }*/
/// Constructor.
RMCUpdatePbyPWithDrift::RMCUpdatePbyPWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg,std::vector<int> act, std::vector<int> tp):  QMCUpdateBase(w,psi,h,rg), Action(act), TransProb(tp)
{
  //add_rmc_timers(myTimers);
}

RMCUpdatePbyPWithDrift::~RMCUpdatePbyPWithDrift()
{
}


void RMCUpdatePbyPWithDrift::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
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
  }
}

void RMCUpdatePbyPWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  IndexType direction = W.reptile->direction;
  IndexType forward =(1-direction)/2;
  IndexType backward=(1+direction)/2;
  Walker_t& curhead=W.reptile->getHead();
  Walker_t::Buffer_t& w_buffer(curhead.DataSet);
  W.loadWalker(curhead, true);
  W.R=curhead.R;
  W.update();
  //W.loadWalker(awalker,UpdatePbyP);
  if (curhead.DataSet.size())
    curhead.DataSet.clear();
  curhead.DataSet.rewind();
  RealType logpsi=Psi.registerData(W,curhead.DataSet);
  RealType logpsi2=Psi.updateBuffer(W,curhead.DataSet,false);
  // curhead.G=W.G;
  //  curhead.L=W.L;
  // Walker_t& thisWalker(**it);
  Psi.copyFromBuffer(W,w_buffer);
  //  RealType* restrict old_headProp ((*it)->getPropertyBase());
  //  old_headProp[TransProb[forward]]=0;
  //old_headProp[Action[forward]]=0;
  // RealType* restrict new_headProp (W.getPropertyBase());
  // new_headProp[Action[backward]]= 0;
  // new_headProp[TransProb[backward]]= 0;
  //create a 3N-Dimensional Gaussian with variance=1
  /* makeGaussRandomWithEngine(deltaR,RandomGen);
   int nAcceptTemp(0);
   int nRejectTemp(0);
   //copy the old energy and scale factor of drift
   RealType eold(curhead.Properties(LOCALENERGY));
   RealType vqold(curhead.Properties(DRIFTSCALE));
   RealType enew(eold);
   RealType rr_proposed=0.0;
   RealType rr_accepted=0.0;
   RealType gf_acc=1.0;
   RealType ratio_acc=1.0;
   bool rejected(false);
  RealType logpsiold = curhead.Properties(LOGPSI);
   for(int iat=0; iat<NumPtcl; ++iat)
   {
     GradType grad_iat=Psi.evalGrad(W,iat);
     PosType dr;
    // PosType oldpos=W.R[iat];
     getScaledDrift(m_tauovermass, grad_iat, dr);
     dr += m_sqrttau * deltaR[iat];
   // dr = m_tauovermass*grad_iat + m_sqrttau * deltaR[iat];
    // app_log()<<"iat "<<iat<<"  dr="<<dr<<endl;
     RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
     rr_proposed+=rr;
     if (rr>m_r2max)
     {
    rejected=true;
       W.reptile->flip();
    }

     if(!W.makeMoveAndCheck(iat,dr))
     {
       rejected=true;
       W.reptile->flip();
       //return;
     }
     PosType newpos(W.R[iat]);
     RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
     if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
     {
       rejected=true;
       ++nRejectTemp;
       ++nNodeCrossing;
       W.rejectMove(iat);
       Psi.rejectMove(iat);
       W.reptile->flip();
       //return;
     }


     else
     {
   // G = W.G+dG;
         RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);

         //Use the force of the particle iat
         //RealType scale=getDriftScale(m_tauovermass,G[iat]);
         //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
        // getScaledDrift(m_tauovermass, G[iat], dr);
         getScaledDrift(m_tauovermass, grad_iat, dr);
         dr = curhead.R[iat] - newpos - dr;

         RealType logGb = -m_oneover2tau*dot(dr,dr);
         RealType prob = ratio*ratio*std::exp(logGb-logGf);

       rr_accepted+=rr;
       ratio_acc*=ratio;

       if(RandomGen() < prob)
         {
           valid_move=true;
           ++nAcceptTemp;
           W.acceptMove(iat);
           Psi.acceptMove(W,iat);
           W.G = G;
           W.L += dL;
           rr_accepted+=rr;
           gf_acc *=prob;//accumulate the ratio
         }
         else
         {
           ++nRejectTemp;
           W.rejectMove(iat); Psi.rejectMove(iat);
         }

        //dr = m_tauovermass*grad_iat;
    //   dr = curhead.R[iat] - newpos - dr;
       W.acceptMove(iat);
       Psi.acceptMove(W,iat);
     }
   }
   */
  RealType logpsiold = curhead.Properties(LOGPSI);
  makeGaussRandomWithEngine(deltaR,RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  RealType eold(curhead.Properties(LOCALENERGY));
  RealType vqold(curhead.Properties(DRIFTSCALE));
  RealType enew(eold);
  RealType rr_proposed=0.0;
  RealType rr_accepted=0.0;
  RealType gf_acc=1.0;
  RealType logGf=0;
  RealType logGb=0;
  RealType ratio_acc=1;
  //   myTimers[1]->start();
  for(int iat=0; iat<NumPtcl; ++iat)
  {
    PosType dr;
    //get the displacement
    //RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
    //PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
    getScaledDrift(m_tauovermass,W.G[iat],dr);
    dr += m_sqrttau*deltaR[iat];
    //RealType rr=dot(dr,dr);
    RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
    rr_proposed+=rr;
    // if (W.reptile->r2
    if(rr>m_r2max)
    {
      ++nRejectTemp;
      W.reptile->flip();
      return;
    }
    //PosType newpos(W.makeMove(iat,dr));
    if(!W.makeMoveAndCheck(iat,dr))
    {
      ++nRejectTemp;
      W.reptile->flip();
      return;
    }
    PosType newpos(W.R[iat]);
    RealType ratio=Psi.ratio(W,iat,dG,dL);
    bool valid_move=false;
    //node is crossed reject the move
    //if(Psi.getPhase() > numeric_limits<RealType>::epsilon())
    if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      ++nRejectTemp;
      ++nNodeCrossing;
      W.rejectMove(iat);
      Psi.rejectMove(iat);
      W.reptile->flip();
      return;
    }
    else
    {
      G = W.G+dG;
      RealType logGf_pbyp = -0.5*dot(deltaR[iat],deltaR[iat]);
      logGf+=logGf_pbyp;
      //Use the force of the particle iat
      //RealType scale=getDriftScale(m_tauovermass,G[iat]);
      //dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
      getScaledDrift(m_tauovermass, G[iat], dr);
      dr = curhead.R[iat] - newpos - dr;
      RealType logGb_pbyp = -m_oneover2tau*dot(dr,dr);
      logGb+=logGb_pbyp;
      RealType prob = ratio*ratio*std::exp(logGb_pbyp-logGf_pbyp);
      //this is useless
      //RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
      if(RandomGen() < prob)
      {
        valid_move=true;
        ++nAcceptTemp;
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        W.G = G;
        W.L += dL;
        rr_accepted+=rr;
        gf_acc *=prob;//accumulate the ratio
        ratio_acc*=ratio;
      }
      else
      {
        ++nRejectTemp;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
        // W.reptile->flip();
        // return;
      }
    }
  }
//   RealType driftscaleold=getDriftScale(m_tauovermass,curhead.G);
//   assignDrift(driftscaleold, curhead.G, drift);
  //  RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,curhead.G,drift);
  //RealType nodecorr=1;
  //assignDrift(m_tauovermass,curhead.G,drift);
  //app_log()<<"Sign head = "<<curhead.Properties(SIGN)<<endl;
  //app_log()<<"Old phase = "<<Psi.getPhase()<<endl;
  //  makeGaussRandomWithEngine(deltaR,RandomGen);
  //  RealType r2proposed=rr_proposed;
  //  RealType r2accept=0.0;
//  if (W.reptile->r2prop < 0) app_log()<<"r2prop = "<<W.reptile->r2prop<<endl;
  //if (W.reptile->r2samp < 0) app_log()<<"r2samp = "<<W.reptile->r2samp<<endl;
  //if (W.reptile->r2accept < 0) app_log()<<"r2samp = "<<W.reptile->r2accept<<endl;
  W.reptile->r2prop += rr_proposed;
  W.reptile->r2samp++;
  logpsi = std::log(ratio_acc)+logpsiold;
  // 	 RealType logpsi(Psi.evaluateLog(W));
  // 	 if (logpsi2 != logpsi) app_log()<<"no match. "<<logpsi2<<"  vs exact "<<logpsi<<endl;
  // app_log()<<"Sign newhead = "<<W.Properties(SIGN)<<endl;
  RealType tauscale = W.reptile->tauscale;
  RealType tau_eff = Tau*tauscale;
  //RealType* restrict old_headProp ((*it)->getPropertyBase());
  //old_headProp[TransProb[forward]]= 0.5*Dot(deltaR,deltaR);
//     curhead.Properties(W.reptile->TransProb[forward])=-0.5*Dot(deltaR,deltaR);
//	  curhead.Properties(W.reptile->Action[forward])= 0.5*0.5*Dot(deltaR,deltaR);
//	  Walker_t::ParticlePos_t fromdeltaR(deltaR);
  //  RealType driftscalenew=getDriftScale(m_tauovermass,W.G);
  // assignDrift(driftscalenew, W.G, drift);
  //app_log()<<"Driftscalenew="<<driftscalenew<<endl;
  // RealType nodecorrnew = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  // assignDrift(m_tauovermass,W.G,drift);
  //  RealType nodecorrnew = 1.0;
//	  fromdeltaR = curhead.R - W.R - drift;
//	  RealType* restrict new_headProp (W.getPropertyBase());
//	  W.Properties(W.reptile->TransProb[backward])= -m_oneover2tau*Dot(fromdeltaR,fromdeltaR);
//	  W.Properties(W.reptile->Action[backward])= 0.5*m_oneover2tau*Dot(fromdeltaR,fromdeltaR);
  Walker_t& lastbead(W.reptile->getTail()), nextlastbead(W.reptile->getNext());
  //Implementing the fixed-node approximation.  If phase difference is not a multiple of 2pi, bounce away from node.
  //  RealType newphase = Psi.getPhase();
  // RealType phasediff = newphase - curhead.Properties(SIGN);
  RealType eloc=H.evaluate(W);
  //This is going to calculate how many standard deviations eloc is away from the current estimate.  Lets say 10 and it rejects automatically.
  RealType eest = W.reptile->eest;
  RealType stddev = std::fabs(eloc - W.reptile->eest)/std::sqrt(W.reptile->evar);
  RealType fbet = std::max(eest - curhead.Properties(LOCALENERGY), eest - eloc);
  //  app_log()<<"eval = "<<eest<<" estdev="<<stddev<<endl;
  RealType rawcutoff=10*std::sqrt(W.reptile->evar);
  RealType cutoffmax = 1.5*rawcutoff;
  RealType cutoff=1;
  if (fbet > rawcutoff)
    cutoff = 1-(fbet - rawcutoff)/(rawcutoff*0.5);
  if( fbet > cutoffmax )
    cutoff=0;
  /*  if( stddev > 10 )
    {
  	++nReject;
  	H.rejectedMove(W,curhead);
  	curhead.Age+=1;
  	W.reptile->flip();
  	app_log()<<"Rejecting cause of large change in E\n";
  	app_log()<<"\teloc="<<eloc<<" eavg,var = "<<W.reptile->eest<<" , "<<std::sqrt(W.reptile->evar)<<endl;
  	app_log()<<"\t\tDe="<<std::fabs(eloc - W.reptile->eest)<<"  Drel="<<stddev<<endl;


  	return;
    }*/
  //new_headProp[Action[2]]= 0.5*Tau*eloc;
  // W.Properties(W.reptile->Action[2])= 0.5*eloc*tau_eff*driftscalenew/Tau;
  //app_log()<<nodecorrnew<<endl;
  //   RealType scaledaction=eest + (eloc - eest)*nodecorrnew;
  // W.Properties(W.reptile->Action[2])= 0.5*eloc*tau_eff;
  W.Properties(W.reptile->Action[2])=0.5*eloc*tau_eff*cutoff;
  //W.Properties(W.reptile->Action[2])= 0.5*scaledaction*tau_eff*cutoff;
  //app_log()<<"\tDriftscalenew="<<driftscalenew<<endl;
  //H.evaluate()
  /*		RealType dS = (curhead.Properties(LOGPSI) + lastbead.Properties(LOGPSI) - logpsi - nextlastbead.Properties(LOGPSI))
  					+ curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
  					+ curhead.Properties(W.reptile->Action[forward]) + W.Properties(W.reptile->Action[backward])
  					- (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]))
  					- (lastbead.Properties(W.reptile->Action[forward]) + nextlastbead.Properties(W.reptile->Action[backward]));
  	*/
  //	RealType dS = 0.5*tau_eff*((eloc+curhead.Properties(LOCALENERGY))*cutoff - nextlastbead.Properties(LOCALENERGY) - lastbead.Properties(LOCALENERGY));
//			RealType dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
//						-(lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
//
  //RealType acceptProb=std::min(1.0,std::exp(-dS + (nextlastbead.Properties(W.reptile->TransProb[backward]) - curhead.Properties(W.reptile->TransProb[forward]))));
  RealType dS = curhead.Properties(W.reptile->Action[2]) + W.Properties(W.reptile->Action[2])
                - (lastbead.Properties(W.reptile->Action[2]) + nextlastbead.Properties(W.reptile->Action[2]));
  RealType acceptProb=std::min(1.0,std::exp(-dS ));
  if ((RandomGen() <= acceptProb ) || (curhead.Age>=MaxAge))
  {
    // r2accept=r2proposed;
    W.reptile->r2accept+=rr_accepted;
    MCWalkerConfiguration::Walker_t& overwriteWalker(W.reptile->getNewHead());
    if (curhead.Age>=MaxAge)
      app_log()<<"\tForce Acceptance...\n";
    Walker_t::Buffer_t& o_buffer(overwriteWalker.DataSet);
    overwriteWalker.R = W.R;
    RealType logpsi = Psi.updateBuffer(W,o_buffer,false);
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
    W.reptile->accumulateE(eloc);
    ++nAccept;
  }
  else
  {
    ++nReject;
    H.rejectedMove(W,curhead);
    curhead.Age+=1;
    W.reptile->flip();
    // return;
  }
}



}

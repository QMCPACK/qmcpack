//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/// Constructor.
DMCNonLocalUpdatePbyP::DMCNonLocalUpdatePbyP(MCWalkerConfiguration& w
    , TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
  : QMCUpdateBase(w,psi,h,rg)
{
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::advance")); //timer for the walker loop
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::updateMBO")); //timer for measurements
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::energy")); //timer for measurements
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
}

/// destructor
DMCNonLocalUpdatePbyP::~DMCNonLocalUpdatePbyP() { }

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void DMCNonLocalUpdatePbyP::advanceWalkers(WalkerIter_t it
    , WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  for(; it!=it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    RealType eold(thisWalker.Properties(LOCALENERGY));
    //RealType emixed(eold), enew(eold);
    RealType enew(eold);
    W.R = thisWalker.R;
    w_buffer.rewind();
    W.copyFromBuffer(w_buffer);
    Psi.copyFromBuffer(W,w_buffer);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR, RandomGen);
    bool notcrossed(true);
    int nAcceptTemp(0);
    int nRejectTemp(0);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
      RealType sc=getDriftScale(Tau,W.G[iat]);
      PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
      //RealType rr=dot(dr,dr);
      RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
      rr_proposed+=rr;
      if(rr>m_r2max)//too big move
      {
        ++nRejectTemp;
        continue;
      }
      PosType newpos(W.makeMove(iat,dr));
      RealType ratio=Psi.ratio(W,iat,dG,dL);
      //node is crossed reject the move
      if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nRejectTemp;
        ++nNodeCrossing;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
      else
      {
        G = W.G+dG;
        RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
        //RealType scale=getDriftScale(Tau,G);
        RealType scale=getDriftScale(m_tauovermass,G[iat]);
        dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
        RealType logGb = -m_oneover2tau*dot(dr,dr);
        RealType prob = std::min(1.0,ratio*ratio*std::exp(logGb-logGf));
        if(RandomGen() < prob)
        {
          ++nAcceptTemp;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          W.G = G;
          W.L += dL;
          //assignDrift(scale,G,thisWalker.Drift);
          rr_accepted+=rr;
        }
        else
        {
          ++nRejectTemp;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }//end of drift+diffusion for all the particles of a walker
    myTimers[1]->stop();
    nonLocalOps.reset();
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      myTimers[2]->start();
      thisWalker.R = W.R;
      w_buffer.rewind();
      W.copyToBuffer(w_buffer);
      //RealType psi = Psi.evaluate(W,w_buffer);
      RealType logpsi = Psi.evaluateLog(W,w_buffer);
      myTimers[2]->stop();
      myTimers[3]->start();
      enew= H.evaluate(W,nonLocalOps.Txy);
      myTimers[3]->stop();
      //thisWalker.resetProperty(std::log(std::abs(psi)),psi,enew,rr_accepted,rr_proposed,1.0);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
      thisWalker.Drift=W.G;
    }
    else
    {
      H.rejectedMove(W,thisWalker);
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      ++nAllRejected;
      enew=eold;//copy back old energy
    }
    int ibar = nonLocalOps.selectMove(RandomGen());
    //make a non-local move
    if(ibar)
    {
      myTimers[2]->start();
      int iat=nonLocalOps.id(ibar);
      PosType newpos(W.makeMove(iat, nonLocalOps.delta(ibar)));
      RealType ratio=Psi.ratio(W,iat,dG,dL);
      W.acceptMove(iat);
      Psi.acceptMove(W,iat);
      W.G += dG;
      W.L += dL;
      PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
      thisWalker.R[iat]=W.R[iat];
      w_buffer.rewind();
      W.copyToBuffer(w_buffer);
      //RealType psi = Psi.evaluate(W,w_buffer);
      RealType logpsi = Psi.evaluateLog(W,w_buffer);
      ++NonLocalMoveAccepted;
      myTimers[2]->stop();
    }
    thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}

DMCNonLocalUpdatePbyPFast::DMCNonLocalUpdatePbyPFast(MCWalkerConfiguration& w
    , TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator_t& rg)
  : QMCUpdateBase(w,psi,h,rg)
{
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::advance")); //timer for the walker loop
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::updateMBO")); //timer for measurements
  myTimers.push_back(new NewTimer("DMCNonLocalUpdatePbyP::energy")); //timer for measurements
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
}

/// destructor
DMCNonLocalUpdatePbyPFast::~DMCNonLocalUpdatePbyPFast() { }
void DMCNonLocalUpdatePbyPFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end,
    bool measure)
{
  Timer localTimer;
  myTimers[0]->start();
  for(; it != it_end; ++it)
  {
    //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.R = thisWalker.R;
    w_buffer.rewind();
    W.copyFromBuffer(w_buffer);
    Psi.copyFromBuffer(W,w_buffer);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    int nAcceptTemp(0);
    int nRejectTemp(0);
    //copy the old energy and scale factor of drift
    RealType eold(thisWalker.Properties(LOCALENERGY));
    RealType enew(eold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    RealType gf_acc=1.0;
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      //get the displacement
      //PosType dr(m_sqrttau*deltaR[iat]+thisWalker.Drift[iat]);
      GradType grad_iat=Psi.evalGrad(W,iat);
      RealType sc=getDriftScale(m_tauovermass,grad_iat);
      PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_iat));
      //RealType rr=dot(dr,dr);
      RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
      rr_proposed+=rr;
      if(rr>m_r2max)
      {
        ++nRejectTemp;
        continue;
      }
      //PosType newpos(W.makeMove(iat,dr));
      if(!W.makeMoveAndCheck(iat,dr))
        continue;
      PosType newpos(W.R[iat]);
      RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
      bool valid_move=false;
      //node is crossed reject the move
      if(branchEngine->phaseChanged(Psi.getPhaseDiff()))
      {
        ++nRejectTemp;
        ++nNodeCrossing;
        W.rejectMove(iat);
        Psi.rejectMove(iat);
      }
      else
      {
        RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
        //Use the force of the particle iat
        RealType scale=getDriftScale(m_tauovermass,grad_iat);
        dr = thisWalker.R[iat]-newpos-scale*real(grad_iat);
        RealType logGb = -m_oneover2tau*dot(dr,dr);
        RealType prob = ratio*ratio*std::exp(logGb-logGf);
        if(RandomGen() < prob)
        {
          valid_move=true;
          ++nAcceptTemp;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          rr_accepted+=rr;
          gf_acc *=prob;//accumulate the ratio
        }
        else
        {
          ++nRejectTemp;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
    }
    myTimers[1]->stop();
    //RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
    nonLocalOps.reset();
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      myTimers[2]->start();
      thisWalker.Age=0;
      thisWalker.R = W.R;
      w_buffer.rewind();
      W.updateBuffer(w_buffer);
      RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      myTimers[2]->stop();
      myTimers[3]->start();
      enew= H.evaluate(W,nonLocalOps.Txy);
      myTimers[3]->stop();
      //nodecorr=getNodeCorrection(W.G,thisWalker.Drift);
      //thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,nodecorr);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0 );
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    else
    {
      //all moves are rejected: does not happen normally with reasonable wavefunctions
      H.rejectedMove(W,thisWalker);
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      ++nAllRejected;
      enew=eold;//copy back old energy
      gf_acc=1.0;
      std::cerr << "  Failed to move any particle " << thisWalker.ID << std::endl;
    }
    int ibar = nonLocalOps.selectMove(RandomGen());
    //make a non-local move
    if(ibar)
    {
      int iat=nonLocalOps.id(ibar);
      if(!W.makeMoveAndCheck(iat,nonLocalOps.delta(ibar)))
        continue;
      myTimers[2]->start();
      /////RealType ratio=Psi.ratio(W,iat);
      /////W.acceptMove(iat);
      /////Psi.acceptMove(W,iat);
      /////thisWalker.R[iat]=W.R[iat];
      /////w_buffer.rewind();
      /////W.copyToBuffer(w_buffer);
      /////RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      RealType ratio=Psi.ratio(W,iat,dG,dL);
      W.acceptMove(iat);
      Psi.acceptMove(W,iat);
      W.G += dG;
      W.L += dL;
      thisWalker.R[iat]=W.R[iat];
      w_buffer.rewind();
      W.copyToBuffer(w_buffer);
      RealType logpsi = Psi.evaluateLog(W,w_buffer);
      PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
      ++NonLocalMoveAccepted;
      myTimers[2]->stop();
    }
    //2008-06-26: select any
    //bare green function by setting nodecorr=nodecorr_old=1.0
    thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
    //Filtering extreme energies
    //thisWalker.Weight *= branchEngine->branchWeight(eold,enew);
    //using the corrections: see QMCUpdateBase::getNodeCorrection
    //thisWalker.Weight *= branchEngine->branchWeight(enew,eold,nodecorr,nodecorr_old);
    //using the corrections: see QMCUpdateBase::getNodeCorrection  including gf_acc
    //RealType odd=std::min(gf_acc,1.0)
    //thisWalker.Weight *= branchEngine->branchWeight(enew,eold,nodecorr,nodecorr_oldi,odd);
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}
}


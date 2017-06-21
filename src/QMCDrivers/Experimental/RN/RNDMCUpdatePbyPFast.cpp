//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
//#define TEST_INNERBRANCH


namespace qmcplusplus
{

/// Constructor.
RNDMCUpdatePbyPFast::RNDMCUpdatePbyPFast(MCWalkerConfiguration& w, MCWalkerConfiguration& wg,
    TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,guide,h,rg), maxS(100), W_G(wg), efn(0)
{
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::advance")); //timer for the walker loop
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::movePbyP")); //timer for MC, ratio etc
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::updateMBO")); //timer for measurements
  myTimers.push_back(new NewTimer("RNDMCUpdatePbyP::energy")); //timer for measurements
  TimerManager.addTimer(myTimers[0]);
  TimerManager.addTimer(myTimers[1]);
  TimerManager.addTimer(myTimers[2]);
  TimerManager.addTimer(myTimers[3]);
  myParams.add(maxS,"Smax","int");
  myParams.add(efn,"fnenergy","float");
}

/// destructor
RNDMCUpdatePbyPFast::~RNDMCUpdatePbyPFast() { }

/** advance all the walkers with killnode==no
 * @param nat number of particles to move
 *
 * When killnode==no, any move resulting in node-crossing is treated
 * as a normal rejection.
 */
void RNDMCUpdatePbyPFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end
    , bool measure)
{
  myTimers[0]->start();
  for(; it != it_end; ++it)
  {
    //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.loadWalker(thisWalker,true);
    Psi.copyFromBuffer(W,w_buffer);
    RealType oldphase=Psi.getPhase();
    W_G.loadWalker(thisWalker,true);
    Guide.copyFromBuffer(W_G,w_buffer);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandomWithEngine(deltaR,RandomGen);
    int nAcceptTemp(0);
    int nRejectTemp(0);
    //copy the old energy and scale factor of drift
    RealType eold(thisWalker.Properties(ALTERNATEENERGY));
    RealType feold(thisWalker.Properties(LOCALENERGY));
    RealType enew(eold);
    RealType fenew(feold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    int crossedNode(0);
    myTimers[1]->start();
    for(int iat=0; iat<NumPtcl; ++iat)
    {
      //get the displacement
      GradType grad_iat=Guide.evalGrad(W_G,iat);
      PosType dr,drf;
      //RealType sc=getDriftScale(m_tauovermass,grad_iat);
      //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_iat));
//         getScaledDrift(m_tauovermass, grad_iat, dr);
      drf = dr = m_sqrttau *(deltaR[iat]+grad_iat);
      //RealType rr=dot(dr,dr);
      RealType rr=m_tauovermass*dot(deltaR[iat],deltaR[iat]);
      rr_proposed+=rr;
      if(rr>m_r2max)
      {
        ++nRejectTemp;
        continue;
      }
      //PosType newpos(W.makeMove(iat,dr));
      if(!W_G.makeMoveAndCheck(iat,dr))
        continue;
      PosType newpos(W_G.R[iat]);
      RealType ratio = Guide.ratioGrad(W_G,iat,grad_iat);
      bool valid_move=false;
      RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
      dr = thisWalker.R[iat] - newpos - dr;
      RealType logGb = -m_oneover2tau*dot(dr,dr);
      RealType prob = ratio*ratio*std::exp(logGb-logGf);
      if(RandomGen() < prob)
      {
        valid_move=true;
        ++nAcceptTemp;
        W_G.acceptMove(iat);
        Guide.acceptMove(W_G,iat);
        rr_accepted+=rr;
        W.makeMoveAndCheck(iat,drf);
        RealType psi_ratio = Psi.ratio(W,iat);
        if (psi_ratio<0)
          crossedNode+=1;
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
      }
      else
      {
        ++nRejectTemp;
        W_G.rejectMove(iat);
        Guide.rejectMove(iat);
      }
    }
    myTimers[1]->stop();
    //RealType nodecorr_old=thisWalker.Properties(DRIFTSCALE);
    //RealType nodecorr=nodecorr_old;
    bool advanced=true;
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      myTimers[2]->start();
      thisWalker.Age++;
      //w_buffer.rewind();
      //W.updateBuffer(w_buffer);
      RealType logpsi = Guide.updateBuffer(W_G,w_buffer,false);
      W_G.saveWalker(thisWalker);
      myTimers[2]->stop();
      myTimers[3]->start();
      fenew= H.evaluate(W);
      myTimers[3]->stop();
      enew = -0.5*(Sum(W_G.L)+Dot(W_G.G,W_G.G)) + thisWalker.Properties(LOCALPOTENTIAL);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),fenew,rr_accepted,rr_proposed,1.0 );
      thisWalker.resetReleasedNodeProperty(fenew,enew);
      H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    else
    {
      //all moves are rejected: does not happen normally with reasonable wavefunctions
      advanced=false;
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      H.rejectedMove(W,thisWalker);
      ++nAllRejected;
      enew=eold;//copy back old energy
    }
    thisWalker.ReleasedNodeWeight = std::cos(oldphase-Psi.getPhase())*std::exp(2.0*(Psi.getLogPsi() - Guide.getLogPsi()) );
    if ((crossedNode==0)&&(thisWalker.ReleasedNodeAge==0))
    {
//         Case where walker is still in fixed node population
      thisWalker.Weight *= branchEngine->branchWeight(fenew,feold);
    }
    else
      if (thisWalker.ReleasedNodeAge<maxS)
      {
//         Case where walker is in released node population
//         thisWalker.Weight *= branchEngine->branchWeight(enew,eold);
        thisWalker.Weight *= branchEngine->branchWeightReleasedNode(enew,eold,efn);
        thisWalker.ReleasedNodeAge++;
      }
      else
      {
//         case where walker has exceeded the time limit.
        thisWalker.Weight=0;
        thisWalker.ReleasedNodeWeight=0;
        thisWalker.ReleasedNodeAge++;
      }
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }
  myTimers[0]->stop();
}

void RNDMCUpdatePbyPFast::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=true;
  //   Guide.resizeTempP(W);
  for (; it != it_end; ++it)
  {
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,UpdatePbyP);
    Walker_t::Buffer_t tbuffer;
    RealType logpsi=Psi.registerData(W,tbuffer);
    RealType fene = H.evaluate(W);
    thisWalker.resetProperty(logpsi,Psi.getPhase(),fene);
    H.saveProperty(thisWalker.getPropertyBase());
    W_G.loadWalker(thisWalker,UpdatePbyP);
    RealType logguide=Guide.registerData(W_G,tbuffer);
    thisWalker.DataSet=tbuffer;
    RealType enew = -0.5*(Sum(W_G.L)+Dot(W_G.G,W_G.G)) + thisWalker.Properties(LOCALPOTENTIAL);
    thisWalker.resetReleasedNodeProperty(fene,enew);
    thisWalker.ReleasedNodeAge=0;
    thisWalker.ReleasedNodeWeight=std::exp(2.0*(logpsi-logguide));
    thisWalker.Weight=1.0;
  }
}


/*
void RNDMCUpdatePbyPFast::estimateNormWalkers(std::vector<TrialWaveFunction*>& pclone
    , std::vector<MCWalkerConfiguration*>& wclone
    , std::vector<QMCHamiltonian*>& hclone
    , std::vector<RandomGenerator_t*>& rng
    , std::vector<RealType>& ratio_i_0)
{
  int NumThreads(pclone.size());

  //this can be modified for cache etc
  RealType psi2_i_new[64];
  for (int ip=0; ip<NumThreads; ++ip)
    psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];

//     RealType nn = -std::log(1.0*nSubSteps*W.getTotalNum());
#pragma omp parallel
  {
    int nptcl=W.getTotalNum();
    int ip=omp_get_thread_num();
    RandomGenerator_t& rng_loc(*rng[ip]);

    //copy the new to now
    RealType psi2_i_now=psi2_i_new[ip];
    RealType psi2_0_now=psi2_i_new[0];

    for (int iter=0; iter<nSubSteps; ++iter)
    {
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(deltaR,rng_loc);

      for (int iat=0; iat<nptcl; ++iat)
      {
        PosType dr=m_sqrttau*deltaR[iat];

        bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
        //everyone should skip this; could be a problem with compilers
        if (!movePtcl) continue;

        RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
#pragma omp barrier
        psi2_i_new[ip] = 2.0*logl(std::abs(ratio)) + psi2_i_now;
#pragma omp barrier
#pragma omp critical
        {
          ratio_i_0[ip] += expl( psi2_i_new[ip]-psi2_i_new[0]);
        }
        wclone[ip]->rejectMove(iat);
        pclone[ip]->rejectMove(iat);
      }
    }

    //     Walker_t& thisWalker(*W[ip]);
    //     Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    //     RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
    //     wclone[ip]->saveWalker(*W[ip]);
    //     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
    //     //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
    //     thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);

    //     hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
    //     hclone[ip]->saveProperty(thisWalker.getPropertyBase());
  }
  RealType nn = 1.0/ratio_i_0[0];
  for (int ip=0; ip<NumThreads; ++ip) ratio_i_0[ip]*=nn;
  myTimers[0]->stop();
}*/

}


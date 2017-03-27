//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCDrivers/VMCParticleByParticle.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/CommCreate.h"

namespace qmcplusplus
{

/// Constructor.
VMCParticleByParticle::VMCParticleByParticle(MCWalkerConfiguration& w,
    TrialWaveFunction& psi,
    QMCHamiltonian& h):
  QMCDriver(w,psi,h), nSubSteps(1), UseDrift("yes")
{
  RootName = "vmc";
  QMCType ="VMCParticleByParticle";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(nSubSteps,"subSteps","int");
  m_param.add(nSubSteps,"substeps","int");
}

bool VMCParticleByParticle::run()
{
  bool useBlockWithDrift = (UseDrift == "yes");
  Estimators->reportHeader(AppendRun);
  Estimators->reset();
  IndexType block = 0;
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = sqrt(Tau);
  G.resize(W.getTotalNum());
  dG.resize(W.getTotalNum());
  L.resize(W.getTotalNum());
  dL.resize(W.getTotalNum());
  nAcceptTot = 0;
  nRejectTot = 0;
  do
  {
    nAccept = 0;
    nReject=0;
    nAllRejected = 0;
    Estimators->startBlock();
    if(useBlockWithDrift)
      runBlockWithDrift();
    else
      runBlockWithoutDrift();
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    ++block;
    //record the current configuration
    recordBlock(block);
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "Ratio = "
            << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
            << std::endl;
  //finalize a qmc section
  return finalize(block);
}

void VMCParticleByParticle::runBlockWithDrift()
{
  IndexType step = 0;
  MCWalkerConfiguration::iterator it_end(W.end());
  do
  {
    //advanceWalkerByWalker();
    MCWalkerConfiguration::iterator it(W.begin());
    while(it != it_end)
    {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);
      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);
      RealType psi_old = thisWalker.Properties(SIGN);
      RealType psi = psi_old;
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      bool moved = false;
      for(int iat=0; iat<W.getTotalNum(); iat++)
      {
        PosType dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
        PosType newpos = W.makeMove(iat,dr);
        //RealType ratio = Psi.ratio(W,iat);
        RealType ratio = Psi.ratio(W,iat,dG,dL);
        G = W.G+dG;
        //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
        //dr = (*it)->R[iat]-newpos-Tau*G[iat];
        //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
        RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
        RealType scale=getDriftScale(Tau,G);
        //COMPLEX WARNING
        //dr = thisWalker.R[iat]-newpos-scale*G[iat];
        dr = thisWalker.R[iat]-newpos-scale*real(G[iat]);
        RealType logGb = -m_oneover2tau*dot(dr,dr);
        RealType prob = std::min(1.0e0,ratio*ratio*exp(logGb-logGf));
        //alternatively
        if(Random() < prob)
        {
          moved = true;
          ++nAccept;
          W.acceptMove(iat);
          Psi.acceptMove(W,iat);
          W.G = G;
          W.L += dL;
          //thisWalker.Drift = scale*G;
          assignDrift(scale,G,thisWalker.Drift);
        }
        else
        {
          ++nReject;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
      }
      if(moved)
      {
        w_buffer.rewind();
        W.copyToBuffer(w_buffer);
        psi = Psi.evaluate(W,w_buffer);
        thisWalker.R = W.R;
        RealType eloc=H.evaluate(W);
        thisWalker.resetProperty(log(std::abs(psi)), psi,eloc);
        H.saveProperty(thisWalker.getPropertyBase());
      }
      else
      {
        ++nAllRejected;
      }
      ++it;
    }
    ++step;
    ++CurrentStep;
    Estimators->accumulate(W);
    //if(CurrentStep%100 == 0) updateWalkers();
  }
  while(step<nSteps);
}

void VMCParticleByParticle::runBlockWithoutDrift()
{
  IndexType step = 0;
  MCWalkerConfiguration::iterator it_end(W.end());
  do
  {
    //advanceWalkerByWalker();
    MCWalkerConfiguration::iterator it(W.begin());
    while(it != it_end)
    {
      //MCWalkerConfiguration::WalkerData_t& w_buffer = *(W.DataSet[iwalker]);
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);
      W.R = thisWalker.R;
      w_buffer.rewind();
      W.copyFromBuffer(w_buffer);
      Psi.copyFromBuffer(W,w_buffer);
      RealType psi_old = thisWalker.Properties(SIGN);
      RealType psi = psi_old;
      for(int iter=0; iter<nSubSteps; iter++)
      {
        makeGaussRandom(deltaR);
        bool stucked=true;
        for(int iat=0; iat<W.getTotalNum(); iat++)
        {
          PosType dr = m_sqrttau*deltaR[iat];
          PosType newpos = W.makeMove(iat,dr);
          RealType ratio = Psi.ratio(W,iat);
          //RealType ratio = Psi.ratio(W,iat,dG,dL);
          RealType prob = std::min(1.0e0,ratio*ratio);
          //alternatively
          if(Random() < prob)
          {
            stucked=false;
            ++nAccept;
            W.acceptMove(iat);
            Psi.acceptMove(W,iat);
            //W.G+=dG;
            //W.L+=dL;
          }
          else
          {
            ++nReject;
            W.rejectMove(iat);
            Psi.rejectMove(iat);
          }
        }
        if(stucked)
        {
          ++nAllRejected;
        }
      }
      thisWalker.R = W.R;
      w_buffer.rewind();
      W.updateBuffer(w_buffer);
      RealType logpsi = Psi.updateBuffer(W,w_buffer);
      //W.copyToBuffer(w_buffer);
      //RealType logpsi = Psi.evaluate(W,w_buffer);
      RealType eloc=H.evaluate(W);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
      H.saveProperty(thisWalker.getPropertyBase());
      ++it;
    }
    ++step;
    ++CurrentStep;
    Estimators->accumulate(W);
  }
  while(step<nSteps);
}

bool
VMCParticleByParticle::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}


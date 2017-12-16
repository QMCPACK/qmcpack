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
    
    





#include "QMCDrivers/CorrelatedSampling/CSVMCUpdatePbyP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
//#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
CSVMCUpdatePbyP::CSVMCUpdatePbyP(MCWalkerConfiguration& w,
                                 TrialWaveFunction& psi,
                                 QMCHamiltonian& h,
                                 RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{
}

CSVMCUpdatePbyP::~CSVMCUpdatePbyP() { }

void CSVMCUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  int iwalker=0;
  //only used locally
  std::vector<RealType> ratio(nPsi), uw(nPsi);
  while(it != it_end)
  {
    //Walkers loop
    Walker_t& thisWalker(**it);
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    W.R = thisWalker.R;
    w_buffer.rewind();
    // Copy walker info in W
    W.copyFromBuffer(w_buffer);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      // Copy wave function info in W and Psi1
      Psi1[ipsi]->copyFromBuffer(W,w_buffer);
      Psi1[ipsi]->G=W.G;
      Psi1[ipsi]->L=W.L;
    }
    makeGaussRandomWithEngine(deltaR,RandomGen);
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      uw[ipsi]= thisWalker.Properties(ipsi,UMBRELLAWEIGHT);
    }
    // Point to the correct walker in the ratioij buffer
    RealType* restrict ratioijPtr=ratioIJ[iwalker];
    bool moved = false;
    for(int iat=0; iat<W.getTotalNum(); iat++)
      //Particles loop
    {
      PosType dr = m_sqrttau*deltaR[iat];
      PosType newpos = W.makeMove(iat,dr);
      RealType ratio_check=1.0;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        // Compute ratios before and after the move
        ratio_check *= ratio[ipsi] = Psi1[ipsi]->ratio(W,iat,*G1[ipsi],*L1[ipsi]);
        logpsi[ipsi]=std::log(ratio[ipsi]*ratio[ipsi]);
        // Compute Gradient in new position
        //*G1[ipsi]=Psi1[ipsi]->G + dG;
        // Initialize: sumratio[i]=(Psi[i]/Psi[i])^2=1.0
        sumratio[ipsi]=1.0;
      }
      bool accept_move=false;
      if(ratio_check>1e-12)//if any ratio is too small, reject the move automatically
      {
        int indexij(0);
        // Compute new (Psi[i]/Psi[j])^2 and their sum
        for(int ipsi=0; ipsi< nPsi-1; ipsi++)
        {
          for(int jpsi=ipsi+1; jpsi < nPsi; jpsi++, indexij++)
          {
            // Ratio between norms is already included in ratioijPtr from initialize.
            RealType rji=std::exp(logpsi[jpsi]-logpsi[ipsi])*ratioijPtr[indexij];
            instRij[indexij]=rji;
            //ratioij[indexij]=rji;
            sumratio[ipsi] += rji;
            sumratio[jpsi] += 1.0/rji;
          }
        }
        // Evaluate new Umbrella Weight
        for(int ipsi=0; ipsi< nPsi ; ipsi++)
          invsumratio[ipsi]=1.0/sumratio[ipsi];
        RealType td=ratio[0]*ratio[0]*sumratio[0]/(*it)->Multiplicity;
        accept_move=Random()<std::min(1.0,td);
      }
      //RealType prob = std::min(1.0,td);
      //if(Random() < prob)
      if(accept_move)
      {
        /* Electron move is accepted. Update:
           -ratio (Psi[i]/Psi[j])^2 for this walker
           -Gradient and laplacian for each Psi1[i]
           -Drift
           -buffered info for each Psi1[i]*/
        moved = true;
        ++nAccept;
        W.acceptMove(iat);
        // Update Buffer for (Psi[i]/Psi[j])^2
        copy(instRij.begin(),instRij.end(),ratioijPtr);
        // copy new Umbrella weight for averages
        uw=invsumratio;
        // Store sumratio for next Accept/Reject step to Multiplicity
        //thisWalker.Properties(SUMRATIO)=sumratio[0];
        thisWalker.Multiplicity=sumratio[0];
        for(int ipsi=0; ipsi< nPsi; ipsi++)
        {
          //Update local Psi1[i] buffer for the next move
          Psi1[ipsi]->acceptMove(W,iat);
          //Update G and L in Psi1[i]
          //Psi1[ipsi]->G = *G1[ipsi];
          Psi1[ipsi]->G += *G1[ipsi];
          Psi1[ipsi]->L += *L1[ipsi];
          thisWalker.Properties(ipsi,LOGPSI)+=std::log(std::abs(ratio[ipsi]));
        }
      }
      else
      {
        ++nReject;
        W.rejectMove(iat);
        for(int ipsi=0; ipsi< nPsi; ipsi++)
          Psi1[ipsi]->rejectMove(iat);
      }
    }
    if(moved)
    {
      /* The walker moved: Info are copied back to buffers:
         -copy (Psi[i]/Psi[j])^2 to ratioijBuffer
         -Gradient and laplacian for each Psi1[i]
         -Drift
         -buffered info for each Psi1[i]
         Physical properties are updated */
      (*it)->Age=0;
      (*it)->R = W.R;
      w_buffer.rewind();
      W.copyToBuffer(w_buffer);
      for(int ipsi=0; ipsi< nPsi; ipsi++)
      {
        W.G=Psi1[ipsi]->G;
        W.L=Psi1[ipsi]->L;
        //ValueType psi = Psi1[ipsi]->evaluate(W,w_buffer);
        ValueType logpsi = Psi1[ipsi]->evaluateLog(W,w_buffer);
        RealType et = H1[ipsi]->evaluate(W);
        //multiEstimator->updateSample(iwalker,ipsi,et,UmbrellaWeight[ipsi]);
        //Properties is used for UmbrellaWeight and UmbrellaEnergy
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=uw[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
        H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
    }
    else
    {
      ++nAllRejected;
    }
    ++it;
    ++iwalker;
  }
}
}


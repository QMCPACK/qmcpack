//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
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
                                 std::vector<TrialWaveFunction*>& psi,
                                 std::vector<QMCHamiltonian*>& h,
                                 RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{
}

CSVMCUpdatePbyP::~CSVMCUpdatePbyP() { }

void CSVMCUpdatePbyP::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  W.loadWalker(thisWalker,true);

  //First step, we initialize all Psis, and read up the value of logpsi
  //from the last run.
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    Psi1[ipsi]->copyFromBuffer(W,thisWalker.DataSet);
    logpsi[ipsi] = thisWalker.Properties(ipsi,LOGPSI);
  }
 
  //Now we compute sumratio and more importantly, ratioij.
  computeSumRatio(logpsi,avgNorm,RatioIJ,sumratio);
  RealType r(1.0); //a temporary variable for storing ratio^2.
 // myTimers[1]->start();
  for (int iter=0; iter<nSubSteps; ++iter)
  {
    makeGaussRandomWithEngine(deltaR,RandomGen);
    bool stucked=true;
    for(int ig=0; ig<W.groups(); ++ig) //loop over species
    {
      RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);
      for (int iat=W.first(ig); iat<W.last(ig); ++iat)
      {
        W.setActive(iat);
        mPosType dr = sqrttau*deltaR[iat];
        //The move proposal for particle iat. 
        if (W.makeMoveAndCheck(iat,dr))
        {
          for(int ipsi=0; ipsi<nPsi; ipsi++)
          {
            r=Psi1[ipsi]->ratio(W,iat);  
            ratio[ipsi]=r*r;
          }
          //Compute the ratio and acceptance probability.
          RealType prob=0;
          for(int ipsi=0; ipsi<nPsi; ipsi++)
            prob+=ratio[ipsi]/sumratio[ipsi];
   
          if (RandomGen() < prob)
          {
            stucked=false;
            ++nAccept;
            for(int ipsi=0; ipsi<nPsi; ipsi++)
              Psi1[ipsi]->acceptMove(W,iat);
           
            W.acceptMove(iat);
            //Now we update ratioIJ.
            updateRatioMatrix(ratio,RatioIJ);
            computeSumRatio(RatioIJ,sumratio);
            
          }
          else
          {
            ++nReject;
            W.rejectMove(iat);
            for(int ipsi=0; ipsi<nPsi; ipsi++)
              Psi1[ipsi]->rejectMove(iat);
          }
        }
        else //reject illegal moves
          ++nReject;
      } //iat
    }//ig for the species
    if (stucked)
    {
      ++nAllRejected;
    }
  }
//  myTimers[1]->stop();
//  myTimers[2]->start();

  W.donePbyP();

  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    //now recompute logpsi, gradients, and laplacians of the new R. Save them.
    logpsi[ipsi]=Psi1[ipsi]->updateBuffer(W,thisWalker.DataSet,recompute);
    W.saveWalker(thisWalker);
    //Save G and L for this wavefunction in the working G1, L1 arrays.
    *G1[ipsi]=thisWalker.G;
    *L1[ipsi]=thisWalker.L;
    //Set G and L for this wavefunction.
    Psi1[ipsi]->G = W.G;
    Psi1[ipsi]->L = W.L;
  }

  computeSumRatio(logpsi,avgNorm,sumratio);
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      invsumratio[ipsi]=sumratio[ipsi];
      cumNorm[ipsi]+=1.0/invsumratio[ipsi];
    }
  
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    //Now reload the G and L associated with ipsi into the walker and particle set.
    //This is required for correct calculation of kinetic energy.
    thisWalker.G=*G1[ipsi];
    thisWalker.L=*L1[ipsi];
    W.L=thisWalker.L;
    W.G=thisWalker.G;
    thisWalker.Properties(ipsi,LOCALENERGY)=H1[ipsi]->evaluate(W);
    thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
    thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
    thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=1.0/sumratio[ipsi];
    //Use Multiplicity as a temporary container for sumratio.
    thisWalker.Multiplicity=sumratio[0];
    H1[ipsi]->auxHevaluate(W,thisWalker);
    H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
  }
}


/// UpdatePbyP With Drift Fast.
CSVMCUpdatePbyPWithDriftFast::CSVMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w,
                                 std::vector<TrialWaveFunction*>& psi,
                                 std::vector<QMCHamiltonian*>& h,
                                 RandomGenerator_t& rg):
  CSUpdateBase(w,psi,h,rg)
{
  APP_ABORT("CSVMCUpdatePbyPWithDriftFast currently not working.  Please eliminate \
             drift option, or choose all electron moves instead.")
}

CSVMCUpdatePbyPWithDriftFast::~CSVMCUpdatePbyPWithDriftFast() { }

void CSVMCUpdatePbyPWithDriftFast::advanceWalker(Walker_t& thisWalker, bool recompute)
{
}


}


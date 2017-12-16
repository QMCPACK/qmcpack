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
    
    


#include "QMCDrivers/VMC/VMCPbyPMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommCreate.h"
#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
VMCPbyPMultiple::VMCPbyPMultiple(MCWalkerConfiguration& w,
                                 TrialWaveFunction& psi,
                                 QMCHamiltonian& h):
  QMCDriver(w,psi,h)
{
  RootName = "vmc";
  QMCType ="vmc";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_MULTIPLE,1);
  equilBlocks=-1;
  useDriftOpt="no";
  useDrift=false;
  m_param.add(equilBlocks,"equilBlocks","int");
  m_param.add(useDriftOpt,"useDrift","string");
  m_param.add(useDriftOpt,"usedrift","string");
  add_H_and_Psi(&h,&psi);
}

VMCPbyPMultiple::~VMCPbyPMultiple()
{
  for(int i=0; i<G.size(); i++)
    delete G[i];
  for(int i=0; i<dL.size(); i++)
    delete dL[i];
}

void VMCPbyPMultiple::resize(int ncopy, int nptcls)
{
  int m=ncopy*(ncopy-1)/2;
  ratio.resize(ncopy);
  logpsi2.resize(ncopy);
  UmbrellaWeight.resize(ncopy);
  invsumratio.resize(ncopy);
  sumratio.resize(ncopy);
  ratioij.resize(m);
  for(int i=0; i<ncopy; i++)
  {
    G.push_back(new ParticleGradient_t(nptcls));
    dL.push_back(new ParticleLaplacian_t(nptcls));
  }
}
bool VMCPbyPMultiple::run()
{
  useDrift = (useDriftOpt=="yes");
  if(useDrift)
    app_log() << "  VMCPbyPMultiple::run useDrift=yes" << std::endl;
  else
    app_log() << "  VMCPbyPMultiple::run useDrift=no" << std::endl;
  //TEST CACHE
  //Estimators->reportHeader(AppendRun);
  //going to add routines to calculate how much we need
  bool require_register =  W.createAuxDataSet();
  std::vector<RealType>Norm(nPsi),tmpNorm(nPsi);
  if(equilBlocks > 0)
  {
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      Norm[ipsi]=1.0;
      tmpNorm[ipsi]=0.0;
    }
  }
  else
  {
    for(int ipsi=0; ipsi< nPsi; ipsi++)
      Norm[ipsi]=std::exp(branchEngine->LogNorm[ipsi]);
  }
  multiEstimator->initialize(W,H1,Psi1,Tau,Norm,require_register);
  //TEST CACHE
  //Estimators->reset();
  Estimators->start(nBlocks);
  //TEST CACHE
  IndexType block = 0;
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = std::sqrt(Tau);
  RealType nPsi_minus_one = nPsi-1;
  ParticleSet::ParticleGradient_t dG(W.getTotalNum());
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  MCWalkerConfiguration::iterator it;
  MCWalkerConfiguration::iterator it_end(W.end());
  do
    //Blocks loop
  {
    IndexType step = 0;
    nAccept = 0;
    nReject=0;
    IndexType nAllRejected = 0;
    Estimators->startBlock(nSteps);
    do
      //Steps loop
    {
      it = W.begin();
      int iwalker=0;
      while(it != it_end)
        //Walkers loop
      {
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
        // Point to the correct walker in the ratioij buffer
        RealType *ratioijPtr=multiEstimator->RatioIJ[iwalker];
        //This is not used
        //ValueType psi_old = thisWalker.Properties(SIGN);
        //ValueType psi = psi_old;
        //create a 3N-Dimensional Gaussian with variance=1
        makeGaussRandom(deltaR);
        bool moved = false;
        for(int iat=0; iat<W.getTotalNum(); iat++)
          //Particles loop
        {
          PosType dr = m_sqrttau*deltaR[iat];
          if(useDrift)
            dr += thisWalker.Drift[iat];
          PosType newpos = W.makeMove(iat,dr);
          for(int ipsi=0; ipsi<nPsi; ipsi++)
          {
            // Compute ratios before and after the move
            ratio[ipsi] = Psi1[ipsi]->ratio(W,iat,dG,*dL[ipsi]);
            logpsi2[ipsi]=std::log(ratio[ipsi]*ratio[ipsi]);
            // Compute Gradient in new position
            *G[ipsi]=Psi1[ipsi]->G + dG;
            // Initialize: sumratio[i]=(Psi[i]/Psi[i])^2=1.0
            sumratio[ipsi]=1.0;
          }
          // Compute new (Psi[i]/Psi[j])^2 and their sum
          int indexij(0);
          for(int ipsi=0; ipsi< nPsi_minus_one; ipsi++)
          {
            for(int jpsi=ipsi+1; jpsi < nPsi; jpsi++, indexij++)
            {
              // Ratio between norms is already included in ratioijPtr from initialize.
              RealType rji=std::exp(logpsi2[jpsi]-logpsi2[ipsi])*ratioijPtr[indexij];
              ratioij[indexij]=rji;
              sumratio[ipsi] += rji;
              sumratio[jpsi] += 1.0/rji;
            }
          }
          // Evaluate new Umbrella Weight
          for(int ipsi=0; ipsi< nPsi ; ipsi++)
          {
            invsumratio[ipsi]=1.0/sumratio[ipsi];
          }
          RealType td=ratio[0]*ratio[0]*sumratio[0]/(*it)->Multiplicity;
          if(useDrift)
          {
            // Evaluate new drift
            PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
            for(int ipsi=1; ipsi< nPsi ; ipsi++)
            {
              PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
            }
            setScaledDrift(Tau,drift);
            RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
            dr = thisWalker.R[iat]-newpos-drift[iat];
            RealType logGb = -m_oneover2tau*dot(dr,dr);
            td *=std::exp(logGb-logGf);
          }
          // td = Target Density ratio
          //RealType td=pow(ratio[0],2)*sumratio[0]/(*it)->Properties(SUMRATIO);
          //td=ratio[0]*ratio[0]*sumratio[0]/(*it)->Multiplicity;
          //RealType prob = std::min(1.0,td*exp(logGb-logGf));
          RealType prob = std::min(1.0,td);
          if(Random() < prob)
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
            copy(ratioij.begin(),ratioij.end(),ratioijPtr);
            // Update Umbrella weight
            UmbrellaWeight=invsumratio;
            // Store sumratio for next Accept/Reject step to Multiplicity
            //thisWalker.Properties(SUMRATIO)=sumratio[0];
            thisWalker.Multiplicity=sumratio[0];
            for(int ipsi=0; ipsi< nPsi; ipsi++)
            {
              //Update local Psi1[i] buffer for the next move
              Psi1[ipsi]->acceptMove(W,iat);
              //Update G and L in Psi1[i]
              Psi1[ipsi]->G = *G[ipsi];
              Psi1[ipsi]->L += *dL[ipsi];
              thisWalker.Properties(ipsi,LOGPSI)+=std::log(std::abs(ratio[ipsi]));
            }
            // Update Drift
            if(useDrift)
              (*it)->Drift = drift;
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
            ValueType psi = Psi1[ipsi]->evaluate(W,w_buffer);
            RealType et = H1[ipsi]->evaluate(W);
            //multiEstimator->updateSample(iwalker,ipsi,et,UmbrellaWeight[ipsi]);
            //Properties is used for UmbrellaWeight and UmbrellaEnergy
            thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=UmbrellaWeight[ipsi];
            thisWalker.Properties(ipsi,LOCALENERGY)=et;
            H1[ipsi]->auxHevaluate(W);
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
      ++step;
      ++CurrentStep;
      Estimators->accumulate(W);
    }
    while(step<nSteps);
    //Modify Norm.
    if(block < equilBlocks)
    {
      for(int ipsi=0; ipsi< nPsi; ipsi++)
      {
        tmpNorm[ipsi]+=multiEstimator->esum(ipsi,MultipleEnergyEstimator::WEIGHT_INDEX);
      }
      if(block==(equilBlocks-1) || block==(nBlocks-1))
      {
        RealType SumNorm(0.e0);
        for(int ipsi=0; ipsi< nPsi; ipsi++)
          SumNorm+=tmpNorm[ipsi];
        for(int ipsi=0; ipsi< nPsi; ipsi++)
        {
          Norm[ipsi]=tmpNorm[ipsi]/SumNorm;
          branchEngine->LogNorm[ipsi]=std::log(Norm[ipsi]);
        }
      }
    }
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    ++block;
    //record the current configuration
    recordBlock(block);
    //re-evaluate the ratio and update the Norm
    multiEstimator->initialize(W,H1,Psi1,Tau,Norm,false);
  }
  while(block<nBlocks);
  ////Need MPI-IO
  //app_log()
  //  << "Ratio = "
  //  << static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot)
  //  << std::endl;
  return finalize(block);
}


bool
VMCPbyPMultiple::put(xmlNodePtr q)
{
  nPsi=Psi1.size();
  resize(nPsi,W.getTotalNum());
  if(branchEngine->LogNorm.size()==0)
    branchEngine->LogNorm.resize(nPsi);
  if(equilBlocks>0)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      branchEngine->LogNorm[ipsi]=0.e0;
  }
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    H1[ipsi]->add2WalkerProperty(W);
  if(Estimators == 0)
  {
    Estimators = new EstimatorManagerBase(H);
    multiEstimator = new MultipleEnergyEstimator(H,nPsi);
    Estimators->add(multiEstimator,"elocal");
  }
  H1[0]->setPrimary(true);
  for(int ipsi=1; ipsi<nPsi; ipsi++)
  {
    H1[ipsi]->setPrimary(false);
  }
  return true;
}
}


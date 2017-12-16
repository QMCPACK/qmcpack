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
    
    
#include "QMCDrivers/VMC/VMCMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/Communicate.h"
#include "Estimators/MultipleEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{

/// Constructor.
VMCMultiple::VMCMultiple(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h), multiEstimator(0)
{
  RootName = "vmc";
  QMCType ="VMCMultiple";
  equilBlocks=-1;
  m_param.add(equilBlocks,"equilBlocks","int");
  std::cout << "EquilBlocks " << equilBlocks << std::endl;
  QMCDriverMode.set(QMC_MULTIPLE,1);
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  add_H_and_Psi(&h,&psi);
}

/** allocate internal data here before run() is called
 * @author SIMONE
 *
 * See QMCDriver::process
 */
bool VMCMultiple::put(xmlNodePtr q)
{
  nPsi=Psi1.size();
  logpsi.resize(nPsi);
  sumratio.resize(nPsi);
  invsumratio.resize(nPsi);
  Norm.resize(nPsi);
  if(branchEngine->LogNorm.size()==0)
    branchEngine->LogNorm.resize(nPsi);
  if(equilBlocks>0)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      branchEngine->LogNorm[ipsi]=0.e0;
  }
  for(int ipsi=0; ipsi<nPsi; ipsi++)
    H1[ipsi]->add2WalkerProperty(W);
  Estimators = branchEngine->getEstimatorManager();
  if(Estimators == 0)
  {
    Estimators = new EstimatorManagerBase(H);
    multiEstimator = new MultipleEnergyEstimator(H,nPsi);
    Estimators->add(multiEstimator,Estimators->MainEstimatorName);
    branchEngine->setEstimatorManager(Estimators);
  }
  app_log() << "Number of H and Psi " << nPsi << std::endl;
  H1[0]->setPrimary(true);
  for(int ipsi=1; ipsi<nPsi; ipsi++)
  {
    H1[ipsi]->setPrimary(false);
  }
  return true;
}

/** Run the VMCMultiple algorithm.
 *
 * Similar to VMC::run
 */
bool VMCMultiple::run()
{
  //TEST CACHE
  //Estimators->reportHeader(AppendRun);
  bool require_register=false;
  //Check if we need to update the norm of the wave functions
  std::vector<RealType> tmpNorm(nPsi);
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
  //this is where the first values are evaulated
  multiEstimator->initialize(W,H1,Psi1,Tau,Norm,require_register);
  //TEST CACHE
  //Estimators->reset();
  Estimators->start(nBlocks);
  //TEST CACHE
  IndexType block = 0;
  double wh=0.0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  do
  {
    IndexType step = 0;
    nAccept = 0;
    nReject=0;
    Estimators->startBlock(nSteps);
    do
    {
      advanceWalkerByWalker();
      step++;
      CurrentStep++;
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
    Estimators->stopBlock(nAccept/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    nAccept = 0;
    nReject = 0;
    block++;
    //record the current configuration
    recordBlock(block);
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "Ratio = "
            << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
            << std::endl;
  //finalize a qmc section
  return finalize(block);
}

/**  Advance all the walkers one timstep.
 */
void
VMCMultiple::advanceWalkerByWalker()
{
  m_oneover2tau = 0.5/Tau;
  m_sqrttau = std::sqrt(Tau);
  //MCWalkerConfiguration::PropertyContainer_t Properties;
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int iwlk(0);
  int nPsi_minus_one(nPsi-1);
  while(it != it_end)
  {
    MCWalkerConfiguration::Walker_t &thisWalker(**it);
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(deltaR);
    W.R = m_sqrttau*deltaR + thisWalker.R + thisWalker.Drift;
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
        RealType ratioij=Norm[ipsi]/Norm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
    }
    for(int ipsi=0; ipsi< nPsi; ipsi++)
      invsumratio[ipsi]=1.0/sumratio[ipsi];
    // Only these properties need to be updated
    // Using the sum of the ratio Psi^2[j]/Psi^2[iwref]
    // because these are number of order 1. Potentially
    // the sum of Psi^2[j] can get very big
    //Properties(LOGPSI) =logpsi[0];
    //Properties(SUMRATIO) = sumratio[0];
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    //RealType scale = Tau; // ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
    //accumulate the weighted drift: using operators in ParticleBase/ParticleAttribOps.h
    //to handle complex-to-real assignment and axpy operations.
    //drift = invsumratio[0]*Psi1[0]->G;
    //for(int ipsi=1; ipsi< nPsi ;ipsi++) {
    //  drift += invsumratio[ipsi]*Psi1[ipsi]->G;
    //}
    PAOps<RealType,DIM>::scale(invsumratio[0],Psi1[0]->G,drift);
    for(int ipsi=1; ipsi< nPsi ; ipsi++)
    {
      PAOps<RealType,DIM>::axpy(invsumratio[ipsi],Psi1[ipsi]->G,drift);
    }
    //drift *= scale;
    setScaledDrift(Tau,drift);
    deltaR = thisWalker.R - W.R - drift;
    RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
    //Original
    //RealType g = Properties(SUMRATIO)/thisWalker.Properties(SUMRATIO)*
    //	exp(logGb-logGf+2.0*(Properties(LOGPSI)-thisWalker.Properties(LOGPSI)));
    //Reuse Multiplicity to store the sumratio[0]
    RealType g = sumratio[0]/thisWalker.Multiplicity*
                 std::exp(logGb-logGf+2.0*(logpsi[0]-thisWalker.Properties(LOGPSI)));
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
}


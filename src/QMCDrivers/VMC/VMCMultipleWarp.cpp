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
    
    


#include "QMCDrivers/VMC/VMCMultipleWarp.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"
#include "Estimators/MultipleEnergyEstimator.h"
#include "Particle/DistanceTable.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus
{


/// Constructor.
VMCMultipleWarp::VMCMultipleWarp(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, ParticleSetPool& ptclPool):
  QMCDriver(w,psi,h), PtclPool(ptclPool), multiEstimator(0)
{
  RootName = "vmc-warp";
  QMCType ="vmc-warp";
  equilBlocks=-1;
  m_param.add(equilBlocks,"equilBlocks","int");
  refSetName="invalid";
  m_param.add(refSetName,"reference","str");
  QMCDriverMode.set(QMC_MULTIPLE,1);
  JACOBIAN=w.addProperty("Jacobian");
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  add_H_and_Psi(&h,&psi);
}

/** allocate internal data here before run() is called
 * @author SIMONE
 *
 * See QMCDriver::process
 */
bool VMCMultipleWarp::put(xmlNodePtr q)
{
  if(WW.empty())
  {
    W.clearDistanceTables();
  }
  //qmcsystem
  std::vector<ParticleSet*> ionSets;
  DistanceTableData* dtableReference;
  xmlNodePtr cur=q->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "qmcsystem")
    {
      std::string source_name((const char*)xmlGetProp(cur,(const xmlChar*)"source"));
      ionSets.push_back(PtclPool.getParticleSet(source_name));
    }
    cur=cur->next;
  }
  ParticleSet* p(0);
  if(refSetName!="invalid")
  {
    p=PtclPool.getParticleSet(refSetName);
    if(p==0)
    {
      std::cout << "The specified reference cannot be found. Stop." << std::endl;
      abort();
    }
  }
  else
  {
    refSetName=ionSets[0]->getName().c_str();
    p=PtclPool.getParticleSet(refSetName);
  }
  dtableReference=DistanceTable::add(*p,W);
  nptcl=W.R.size();
  nPsi=Psi1.size();
  PtclWarp.initialize(ionSets,dtableReference);
  logpsi.resize(nPsi);
  sumratio.resize(nPsi);
  invsumratio.resize(nPsi);
  Jacobian.resize(nPsi);
  Norm.resize(nPsi);
  if(branchEngine->LogNorm.size()==0)
    branchEngine->LogNorm.resize(nPsi);
  if(equilBlocks>0)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      branchEngine->LogNorm[ipsi]=0.e0;
  }
  //for(int ipsi=0; ipsi<nPsi; ipsi++)
  //  H1[ipsi]->add2WalkerProperty(W);
  if(Estimators == 0)
  {
    Estimators = new EstimatorManagerBase(myComm);
    multiEstimator = new MultipleEnergyEstimator(H,nPsi);
    Estimators->add(multiEstimator,"elocal");
  }
  LOGMSG("Number of H and Psi " << nPsi)
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    H1[ipsi]->setPrimary(true);
  }
  if(WW.empty())
  {
    //WW.push_back(&W);
    char newname[128];
    //for(int ipsi=1; ipsi<nPsi; ipsi++){
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      sprintf(newname,"%s%d", W.getName().c_str(),ipsi);
      ParticleSet* pclone=PtclPool.getParticleSet(newname);
      if(pclone == 0)
      {
        app_log() << "  Cloning particle set in VMCMultipleWarp " << newname << std::endl;
        pclone=new ParticleSet(W);
        pclone->setName(newname);
        pclone->clearDistanceTables(); // NECESSARY???????
        PtclPool.addParticleSet(pclone);
      }
      else
      {
        app_log() << "  Cloned particle exists " << newname << std::endl;
      }
      WW.push_back(pclone);
      Psi1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
      H1[ipsi]->resetTargetParticleSet(*WW[ipsi]);
    }
  }
  return true;
}

/** Run the VMCMultipleWarp algorithm.
 *
 * Similar to VMC::run
 */
bool VMCMultipleWarp::run()
{
  int JACCOL=Estimators->addProperty("LogJacob");
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
  // Have to generalize to WW
  multiEstimator->initialize(W,WW,PtclWarp,H1,Psi1,Tau,Norm,require_register);
  //TEST CACHE
  //Estimators->reset();
  Estimators->start(nBlocks);
  //TEST CACHE
  IndexType block = 0;
  Pooma::Clock timer;
  double wh=0.0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  //bool appendwalker=pStride>0;
  do
  {
    IndexType step = 0;
    nAccept = 0;
    nReject=0;
    Estimators->startBlock(nSteps);
    RealType Jacblk=0.e0;
    do
    {
      advanceWalkerByWalker();
      step++;
      CurrentStep++;
      Estimators->accumulate(W);
      Jacblk+=std::log(std::abs((*W.begin())->Properties(1,JACOBIAN)));
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
        std::cout << "BranchEngine is updated " << branchEngine->LogNorm[0] << " " << branchEngine->LogNorm[1] << std::endl;
      }
    }
    Estimators->stopBlock(nAccept/static_cast<RealType>(nAccept+nReject));
    RealType AveJacobLog=Jacblk/static_cast<RealType>(nSteps);
    Estimators->setProperty(JACCOL,AveJacobLog);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    branchEngine->LogJacobRef+=AveJacobLog;
    nAccept = 0;
    nReject = 0;
    block++;
    //record the current configuration
    recordBlock(block);
  }
  while(block<nBlocks);
  branchEngine->LogJacobRef/=static_cast<RealType>(nBlocks);
  // do {
  //   IndexType step = 0;
  //   timer.start();
  //   nAccept = 0; nReject=0;
  //   do {
  //     advanceWalkerByWalker();
  //     step++;CurrentStep++;
  //     Estimators->accumulate(W);
  //   } while(step<nSteps);
  //   timer.stop();
  //   nAcceptTot += nAccept; nRejectTot += nReject;
  //   Estimators->flush();
  //   RealType TotalConfig=static_cast<RealType>(nAccept+nReject);
  //   Estimators->setColumn(AcceptIndex,nAccept/TotalConfig);
  //   Estimators->report(CurrentStep);
  //
  //   LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << std::endl;
  //   //HDFWalkerOutput WO(RootName,block&&appendwalker, block);
  //   //WO.get(W);
  //   branchEngine->accumulate(Estimators->average(0),1.0);
  //   nAccept = 0; nReject = 0;
  //   block++;
  //   //record the current configuration
  //   recordWalkerConfigurations(block);
  // } while(block<nBlocks);
  app_log()
      << "Ratio = "
      << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
      << std::endl;
  //finalize a qmc section
  return finalize(block);
}

/**  Advance all the walkers one timstep.
 */
void
VMCMultipleWarp::advanceWalkerByWalker()
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
    W.loadWalker(thisWalker,false);
    //CHEAT!!! Tau should tau/mass
    setScaledDrift(Tau,W.G,drift);
    makeGaussRandom(deltaR);
    W.R = m_sqrttau*deltaR + thisWalker.R + drift;
    W.update();
    for(int ipsi=0; ipsi<nPsi; ipsi++)
      Jacobian[ipsi]=1.e0;
    for(int iptcl=0; iptcl< nptcl; iptcl++)
    {
      PtclWarp.warp_one(iptcl,nptcl);
      //Save particle position
      for(int ipsi=1; ipsi<nPsi; ipsi++)
      {
        WW[ipsi]->R[iptcl]=W.R[iptcl]+PtclWarp.get_displacement(iptcl,ipsi);
        Jacobian[ipsi]*=PtclWarp.get_Jacobian(iptcl,ipsi);
      }
    }
    for(int ipsi=0; ipsi< nPsi; ipsi++)
    {
      WW[ipsi]->update();
      logpsi[ipsi]=Psi1[ipsi]->evaluateLog(*WW[ipsi]);
      //Redundant???
      Psi1[ipsi]->L=WW[ipsi]->L;
      Psi1[ipsi]->G=WW[ipsi]->G;
      sumratio[ipsi]=1.0;
    }
    // Compute the sum over j of Psi^2[j]/Psi^2[i] for each i
    for(int ipsi=0; ipsi< nPsi_minus_one; ipsi++)
    {
      for(int jpsi=ipsi+1; jpsi< nPsi; jpsi++)
      {
        RealType ratioij=Norm[ipsi]/Norm[jpsi]*std::exp(2.0*(logpsi[jpsi]-logpsi[ipsi]));
        ratioij*=(Jacobian[jpsi]/Jacobian[ipsi]);
        sumratio[ipsi] += ratioij;
        sumratio[jpsi] += 1.0/ratioij;
      }
    }
    for(int ipsi=0; ipsi< nPsi; ipsi++)
      invsumratio[ipsi]=1.0/sumratio[ipsi];
    // std::cout << invsumratio[0] << " " << invsumratio[1] << std::endl;
    // Only these properties need to be updated
    // Using the sum of the ratio Psi^2[j]/Psi^2[iwref]
    // because these are number of order 1. Potentially
    // the sum of Psi^2[j] can get very big
    //Properties(LOGPSI) =logpsi[0];
    //Properties(SUMRATIO) = sumratio[0];
    RealType logGf = -0.5*Dot(deltaR,deltaR);
    RealType scale = Tau; // ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
    //accumulate the weighted drift
    //drift = Psi1[0]->G;
    /*drift = invsumratio[0]*Psi1[0]->G;
      for(int ipsi=1; ipsi< nPsi ;ipsi++) {
      drift += invsumratio[ipsi]*Warp.G[ipsi];
      }*/
    //drift *= scale;
    setScaledDrift(Tau,Psi1[0]->G,drift);
    drift*=invsumratio[0];
    for(int ipsi=1; ipsi< nPsi ; ipsi++)
    {
      setScaledDrift(Tau,Psi1[ipsi]->G,deltaR);
      drift+= (invsumratio[ipsi]*deltaR);
    }
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
      //cout << "REJECTED" << std::endl << std::endl;
      thisWalker.Age++;
      ++nReject;
    }
    else
    {
      thisWalker.Age=0;
      thisWalker.Multiplicity=sumratio[0];
      W.saveWalker(thisWalker);
      for(int ipsi=0; ipsi<nPsi; ipsi++)
      {
        WW[ipsi]->L=Psi1[ipsi]->L; //NECESSARY???????
        WW[ipsi]->G=Psi1[ipsi]->G;
        RealType et = H1[ipsi]->evaluate(*WW[ipsi]);
        thisWalker.Properties(ipsi,LOGPSI)=logpsi[ipsi];
        thisWalker.Properties(ipsi,SIGN) =Psi1[ipsi]->getPhase();
        thisWalker.Properties(ipsi,UMBRELLAWEIGHT)=invsumratio[ipsi];
        thisWalker.Properties(ipsi,LOCALENERGY)=et;
        thisWalker.Properties(ipsi,JACOBIAN)=Jacobian[ipsi];
        H1[ipsi]->auxHevaluate(W);
        H1[ipsi]->saveProperty(thisWalker.getPropertyBase(ipsi));
      }
      //cout << "ACCEPTED" << std::endl << std::endl;
      ++nAccept;
    }
    ++it;
    ++iwlk;
  }
  //cout << std::endl;
}
}


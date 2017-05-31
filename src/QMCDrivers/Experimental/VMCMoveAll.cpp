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
    
    
#include "QMC/VMCMoveAll.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Utilities/Clock.h"

namespace qmcplusplus
{

/// Constructor.
VMCMoveAll::VMCMoveAll(MCWalkerConfiguration& w,
                       TrialWaveFunction& psi,
                       QMCHamiltonian& h,
                       xmlNodePtr q):
  QMCDriver(w,psi,h,q)
{
  RootName = "vmc";
  QMCType ="vmc";
}

/** Run the VMCMoveAll algorithm.
 *
 * Advance the walkers nblocks*nsteps timesteps.
 * For each timestep:
 * <ul>
 * <li> Move all the particles of a walker.
 * <li> Calculate the properties for the new walker configuration.
 * <li> Accept/reject the new configuration.
 * <li> Accumulate the estimators.
 * </ul>
 * For each block:
 * <ul>
 * <li> Flush the estimators and print to file.
 * <li> (Optional) Print the ensemble of walker configurations.
 * </ul>
 *
 * Default mode: Print the ensemble of walker configurations
 * at the end of the run.
 */
bool VMCMoveAll::run()
{
  DistanceTable::create(W.getActiveWalkers());
  Psi.resizeByWalkers(W.getActiveWalkers());
  if(put(qmc_node))
  {
    H.setTau(0.0);
    //set the data members to start a new run
    getReady();
    //probably unnecessary
    MCWalkerConfiguration::iterator it = W.begin();
    MCWalkerConfiguration::iterator it_end = W.end();
    while(it != it_end)
    {
      (*it)->Properties(WEIGHT) = 1.0;
      ++it;
    }
    deltaR.resize(W.getTotalNum());
    drift.resize(W.getTotalNum());
    Estimators.reset();
    //create an output engine
    HDFWalkerOutput WO(RootName);
    IndexType block = 0;
    Pooma::Clock timer;
    double wh=0.0;
    IndexType nAcceptTot = 0;
    IndexType nRejectTot = 0;
    do
    {
      IndexType step = 0;
      timer.start();
      nAccept = 0;
      nReject=0;
      do
      {
        advanceAllWalkers();
        step++;
        CurrentStep++;
        Estimators.accumulate(W);
      }
      while(step<nSteps);
      timer.stop();
      nAcceptTot += nAccept;
      nRejectTot += nReject;
      Estimators.flush();
      Estimators.setColumn(AcceptIndex,
                           static_cast<double>(nAccept)/static_cast<double>(nAccept+nReject));
      Estimators.report(CurrentStep);
      LogOut->getStream() << "Block " << block << " " << timer.cpu_time() << std::endl;
      if(pStride)
        WO.get(W);
      nAccept = 0;
      nReject = 0;
      block++;
    }
    while(block<nBlocks);
    LogOut->getStream()
        << "Ratio = "
        << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
        << std::endl;
    if(!pStride)
      WO.get(W);
    Estimators.finalize();
    return true;
  }
  else
  {
    ERRORMSG("Error with Input")
    return false;
  }
}


bool
VMCMoveAll::put(xmlNodePtr q)
{
  xmlNodePtr qsave=q;
  bool success = putQMCInfo(q);
  success = Estimators.put(qsave);
  return success;
}

/**  Advance all the walkers simultaneously.
 *
 * Broken and does not help us at all
 */
void VMCMoveAll::advanceAllWalkers()
{
  // deltaR.resize(W.getTotalNum());
  // WalkerSetRef Wref(W);
  // Wref.resize(W.getActiveWalkers(),W.getTotalNum());
  //
  // //Pooma::Clock timer;
  // RealType oneovertau = 1.0/Tau;
  // RealType oneover2tau = 0.5*oneovertau;
  // RealType g = sqrt(Tau);
  //
  // MCWalkerConfiguration::PropertyContainer_t Properties;
  // makeGaussRandom(Wref.R);
  //
  // Wref.R *= g;
  //
  // int nptcl = W.getTotalNum();
  // int iw = 0;
  // MCWalkerConfiguration::iterator it = W.begin();
  // while(it !=  W.end()) {
  //   const ParticleSet::ParticlePos_t& r = (*it)->R;
  //   for(int jat=0; jat<nptcl; jat++) {
  //     Wref.R(iw,jat) += r(jat) + (*it)->Drift(jat);
  //   }
  //   iw++; it++;
  // }
  //
  // DistanceTable::update(Wref);
  //
  // OrbitalBase::ValueVectorType   logpsi(iw), energy(iw);
  //
  // //THIS IS TOTALLY BROKEN
  // Psi.evaluate(Wref,logpsi);
  //
  // H.evaluate(Wref,energy);
  //
  // //multiply tau to convert gradient to drift term
  // Wref.G *= Tau;
  //
  // iw = 0;
  // it = W.begin();
  // while(it !=  W.end()) {
  //
  //   ValueType eold = Properties(LOCALENERGY);
  //
  //   for(int iat=0; iat<nptcl; iat++)
  //     deltaR(iat) = Wref.R(iw,iat) - (*it)->R(iat) - (*it)->Drift(iat);
  //   //RealType forwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
  //   RealType logforwardGF = -oneover2tau*Dot(deltaR,deltaR);
  //
  //   for(int iat=0; iat<nptcl; iat++)
  //     deltaR(iat) = (*it)->R(iat) - Wref.R(iw,iat) - Wref.G(iw,iat);
  //
  //   ////RealType backwardGF = exp(-oneover2tau*Dot(deltaR,deltaR));
  //   RealType logbackwardGF = -oneover2tau*Dot(deltaR,deltaR);
  //
  //   RealType logpsi=psi(iw);
  //   RealType g=exp(logbackwardGF-logforwardGF+2.0*(logpsi-(*it)->Properties(LOGPSI)));
  //   //ValueType psisq = psi(iw)*psi(iw);
  //   if(Random() > g) {
  //     ++nReject;
  //     (*it)->Properties(AGE) += 1;
  //   } else {
  //     (*it)->Properties(AGE) = 0;
  //     for(int iat=0; iat<nptcl; iat++) (*it)->R(iat) = Wref.R(iw,iat);
  //     for(int iat=0; iat<nptcl; iat++) (*it)->Drift(iat) = Wref.G(iw,iat);
  //     (*it)->Properties(PSI) = psi(iw);
  //     (*it)->Properties(LOGPSI) = logpsi;//log(std::abs(psi));
  //     (*it)->Properties(LOCALENERGY) = energy(iw);
  //     ++nAccept;
  //   }
  //   iw++;it++;
  // }
}
}


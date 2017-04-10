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
    
    


#include "QMCDrivers/ReptationMC.h"
#include "QMCDrivers/PolymerChain.h"
#include "QMCDrivers/PolymerEstimator.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
namespace qmcplusplus
{

ReptationMC::ReptationMC(MCWalkerConfiguration& w,
                         TrialWaveFunction& psi,
                         QMCHamiltonian& h):
  QMCDriver(w,psi,h),
  UseBounce(false),
  ClonePolymer(true),
  PolymerLength(21),
  NumCuts(1),
  NumTurns(0),
  Reptile(0)
{
  RootName = "rmc";
  QMCType ="rmc";
  m_param.add(PolymerLength,"chains","int");
  m_param.add(NumCuts,"cuts","int");
  m_param.add(UseBounce,"bounce","int");
  m_param.add(ClonePolymer,"clone","int");
}

ReptationMC::~ReptationMC()
{
  if(Reptile)
    delete Reptile;
}

///initialize polymers
void ReptationMC::initReptile()
{
  //overwrite the number of cuts for Bounce algorithm
  if(UseBounce)
    NumCuts = 1;
  app_log() << "Moving " << NumCuts << " for each reptation step" << std::endl;
  //Reptile is NOT allocated. Create one.
  if(Reptile == 0)
  {
    MCWalkerConfiguration::iterator it(W.begin());
    Walker_t* cur(*it);
    cur->Weight=1.0;
    W.R = cur->R;
    //DistanceTable::update(W);
    W.update();
    RealType logpsi(Psi.evaluateLog(W));
    RealType  eloc_cur = H.evaluate(W);
    cur->resetProperty(logpsi,Psi.getPhase(),eloc_cur);
    H.saveProperty(cur->getPropertyBase());
    cur->Drift = W.G;
    Reptile  = new PolymerChain(cur,PolymerLength,NumCuts);
    Reptile->resizeArrays(1);
  }
  //If ClonePolyer==false, generate configuration using diffusion
  //Not so useful
  if(!ClonePolymer)
  {
    Walker_t* cur((*Reptile)[0]);
    RealType g = std::sqrt(Tau);
    for(int i=0; i<NumCuts-1; i++ )
    {
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandom(deltaR);
      W.R = cur->R + g*deltaR + Tau*cur->Drift;
      //update the distance table associated with W
      //DistanceTable::update(W);
      W.update();
      //evaluate wave function
      RealType logpsic(Psi.evaluateLog(W));
      RealType e0=H.evaluate(W);
      cur = (*Reptile)[i+1];
      cur->resetProperty(logpsic,Psi.getPhase(),e0);
      H.saveProperty(cur->getPropertyBase());
      cur->R = W.R;
      cur->Drift = W.G;
    }
  }
}

bool ReptationMC::run()
{
  Estimators->reportHeader(AppendRun);
  initReptile();
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  PolymerEstimator pe(*Reptile);
  pe.resetReportSettings(RootName);
  //accumulate configuration: probably need to reorder
  HDFWalkerOutput WO(RootName);
  do
  {
    IndexType step = 0;
    NumTurns = 0;
    Estimators->startBlock();
    do
    {
      moveReptile();
      step++;
      CurrentStep++;
      W.copyWalkerRefs(Reptile->front(),Reptile->back());
      Estimators->accumulate(W);
      pe.accumulate();
    }
    while(step<nSteps);
    Estimators->stopBlock(static_cast<RealType>(nAccept)/static_cast<RealType>(nAccept+nReject));
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    Estimators->report(CurrentStep);
    pe.report(CurrentStep);
    nAccept = 0;
    nReject = 0;
    block++;
  }
  while(block<nBlocks);
  //Need MPI-IO
  app_log() << "ratio = "
            << static_cast<double>(nAcceptTot)/static_cast<double>(nAcceptTot+nRejectTot)
            << std::endl;
  return finalize(block);
}

bool
ReptationMC::put(xmlNodePtr q)
{
  //nothing to do yet
  return true;
}

void
ReptationMC::moveReptile()
{
  //RealType oneovertau = 1.0/Tau;
  //RealType oneover2tau = 0.5*oneovertau;
  RealType tauover2 = 0.5*Tau;
  RealType g = std::sqrt(Tau);
  typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;
  if(!UseBounce && Random()<0.5)
  {
    Reptile->flip();
    NumTurns++;
  }
  Walker_t* anchor = Reptile->makeEnds();
  //save the local energies of the anchor and tails
  //eloc_xp = the energy of the front
  //eloc_yp = the energy of the proposed move
  //eloc_x = the energy of the tail
  //eloc_y = the energy of the tail-1
  RealType eloc_xp = anchor->Properties(LOCALENERGY);
  RealType eloc_x =  Reptile->tails[0]->Properties(LOCALENERGY);
  RealType eloc_y =  Reptile->tails[1]->Properties(LOCALENERGY);
  NumCuts = Reptile->NumCuts;
  RealType Wpolymer=0.0;
  for(int i=0; i<NumCuts; )
  {
    Walker_t* head=Reptile->heads[i];
    //create a 3N-Dimensional Gaussian with variance=1
    makeGaussRandom(deltaR);
    W.R = anchor->R + g*deltaR + Tau* anchor->Drift;
    //update the distance table associated with W
    //DistanceTable::update(W);
    W.update();
    //evaluate wave function
    RealType logpsi(Psi.evaluateLog(W));
    //update the properties of the front chain
    //RealType eloc_yp = head->Properties(LOCALENERGY) = H.evaluate(W);
    //H.copy(head->getEnergyBase());
    //head->Properties(LOCALPOTENTIAL) = H.getLocalPotential();
    RealType eloc_yp = H.evaluate(W);
    head->resetProperty(logpsi,Psi.getPhase(),eloc_yp);
    H.saveProperty(head->getPropertyBase());
    head->R = W.R;
    //ValueType vsq = Dot(W.G,W.G);
    //ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
    //head->Drift = scale*W.G;
    head->Drift = W.G;
    //\f${x-y-\tau\nabla \ln \Psi_{T}(y))\f$
    //deltaR = anchor->R - W.R - heads[i]->Drift;
    //Gdrift *= exp(-oneover2tau*Dot(deltaR,deltaR));
    /*
       \f$ X= \{R_0, R_1, ... , R_M\}\f$
       \f$ X' = \{R_1, .., R_M, R_{M+1}\}\f$
       \f[ G_B(R_{M+1}\leftarrow R_{M}, \tau)/G_B(R_{0}\leftarrow R_{1}, \tau)
       = exp\(-\tau/2[E_L(R_{M+1})+E_L(R_M)-E_L(R_1)-E_L(R_0)]\)\f]
       *
       -  eloc_yp = \f$E_L(R_{M+1})\f$
       -  eloc_xp = \f$E_L(R_{M})\f$
       -  eloc_y = \f$E_L(R_{1})\f$
       -  eloc_x = \f$E_L(R_{0})\f$
    */
    //Wpolymer *= exp(-oneover2tau*(eloc_yp+eloc_xp-eloc_x-eloc_y));
    Wpolymer +=(eloc_yp+eloc_xp-eloc_x-eloc_y);
    //move the anchor and swap the local energies for Wpolymer
    anchor=head;
    //increment the index
    i++;
    if(i<NumCuts)
    {
      eloc_xp  = eloc_yp;
      eloc_x = eloc_y;
      eloc_y = Reptile->tails[i+1]->Properties(LOCALENERGY);
    }
  }
  Wpolymer = std::exp(-tauover2*Wpolymer);
  double accept = std::min(1.0,Wpolymer);
  if(Random() < accept)
    //move accepted
  {
    Reptile->updateEnds();
    ++nAccept;
  }
  else
  {
    ++nReject;
    if(UseBounce)
    {
      NumTurns++;
      Reptile->flip();
    }
  }
  //RealType Bounce =  UseBounce ? 1.0-accept: 0.5;
  //if(Random()<Bounce) {
  //  Reptile->flip();
  //  LogOut->getStream() << "Bounce = " << Bounce << " " << NumTurns << " " << polymer.MoveHead << std::endl;
  //  NumTurns++;//increase the number of turns
  //}
}
}

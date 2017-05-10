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
    
    


#include "QMCDrivers/ZeroVarianceOptimize.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
//#include "QMCWaveFunctions/Jastrow/DerivPadeFunctors.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"

namespace qmcplusplus
{

/// Constructor.
ZeroVarianceOptimize::ZeroVarianceOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  QMCDriver(w,psi,h), WarmupBlocks(1), Mover(0), UseDrift("yes")
{
  RootName = "opt";
  QMCType ="ZeroVarianceOptimize";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
}

///destructor
ZeroVarianceOptimize::~ZeroVarianceOptimize()
{
  //nothing to do yet
}

bool ZeroVarianceOptimize::run()
{
  resetRun();
  Mover->startRun(nBlocks,true);
  //do the warmup blocks
  for(int block=0; block<WarmupBlocks; ++block)
  {
    for(int step=0; step<nSteps; ++step)
      Mover->advanceWalkers(W.begin(),W.end(),true);
  }
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  dLogPsi=0.0;
  dHPsi=0.0;
  Hessian=0.0;
  Overlap=0.0;
  for(int block=0; block<nBlocks; ++block)
  {
    Mover->startBlock(nSteps);
    for(int step=0; step<nSteps; ++step)
    {
      Mover->advanceWalkers(W.begin(),W.end(),true); //step==nSteps);
      Estimators->accumulate(W);
    }
    Mover->stopBlock();
    //accumulate Hessian and Overlap matrices
    std::accumulate(W.begin(),W.end());
    nAcceptTot += Mover->nAccept;
    nRejectTot += Mover->nReject;
    //periodically re-evaluate everything for pbyp
    if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0)
      Mover->updateWalkers(W.begin(),W.end());
  }
  Mover->stopRun();
  //cout << "Hessian matrix " << std::endl;
  //for(int i=0; i<NoptPlusOne; ++i)
  //{
  //for(int j=0; j<NoptPlusOne; ++j)
  //{
  //cout << std::setw(22) << Hessian(i,j);
  //}
  //cout << std::endl;
  //}
  //cout << "Overlap matrix " << std::endl;
  //for(int i=0; i<NoptPlusOne; ++i)
  //{
  //for(int j=0; j<NoptPlusOne; ++j)
  //{
  //cout << std::setw(22) << Overlap(i,j);
  //}
  //cout << std::endl;
  //}
  //finalize a qmc section
  return finalize(nBlocks);
}

void ZeroVarianceOptimize::accumulate(WalkerIter_t it, WalkerIter_t it_end)
{
  for(; it != it_end; ++it)
  {
    //Walker_t& thisWalker(**it);
    ////only need to update distance tables and structure factor
    //W.R = thisWalker.R;
    //if(QMCDriverMode[QMC_UPDATE_MODE])
    //{
    //Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    //w_buffer.rewind();
    //W.copyFromBuffer(w_buffer);
    //}
    //else
    //W.update();
    //dPsi[0]->dHPsi=thisWalker.Properties(LOCALENERGY);
    //RealType ke0=thisWalker.Properties(NUMPROPERTIES);
    //for(int i=1; i<NoptPlusOne; ++i)
    //dPsi[i]->evaluateDerivatives(W,ke0);
    //for(int i=0; i<NoptPlusOne; ++i)
    //{
    //RealType dlogpsi= dPsi[i]->dLogPsi;
    //for(int j=0; j<NoptPlusOne; ++j)
    //{
    //Hessian(i,j) += dlogpsi*dPsi[j]->dHPsi;
    //Overlap(i,j) += dlogpsi*dPsi[j]->dLogPsi;
    //}
    //}
  }
}

void ZeroVarianceOptimize::resetRun()
{
  if(Mover ==0)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      app_log() << "  Update particle by particle " << std::endl;
      if(UseDrift == "yes")
        Mover=new VMCUpdatePbyPWithDrift(W,Psi,H,Random);
      else
        Mover=new VMCUpdatePbyP(W,Psi,H,Random);
      Mover->resetRun(branchEngine,Estimators);
    }
    else
    {
      app_log() << "  Update walker by walker " << std::endl;
      if(UseDrift == "yes")
        Mover=new VMCUpdateAllWithDrift(W,Psi,H,Random);
      else
        Mover=new VMCUpdateAll(W,Psi,H,Random);
      Mover->resetRun(branchEngine,Estimators);
    }
    //push an dummy AnalyticDiffOrbital
    //dPsi.push_back(new AnalyticDiffOrbital(0));
    //RealType p0=Psi.VarList["jee_b"];
    //DPadeDBFunctor<RealType>* dpade = new DPadeDBFunctor<RealType>(-0.5,p0);
    //TwoBodyJastrowOrbital<DPadeDBFunctor<RealType> > *J2=new TwoBodyJastrowOrbital<DPadeDBFunctor<RealType> >(W);
    //dpade->ID_B="jee_b";
    //J2->insert("j2",dpade);
    //for(int i=0; i<4; i++) J2->addFunc(dpade);
    //DiffOrbitalBase *o= new AnalyticDiffOrbital(J2);
    //o->resize(W.getTotalNum());
    //o->setParameter("jee_b",p0);
    //dPsi.push_back(o);
    //NoptPlusOne=dPsi.size();
    //dLogPsi.resize(NoptPlusOne);
    //dHPsi.resize(NoptPlusOne);
    //Hessian.resize(NoptPlusOne,NoptPlusOne);
    //Overlap.resize(NoptPlusOne,NoptPlusOne);
  }
  if(QMCDriverMode[QMC_UPDATE_MODE])
    Mover->initWalkersForPbyP(W.begin(),W.end());
  else
    Mover->initWalkers(W.begin(),W.end());
  Mover->put(qmcNode);
}

bool
ZeroVarianceOptimize::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}

//Debug
//{
//  IndexType nskipped = 0;
//  RealType sig2Enloc=0, sig2Drift=0;
//  RealType delta = 0.0001;
//  RealType delta2 = 2*delta;
//  ValueType c1 = 1.0/delta/2.0;
//  ValueType c2 = 1.0/delta/delta;

//  int nat = W.getTotalNum();

//  ParticleSet::ParticlePos_t deltaR(nat);
//  MCWalkerConfiguration::PropertyContainer_t Properties;
//  //pick the first walker
//  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

//  //copy the properties of the working walker
//  Properties = awalker->Properties;
//
//  W.R = awalker->R;

//  W.update();
//  ValueType psi = Psi.evaluate(W);
//  RealType eloc=H.evaluate(W);

//  app_log() << "  HamTest " << "  Total " <<  eloc << std::endl;
//  for(int i=0; i<H.size(); i++)
//    app_log() << "  HamTest " << H.getName(i) << " " << H[i] << std::endl;


//  RealType p0=Psi.VarList["jee_b"];
//  ///testing with pade
//  typedef PadeFunctor<RealType> FuncType;
//  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
//  JeeType *J2 = new JeeType(W);
//  FuncType *func=new FuncType(-0.5,p0);
//  func->ID_B="jee_b";
//  J2->insert("j2",func);
//  for(int i=0; i<4; i++) J2->addFunc(func);

//  NumericalDiffOrbital dpsi(J2);
//  dpsi.resize(nat);
//  dpsi.setParameter("jee_b",p0);
//  dpsi.evaluateDerivatives(W,H[0]);

//  DPadeDBFunctor<RealType> dpade(-0.5,p0);
//  TwoBodyJastrowOrbital<DPadeDBFunctor<RealType> > J2_ana(W);
//  dpade.ID_B="jee_b";
//  J2_ana.insert("j2",&dpade);
//  for(int i=0; i<4; i++) J2_ana.addFunc(&dpade);

//  AnalyticDiffOrbital dpsiA(&J2_ana);
//  dpsiA.resize(nat);
//  dpsiA.setParameter("jee_b",p0);
//  dpsiA.evaluateDerivatives(W,H[0]);

//  OrbitalBase::OptimizableSetType v;
//  delta=0.0001;

//  v["jee_b"]=p0+delta;
//  Psi.resetParameters(v);
//  ValueType logpsi_p = Psi.evaluateLog(W);
//  RealType eloc_p=H.evaluate(W);
//  std::cout << "### eloc_p " << eloc_p << " " << H[0] << std::endl;
//  eloc_p=H[0];

//  v["jee_b"]=p0-delta;
//  Psi.resetParameters(v);
//  ValueType logpsi_m = Psi.evaluateLog(W);
//  RealType eloc_m=H.evaluate(W);
//  std::cout << "### eloc_m " << eloc_m << " " << H[0] << std::endl;
//  eloc_m=H[0];


//  std::cout << "### Direct evaluation " << std::endl;
//  std::cout << "### logpsi_p " << logpsi_p << std::endl;
//  std::cout << "### logpsi_m " << logpsi_m << std::endl;
//  std::cout << "### d logpsi " << (logpsi_p-logpsi_m)/delta/2.0 << std::endl;
//  std::cout << "### d h logpsi " << (eloc_p-eloc_m)/delta/2.0/psi << std::endl;
//  std::cout << "### d h logpsi " << (eloc_p-eloc_m)/delta/2.0 << std::endl;
//}
}


//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMCLinearOptOMP.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Optimize/VarList.h"
#include "Numerics/LinearFit.h"
//#define ENABLE_VMC_OMP_MASTER
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

/// Constructor.
VMCLinearOptOMP::VMCLinearOptOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                                 HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
  UseDrift("yes"), NumOptimizables(0), w_beta(0.0), GEVtype("mixed"),
  w_alpha(0.0),printderivs("no")
//     myRNWarmupSteps(0), logoffset(2.0), logepsilon(0), beta_errorbars(0), alpha_errorbars(0),
{
  RootName = "vmc";
  QMCType ="VMCLinearOptOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  DumpConfig=false;
  //default is 10
  nWarmupSteps=10;
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(UseDrift,"use_drift","string");
  m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(nTargetSamples,"targetwalkers","int");
  m_param.add(nTargetSamples,"target_walkers","int");
//     m_param.add(beta_errorbars,"beta_error","double");
//     m_param.add(alpha_errorbars,"alpha_error","double");
  m_param.add(w_beta,"beta","double");
  m_param.add(w_alpha,"alpha","double");
//     m_param.add(logepsilon,"logepsilon","double");
//     m_param.add(logoffset,"logoffset","double");
  m_param.add(printderivs,"printderivs","string");
  m_param.add(GEVtype,"GEVMethod","string");
//     m_param.add(myRNWarmupSteps,"rnwarmupsteps","int");
//     m_param.add(myRNWarmupSteps,"cswarmupsteps","int");
}

bool VMCLinearOptOMP::run()
{
  RngSaved.resize(NumThreads);
  resetRun();
  std::vector<opt_variables_type> dummyOptVars;
  for (int ip=0; ip<NumThreads; ++ip)
  {
    opt_variables_type dummy;
    psiClones[ip]->checkInVariables(dummy);
    dummy.resetIndex();
    psiClones[ip]->checkOutVariables(dummy);
    dummyOptVars.push_back(dummy);
  }
  NumOptimizables=dummyOptVars[0].size();
  resizeForOpt(NumOptimizables);
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip=0; ip<NumThreads; ++ip)
    Movers[ip]->startRun(nBlocks,false);
//     RealType target_errorbars;
//     target_errorbars = beta_errorbars;
//     RealType errorbars = target_errorbars+1;
  CurrentStep=0;
  int CurrentBlock=0;
//     int minBlocks=4;
  while (CurrentBlock<nBlocks)
  {
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
    {
      Movers[ip]->startBlock(nSteps);
      int now_loc=CurrentStep;
      //rest the collectables and keep adding
      wClones[ip]->resetCollectables();
      //rest the collectables and keep adding
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      for (int step=0; step<nSteps; ++step)
      {
        Movers[ip]->advanceWalkers(wit,wit_end,false);
        Movers[ip]->accumulate(wit,wit_end);
        ++now_loc;
        if (Period4WalkerDump&& now_loc%myPeriod4WalkerDump==0)
          wClones[ip]->saveEnsemble(wit,wit_end);
      }
      Movers[ip]->stopBlock(false);
    }//end-of-parallel for
    CurrentStep+=nSteps;
//       Estimators->accumulateCollectables(wClones,nSteps);
    Estimators->stopBlock(estimatorClones);
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
    {
      std::vector<RealType> Dsaved(NumOptimizables);
      std::vector<RealType> HDsaved(NumOptimizables);
      psiClones[ip]->evaluateDerivatives(*wClones[ip],dummyOptVars[ip],Dsaved,HDsaved);
      #pragma omp critical
      {
        copy(Dsaved.begin(),Dsaved.end(),&DerivRecords(ip,0));
        copy(HDsaved.begin(),HDsaved.end(),&HDerivRecords(ip,0));
      }
    }
    fillComponentMatrices();
    CurrentBlock++;
  }//block
//     app_log()<<" Blocks used   : "<<CurrentBlock<< std::endl;
//     app_log()<<" Errorbars are : "<<errorbars<< std::endl;
  Estimators->stop(estimatorClones);
  //copy back the random states
  for (int ip=0; ip<NumThreads; ++ip)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  //finalize a qmc section
  return finalize(nBlocks);
}



//   void VMCLinearOptOMP::initCS()
//   {
//     firstWalker=(*W[0]);
// #pragma omp parallel
//     {
//       int ip=omp_get_thread_num();
//       if (QMCDriverMode[QMC_UPDATE_MODE])
//         CSMovers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       else
//         CSMovers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//     }
//     //resetting containers
//     clearCSEstimators();
//     w_i.resize(NumThreads);
//     for (int ip=0; ip<NumThreads; ++ip) w_i[ip]=0;
//
//     //  set all walker positions to the same place
//     setWalkersEqual(firstWalker);
//     for(int ip=0; ip<NumThreads; ++ip)
//     {
//       //    synchronize the random number generator with the node
//       *Rng[ip]=*Rng[0];
//       hClones[ip]->setRandomGenerator(Rng[ip]);
//     }
//     if(myRNWarmupSteps>0)
//     {
//       for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//       {
//         CSMovers[0]->estimateNormWalkers(psiClones, wClones, hClones, Rng, w_i);
//         //         for (int ip=0; ip<NumThreads; ip++) app_log()<<"  w_i:"<<w_i[ip]<< std::endl;
//       }
//       myComm->allreduce(w_i);
//       RealType w_0=w_i[0];
//       for (int ip=0; ip<NumThreads; ++ip) w_i[ip] = -std::log(w_i[ip]/w_0);
//       for(int ip=0; ip<NumThreads; ++ip)
//       {
//         //    synchronize the random number generator with the node
//         *Rng[ip]=*Rng[0];
//         hClones[ip]->setRandomGenerator(Rng[ip]);
//       }
//     }
//
//     RealType overNT= 1.0/NumThreads;
//     for (int step=0; step<nWarmupSteps; ++step)
//     {
//       CSMovers[0]->advanceCSWalkers(psiClones, wClones, hClones, Rng, w_i);
//       estimateCS();
//       int max_i(0);
//       int min_i(0);
//       for (int ip=1; ip<NumThreads; ip++) if(Norms[ip]>Norms[max_i]) max_i=ip;
//       for (int ip=1; ip<NumThreads; ip++) if(Norms[ip]<Norms[min_i]) min_i=ip;
//       if ((Norms[max_i]-Norms[min_i])< 0.1*overNT)
//       {
//         step=nWarmupSteps;
//         clearCSEstimators();
//         continue;
//       }
//       //   rebalance weights
//       for (int ip=0; ip<NumThreads; ip++)
//       {
//         w_i[ip] += overNT*std::log(gNorms[0]/gNorms[ip]);
// //         app_log()<<"Norm["<<ip<<"]: "<<Norms[ip]<<"  w_i:"<<w_i[ip]<< std::endl;
//       }
//       clearCSEstimators();
//     }
//     setWalkersEqual(firstWalker);
//     clearCSEstimators();
//   }

//   int VMCLinearOptOMP::runCS(std::vector<std::vector<RealType> >& bestParams, RealType& errorbars)
//   {
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       opt_variables_type dummy;
//       psiClones[ip]->checkInVariables(dummy);
//       dummy.resetIndex();
//       psiClones[ip]->checkOutVariables(dummy);
//       for (int i=0;i<bestParams[0].size();i++)  dummy[i] = bestParams[ip][i];
//       psiClones[ip]->resetParameters(dummy);
//       //       app_log()<<ip<< std::endl;
//       //       psiClones[ip]->reportStatus(app_log());
//     }
//
//
//     // save the state of current generators
//
//     for(int ip=0; ip<NumThreads; ++ip) RngSaved[ip]=*Rng[ip];
//     initCS();
//
//     errorbars=alpha_errorbars+1;
//     CurrentStep=0;
//     CSBlock=0;
//     int minCSBlocks(4);
//     //     run long enough to get accurate errorbars ~4 blocks.
//     //     run until errorbars are small enough or when the energy difference is greater than 3 errorbars.
//     //     max run is defined by nBlocks
//     while ((CSBlock<minCSBlocks)||((errorbars>alpha_errorbars)&&(CSBlock<nBlocks)))
//     {
//       for (int step=0; step<nSteps; ++step)
//         CSMovers[0]->advanceCSWalkers(psiClones, wClones, hClones, Rng, w_i);
//       //         app_log()<<CSBlock<< std::endl;
//       CurrentStep+=nSteps;
//       errorbars = estimateCS();
//       CSBlock++;
//     }//block
//     app_log()<<" Blocks used   : "<<CSBlock<< std::endl;
//     app_log()<<" Errorbars are : "<<errorbars<< std::endl;
// //         app_log()<<" Min E["<<minE<<"] estimate: "<<NE_i[minE]<< std::endl;
//
//     ///restore the state
//     for(int ip=1; ip<NumThreads; ++ip)
//     {
//       *Rng[ip]=RngSaved[ip];
//       hClones[ip]->setRandomGenerator(Rng[ip]);
//     }
//     //copy back the random states
//     for (int ip=0; ip<NumThreads; ++ip) *(RandomNumberControl::Children[ip])=*(Rng[ip]);
//     if (std::abs(NE_i[minE])<1e6)  return minE;
//     else return -1;
//   }

//   bool VMCLinearOptOMP::bracketing(std::vector<RealType>& lambdas, RealType errorbars)
//   {
//     //Do some bracketing and line searching if we need to
//     RealType dl = std::abs(lambdas[1]-lambdas[0]);
//     RealType mL= lambdas[minE];
//     RealType DE = NE_i[nE] - NE_i[minE];
//
//     if (minE==(NumThreads-1))
//     {
//       if(moved_left)
//       {
//         app_log()<<" Bracketed minimum between CS runs"<< std::endl;
//         if (lambdas[minE]==0)
//         {
//           moved_right=false;
//           moved_left=false;
//           dl = std::abs(lambdas[1]-lambdas[0]);
//           mL = -0.5*(NumThreads-1)*dl;
//           for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
//         }
//         else return false;
//       }
//       else
//       {
//         app_log()<<" Move Right"<< std::endl;
//         moved_right=true;
//         //minE is an extreme value on the line search, move over and repeat.
//         for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = ip*dl + mL;
//       }
//     }
//     else if(minE==0)
//     {
//       if (moved_right)
//       {
//         app_log()<<" Bracketed minimum between CS runs"<< std::endl;
//         if (lambdas[minE]==0)
//         {
//           moved_right=false;
//           moved_left=false;
//           dl = std::abs(lambdas[1]-lambdas[0]);
//           mL = -0.5*(NumThreads-1)*dl;
//           for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
//         }
//         else return false;
//       }
//       else
//       {
//         app_log()<<" Move Left"<< std::endl;
//         moved_left=true;
//         //minE is an extreme value on the line search, move over and repeat.
//         for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = (ip-NumThreads+1.0)*dl + mL;
//       }
//     }
//     else
//     {
//       //         minimum is bracketed
//       // if energy difference is smaller than the errorbars and deltaP<1e-4 we computed then we are done
//       //         if (DE<errorbars)
//       //         {
//       int nms=3;
//       std::vector<RealType>  Y(nms), Coefs(3);
//       Matrix<RealType> X(nms,3);
//       for (int i=0; i<nms; i++) X(i,0)=1.0;
//       for (int i=0; i<nms; i++) X(i,1)=lambdas[i+minE-1];
//       for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
//       for (int i=0; i<nms; i++) Y[i]=NE_i[i+minE-1];
//       LinearFit(Y,X,Coefs);
//
//       RealType quadraticMinimum(-0.5*Coefs[1]/Coefs[2]);
//       lambdas[minE]=quadraticMinimum;
//       app_log()<<"Min predicted at: "<<quadraticMinimum<< std::endl;
//       return false;
//       //         }
//       //         else
//       //         {
//       //           app_log()<<" Bracketed minimum, refine"<< std::endl;
//       //           moved_right=false;
//       //           moved_left=false;
//       // // energy difference between the points is still larger than the error bars we require
//       // // need to "zoom" into find minimum more precisely
//       //             dl = 2.0*std::abs(lambdas[1]-lambdas[0])/(NumThreads-1.0);
//       //             mL = std::min(lambdas[minE],lambdas[nE])-0.5*dl;
//       //             for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + dl*ip;
//       //         }
//     }
//     return true;
//   }

//   VMCLinearOptOMP::RealType VMCLinearOptOMP::runCS(std::vector<RealType>& curParams, std::vector<RealType>& curDir, std::vector<RealType>& lambdas)
//   {
//     bool notConverged(true);
//
//     moved_right=false;
//     moved_left=false;
//     //   setWalkersEqual(firstWalker);
//     while (notConverged)
//     {
//       std::vector<std::vector<RealType> > dummy(NumThreads,std::vector<RealType>(curParams.size()));
//       for (int ip=0; ip<NumThreads; ++ip) for (int i=0;i<curParams.size();i++)  dummy[ip][i] = curParams[i] + lambdas[ip]*curDir[i+1];
//
//       RealType errorbars;
//       int bombed = runCS(dummy, errorbars);
//       for (int ip=0; ip<NumThreads; ++ip) app_log()<<"E["<<lambdas[ip]<<"] estimate: "<<NE_i[ip]<< std::endl;
//       if (bombed<0){
//         lambdas[0]=0;
//         return 0;
//       }
//
//       int maxI(0);
//       RealType maxV(0);
//       for (int i=0;i<curParams.size();i++) if (maxV<std::abs(curDir[i+1])){ maxI=i; maxV=std::abs(curDir[i+1]);};
//       RealType maxPChange(maxV*(lambdas[1]-lambdas[0]));
//       app_log()<<" Parameter diffs: "<<maxPChange<< std::endl;
//       if (maxPChange<1e-6)
//       {
//         notConverged = false;
//       }
//       else
//         notConverged = bracketing(lambdas, alpha_errorbars);
//
//     }
//     setWalkersEqual(firstWalker);
//     lambdas[0]=lambdas[minE];
//     return NE_i[minE];
//   }

//   VMCLinearOptOMP::RealType VMCLinearOptOMP::estimateCS()
//   {
//     std::vector<long double> e_i(NumThreads), psi2_i(NumThreads);
//     long double psi2(0);
//
//     for (int ip=0; ip<NumThreads; ip++)
//     {
//       e_i[ip]    = (W[ip])->getPropertyBase()[LOCALENERGY];
//       psi2_i[ip] = expl(2.0*(W[ip])->getPropertyBase()[LOGPSI] + w_i[ip] - logpsi2_0_0);
//       psi2       += psi2_i[ip];
//     }
//     for (int ip=0; ip<NumThreads; ip++) Norms[ip]  += psi2_i[ip]/psi2;
//     Norms[NumThreads] += psi2;
//     Norm2s(NumThreads,NumThreads) += psi2*psi2;
//
//
//     for (int ip=0; ip<NumThreads; ip++)
//       Energies[ip] += e_i[ip]*(psi2_i[ip]/psi2);
//
//     for (int ip=0; ip<NumThreads; ip++)
//       for (int ip2=0; ip2<NumThreads; ip2++)
//       {
//         Norm2s(ip,ip2)      += psi2_i[ip]*(psi2_i[ip2]/psi2);
//         CorrelatedH(ip,ip2) += psi2_i[ip]*(psi2_i[ip2]/psi2)*e_i[ip]*e_i[ip2];
//       }
//
//     //   // global quantities for mpi collection
//     //   std::vector<RealType> gEnergies(Energies), gNorms(Norms);
//     //   Matrix<RealType> gNorm2s(Norm2s), gCorrelatedH(CorrelatedH);
//     gNorm2s=Norm2s;
//     gCorrelatedH=CorrelatedH;
//     for (int ip=0; ip<NumThreads; ip++) gEnergies[ip]=Energies[ip];
//     for (int ip=0; ip<NumThreads+1; ip++) gNorms[ip]=Norms[ip];
//
//     myComm->allreduce(gEnergies);
//     myComm->allreduce(gNorms);
//     myComm->allreduce(gNorm2s);
//     myComm->allreduce(gCorrelatedH);
//
//     //   Here are the global energy estimates
//     for (int ip=0; ip<NumThreads; ip++) NE_i[ip] = gEnergies[ip]/gNorms[ip];
//     //   for (int ip=0; ip<NumThreads; ip++) app_log()<<ip<<": "<<gEnergies[ip]<<"  "<<gNorms[ip]<<"  "<<gNorm2s(ip,ip)<< std::endl;
//     //   app_log()<<NumThreads<<": "<<gNorms[NumThreads]<<"  "<<gNorm2s(NumThreads,NumThreads)<< std::endl;
//     //   app_log()<< std::endl;
//     //   find lowest energy
//     minE=0;
//     for (int ip=1; ip<NumThreads; ip++) if (NE_i[ip]<NE_i[minE]) minE=ip;
//
//     //   nE is the next lowest energy
//     nE=minE;
//     if (minE==0) nE=1;
//     else if (minE==NumThreads-1) nE=NumThreads-2;
//     else nE=(NE_i[minE+1]>NE_i[minE-1])?minE-1:minE+1;
//
// //     return the error in the energy differences between lowest two. not quite right.
//     long double rval = (gCorrelatedH(minE,minE)+gCorrelatedH(nE,nE)-2.0*gCorrelatedH(minE,nE))/Norm2s(minE,nE);
//
//     //rval = ((rval<0)?-1.0:(std::sqrt(rval/(CSBlock+1))));
//     rval = ((rval<0)?1.0:(std::sqrt(rval/(CSBlock+1.0))));
//     return rval;
//   }

void VMCLinearOptOMP::resetRun()
{
  //only VMC can overwrite this
  if(nTargetPopulation>0)
    branchEngine->iParam[SimpleFixedNodeBranch::B_TARGETWALKERS]=static_cast<int>(std::ceil(nTargetPopulation));
  //     firstWalker=(*W[0]);
  makeClones(W,Psi,H);
//     clearCSEstimators();
  std::vector<IndexType> samples_th(omp_get_max_threads(),0);
  myPeriod4WalkerDump=(Period4WalkerDump>0)?Period4WalkerDump:(nBlocks+1)*nSteps;
  int samples_this_node = nTargetSamples/myComm->size();
  if (nTargetSamples%myComm->size() > myComm->rank())
    samples_this_node+=1;
  int samples_each_thread = samples_this_node/omp_get_max_threads();
  for (int ip=0; ip<omp_get_max_threads(); ++ip)
    samples_th[ip]=samples_each_thread;
  if(samples_this_node%omp_get_max_threads())
    for (int ip=0; ip < samples_this_node%omp_get_max_threads(); ++ip)
      samples_th[ip] +=1;
  app_log() << "  Samples are dumped every " << myPeriod4WalkerDump << " steps " << std::endl;
  app_log() << "  Total Sample Size =" << nTargetSamples << std::endl;
  app_log() << "  Nodes Sample Size =" << samples_this_node << std::endl;
  for (int ip=0; ip<NumThreads; ++ip)
    app_log()  << "    Sample size for thread " <<ip<<" = " << samples_th[ip] << std::endl;
  app_log() << "  Warmup Steps " << nWarmupSteps << std::endl;
//     if (UseDrift == "rn") makeClones( *(psiPool.getWaveFunction("guide")) );
//    app_log() << "  Warmup Steps " << nWarmupSteps << std::endl;
  if (Movers.empty())
  {
    Movers.resize(NumThreads,0);
//       CSMovers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    traceClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
    FairDivideLow(nwtot,NumThreads,wPerNode);
    app_log() << "  Initial partition of walkers ";
    copy(wPerNode.begin(),wPerNode.end(),std::ostream_iterator<int>(app_log()," "));
    app_log() << std::endl;
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
    {
      std::ostringstream os;
      estimatorClones[ip]= new EstimatorManagerBase(*Estimators);//,*hClones[ip]);
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
#if !defined(REMOVE_TRACEMANAGER)
      traceClones[ip] = Traces->makeClone();
#endif
      Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if (QMCDriverMode[QMC_UPDATE_MODE])
      {
//           if (UseDrift == "rn")
//           {
//             os <<"  PbyP moves with RN, using VMCUpdatePbyPSampleRN"<< std::endl;
//             Movers[ip]=new VMCUpdatePbyPSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
//             Movers[ip]->setLogEpsilon(logepsilon);
//
//             CSMovers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//             //               Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
//           else if (UseDrift == "yes")
//           {
//             os <<"  PbyP moves with drift, using VMCUpdatePbyPWithDriftFast"<< std::endl;
//             CSMovers[ip]=Movers[ip]=new VMCUpdatePbyPWithDriftFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
// //             CSMovers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
//           else
//           {
        os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<< std::endl;
//             CSMovers[ip]=
        Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      else
      {
//           if (UseDrift == "rn")
//           {
//             os <<"  walker moves with RN, using VMCUpdateAllSampleRN"<< std::endl;
//             Movers[ip] =new VMCUpdateAllSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
//             Movers[ip]->setLogEpsilon(logepsilon);
//
//             CSMovers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
        //             else if (UseDrift == "yes")
        //             {
        //               os <<"  walker moves with drift, using VMCUpdateAllWithDriftFast"<< std::endl;
        //               Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        //             }
//           else
//           {
        os <<"  walker moves with |psi|^2, using VMCUpdateAll"<< std::endl;
//             CSMovers[ip]=
        Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      if (ip==0)
        app_log() << os.str() << std::endl;
    }
  }
#if !defined(REMOVE_TRACEMANAGER)
  else
  {
    #pragma omp parallel for
    for(int ip=0; ip<NumThreads; ++ip)
    {
      traceClones[ip]->transfer_state_from(*Traces);
    }
  }
#endif
  #pragma omp parallel
  {
    int ip=omp_get_thread_num();
    Movers[ip]->put(qmcNode);
//       CSMovers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip],traceClones[ip]);
//       CSMovers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
    if (QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    else
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       if (UseDrift != "rn")
//       {
    for (int prestep=0; prestep<nWarmupSteps; ++prestep)
      Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
    if (nWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    #pragma omp critical
    {
      wClones[ip]->clearEnsemble();
      wClones[ip]->setNumSamples(samples_th[ip]);
    }
//       }
  }
//     if (UseDrift == "rn")
//     {
//       RealType avg_w(0);
//       RealType n_w(0);
// #pragma omp parallel
//       {
//         int ip=omp_get_thread_num();
//
//         for (int step=0; step<nWarmupSteps; ++step)
//         {
//           avg_w=0;
//           n_w=0;
//           for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//           {
//             Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
// #pragma omp single
//             {
//               MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//               while (wit!=wit_end)
//               {
//                 avg_w += (*wit)->Weight;
//                 n_w +=1;
//                 wit++;
//               }
//             }
// #pragma omp barrier
//           }
// #pragma omp single
//           {
//             avg_w *= 1.0/n_w;
//             RealType w_m = avg_w/(1.0-avg_w);
//             w_m = std::log(0.5+0.5*w_m);
//             if (std::abs(w_m)>0.01)
//               logepsilon += w_m;
//           }
// #pragma omp barrier
//           Movers[ip]->setLogEpsilon(logepsilon);
//         }
//
//         for (int prestep=0; prestep<nWarmupSteps; ++prestep)
//           Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//
//         if (nWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
//           Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//
// #pragma omp critical
//         {
//             wClones[ip]->clearEnsemble();
//             wClones[ip]->setNumSamples(samples_th[ip]);
//         }
//       }
//     }
}

//   void VMCLinearOptOMP::fillMatrices(Matrix<RealType>& H2, Matrix<RealType>& Hamiltonian, Matrix<RealType>& Variance, Matrix<RealType>& Overlap)
//   {
//     RealType nrm = 1.0/sW;
//     //     RealType nrm2 = nrm*nrm;
//     for (int i=0; i<NumOptimizables; i++)
//     {
//       HDiE[i]*= nrm;
//       HDi[i] *= nrm;
//       DiE2[i]*= nrm;
//       DiE[i] *= nrm;
//       Di[i]  *= nrm;
//     }
//     HDiHDj*= nrm;
//     DiHDjE*= nrm;
//     DiHDj *= nrm;
//     DiDjE2*= nrm;
//     DiDjE *= nrm;
//     DiDj  *= nrm;
//
//     RealType H2_avg = sE2*nrm;
//     E_avg = sE*nrm;
//     V_avg = H2_avg - E_avg*E_avg;
//
//
//     for (int pm=0; pm<NumOptimizables; pm++)
//     {
//       RealType wfe = HDi[pm] + DiE[pm]-Di[pm]*E_avg;
//       RealType wfm = HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg;
//       //         Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
//
//       H2(0,pm+1) = HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg;
//       H2(pm+1,0) = HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg;
//
//       //         HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
//       RealType vterm = HDiE[pm]-HDi[pm]*E_avg + DiE2[pm]-Di[pm]*H2_avg -2.0*E_avg*(DiE[pm]-Di[pm]*E_avg);
//       Variance(0,pm+1) = vterm;
//       Variance(pm+1,0) = vterm;
//
//       Hamiltonian(0,pm+1) = wfe;
//       Hamiltonian(pm+1,0) = DiE[pm]-Di[pm]*E_avg;
//
//       for (int pm2=0; pm2<NumOptimizables; pm2++)
//       {
//         //           H2(pm+1,pm2+1) += wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new - curAvg_w));
//         //           Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
//         //           Variance(pm+1,pm2+1) += wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
//         //           Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);
//
//         //        Symmetric  (HDi[pm] + DiE[pm]-Di[pm]*E_avg)(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg)
//         H2(pm+1,pm2+1) = HDiHDj(pm,pm2) + DiHDjE(pm2,pm) - DiHDj(pm2,pm)*E_avg
//           + DiHDjE(pm,pm2) + DiDjE2(pm,pm2) - DiDjE(pm,pm2)*E_avg
//           + E_avg*(DiHDj(pm,pm2) + DiDjE(pm,pm2) - DiDj(pm,pm2)*E_avg);
//         //        Non-symmetric    (Dsaved[pm]-D_avg[pm])*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w))
//         Hamiltonian(pm+1,pm2+1) = DiHDj(pm,pm2) + DiDjE(pm,pm2) - Di[pm2]*DiE[pm] - Di[pm]*(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg);
//         //        Symmetric  (HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg)*( HDi[pm2] - 2.0* DiE[pm2]+2.0*Di[pm2]*E_avg)
//         Variance(pm+1,pm2+1) = HDiHDj(pm,pm2) -2.0*DiHDjE(pm2,pm) +2.0*DiHDj(pm,pm2)*E_avg
//           -2.0*( DiHDjE(pm,pm2) - 2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2))
//           +2.0*E_avg*(DiHDj(pm,pm2) -2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2));
//         //        Symmetric
//         Overlap(pm+1,pm2+1) = DiDj(pm,pm2)-Di[pm]*Di[pm2];
//
//       }
//     }
//
//     Hamiltonian(0,0) = E_avg;
//     Overlap(0,0) = 1.0;
//     H2(0,0) = H2_avg;
//     Variance(0,0) = V_avg;
//
//     for (int pm=1; pm<NumOptimizables+1; pm++)
//       for (int pm2=1; pm2<NumOptimizables+1; pm2++)
//         Variance(pm,pm2) += V_avg*Overlap(pm,pm2);
//
//     //     app_log()<<V_avg<<"  "<<E_avg<<"  "<<sW<< std::endl;
//   }

//   VMCLinearOptOMP::RealType VMCLinearOptOMP::fillOverlapHamiltonianMatrices(Matrix<RealType>& LeftM, Matrix<RealType>& RightM)
//   {
//     RealType b1,b2;
//     if (GEVtype=="H2")
//     {
//       b1=w_beta; b2=0;
//     }
//     else
//     {
//       b2=w_beta; b1=0;
//     }
//
//     RealType nrm = 1.0/sW;
//     //     RealType nrm2 = nrm*nrm;
//     for (int i=0; i<NumOptimizables; i++)
//     {
//       HDiE[i]*= nrm;
//       HDi[i] *= nrm;
//       DiE2[i]*= nrm;
//       DiE[i] *= nrm;
//       Di[i]  *= nrm;
//     }
//     HDiHDj*= nrm;
//     DiHDjE*= nrm;
//     DiHDj *= nrm;
//     DiDjE2*= nrm;
//     DiDjE *= nrm;
//     DiDj  *= nrm;
//
//     RealType H2_avg = 1.0/(sE2*nrm);
//     E_avg = sE*nrm;
//     V_avg = H2_avg - E_avg*E_avg;
//
//
//     for (int pm=0; pm<NumOptimizables; pm++)
//     {
//       RealType wfe = HDi[pm] + DiE[pm]-Di[pm]*E_avg;
//       RealType wfm = HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg;
//       //         Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
//
// //       H2
//       RightM(0,pm+1) += b1*H2_avg*(HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg);
//       RightM(pm+1,0) += b1*H2_avg*( HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg);
//
//       //         HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
//       RealType vterm = HDiE[pm]-HDi[pm]*E_avg + DiE2[pm]-Di[pm]*H2_avg -2.0*E_avg*(DiE[pm]-Di[pm]*E_avg);
// //       variance
//       LeftM(0,pm+1) += b2*vterm;
//       LeftM(pm+1,0) += b2*vterm;
//
// //       hamiltonian
//       LeftM(0,pm+1) += (1-b2)*wfe;
//       LeftM(pm+1,0) += (1-b2)*(DiE[pm]-Di[pm]*E_avg);
//
//       for (int pm2=0; pm2<NumOptimizables; pm2++)
//       {
//         //           H2(pm+1,pm2+1) += wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new - curAvg_w));
//         //           Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
//         //           Variance(pm+1,pm2+1) += wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
//         //           Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);
//
//         //        Symmetric  (HDi[pm] + DiE[pm]-Di[pm]*E_avg)(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg)
// //         H2
//         RightM(pm+1,pm2+1) += (b1*H2_avg)*(HDiHDj(pm,pm2) + DiHDjE(pm2,pm) - DiHDj(pm2,pm)*E_avg
//           + DiHDjE(pm,pm2) + DiDjE2(pm,pm2) - DiDjE(pm,pm2)*E_avg
//           + E_avg*(DiHDj(pm,pm2) + DiDjE(pm,pm2) - DiDj(pm,pm2)*E_avg));
//         //        Non-symmetric    (Dsaved[pm]-D_avg[pm])*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w))
// //         Hamiltonian
//         LeftM(pm+1,pm2+1) += (1-b2)*(DiHDj(pm,pm2) + DiDjE(pm,pm2) - Di[pm2]*DiE[pm] - Di[pm]*(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg));
//         //        Symmetric  (HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg)*( HDi[pm2] - 2.0* DiE[pm2]+2.0*Di[pm2]*E_avg)
// //         Variance
//         LeftM(pm+1,pm2+1) += b2*(HDiHDj(pm,pm2) -2.0*DiHDjE(pm2,pm) +2.0*DiHDj(pm,pm2)*E_avg
//           -2.0*( DiHDjE(pm,pm2) - 2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2))
//           +2.0*E_avg*(DiHDj(pm,pm2) -2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2)));
//         //        Symmetric
// //         Overlap
//         RightM(pm+1,pm2+1) += (1-b1)*(DiDj(pm,pm2)-Di[pm]*Di[pm2]);
// //         LeftM(pm+1,pm2+1) += b2*V_avg*(DiDj(pm,pm2)-Di[pm]*Di[pm2]);
//       }
//     }
//
//     LeftM(0,0) += (1-b2)*E_avg;
//     RightM(0,0) = 1.0;
//     LeftM(0,0) += b2*V_avg;
//
// //     for (int pm=0; pm<NumOptimizables; pm++)
// //       for (int pm2=0; pm2<NumOptimizables; pm2++)
// //         LeftM(pm+1,pm2+1) += b2*V_avg*RightM(pm+1,pm2+1);
//
//     //     app_log()<<V_avg<<"  "<<E_avg<<"  "<<sW<< std::endl;
//     if (GEVtype=="H2")
//       return 1.0/H2_avg;
//     return 1.0;
//
//   }

VMCLinearOptOMP::RealType VMCLinearOptOMP::fillOverlapHamiltonianMatrices(Matrix<RealType>& LeftM, Matrix<RealType>& RightM)
{
  RealType b1,b2;
  if (GEVtype=="H2")
  {
    b1=w_beta;
    b2=0;
  }
  else
  {
    b2=w_beta;
    b1=0;
  }
  std::vector<RealType> g_stats(5,0);
  g_stats[0]=s_vec[0];
  g_stats[1]=s_vec[1];
  g_stats[2]=s_vec[2];
  g_stats[3]=s_vec[3];
  g_stats[4]=s_vec[4];
  myComm->allreduce(g_stats);
  RealType g_nrm = 1.0/g_stats[3];
  E_avg = g_nrm*g_stats[0];
  RealType E_avg2=E_avg*E_avg;
  RealType E2_avg = g_nrm*g_stats[1];
  V_avg = E2_avg-E_avg2;
//     app_log()<<E_avg<<" "<<V_avg<<" "<<E2_avg<< std::endl;
  myComm->allreduce(Ham2);
  Ham2*=g_nrm;
  myComm->allreduce(Ham);
  Ham *=g_nrm;
  myComm->allreduce(Olp);
  Olp *=g_nrm;
  myComm->allreduce(m_vec);
  m_vec*=g_nrm;
  if ((printderivs=="yes")&&(myComm->rank()==0))
  {
    std::stringstream fn;
    fn<<RootName.c_str()<<".derivs";
    std::ofstream d_out(fn.str().c_str());
    d_out.precision(6);
    d_out<<"#csf    D        HD"<< std::endl;
    for (int i=0; i<NumOptimizables; i++)
      d_out<<i+1<<" "<<m_vec(0,i)<<"  "<<m_vec(1,i)<< std::endl;
  }
  for (int i=0; i<NumOptimizables; i++)
    for (int j=0; j<NumOptimizables; j++)
      Ham(i,j) += -m_vec(0,i)*(m_vec(1,j)+ m_vec(2,j) - m_vec(0,j)*E_avg)  -m_vec(0,j)*m_vec(2,i);
  for (int i=0; i<NumOptimizables; i++)
    for (int j=0; j<NumOptimizables; j++)
      Olp(i,j) -= m_vec(0,i)*m_vec(0,j);
  for (int i=0; i<NumOptimizables; i++)
    for (int j=0; j<NumOptimizables; j++)
      Ham2(i,j) += 2*m_vec(0,j)*(m_vec(3,i)-2.0*m_vec(4,i)) + 2*m_vec(0,i)*(m_vec(3,j)-2.0*m_vec(4,j)) +4*m_vec(0,i)*m_vec(0,j)*E2_avg;
  RealType b1_rat = b1/E_avg2;
  for (int i=1; i<NumOptimizables+1; i++)
    for (int j=1; j<NumOptimizables+1; j++)
    {
      LeftM(i,j) = (1-b2)*Ham(i-1,j-1) + b2*(Ham2(i-1,j-1) + V_avg*Olp(i-1,j-1)) ;
      RightM(i,j) = Olp(i-1,j-1) + b1_rat*Ham2(i-1,j-1);
    }
  RightM(0,0)= 1.0;
  LeftM(0,0)=(1-b2)*E_avg+b2*V_avg;
  for (int i=1; i<NumOptimizables+1; i++)
  {
    RealType vterm=m_vec(3,i-1)-m_vec(1,i-1)*E_avg + m_vec(4,i-1)-m_vec(0,i-1)*E2_avg -2.0*(m_vec(2,i-1)*E_avg-m_vec(0,i-1)*E_avg2);
    RightM(0,i)= RightM(i,0) = b1_rat*vterm;
    LeftM(i,0) = (1-b2)*(m_vec(2,i-1)-E_avg*m_vec(0,i-1))
                 +b2*vterm;
    LeftM(0,i) = (1-b2)*(m_vec(1,i-1)+m_vec(2,i-1)-E_avg*m_vec(0,i-1))
                 +b2*vterm;
  }
  return 1.0;
}

VMCLinearOptOMP::RealType VMCLinearOptOMP::fillComponentMatrices()
{
  std::vector<RealType> g_stats(5,0);
  for (int ip=0; ip<NumThreads; ip++)
  {
    RealType E_L = W[ip]->getPropertyBase()[LOCALENERGY];
    RealType E_L2= E_L*E_L;
    RealType wW  = W[ip]->Weight;
    if(std::isnan(wW)||std::isinf(wW))
      wW=0;
    s_vec[0]+=E_L*wW;
    s_vec[1]+=E_L2*wW;
    s_vec[2]+=E_L2*E_L2*wW;
    s_vec[3]+=wW;
    s_vec[4]+=1;
    for (int i=0; i<NumOptimizables; i++)
    {
      RealType di  = DerivRecords(ip,i);
      RealType hdi = HDerivRecords(ip,i);
      //             vectors
      m_vec(0,i)+= wW*di;
      m_vec(1,i)+= wW*hdi;
      m_vec(2,i)+= wW*di*E_L;
      m_vec(3,i)+= wW*hdi*E_L;
      m_vec(4,i)+= wW*di*E_L2;
      m_vec(5,i)+= wW*E_L*(hdi+di*E_L);
      for (int j=0; j<NumOptimizables; j++)
      {
        RealType dj  = DerivRecords(ip,j);
        RealType hdj = HDerivRecords(ip,j);
        Ham(i,j) += wW*di*(hdj+dj*E_L);
        Olp(i,j) += wW*di*dj;
        Ham2(i,j)+= wW*(hdj-2*dj*E_L)*(hdi-2*di*E_L);
      }
    }
  }
  g_stats[0]=s_vec[0];//     sE
  g_stats[1]=s_vec[1];//     sE2
  g_stats[2]=s_vec[2];//     sE4
  g_stats[3]=s_vec[3];//     sW
  g_stats[4]=s_vec[4];//     sN
  myComm->allreduce(g_stats);
  RealType nrm = 1.0/g_stats[3];
  E_avg = nrm*g_stats[0];
  V_avg = nrm*g_stats[1]-E_avg*E_avg;
  RealType g_nrm = 1.0/g_stats[4];
  RealType err_E(std::sqrt( ((V_avg<0.0)?(1.0):(V_avg*g_nrm)) ));
  RealType err_E2(nrm*g_stats[2]-nrm*nrm*s_vec[1]*s_vec[1]);
  err_E2 = std::sqrt( ((err_E2<0.0)?(1.0):(err_E2*g_nrm)) );
  return w_beta*err_E2+(1.0-w_beta)*err_E;
}


bool
VMCLinearOptOMP::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}


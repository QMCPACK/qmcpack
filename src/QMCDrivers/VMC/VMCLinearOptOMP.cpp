//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/VMC/VMCLinearOptOMP.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Optimize/VarList.h"
#include "Numerics/LinearFit.h"
//#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus
{

  /// Constructor.
  VMCLinearOptOMP::VMCLinearOptOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
      HamiltonianPool& hpool, WaveFunctionPool& ppool):
    QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
    myRNWarmupSteps(100), myWarmupSteps(10),UseDrift("yes"), NumOptimizables(0), w_beta(0.0), GEVtype("mixed"), logoffset(2.0), logepsilon(0)
  {
    RootName = "vmc";
    QMCType ="VMCLinearOptOMP";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_WARMUP,0);
    DumpConfig=false;
    m_param.add(UseDrift,"useDrift","string");
    m_param.add(UseDrift,"usedrift","string");
    m_param.add(UseDrift,"use_drift","string");
    m_param.add(myWarmupSteps,"warmupSteps","int");
    m_param.add(myWarmupSteps,"warmupsteps","int");
    m_param.add(myWarmupSteps,"warmup_steps","int");
    m_param.add(nTargetSamples,"targetWalkers","int");
    m_param.add(nTargetSamples,"targetwalkers","int");
    m_param.add(nTargetSamples,"target_walkers","int");
    m_param.add(beta_errorbars,"beta_error","double");
    m_param.add(alpha_errorbars,"alpha_error","double");
    m_param.add(w_beta,"beta","double");
    m_param.add(logepsilon,"logepsilon","double");
    m_param.add(logoffset,"logoffset","double");
    m_param.add(GEVtype,"GEVMethod","string");
    m_param.add(myRNWarmupSteps,"rnwarmupsteps","int");
    m_param.add(myRNWarmupSteps,"cswarmupsteps","int");
  }

  bool VMCLinearOptOMP::run()
  {
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
    for (int ip=0; ip<NumThreads; ++ip) Movers[ip]->startRun(nBlocks,false);


    RealType target_errorbars;
    target_errorbars = beta_errorbars;
    RealType errorbars = target_errorbars+1;
    CurrentStep=0;
    int CurrentBlock=0;
    int minBlocks=4;
    while (((errorbars>target_errorbars)&&(CurrentBlock<nBlocks))||(CurrentBlock<minBlocks))
    {
#pragma omp parallel for
      for (int ip=0; ip<NumThreads; ++ip)
      { 
        Movers[ip]->startBlock(nSteps);
        //rest the collectables and keep adding
        wClones[ip]->resetCollectables();
        //rest the collectables and keep adding

        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        for (int step=0; step<nSteps; ++step)
        {
          Movers[ip]->advanceWalkers(wit,wit_end,false);
          Movers[ip]->accumulate(wit,wit_end);
        }
        Movers[ip]->stopBlock(false);
      }//end-of-parallel for
      CurrentStep+=nSteps;

      Estimators->accumulateCollectables(wClones,nSteps);
      Estimators->stopBlock(estimatorClones);
#pragma omp parallel for
      for (int ip=0; ip<NumThreads; ++ip)
      {
        vector<RealType> Dsaved(NumOptimizables);
        vector<RealType> HDsaved(NumOptimizables);
        psiClones[ip]->evaluateDerivatives(*wClones[ip],dummyOptVars[ip],Dsaved,HDsaved);
#pragma omp critical
        {
          std::copy(Dsaved.begin(),Dsaved.end(),&DerivRecords(ip,0));
          std::copy(HDsaved.begin(),HDsaved.end(),&HDerivRecords(ip,0));
        }
      }
      errorbars = fillComponentMatrices();
      CurrentBlock++;

    }//block
    app_log()<<" Blocks used   : "<<CurrentBlock<<endl;
    app_log()<<" Errorbars are : "<<errorbars<<endl;

    Estimators->stop(estimatorClones);
    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip)
      *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    //finalize a qmc section
    return finalize(nBlocks);
  }



  void VMCLinearOptOMP::initCS()
  {
    firstWalker=(*W[0]);
#pragma omp parallel
    {
      int ip=omp_get_thread_num();
      if (QMCDriverMode[QMC_UPDATE_MODE])
        CSMovers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      else
        CSMovers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
    //resetting containers
    clearCSEstimators();
    w_i.resize(NumThreads);
    for (int ip=0; ip<NumThreads; ++ip) w_i[ip]=0;    

    //  set all walker positions to the same place
    setWalkersEqual(firstWalker);
    if(myRNWarmupSteps>0)
    {
      for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
      {
        CSMovers[0]->estimateNormWalkers(psiClones, wClones, hClones, Rng, w_i);
        //         for (int ip=0; ip<NumThreads; ip++) app_log()<<"  w_i:"<<w_i[ip]<<endl;
      }
      myComm->allreduce(w_i);
      RealType w_0=w_i[0];
      for (int ip=0; ip<NumThreads; ++ip) w_i[ip] = -std::log(w_i[ip]/w_0);
    }
    else
    {
      for (int ip=0; ip<NumThreads; ++ip) w_i[ip]=1.0;
    }
    RealType overNT= 1.0/NumThreads;
    for (int step=0; step<myWarmupSteps; ++step)
    {
      CSMovers[0]->advanceCSWalkers(psiClones, wClones, hClones, Rng, w_i);
      estimateCS();
      int max_i(0);
      int min_i(0);
      for (int ip=1; ip<NumThreads; ip++) if(Norms[ip]>Norms[max_i]) max_i=ip;
      for (int ip=1; ip<NumThreads; ip++) if(Norms[ip]<Norms[min_i]) min_i=ip;
      if ((Norms[max_i]-Norms[min_i])< 0.1*overNT)
      {
        step=myWarmupSteps;
        clearCSEstimators();
        continue;
      }
      //   rebalance weights
      for (int ip=0; ip<NumThreads; ip++)
      {
        w_i[ip] += overNT*std::log(gNorms[0]/gNorms[ip]);
//         app_log()<<"Norm["<<ip<<"]: "<<Norms[ip]<<"  w_i:"<<w_i[ip]<<endl;
      }
      clearCSEstimators();
    }
    setWalkersEqual(firstWalker);
    clearCSEstimators();
  }

  int VMCLinearOptOMP::runCS(vector<vector<RealType> >& bestParams, RealType& errorbars)
  {
    for (int ip=0; ip<NumThreads; ++ip)
    {
      opt_variables_type dummy;
      psiClones[ip]->checkInVariables(dummy);
      dummy.resetIndex();
      psiClones[ip]->checkOutVariables(dummy);
      for (int i=0;i<bestParams[0].size();i++)  dummy[i] = bestParams[ip][i];
      psiClones[ip]->resetParameters(dummy);
      //       app_log()<<ip<<endl;
      //       psiClones[ip]->reportStatus(app_log());
    }

    // save the state of current generators
    vector<RandomGenerator_t> RngSaved(NumThreads);
    for(int ip=1; ip<NumThreads; ++ip) RngSaved[ip]=*Rng[ip];
    for(int ip=0; ip<NumThreads; ++ip)
    {
      //    synchronize the random number generator with the node
      *Rng[ip]=*Rng[0];
      hClones[ip]->setRandomGenerator(Rng[ip]);
    }
    initCS();

    errorbars=alpha_errorbars+1;
    CurrentStep=0;
    CSBlock=0;
    int minCSBlocks(4);
    //     run long enough to get accurate errorbars ~4 blocks.
    //     run until errorbars are small enough or when the energy difference is greater than 3 errorbars.
    //     max run is defined by nBlocks
    while ((CSBlock<minCSBlocks)||((errorbars>alpha_errorbars)&&(CSBlock<nBlocks)))
    {
      for (int step=0; step<nSteps; ++step)
        CSMovers[0]->advanceCSWalkers(psiClones, wClones, hClones, Rng, w_i);
      //         app_log()<<CSBlock<<endl;
      CurrentStep+=nSteps;
      errorbars = estimateCS();
      CSBlock++;
    }//block
    app_log()<<" Blocks used   : "<<CSBlock<<endl;
    app_log()<<" Errorbars are : "<<errorbars<<endl;
//         app_log()<<" Min E["<<minE<<"] estimate: "<<NE_i[minE]<<endl;

    ///restore the state
    for(int ip=1; ip<NumThreads; ++ip)
    {
      *Rng[ip]=RngSaved[ip];
      hClones[ip]->setRandomGenerator(Rng[ip]);
    }
    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip) *(RandomNumberControl::Children[ip])=*(Rng[ip]);
    if (std::abs(NE_i[minE])<1e6)  return minE;
    else return -1;
  }

  bool VMCLinearOptOMP::bracketing(vector<RealType>& lambdas, RealType errorbars)
  { 
    //Do some bracketing and line searching if we need to
    RealType dl = std::abs(lambdas[1]-lambdas[0]);
    RealType mL= lambdas[minE];
    RealType DE = NE_i[nE] - NE_i[minE];

    if (minE==(NumThreads-1))
    {
      if(moved_left)
      {
        app_log()<<" Bracketed minimum between CS runs"<<endl;
        if (lambdas[minE]==0)
        {
          moved_right=false;
          moved_left=false;
          dl = std::abs(lambdas[1]-lambdas[0]);
          mL = -0.5*(NumThreads-1)*dl;
          for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
        }
        else return false;
      }
      else
      {
        app_log()<<" Move Right"<<endl;
        moved_right=true;
        //minE is an extreme value on the line search, move over and repeat.
        for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = ip*dl + mL;
      }
    }
    else if(minE==0)
    {
      if (moved_right)
      {
        app_log()<<" Bracketed minimum between CS runs"<<endl;
        if (lambdas[minE]==0)
        {
          moved_right=false;
          moved_left=false;
          dl = std::abs(lambdas[1]-lambdas[0]);
          mL = -0.5*(NumThreads-1)*dl;
          for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
        }
        else return false;
      }
      else
      {
        app_log()<<" Move Left"<<endl;
        moved_left=true;
        //minE is an extreme value on the line search, move over and repeat.
        for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = (ip-NumThreads+1.0)*dl + mL;
      }
    }
    else
    {
      //         minimum is bracketed
      // if energy difference is smaller than the errorbars and deltaP<1e-4 we computed then we are done
      //         if (DE<errorbars)
      //         {
      int nms=3;
      vector<RealType>  Y(nms), Coefs(3);
      Matrix<RealType> X(nms,3);
      for (int i=0; i<nms; i++) X(i,0)=1.0;
      for (int i=0; i<nms; i++) X(i,1)=lambdas[i+minE-1];
      for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
      for (int i=0; i<nms; i++) Y[i]=NE_i[i+minE-1];
      LinearFit(Y,X,Coefs);

      RealType quadraticMinimum(-0.5*Coefs[1]/Coefs[2]);
      lambdas[minE]=quadraticMinimum;
      app_log()<<"Min predicted at: "<<quadraticMinimum<<endl;
      return false;
      //         }
      //         else
      //         {
      //           app_log()<<" Bracketed minimum, refine"<<endl;
      //           moved_right=false;
      //           moved_left=false;
      // // energy difference between the points is still larger than the error bars we require
      // // need to "zoom" into find minimum more precisely
      //             dl = 2.0*std::abs(lambdas[1]-lambdas[0])/(NumThreads-1.0);
      //             mL = std::min(lambdas[minE],lambdas[nE])-0.5*dl;
      //             for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + dl*ip;
      //         }
    }
    return true;
  }

  VMCLinearOptOMP::RealType VMCLinearOptOMP::runCS(vector<RealType>& curParams, vector<RealType>& curDir, vector<RealType>& lambdas)
  {
    bool notConverged(true);

    moved_right=false;
    moved_left=false;
    //   setWalkersEqual(firstWalker);
    while (notConverged)
    {
      vector<vector<RealType> > dummy(NumThreads,std::vector<RealType>(curParams.size()));
      for (int ip=0; ip<NumThreads; ++ip) for (int i=0;i<curParams.size();i++)  dummy[ip][i] = curParams[i] + lambdas[ip]*curDir[i+1];

      RealType errorbars;
      int bombed = runCS(dummy, errorbars);
      for (int ip=0; ip<NumThreads; ++ip) app_log()<<"E["<<lambdas[ip]<<"] estimate: "<<NE_i[ip]<<endl;
      if (bombed<0){
        lambdas[0]=0;
        return 0;
      }    

      int maxI(0);
      RealType maxV(0);
      for (int i=0;i<curParams.size();i++) if (maxV<abs(curDir[i+1])){ maxI=i; maxV=abs(curDir[i+1]);};
      RealType maxPChange(maxV*(lambdas[1]-lambdas[0]));
      app_log()<<" Parameter diffs: "<<maxPChange<<endl;
      if (maxPChange<1e-6)
      {
        notConverged = false;
      }
      else
        notConverged = bracketing(lambdas, alpha_errorbars);

    }
    setWalkersEqual(firstWalker);
    lambdas[0]=lambdas[minE];
    return NE_i[minE];
  }

  VMCLinearOptOMP::RealType VMCLinearOptOMP::estimateCS()
  {
    vector<long double> e_i(NumThreads), psi2_i(NumThreads);
    long double psi2(0);

    for (int ip=0; ip<NumThreads; ip++)
    {
      e_i[ip]    = (W[ip])->getPropertyBase()[LOCALENERGY];
      psi2_i[ip] = expl(2.0*(W[ip])->getPropertyBase()[LOGPSI] + w_i[ip] - logpsi2_0_0);
      psi2       += psi2_i[ip];
    }
    for (int ip=0; ip<NumThreads; ip++) Norms[ip]  += psi2_i[ip]/psi2;
    Norms[NumThreads] += psi2;
    Norm2s(NumThreads,NumThreads) += psi2*psi2;


    for (int ip=0; ip<NumThreads; ip++)
      Energies[ip] += e_i[ip]*(psi2_i[ip]/psi2);

    for (int ip=0; ip<NumThreads; ip++)
      for (int ip2=0; ip2<NumThreads; ip2++)
      {
        Norm2s(ip,ip2)      += (psi2_i[ip]/psi2)*(psi2_i[ip2]/psi2) ;
        CorrelatedH(ip,ip2) += (psi2_i[ip]/psi2)*(psi2_i[ip2]/psi2)*e_i[ip]*e_i[ip2];
      }

    //   // global quantities for mpi collection
    //   std::vector<RealType> gEnergies(Energies), gNorms(Norms);
    //   Matrix<RealType> gNorm2s(Norm2s), gCorrelatedH(CorrelatedH);
    gNorm2s=Norm2s;
    gCorrelatedH=CorrelatedH;
    for (int ip=0; ip<NumThreads; ip++) gEnergies[ip]=Energies[ip];
    for (int ip=0; ip<NumThreads+1; ip++) gNorms[ip]=Norms[ip];

    myComm->allreduce(gEnergies);
    myComm->allreduce(gNorms);
    myComm->allreduce(gNorm2s);
    myComm->allreduce(gCorrelatedH);

    //   Here are the global energy estimates
    for (int ip=0; ip<NumThreads; ip++) NE_i[ip] = gEnergies[ip]/gNorms[ip];
    //   for (int ip=0; ip<NumThreads; ip++) app_log()<<ip<<": "<<gEnergies[ip]<<"  "<<gNorms[ip]<<"  "<<gNorm2s(ip,ip)<<endl;
    //   app_log()<<NumThreads<<": "<<gNorms[NumThreads]<<"  "<<gNorm2s(NumThreads,NumThreads)<<endl;
    //   app_log()<<endl;
    //   find lowest energy
    minE=0;
    for (int ip=1; ip<NumThreads; ip++) if (NE_i[ip]<NE_i[minE]) minE=ip;

    //   nE is the next lowest energy
    nE=minE;
    if (minE==0) nE=1;
    else if (minE==NumThreads-1) nE=NumThreads-2;
    else nE=(NE_i[minE+1]>NE_i[minE-1])?minE-1:minE+1;

    //return the error in the energy differences between lowest two
    long double rval = (gCorrelatedH(minE,minE)/(gNorms[minE]*gNorms[minE]) + gCorrelatedH(nE,nE)/(gNorms[nE]*gNorms[nE]) - 2.0*gCorrelatedH(minE,nE)/(gNorms[minE]*gNorms[nE]))-(NE_i[minE]-NE_i[nE])*(NE_i[minE]-NE_i[nE]);
//     long double rval = (gCorrelatedH(minE,minE)/(gNorm2s(minE,minE)) + gCorrelatedH(nE,nE)/(gNorm2s(nE,nE)) - 2.0*gCorrelatedH(minE,nE)/gNorm2s(minE,nE))-(NE_i[minE]-NE_i[nE])*(NE_i[minE]-NE_i[nE]);

    
// //     app_log()<<"CS_new_ED: "<<(gCorrelatedH(minE,minE)/(gNorms[minE]*gNorms[minE]) + gCorrelatedH(nE,nE)/(gNorms[nE]*gNorms[nE]) - 2.0*gCorrelatedH(minE,nE)/(gNorms[minE]*gNorms[nE]))<<endl;
//     app_log()<<"CS_old_ED: "<<(gCorrelatedH(minE,minE)/(gNorm2s(minE,minE)) + gCorrelatedH(nE,nE)/(gNorm2s(nE,nE)) - 2.0*gCorrelatedH(minE,nE)/gNorm2s(minE,nE))<<endl;
//     app_log()<<"AV_ED: "<<(NE_i[minE]-NE_i[nE])*(NE_i[minE]-NE_i[nE])<<endl;
    
    //rval = ((rval<0)?-1.0:(std::sqrt(rval/(CSBlock+1))));
    rval = ((rval<0)?1.0:(std::sqrt(rval/(CSBlock+1.0))));
    return rval;
  }

  void VMCLinearOptOMP::resetRun()
  {
    //     firstWalker=(*W[0]);
    makeClones(W,Psi,H);
    clearCSEstimators();
    if (UseDrift == "rn") makeClones( *(psiPool.getWaveFunction("guide")) );
    app_log() << "  Warmup Steps " << myWarmupSteps << endl;


    if (Movers.empty())
    {
      Movers.resize(NumThreads,0);
      CSMovers.resize(NumThreads,0);
      branchClones.resize(NumThreads,0);
      estimatorClones.resize(NumThreads,0);
      Rng.resize(NumThreads,0);
      int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
      FairDivideLow(nwtot,NumThreads,wPerNode);

      app_log() << "  Initial partition of walkers ";
      std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
      app_log() << endl;

#pragma omp parallel for
      for (int ip=0; ip<NumThreads; ++ip)
      {
        ostringstream os;
        estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
        estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
        estimatorClones[ip]->setCollectionMode(false);
        Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
        hClones[ip]->setRandomGenerator(Rng[ip]);

        branchClones[ip] = new BranchEngineType(*branchEngine);

        if (QMCDriverMode[QMC_UPDATE_MODE])
        {
          if (UseDrift == "rn")
          {
            os <<"  PbyP moves with RN, using VMCUpdatePbyPSampleRN"<<endl;
            Movers[ip]=new VMCUpdatePbyPSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
            Movers[ip]->setLogEpsilon(logepsilon);

            CSMovers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            //               Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          }
//           else if (UseDrift == "yes")
//           {
//             os <<"  PbyP moves with drift, using VMCUpdatePbyPWithDriftFast"<<endl;
//             CSMovers[ip]=Movers[ip]=new VMCUpdatePbyPWithDriftFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
// //             CSMovers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           }
          else
          {
            os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<<endl;
            CSMovers[ip]=Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          }
          //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }
        else
        {
          if (UseDrift == "rn")
          {
            os <<"  walker moves with RN, using VMCUpdateAllSampleRN"<<endl;
            Movers[ip] =new VMCUpdateAllSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
            Movers[ip]->setLogEpsilon(logepsilon);

            CSMovers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          }
          //             else if (UseDrift == "yes")
          //             {
          //               os <<"  walker moves with drift, using VMCUpdateAllWithDriftFast"<<endl;
          //               Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          //             }
          else
          {
            os <<"  walker moves with |psi|^2, using VMCUpdateAll"<<endl;
            CSMovers[ip]=Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          }
          //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
        }

        if (ip==0) app_log() << os.str() << endl;
      }
    }

#pragma omp parallel
    {
      int ip=omp_get_thread_num();
      Movers[ip]->put(qmcNode);
      CSMovers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      CSMovers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);

      if (QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      else
        Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);

      if (UseDrift != "rn")
      {
        for (int prestep=0; prestep<myWarmupSteps; ++prestep)
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);

        if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
          Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);



#pragma omp critical
        {
          wClones[ip]->clearEnsemble();
          wClones[ip]->setNumSamples(0);
        }
      }
    }


    if (UseDrift == "rn")
    {
      RealType avg_w(0);
      RealType n_w(0);
#pragma omp parallel
      {
        int ip=omp_get_thread_num();

        for (int step=0; step<myWarmupSteps; ++step)
        {
          avg_w=0;
          n_w=0;
          for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
          {
            Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
#pragma omp single
            {
              MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
              while (wit!=wit_end)
              {
                avg_w += (*wit)->Weight;
                n_w +=1;
                wit++;
              }
            }
#pragma omp barrier
          }
#pragma omp single
          {
            avg_w *= 1.0/n_w;
            RealType w_m = avg_w/(1.0-avg_w);
            w_m = std::log(0.5+0.5*w_m);
            if (std::abs(w_m)>0.01)
              logepsilon += w_m;
          }
#pragma omp barrier
          Movers[ip]->setLogEpsilon(logepsilon);
        }

        for (int prestep=0; prestep<myWarmupSteps; ++prestep)
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);

        if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
          Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);

#pragma omp critical
        {
          wClones[ip]->clearEnsemble();
          wClones[ip]->setNumSamples(0);
        }
      }
    }

  }

  void VMCLinearOptOMP::fillMatrices(Matrix<RealType>& H2, Matrix<RealType>& Hamiltonian, Matrix<RealType>& Variance, Matrix<RealType>& Overlap)
  {
    RealType nrm = 1.0/sW;
    //     RealType nrm2 = nrm*nrm;
    for (int i=0; i<NumOptimizables; i++)
    {
      HDiE[i]*= nrm;
      HDi[i] *= nrm;
      DiE2[i]*= nrm;
      DiE[i] *= nrm;
      Di[i]  *= nrm;
    }
    HDiHDj*= nrm;
    DiHDjE*= nrm;
    DiHDj *= nrm;
    DiDjE2*= nrm;
    DiDjE *= nrm;
    DiDj  *= nrm;

    RealType H2_avg = sE2*nrm;
    E_avg = sE*nrm;
    V_avg = H2_avg - E_avg*E_avg;


    for (int pm=0; pm<NumOptimizables; pm++)
    {
      RealType wfe = HDi[pm] + DiE[pm]-Di[pm]*E_avg;
      RealType wfm = HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg;
      //         Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;

      H2(0,pm+1) = HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg;
      H2(pm+1,0) = HDiE[pm] + DiE2[pm]-DiE[pm]*E_avg;

      //         HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
      RealType vterm = HDiE[pm]-HDi[pm]*E_avg + DiE2[pm]-Di[pm]*H2_avg -2.0*E_avg*(DiE[pm]-Di[pm]*E_avg);
      Variance(0,pm+1) = vterm;
      Variance(pm+1,0) = vterm;

      Hamiltonian(0,pm+1) = wfe;
      Hamiltonian(pm+1,0) = DiE[pm]-Di[pm]*E_avg;

      for (int pm2=0; pm2<NumOptimizables; pm2++)
      {
        //           H2(pm+1,pm2+1) += wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new - curAvg_w));
        //           Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
        //           Variance(pm+1,pm2+1) += wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
        //           Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);

        //        Symmetric  (HDi[pm] + DiE[pm]-Di[pm]*E_avg)(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg)
        H2(pm+1,pm2+1) = HDiHDj(pm,pm2) + DiHDjE(pm2,pm) - DiHDj(pm2,pm)*E_avg
          + DiHDjE(pm,pm2) + DiDjE2(pm,pm2) - DiDjE(pm,pm2)*E_avg
          + E_avg*(DiHDj(pm,pm2) + DiDjE(pm,pm2) - DiDj(pm,pm2)*E_avg);
        //        Non-symmetric    (Dsaved[pm]-D_avg[pm])*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w))
        Hamiltonian(pm+1,pm2+1) = DiHDj(pm,pm2) + DiDjE(pm,pm2) - Di[pm2]*DiE[pm] - Di[pm]*(HDi[pm2] + DiE[pm2]-Di[pm2]*E_avg);
        //        Symmetric  (HDi[pm] - 2.0*DiE[pm] + 2.0*Di[pm]*E_avg)*( HDi[pm2] - 2.0* DiE[pm2]+2.0*Di[pm2]*E_avg)
        Variance(pm+1,pm2+1) = HDiHDj(pm,pm2) -2.0*DiHDjE(pm2,pm) +2.0*DiHDj(pm,pm2)*E_avg
          -2.0*( DiHDjE(pm,pm2) - 2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2))
          +2.0*E_avg*(DiHDj(pm,pm2) -2.0*DiDjE(pm,pm2) +2.0*E_avg*DiDj(pm,pm2));
        //        Symmetric
        Overlap(pm+1,pm2+1) = DiDj(pm,pm2)-Di[pm]*Di[pm2];

      }
    }

    Hamiltonian(0,0) = E_avg;
    Overlap(0,0) = 1.0;
    H2(0,0) = H2_avg;
    Variance(0,0) = V_avg;

    for (int pm=1; pm<NumOptimizables+1; pm++)
      for (int pm2=1; pm2<NumOptimizables+1; pm2++)
        Variance(pm,pm2) += V_avg*Overlap(pm,pm2);

    //     app_log()<<V_avg<<"  "<<E_avg<<"  "<<sW<<endl;
  }

  VMCLinearOptOMP::RealType VMCLinearOptOMP::fillComponentMatrices()
  {
    int n(NumOptimizables);
    ///These are the values we collect to build the Matrices LOCAL
    Matrix<RealType> lHDiHDj(n,n), lDiHDj(n,n), lDiHDjE(n,n), lDiDj(n,n), lDiDjE(n,n), lDiDjE2(n,n);
    lHDiHDj=0; lDiHDj=0; lDiHDjE=0; lDiDj=0; lDiDjE=0; lDiDjE2=0;


    std::vector<RealType> lHDi(n), lHDiE(n), lDi(n), lDiE(n), lDiE2(n);
    for (int i=0; i<NumOptimizables; i++)
    {
      lHDi[i]=0;
      lHDiE[i]=0;
      lDi[i]=0;
      lDiE[i]=0;
      lDiE2[i]=0;
    }
    RealType lsE(0),lsE2(0),lsE4(0),lsW(0),lsN(0);

    for (int ip=0; ip<NumThreads; ip++)
    {
      RealType E_L = W[ip]->getPropertyBase()[LOCALENERGY];
      RealType E_L2= E_L*E_L;
      RealType wW  = W[ip]->Weight;
      lsE +=E_L*wW;
      lsE2+=E_L2*wW;
      lsE4+=E_L2*E_L2*wW;
      lsW +=wW;
      lsN+=1;
    }



    for (int ip=0; ip<NumThreads; ip++)
    {
      RealType E_L = (*W[ip]).getPropertyBase()[LOCALENERGY];
      RealType E_L2= E_L*E_L;
      RealType wW  = (*W[ip]).Weight;
      for (int i=0; i<NumOptimizables; i++)
      {
        RealType di  = DerivRecords(ip,i);
        RealType hdi = HDerivRecords(ip,i);
        //             vectors
        lHDiE[i]+= wW*E_L* hdi;
        lHDi[i] += wW*     hdi;
        lDiE2[i]+= wW*E_L2*di;
        lDiE[i] += wW*E_L* di;
        lDi[i]  += wW*     di;

        for (int j=0; j<NumOptimizables; j++)
        {
          RealType dj  = DerivRecords(ip,j);
          RealType hdj = HDerivRecords(ip,j);

          lHDiHDj(i,j) += wW*    hdi*hdj;
          lDiHDjE(i,j) += wW* E_L*di*hdj;
          lDiHDj(i,j)  += wW*     di*hdj;
          lDiDjE2(i,j) += wW*E_L2*di*dj;
          lDiDjE(i,j)  += wW* E_L*di*dj;
          lDiDj(i,j)   += wW*     di*dj;
        }
      }
    }

    //         //Lazy. Pack these for better performance.
    //         myComm->allreduce(lsE);
    //         myComm->allreduce(lsE2);
    //         myComm->allreduce(lsE4);
    //         myComm->allreduce(lsW);
    //         myComm->allreduce(lHDiE);
    //         myComm->allreduce(lHDi);
    //         myComm->allreduce(lDiE2);
    //         myComm->allreduce(lDiE);
    //         myComm->allreduce(lDi);
    //         myComm->allreduce(lHDiHDj);
    //         myComm->allreduce(lDiHDjE);
    //         myComm->allreduce(lDiHDj);
    //         myComm->allreduce(lDiDjE2);
    //         myComm->allreduce(lDiDjE);
    //         myComm->allreduce(lDiDj);
    if(myComm->size() > 1)
    {
      Walker_t::Buffer_t tmpBuffer;
      tmpBuffer.rewind();
      tmpBuffer.add(lsE);
      tmpBuffer.add(lsE2);
      tmpBuffer.add(lsE4);
      tmpBuffer.add(lsW);
      tmpBuffer.add(lsN);

      tmpBuffer.add(lHDiE.begin(),lHDiE.end());
      tmpBuffer.add(lHDi.begin(),lHDi.end());
      tmpBuffer.add(lDiE2.begin(),lDiE2.end());
      tmpBuffer.add(lDiE.begin(),lDiE.end());
      tmpBuffer.add(lDi.begin(),lDi.end());

      tmpBuffer.add(lHDiHDj.begin(),lHDiHDj.end());
      tmpBuffer.add(lDiHDjE.begin(),lDiHDjE.end());
      tmpBuffer.add(lDiHDj.begin(),lDiHDj.end());
      tmpBuffer.add(lDiDjE2.begin(),lDiDjE2.end());
      tmpBuffer.add(lDiDjE.begin(),lDiDjE.end());
      tmpBuffer.add(lDiDj.begin(),lDiDj.end());

      myComm->allreduce(tmpBuffer);
      tmpBuffer.rewind();

      tmpBuffer.get(lsE);
      tmpBuffer.get(lsE2);
      tmpBuffer.get(lsE4);
      tmpBuffer.get(lsW);
      tmpBuffer.get(lsN);

      tmpBuffer.get(lHDiE.begin(),lHDiE.end());
      tmpBuffer.get(lHDi.begin(),lHDi.end());
      tmpBuffer.get(lDiE2.begin(),lDiE2.end());
      tmpBuffer.get(lDiE.begin(),lDiE.end());
      tmpBuffer.get(lDi.begin(),lDi.end());

      tmpBuffer.get(lHDiHDj.begin(),lHDiHDj.end());
      tmpBuffer.get(lDiHDjE.begin(),lDiHDjE.end());
      tmpBuffer.get(lDiHDj.begin(),lDiHDj.end());
      tmpBuffer.get(lDiDjE2.begin(),lDiDjE2.end());
      tmpBuffer.get(lDiDjE.begin(),lDiDjE.end());
      tmpBuffer.get(lDiDj.begin(),lDiDj.end());    
    }

    //add locals to globals
    sE +=lsE;
    sE2+=lsE2;
    sE4+=lsE4;
    sW +=lsW;
    sN +=lsN;
    for (int j=0; j<NumOptimizables; j++)
    {
      HDiE[j]+=lHDiE[j];
      HDi[j] +=lHDi[j] ;
      DiE2[j]+=lDiE2[j];
      DiE[j] +=lDiE[j] ;
      Di[j]  +=lDi[j]  ;
    }

    HDiHDj += lHDiHDj;
    DiHDjE += lDiHDjE;
    DiHDj  += lDiHDj ;
    DiDjE2 += lDiDjE2;
    DiDjE  += lDiDjE ;
    DiDj   += lDiDj  ;

    RealType nrm = 1.0/sW;
    E_avg = nrm*sE;
    V_avg = nrm*sE2-E_avg*E_avg;
    //         app_log()<<V_avg<<"  "<<E_avg<<"  "<<std::log(sW)<<endl;
    RealType g_nrm = 1.0/sN;
    RealType err_E(std::sqrt( ((V_avg<0.0)?(1.0):(V_avg*g_nrm)) ));
    RealType err_E2(nrm*sE4-nrm*nrm*sE2*sE2);
    err_E2 *= g_nrm;
    err_E2 = std::sqrt( ((err_E2<0.0)?(1.0):(err_E2)) );

    return w_beta*err_E2+(1.0-w_beta)*err_E;
  }


  bool
    VMCLinearOptOMP::put(xmlNodePtr q)
    {
      //nothing to add
      return true;
    }
}

/***************************************************************************
 * $RCSfile: VMCLinearOptOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCLinearOptOMP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/

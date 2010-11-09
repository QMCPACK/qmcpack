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
                                 HamiltonianPool& hpool):
        QMCDriver(w,psi,h),  CloneManager(hpool),
        myWarmupSteps(0),UseDrift("yes"), NumOptimizables(0), w_beta(0.0), GEVtype("mixed"), logoffset(2.0), logepsilon(0)
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
}

bool VMCLinearOptOMP::run(bool needMatrix)
{
    resetRun();
    std::vector<opt_variables_type> dummyOptVars;
    if (needMatrix)
    {
        
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
        
    }
    else
    {
        for (int ip=0; ip<NumThreads; ++ip)
            for (int prestep=0; prestep<myWarmupSteps; ++prestep)
                Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
    }

    RealType target_errorbars;
    if (needMatrix) target_errorbars = beta_errorbars;
    else target_errorbars = alpha_errorbars;
    RealType errorbars = target_errorbars+1;
    CurrentStep=0;
    int CurrentBlock=0;
    while ((errorbars>target_errorbars)&&(CurrentBlock<nBlocks))
    {
#pragma omp parallel for
        for (int ip=0; ip<NumThreads; ++ip)
        {
            //assign the iterators and resuse them
            MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
            if (needMatrix) Movers[ip]->startBlock(nSteps);
            int now_loc=CurrentStep;
            //rest the collectables and keep adding
            if (needMatrix) wClones[ip]->resetCollectables();
            //rest the collectables and keep adding
            for (int step=0; step<nSteps; ++step)
            {
                Movers[ip]->advanceWalkers(wit,wit_end,false);
                if (needMatrix) Movers[ip]->accumulate(wit,wit_end);
                ++now_loc;
            }
            if (needMatrix) Movers[ip]->stopBlock(false);
        }//end-of-parallel for
        CurrentStep+=nSteps;

        if (needMatrix)
        {
          Estimators->accumulateCollectables(wClones,nSteps);
        
          Estimators->stopBlock(estimatorClones);
#pragma omp parallel for
            for (int ip=0; ip<NumThreads; ++ip)
            {
                vector<RealType> Dsaved(NumOptimizables);
                vector<RealType> HDsaved(NumOptimizables);
                psiClones[ip]->evaluateDerivatives(*wClones[ip],dummyOptVars[ip],Dsaved,HDsaved);
                std::copy(Dsaved.begin(),Dsaved.end(),&DerivRecords(ip,0));
                std::copy(HDsaved.begin(),HDsaved.end(),&HDerivRecords(ip,0));
            }
        }
        errorbars = fillComponentMatrices(needMatrix);
        CurrentBlock++;

    }//block
    app_log()<<" Blocks used   : "<<CurrentBlock<<endl;
    app_log()<<" Errorbars are : "<<errorbars<<endl;

    if (needMatrix) Estimators->stop(estimatorClones);
    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip)
        *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    //finalize a qmc section
    if (needMatrix) return finalize(nBlocks);
    else return true;
}

void VMCLinearOptOMP::initCS()
{
// #pragma omp parallel for
  
  //resetting containersx
    CorrelatedH.resize(NumThreads,NumThreads);
    Norms.resize(NumThreads+1); Norm2s.resize(NumThreads+1,NumThreads+1);
    Energies.resize(NumThreads);
    CorrelatedH=0;
    Norm2s=0;
    NE_i.resize(NumThreads);
    for (int ip=0; ip<NumThreads+1; ++ip) Norms[ip]=0;
    for (int ip=0; ip<NumThreads; ++ip) Energies[ip]=0;
    for (int ip=0; ip<NumThreads; ++ip) NE_i[ip]=0;
    
//         set all walker positions to the same place
    Walker_t& firstWalker(*W[0]);
    for (int ip=1; ip<NumThreads; ++ip)
    {
      (*W[ip]).makeCopy(firstWalker);
    }

    for (int ip=0; ip<NumThreads; ++ip)
      {
        Walker_t& thisWalker(*W[ip]);
        wClones[ip]->loadWalker(thisWalker,true);

        Walker_t::Buffer_t tbuffer;
        RealType logpsi=psiClones[ip]->evaluateLog(*wClones[ip]);
        logpsi=psiClones[ip]->registerData(*wClones[ip],tbuffer);
//               logpsi=psiClones[ip]->updateBuffer(*wClones[ip],tbuffer,true);
        thisWalker.DataSet=tbuffer;
        thisWalker.Weight = 1.0;
        RealType ene = hClones[ip]->evaluate( *wClones[ip]);
//         app_log()<<ene<<" "<<logpsi<<endl;
        thisWalker.resetProperty(logpsi,psiClones[ip]->getPhase(),ene);
        hClones[ip]->saveProperty(thisWalker.getPropertyBase());
        wClones[ip]->saveWalker(thisWalker);
      }
    
    for (int prestep=0; prestep<myWarmupSteps; ++prestep)
    {
      Movers[0]->advanceCSWalkers(psiClones, wClones, hClones);
    }
}

int VMCLinearOptOMP::runCS(vector<vector<RealType> > bestParams, RealType& errorbars)
{
    for (int ip=0; ip<NumThreads; ++ip)
    {
      opt_variables_type dummy;
      psiClones[ip]->checkInVariables(dummy);
      dummy.resetIndex();
      psiClones[ip]->checkOutVariables(dummy);
      for (int i=0;i<bestParams[0].size();i++)  dummy[i] = bestParams[ip][i];
      psiClones[ip]->resetParameters(dummy);
    }
    initCS();
  
    errorbars=alpha_errorbars+1;
    CurrentStep=0;
    CSBlock=0;
//     run long enough to get accurate errorbars ~4 blocks.
//     run until errorbars are small enough or when the energy difference is resolved to with 2 errorbars.
//     max run is defined by nBlocks
    while ((CSBlock<4)||
      (((errorbars>alpha_errorbars)&&((NE_i[nE]-NE_i[minE])/2.0>errorbars))&&(CSBlock<nBlocks) ))
    {
        int now_loc=CurrentStep;
        for (int step=0; step<nSteps; ++step)
        {
            Movers[0]->advanceCSWalkers(psiClones, wClones, hClones);
            ++now_loc;
        }
        CurrentStep+=nSteps;
        errorbars = estimateCS();
        CSBlock++;

    }//block
    app_log()<<" Blocks used   : "<<CSBlock<<endl;
    app_log()<<" Errorbars are : "<<errorbars<<endl;
//     app_log()<<" Min E["<<minE<<"] estimate: "<<NE_i[minE]<<endl;
    
    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip)
        *(RandomNumberControl::Children[ip])=*(Rng[ip]);
    if (std::abs(NE_i[minE])<1e6)  return minE;
    else return -1;
}

bool VMCLinearOptOMP::bracketing(vector<RealType>& lambdas, RealType errorbars)
{
      //Do some bracketing and line searching if we need to
    RealType dl = std::abs(lambdas[1]-lambdas[0]);
    RealType mL= lambdas[minE];
    RealType DE = NE_i[nE] - NE_i[minE];
    
    
    
    if (moved_left&&moved_right&&(DE<errorbars))
    {
      app_log()<<"   Decrease target error bars to resolve energy difference."<<endl;
      return false;
    }
    else if (minE==(NumThreads-1))
    {
      if(moved_left)
      {
        app_log()<<" Bracketed minimum between CS runs"<<endl;
        moved_right=true;
        mL=lambdas[minE]-dl;
        dl = 2.0*std::abs(lambdas[0]-lambdas[1])/(NumThreads-1.0);
        for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
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
        moved_left=true;
        mL = lambdas[minE]-dl;
        dl = 2.0*std::abs(lambdas[1]-lambdas[0])/(NumThreads-1.0);
        for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + ip*dl;
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
// if energy difference is smaller than the errorbars we computed then we are done
        if (DE<errorbars)
        {
          int nms=3;
          vector<RealType>  Y(nms), Coefs(3);
          Matrix<RealType> X(nms,3);
          for (int i=0; i<nms; i++) X(i,0)=1.0;
          for (int i=0; i<nms; i++) X(i,1)=lambdas[i+minE-1];
          for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
          for (int i=0; i<nms; i++) Y[i]=NE_i[i+minE-1];
          LinearFit(Y,X,Coefs);
          
          RealType quadraticMinimum(-1.0*Coefs[1]/Coefs[2]);
          lambdas[minE]=quadraticMinimum;
          
          return false;
        }
        else
        {
          app_log()<<" Bracketed minimum, refine"<<endl;
          moved_left=moved_right=false;
// energy difference between the points is still larger than the error bars we require
// need to "zoom" into find minimum more precisely
            dl = std::abs(lambdas[0]-lambdas[1])/(NumThreads-1.0);
            mL = std::min(lambdas[minE],lambdas[nE]);
            for (int ip=0; ip<NumThreads; ++ip) lambdas[ip] = mL + dl*ip;
        }
    }
    return true;
}

VMCLinearOptOMP::RealType VMCLinearOptOMP::runCS(vector<RealType> curParams, vector<RealType> curDir, vector<RealType>& lambdas)
{
  bool notConverged(true);
  
  moved_right=false;
  moved_left=false;
  
  while (notConverged)
  {
    vector<vector<RealType> > dummy(NumThreads,std::vector<RealType>(curParams.size()));
    for (int ip=0; ip<NumThreads; ++ip) for (int i=0;i<curParams.size();i++)  dummy[ip][i] = curParams[i] + lambdas[ip]*curDir[i+1];
    
    RealType errorbars;
    int bombed = runCS(dummy, errorbars);
    if (bombed<0){
      lambdas[0]=0;
      return 0;
    }
    for (int ip=0; ip<NumThreads; ++ip) app_log()<<"E["<<lambdas[ip]<<"] estimate: "<<NE_i[ip]<<endl;
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
      notConverged = bracketing(lambdas, errorbars); 
  }

    lambdas[0]=lambdas[minE];
    return NE_i[minE];
}

VMCLinearOptOMP::RealType VMCLinearOptOMP::estimateCS()
{
  vector<long double> e_i(NumThreads), psi2_i(NumThreads);
  long double psi2(0);
  if (CSBlock==0)
  {
    if (!myComm->rank())
      logpsi2_0_0 = 2.0*(W[0])->getPropertyBase()[LOGPSI];
    myComm->bcast(logpsi2_0_0);
  }
  
  for (int ip=0; ip<NumThreads; ip++)
  {
    e_i[ip]    = (W[ip])->getPropertyBase()[LOCALENERGY];
    psi2_i[ip] = expl(2.0*(W[ip])->getPropertyBase()[LOGPSI] - logpsi2_0_0);
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
    
  // global quantities for mpi collection
  std::vector<RealType> gEnergies(Energies), gNorms(Norms);
  Matrix<RealType> gNorm2s(Norm2s), gCorrelatedH(CorrelatedH);
  myComm->allreduce(gEnergies);
  myComm->allreduce(gNorms);
  myComm->allreduce(gNorm2s);
  myComm->allreduce(gCorrelatedH);
  
//   Here are the global energy estimates
  for (int ip=0; ip<NumThreads; ip++) NE_i[ip] = gEnergies[ip]/gNorms[ip];

//   find lowest energy
  minE=0;
  for (int ip=1; ip<NumThreads; ip++) if (NE_i[ip]<NE_i[minE]) minE=ip;
  
//   nE is the next lowest energy
  nE=minE;
  if (minE==0) nE=1;
  else if (minE==NumThreads-1) nE=NumThreads-2;
  else nE=(NE_i[minE+1]>NE_i[minE-1])?minE-1:minE+1;

//return the error in the energy differences between lowest two
  RealType rval = (gCorrelatedH(minE,minE)/(gNorms[minE]*gNorms[minE]) + gCorrelatedH(nE,nE)/(gNorms[nE]*gNorms[nE]) - 2.0*gCorrelatedH(minE,nE)/(gNorms[minE]*gNorms[nE])) - (NE_i[minE]-NE_i[nE])*(NE_i[minE]-NE_i[nE]);
  
  rval = std::sqrt(rval/(CSBlock+1));
  return rval;
}

void VMCLinearOptOMP::resetRun()
{
    makeClones(W,Psi,H);
    app_log() << "  Warmup Steps " << myWarmupSteps << endl;


    if (Movers.empty())
    {
        Movers.resize(NumThreads,0);
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
              Movers[ip]=new VMCUpdatePbyPSampleRN(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
              Movers[ip]->setLogEpsilon(logepsilon);
              // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            else if (UseDrift == "yes")
            {
              os <<"  PbyP moves with drift, using VMCUpdatePbyPWithDriftFast"<<endl;
              Movers[ip]=new VMCUpdatePbyPWithDriftFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
              // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            else
            {
              os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<<endl;
              Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
          }
          else
          {
            if (UseDrift == "rn")
            {
              os <<"  walker moves with RN, using VMCUpdateAllSampleRN"<<endl;
              Movers[ip] =new VMCUpdateAllSampleRN(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
              Movers[ip]->setLogEpsilon(logepsilon);
            }
            else if (UseDrift == "yes")
            {
              os <<"  walker moves with drift, using VMCUpdateAllWithDriftFast"<<endl;
              Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            else
            {
              os <<"  walker moves with |psi|^2, using VMCUpdateAll"<<endl;
              Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
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
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);

        if (QMCDriverMode[QMC_UPDATE_MODE])
          Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        else
          Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);

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
    
    if ((UseDrift == "rn")&&(logepsilon==0.0))
    {
      long double psi2_sum(0.0);
      long double psi4_sum(0.0);
      RealType psi2_0_0(0.0);
      RealType nw_local(0.0);
      
      if (myRNWarmupSteps==0) myRNWarmupSteps=myWarmupSteps;
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        Movers[ip]->setLogEpsilon(logepsilon);
        for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
        {
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
          #pragma omp flush
          if((ip==0)&&(prestep>5))
          {
            MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
            if (prestep==0)
            {
              if (myComm->rank()==0) psi2_0_0=2.0* (*wit)->getPropertyBase()[LOGPSI];
              myComm->bcast(psi2_0_0);
            }
            nw_local += wit_end-wit;
            while (wit!=wit_end)
            {
              long double t(expl(2.0*(*wit)->getPropertyBase()[LOGPSI]-psi2_0_0));
              psi2_sum += t;
              psi4_sum += t*t;
              wit++;
            }
//             app_log()<<logl(psi2_sum/nw_local)<<"  "<<logl(sqrtl((psi4_sum/nw_local - (psi2_sum/nw_local * psi2_sum/nw_local))/(1+prestep)) )<<endl;
          }
          #pragma omp barrier
          
        }
      }

//           we need to set the logepsilon to something reasonable
        myComm->allreduce(psi2_sum);
        myComm->allreduce(psi4_sum);
        myComm->allreduce(nw_local);
        nw_local = 1.0/nw_local;

        long double p2(psi2_sum*nw_local);
        RealType eps0 = logl(p2) + psi2_0_0 - logoffset;
        app_log()<<" Using logepsilon "<<eps0<<endl;
        app_log()<<" log<Psi^2>= "<<logl(p2)+ psi2_0_0<<endl;
        long double sigp2(psi4_sum*nw_local-p2*p2);
        app_log()<<" logError  = "<< logl(sqrtl(sigp2*nw_local))<<endl;
        for (int ip=0; ip<NumThreads; ++ip) Movers[ip]->setLogEpsilon(eps0);
        #pragma omp parallel
        {
          int ip=omp_get_thread_num();
          for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
          {
            Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
          }
        }
//         bool badEPS(true);
//         RealType stabilizer=0;
//         RealType eps0 = logl(psi2_sum*nw_local) + psi2_0_0;
// 
//         int Accept(0), Reject(0);
//           for (int ip=0; ip<NumThreads; ++ip) Accept+=Movers[ip]->nAccept;
//           for (int ip=0; ip<NumThreads; ++ip) Reject+=Movers[ip]->nReject;
//           myComm->allreduce(Accept);
//           myComm->allreduce(Reject);
//           RealType accpt=(1.0*Accept)/(Accept+Reject);
//         app_log()<<"Pre: "<<accpt<<endl;
//         while (badEPS)
//         {
//           eps0 = logl(psi2_sum*nw_local) + psi2_0_0 + stabilizer;
//             #pragma omp parallel
//             {
//               int ip=omp_get_thread_num();
//               Movers[ip]->setLogEpsilon(eps0);
//               for (int prestep=0; prestep<myWarmupSteps; ++prestep)
//               {
//                 Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//                 if(ip==0)
//                 {
//                   MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//                   nw_local += wit_end-wit;
//                   while (wit!=wit_end)
//                   {
//                     psi2_sum += expl(2.0*(*wit)->getPropertyBase()[LOGPSI]-psi2_0_0);
//                     wit++;
//                   }
//                 }
//                 #pragma omp barrier
//               }
//             }
//           myComm->allreduce(psi2_sum);
//           myComm->allreduce(nw_local);
//         
//           for (int ip=0; ip<NumThreads; ++ip) 
//             {
//               Movers[ip]->nAccept=0;
//               Movers[ip]->nReject=0;
//               Accept=0;
//               Reject=0;
//             }
//             #pragma omp parallel
//             {
//               int ip=omp_get_thread_num();
//               for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//                 Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//             }
//             for (int ip=0; ip<NumThreads; ++ip) Accept+=Movers[ip]->nAccept;
//             for (int ip=0; ip<NumThreads; ++ip) Reject+=Movers[ip]->nReject;
//             myComm->allreduce(Accept);
//             myComm->allreduce(Reject);
//             RealType accptE=(1.0*Accept)/(Accept+Reject);
//             app_log()<<Accept<<"  "<<Reject<<endl;
//             app_log()<<stabilizer<<"  "<<accptE<<endl;
//             if (accptE<rn_accept_target) badEPS=false;
//             else if (accptE>99) stabilizer -= 8.0;
//             else if (accptE>rn_accept_target+0.3) stabilizer -=1.0;
//             else if (accptE>rn_accept_target+0.2) stabilizer -=0.1;
//             else stabilizer -= 0.01;
//         }
//         app_log()<<"  Using LogEpsilon of: "<<eps0<<endl;
    }
    else if (UseDrift == "rn")
    {
      app_log()<<"  Using LogEpsilon of: "<<logepsilon<<endl;
#pragma omp parallel
      {
        int ip=omp_get_thread_num();
        Movers[ip]->setLogEpsilon(logepsilon);
        for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
        {
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
        }
      }
    }
      
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
    //Used to debug and benchmark opnemp
    //#pragma omp parallel for
    //    for(int ip=0; ip<NumThreads; ip++)
    //    {
    //      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
    //    }
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

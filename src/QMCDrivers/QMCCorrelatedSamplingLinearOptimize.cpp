//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/QMCCorrelatedSamplingLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCApp/HamiltonianPool.h"
#include "Numerics/Blasf.h"
#include <cassert>
#include "Numerics/LinearFit.h"
#include "Optimize/NRCOptimization.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

QMCCSLinearOptimize::QMCCSLinearOptimize(MCWalkerConfiguration& w,
        TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCLinearOptimize(w,psi,h,hpool,ppool), psipool(ppool),
        vmcCSEngine(0), Max_iterations(1), exp0(-16), nstabilizers(10), stabilizerScale(0.5), bigChange(1), w_beta(0.0),
        MinMethod("quartic"), GEVtype("mixed")
{
    //set the optimization flag
    QMCDriverMode.set(QMC_OPTIMIZE,1);
    //read to use vmc output (just in case)
    RootName = "pot";
    QMCType ="QMCCSLinearOptimize";
    m_param.add(Max_iterations,"max_its","int");
    m_param.add(nstabilizers,"nstabilizers","int");
    m_param.add(stabilizerScale,"stabilizerscale","double");
    m_param.add(bigChange,"bigchange","double");
    m_param.add(w_beta,"beta","double");
    stepsize=-1.0;
    m_param.add(stepsize,"stepsize","double");
    m_param.add(exp0,"exp0","double");
    m_param.add(MinMethod,"MinMethod","string");
    m_param.add(GEVtype,"GEVMethod","string");
    //Set parameters for line minimization:
}

/** Clean up the vector */
QMCCSLinearOptimize::~QMCCSLinearOptimize()
{
}

bool QMCCSLinearOptimize::run()
{
    start();
    bool Valid(true);
    int Total_iterations(0);

//size of matrix
    numParams = optTarget->NumParams();
    N = numParams + 1;

//  solve CSFs and other parameters separately then rescale elements accordingly
    int first,last;
    getNonLinearRange(first,last);

//     initialize our parameters
    vector<RealType> currentParameterDirections(N,0);
    vector<RealType> currentParameters(numParams,0);
    for (int i=0; i<numParams; i++) currentParameters[i] = optTarget->Params(i);
    optdir.resize(numParams,0);
    optparm.resize(numParams,0);

    
    Matrix<RealType> Ham(N,N);
    Matrix<RealType> Ham2(N,N);
    Matrix<RealType> Var(N,N);
    Matrix<RealType> S(N,N);
    vmcCSEngine->fillMatrices(Ham2,Ham,Var,S);

    vector<RealType> bestParameters(currentParameters);


//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
    int tooManyTries(200);

    Matrix<RealType> Left(N,N);
    Matrix<RealType> Left_tmp(N,N);
    Matrix<RealType> Right(N,N);
    
    vector<std::pair<RealType,RealType> > mappedStabilizers;
    vector<vector<RealType> > savedCSparameters;
          
    if (GEVtype=="H2")
    {
        Left_tmp=Ham;
        RealType H2rescale=1.0/Ham2(0,0);
        Right=(1-w_beta)*S + w_beta*H2rescale*Ham2;
    }
    else
    {
        Right=S;
        Left_tmp=(1.0-w_beta)*Ham + w_beta*Var;
    }
    
    bool apply_inverse(true);
    if(apply_inverse)
    {
      invert_matrix(Right,false);
      MatrixOperators MO;
      MO.product(Right,Left_tmp,Left);
    }
    else
      Left=Left_tmp;

    //Find largest off-diagonal element compared to diagonal element.
    //This gives us an idea how well conditioned it is and can be used to stabilize.
    RealType od_largest(0);
    for (int i=0; i<N; i++) for (int j=0; j<N; j++)
            od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
    if (od_largest>0) od_largest = std::log(od_largest);
    else od_largest=1.0;

    RealType safe = Left(0,0);
    RealType XS(0);

    if (nstabilizers<=omp_get_max_threads()+1) nstabilizers=omp_get_max_threads()+1;
    else
    {
      nstabilizers -= omp_get_max_threads()+1;
      int Ns(nstabilizers/(omp_get_max_threads()-1));
      nstabilizers = omp_get_max_threads() + 1 + Ns*(omp_get_max_threads()-1);
    }
    
    for (int stability=0; stability<nstabilizers; stability++)
    {
      app_log()<<"Iteration: "<<stability+1<<"/"<<nstabilizers<<endl;
        Matrix<RealType> LeftT(N,N), RightT(N,N);
        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                LeftT(i,j)= Left(j,i);
                RightT(i,j)= Right(j,i);
            }


        
        if(savedCSparameters.size()==omp_get_max_threads())
        {
          app_log()<<"   Choosing best"<<endl;
          RealType error(0);
          int bestP = vmcCSEngine->runCS(savedCSparameters,error);
          if (bestP<0)
          {
            app_log()<<"   Error in CS cost function. Unchanged parameters."<<endl;
            bestP=0;
            for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i];
            finish();
            return false;
          }
          
          vector<RealType> cs(numParams);
          for (int i=0; i<numParams; i++) cs[i] = bestParameters[i] = savedCSparameters[bestP][i];
          savedCSparameters.clear();
          savedCSparameters.push_back(cs);
          std::pair<RealType,RealType> ms;
          ms.first=mappedStabilizers[bestP].first;
          ms.second=mappedStabilizers[bestP].second;
          
          if (stability==nstabilizers-1) continue;
//           if (MinMethod=="rescale")
//           {
            int nms=mappedStabilizers.size();
            std::vector<RealType> csts(nms);
            vmcCSEngine->getDeltaCosts(csts);
            
            if (nms>=5)
            {//Quartic fit the stabilizers we have tried and try to choose the best we can
              vector<RealType>  Y(nms), Coefs(5);
              Matrix<RealType> X(nms,5);
              for (int i=0; i<nms; i++) X(i,0)=1.0;
              for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
              for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
              for (int i=0; i<nms; i++) X(i,3)=X(i,2)*X(i,1);
              for (int i=0; i<nms; i++) X(i,4)=X(i,3)*X(i,1);
              for (int i=0; i<nms; i++) Y[i]=csts[i];
              LinearFit(Y,X,Coefs);
  //lowest we will allow is a little less than the bare base stabilizer
              RealType dltaBest=std::max(stabilityBase-0.1, QuarticMinimum(Coefs));
              XS = dltaBest;
            }
            else
            {//Quadratic fit the stabilizers we have tried and try to choose the best we can
              vector<RealType>  Y(nms), Coefs(3);
              Matrix<RealType> X(nms,3);
              for (int i=0; i<nms; i++) X(i,0)=1.0;
              for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
              for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
              for (int i=0; i<nms; i++) Y[i]=csts[i];
              LinearFit(Y,X,Coefs);
              
              RealType quadraticMinimum(-1.0*Coefs[1]/Coefs[2]);
              RealType dltaBest=std::max(stabilityBase-0.1, quadraticMinimum);
  //               app_log()<<"smallest XS:      "<<X(0,1)<<endl;
  //               app_log()<<"quadraticMinimum: "<<quadraticMinimum<<endl;
              XS = dltaBest;
            }
//           }
//           else
//             XS=0;
          mappedStabilizers.clear();
          mappedStabilizers.push_back(ms);
        }
        else
          XS=0;
            
        if (XS==0)
        {
            XS     = std::exp(stabilityBase +  stability*od_largest/nstabilizers);
            for (int i=1; i<N; i++) LeftT(i,i) += XS;
        }
        else
        {
          //else XS is from the fit
          for (int i=1; i<N; i++) LeftT(i,i) += std::exp(XS);
        }
        
        RealType lowestEV;
        myTimers[2]->start();
        if(apply_inverse)
          lowestEV =getLowestEigenvector(LeftT,currentParameterDirections);
        else
          lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
        myTimers[2]->stop();

//         RealType Lambda_Last(Lambda);
        myTimers[3]->start();
//         if (GEVtype=="H2")
//         {
// //           the rescaling isn't right for this?
//           RealType bigVec(0);
//           for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
//           Lambda=0.5*bigChange/bigVec;
//         }
//           else 
        Lambda = getNonLinearRescale(currentParameterDirections,S);
        
          myTimers[3]->stop();
//         app_log()<<"Computed Lambda is: "<<Lambda<<endl;
        RealType bigVec(0);
        for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
        if (Lambda*bigVec>bigChange)
        {
            app_log()<<"  Failed Step. Largest parameter change: "<<Lambda*bigVec<<endl;
            tooManyTries--;
            if (tooManyTries>0)
            {
                if(stability==0) stabilityBase+=stabilizerScale;
                else stabilityBase-=stabilizerScale;
                stability-=1;
                app_log()<<" Re-run with larger stabilityBase"<<endl;
                continue;
            }
        }

        if (MinMethod=="rescale")
        {
             vector<RealType> cs(numParams);
             for (int i=0; i<numParams; i++) cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
             app_log()<<" Lambda: "<<Lambda<<endl;
             app_log()<<" Largest parameter change: "<<Lambda*bigVec<<endl;
//              optTarget->resetPsi(false);
             savedCSparameters.push_back(cs);
        }
        else
        {
             int nthreads = omp_get_max_threads();
             std::vector<RealType> lambdas(nthreads);
             if (stepsize>0) Lambda = stepsize/bigVec;
             for(int i=0;i<nthreads;i++) lambdas[i] = i/(nthreads-1.0)*Lambda;
             
             vmcCSEngine->runCS(currentParameters,currentParameterDirections,lambdas);
             Lambda=lambdas[0];
             if (Lambda==0)
             {
               app_log()<<"  Failed Step. Lambda=0."<<endl;
               tooManyTries--;
               if (tooManyTries>0)
               {
                if(stability==0) stabilityBase+=stabilizerScale;
                else stabilityBase-=stabilizerScale;
                   stability-=1;
//                    app_log()<<" Re-run with larger stabilityBase"<<endl;
                   continue;
               }
             }
             vector<RealType> cs(numParams);
             for (int i=0; i<numParams; i++) cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
//              optTarget->resetPsi(false);
             savedCSparameters.push_back(cs);
        }
        std::pair<RealType,RealType> ms;
        ms.first=stability;
// the log fit seems to work best
        ms.second=std::log(XS);
        mappedStabilizers.push_back(ms);

    }
    
    for (int i=0; i<numParams; i++) optTarget->Params(i) = bestParameters[i];
//     optTarget->resetPsi(false);
    vmcCSEngine->clearComponentMatrices();


    
    finish();
    if (tooManyTries==0) return false;
    return true;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCCSLinearOptimize::put(xmlNodePtr q)
{
    string useGPU("no");
    string vmcMove("pbyp");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(useGPU,"gpu");
    oAttrib.add(vmcMove,"move");
    oAttrib.put(q);

    xmlNodePtr qsave=q;
    xmlNodePtr cur=qsave->children;


    int pid=OHMMS::Controller->rank();
    //no walkers exist, add 10
    if (W.getActiveWalkers() == 0) addWalkers(omp_get_max_threads());

    NumOfVMCWalkers=W.getActiveWalkers();

#if not defined(ENABLE_OPENMP)
    APP_ABORT(" CS LINEAR OPT: DESIGNED FOR OPENMP-MPI HYBRID. ");
#endif

    //create VMC engine
    if (vmcEngine ==0)
    {
        vmcEngine = vmcCSEngine = new VMCLinearOptOMP(W,Psi,H,hamPool,psipool);
        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
    }
    vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
    vmcEngine->process(qsave);

    bool success=true;
    if (optTarget == 0)
    {
        optTarget = new QMCCSLinearOptimizeWFmanagerOMP(W,Psi,H,hamPool);
        optTarget->setStream(&app_log());
        success=optTarget->put(q);
    }
    return success;
}

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCCSLinearOptimize.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

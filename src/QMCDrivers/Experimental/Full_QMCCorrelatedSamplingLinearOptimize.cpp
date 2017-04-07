//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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

  using MatrixOperators::product;


QMCCSLinearOptimize::QMCCSLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCLinearOptimize(w,psi,h,hpool,ppool),
  vmcCSEngine(0), Max_iterations(1), exp0(-16), nstabilizers(10), stabilizerScale(0.5), bigChange(2), w_beta(0.0),
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
  stepsize=0.75;
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
  std::vector<RealType> currentParameterDirections(N,0);
  std::vector<RealType> currentParameters(numParams,0);
  for (int i=0; i<numParams; i++)
    currentParameters[i] = optTarget->Params(i);
  optdir.resize(numParams,0);
  optparm.resize(numParams,0);
//     Matrix<RealType> Ham(N,N);
//     Matrix<RealType> Ham2(N,N);
//     Matrix<RealType> Var(N,N);
//     Matrix<RealType> S(N,N);
//     vmcCSEngine->fillMatrices(Ham2,Ham,Var,S);
  std::vector<RealType> bestParameters(currentParameters);
//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
  RealType stabilityBase(exp0);
  int tooManyTries(20);
  int failedTries(0);
  Matrix<RealType> Left(N,N);
  Matrix<RealType> LeftT(N,N);//temp and transpose
  Matrix<RealType> Right(N,N);
  RealType H2rescale = vmcCSEngine->fillOverlapHamiltonianMatrices(LeftT,Right);
  std::vector<std::pair<RealType,RealType> > mappedStabilizers;
  std::vector<std::vector<RealType> > savedCSparameters;
//     RealType H2rescale(1.0);
//     if (GEVtype=="H2")
//     {
//         Left_tmp=Ham;
//         H2rescale=1.0/Ham2(0,0);
//         Right=(1-w_beta)*S + w_beta*H2rescale*Ham2;
//     }
//     else
//     {
//         Right=S;
//         Left_tmp=(1.0-w_beta)*Ham + w_beta*Var;
//     }
  bool apply_inverse(true);
  if(apply_inverse)
  {
    Matrix<RealType> Right_tmp(Right);
    invert_matrix(Right,false);
    product(Right,LeftT,Left);
    Right=Right_tmp;
  }
  else
    Left=LeftT;
  //Find largest off-diagonal element compared to diagonal element.
  //This gives us an idea how well conditioned it is and can be used to stabilize.
  RealType od_largest(0);
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
//             app_log()<<"od_largest "<<od_largest<< std::endl;
  if (od_largest>0)
    od_largest = std::log(od_largest)/(nstabilizers-1);
  if (od_largest<=0)
    od_largest = stabilizerScale;
  RealType safe = Left(0,0);
  RealType XS(0);
  RealType lastCost(0);
  nstabilizers=omp_get_max_threads();
  for (int stability=0; stability<nstabilizers; stability++)
  {
    app_log()<<"Iteration: "<<stability+1<<"/"<<nstabilizers<< std::endl;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
      {
        LeftT(i,j)= Left(j,i);
      }
    RealType XS(stabilityBase+od_largest*stability);
    int nms(0);
    for (int i=0; i<mappedStabilizers.size(); i++)
      if (mappedStabilizers[i].second==mappedStabilizers[i].second)
        nms++;
    if (nms>=3)
    {
      RealType estval(0);
      bool SuccessfulFit(fitMappedStabilizers(mappedStabilizers,XS,estval,100));
      if (!SuccessfulFit)
      {
        if (stability==0)
          stabilityBase+=stabilizerScale;
        else
        {
          RealType maxXS(mappedStabilizers[0].second);
          RealType minXS(mappedStabilizers[0].second);
          for (int i=1; i<mappedStabilizers.size(); i++)
            if (mappedStabilizers[i].second==mappedStabilizers[i].second)
            {
              maxXS=std::max(maxXS,mappedStabilizers[i].second);
              minXS=std::min(minXS,mappedStabilizers[i].second);
            }
          //resetting XS range
          od_largest=(maxXS-minXS)/(nstabilizers-stability+1);
          nstabilizers=nstabilizers-stability;
          stability=1;
          stabilityBase=minXS;
          app_log()<<" Resetting exp0 range\n  New stability: "<<stability<<"/"<<nstabilizers;
        }
        XS = stabilityBase+od_largest*stability;
      }
    }
    for (int i=1; i<N; i++)
      LeftT(i,i) += std::exp(XS);
    app_log()<<"Trying exp0: "<<XS<< std::endl;
    RealType lowestEV;
    myTimers[2]->start();
    lowestEV =getLowestEigenvector(LeftT,currentParameterDirections);
    myTimers[2]->stop();
    Lambda = H2rescale*getNonLinearRescale(currentParameterDirections,Right);
    RealType bigVec(0);
    for (int i=0; i<numParams; i++)
      bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
    if (std::abs(Lambda*bigVec)>bigChange)
    {
      app_log()<<"  Failed Step. Largest EV parameter change: "<<Lambda*bigVec<< std::endl;
      if (stability==0)
      {
        failedTries++;
        stability--;
        stabilityBase+=stabilizerScale;
      }
      else
        stability=nstabilizers;
      continue;
    }
    RealType newCost(lowestEV);
    if (MinMethod=="rescale")
    {
      std::vector<RealType> cs(numParams);
      for (int i=0; i<numParams; i++)
        cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
//              optTarget->resetPsi(false);
      savedCSparameters.push_back(cs);
      mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(XS,newCost)));
    }
    else
    {
      //some flavor of linear fit to a "reasonable" linemin search
      int nthreads = omp_get_max_threads();
      std::vector<std::vector<RealType> > params_lambdas(nthreads,std::vector<RealType>(numParams,0));
      std::vector<RealType> csts(nthreads);
      std::vector<RealType> cs(numParams);
      RealType error(0);
      RealType fitCost;
      int retry(0);
      std::vector<std::pair<RealType,RealType> > mappedCosts;
      bool goodTry(true);
      RealType rescaledLambda(Lambda);
      while((goodTry)&&(retry<=1))
      {
        if (!goodTry)
          Lambda*=-1;
        for(int j=0; j<nthreads; j++)
          for (int i=0; i<numParams; i++)
            params_lambdas[j][i] = currentParameters[i] + stepsize*Lambda*(j-1.0)*currentParameterDirections[i+1];
        vmcCSEngine->runCS(params_lambdas,error);
        vmcCSEngine->getDeltaCosts(csts);
        mappedCosts.clear();
        for(int i=0; i<nthreads; i++)
          mappedCosts.push_back(*(new std::pair<RealType,RealType>(stepsize*Lambda*(i-1.0),csts[i]-csts[1])));
        int nmcs(0);
        for (int i=0; i<mappedCosts.size(); i++)
          if (mappedCosts[i].second==mappedCosts[i].second)
            nmcs++;
        if (nmcs<3)
        {
          if (stability==0)
          {
            failedTries++;
            stability--;
            stabilityBase+=stabilizerScale;
          }
          else
            stability=nstabilizers;
          continue;
        }
        if (fitMappedStabilizers(mappedCosts,Lambda,fitCost,bigChange/bigVec))
        {
          //             use fit if it works
          for (int i=0; i<numParams; i++)
            bestParameters[i] = cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
          app_log()<<" : Using fit."<< std::endl;
        }
        else
        {
          //             otherwise use the best value
          int indx(0);
          for(int i=1; i<nthreads; i++)
            if (csts[i]<csts[indx])
              indx=i;
          Lambda=stepsize*Lambda*(indx-1.0);
          fitCost=csts[indx]-csts[1];
          bestParameters=cs=params_lambdas[indx];
          app_log()<<" :\n Using best, bestCost: "<<fitCost<< std::endl;
        }
        if (Lambda*rescaledLambda<0)
        {
          goodTry=true;
          retry++;
        }
        else
        {
          goodTry=false;
          newCost=fitCost;
          mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(XS,newCost)));
          savedCSparameters.push_back(cs);
          app_log()<<"COSTS: ";
          for(int i=0; i<nthreads; i++)
            app_log()<<csts[i]-csts[1]<<" ";
          app_log()<< std::endl;
        }
      }
    }
    if(savedCSparameters.size()==omp_get_max_threads())
    {
      app_log()<<"   Finalizing iteration: Choosing best"<< std::endl;
      RealType error(0);
      int bestP = vmcCSEngine->runCS(savedCSparameters,error);
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = bestParameters[i] = savedCSparameters[bestP][i];
      if (bestP<0)
      {
        app_log()<<"   Error in CS cost function. Unchanged parameters."<< std::endl;
        bestP=0;
        vmcCSEngine->clearComponentMatrices();
        for (int i=0; i<numParams; i++)
          optTarget->Params(i) = currentParameters[i];
        finish();
        return false;
      }
    }
//         else if ((stability>0) && (newCost>lastCost))
//           stability=nstabilizers;
//         else
//           lastCost=newCost;
  }
  for (int i=0; i<numParams; i++)
    optTarget->Params(i) = bestParameters[i];
//     optTarget->resetPsi(false);
  vmcCSEngine->clearComponentMatrices();
  finish();
  if (tooManyTries==0)
    return false;
  return true;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCCSLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("no");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU,"gpu");
  oAttrib.add(vmcMove,"move");
  oAttrib.put(q);
  xmlNodePtr qsave=q;
  xmlNodePtr cur=qsave->children;
  int pid=OHMMS::Controller->rank();
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers=W.getActiveWalkers();
#if not defined(ENABLE_OPENMP)
  APP_ABORT(" CS LINEAR OPT: DESIGNED FOR OPENMP-MPI HYBRID. ");
#endif
  //create VMC engine
  if (vmcEngine ==0)
  {
    vmcEngine = vmcCSEngine = new VMCLinearOptOMP(W,Psi,H,hamPool,psiPool);
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

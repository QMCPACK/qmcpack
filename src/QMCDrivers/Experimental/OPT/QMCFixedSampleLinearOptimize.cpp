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
    
    


#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#endif
#include "QMCDrivers/VMC/VMCSingle.h"
#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCApp/HamiltonianPool.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

QMCFixedSampleLinearOptimize::QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCLinearOptimize(w,psi,h,hpool,ppool),
  Max_iterations(1), exp0(-16), exp1(0),  nstabilizers(3), stabilizerScale(2.0), bigChange(1), eigCG(1), TotalCGSteps(1), w_beta(0.0),
  MinMethod("quartic"), GEVtype("mixed"), StabilizerMethod("best"), GEVSplit("no")
{
  //set the optimization flag
  QMCDriverMode.set(QMC_OPTIMIZE,1);
  //read to use vmc output (just in case)
  RootName = "pot";
  QMCType ="QMCFixedSampleLinearOptimize";
  m_param.add(WarmupBlocks,"warmupBlocks","int");
  m_param.add(Max_iterations,"max_its","int");
  m_param.add(nstabilizers,"nstabilizers","int");
  m_param.add(stabilizerScale,"stabilizerscale","double");
  m_param.add(bigChange,"bigchange","double");
  m_param.add(eigCG,"eigcg","int");
  m_param.add(TotalCGSteps,"cgsteps","int");
  m_param.add(w_beta,"beta","double");
  quadstep=-1.0;
  stepsize=0.75;
  m_param.add(quadstep,"quadstep","double");
  m_param.add(stepsize,"stepsize","double");
  m_param.add(exp0,"exp0","double");
  m_param.add(exp1,"exp1","double");
  m_param.add(MinMethod,"MinMethod","string");
  m_param.add(GEVtype,"GEVMethod","string");
  m_param.add(GEVSplit,"GEVSplit","string");
  m_param.add(StabilizerMethod,"StabilizerMethod","string");
  m_param.add(LambdaMax,"LambdaMax","double");
  //Set parameters for line minimization:
  this->add_timers(myTimers);
}

/** Clean up the vector */
QMCFixedSampleLinearOptimize::~QMCFixedSampleLinearOptimize()
{
}

QMCFixedSampleLinearOptimize::RealType QMCFixedSampleLinearOptimize::Func(RealType dl)
{
  for (int i=0; i<optparm.size(); i++)
    optTarget->Params(i) = optparm[i] + dl*optdir[i];
  QMCLinearOptimize::RealType c = optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
//    if (validFuncVal)
  validFuncVal = optTarget->IsValid;
  return c;
}

bool QMCFixedSampleLinearOptimize::run()
{
  start();
  bool Valid(true);
  int Total_iterations(0);
  savedQuadstep=quadstep;
//size of matrix
  numParams = optTarget->NumParams();
  N = numParams + 1;
//  solve CSFs and other parameters separately then rescale elements accordingly
  int first,last;
  getNonLinearRange(first,last);
//  There is only one type of parameter and all are non-linear, don't split it up
  if (last-first==numParams)
    GEVSplit=="no";
//     initialize our parameters
  std::vector<RealType> currentParameterDirections(N,0);
  std::vector<RealType> currentParameters(numParams,0);
  optdir.resize(numParams,0);
  optparm.resize(numParams,0);
  for (int i=0; i<numParams; i++)
    currentParameters[i] = optTarget->Params(i);
  std::vector<RealType> BestDirection(N,0);
  std::vector<RealType> bestParameters(currentParameters);
  std::vector<RealType> GEVSplitParameters(numParams,0);
  while (Total_iterations < Max_iterations)
  {
    Total_iterations+=1;
    app_log()<<"Iteration: "<<Total_iterations<<"/"<<Max_iterations<< std::endl;
// mmorales
    if (!ValidCostFunction(Valid))
      continue;
//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
//         if (runningStabilityBase>stabilityBase) stabilityBase=runningStabilityBase;
//         app_log()<<"  Starting with Stability Base of: "<<stabilityBase<< std::endl;
    //This is the amount we add to the linear parameters
    RealType linearStabilityBase(exp1);
    std::vector<std::vector<RealType> > LastDirections;
    RealType deltaPrms(-1.0);
    for (int tries=0; tries<TotalCGSteps; tries++)
    {
      bool acceptedOneMove(false);
      int tooManyTries(20);
      int failedTries(0);
      std::vector<std::pair<RealType,RealType> > mappedStabilizers;
      if (nstabilizers<2)
      {
        if (StabilizerMethod=="fit")
          app_log()<<" Need 2 stabilizers minimum for the fit"<< std::endl;
        StabilizerMethod="best";
      }
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = currentParameters[i];
      myTimers[4]->start();
      RealType lastCost(optTarget->Cost(true));
      myTimers[4]->start();
      // mmorales
      Valid=optTarget->IsValid;
      if (!ValidCostFunction(Valid))
        continue;
      RealType newCost(lastCost);
      RealType startCost(lastCost);
      Matrix<RealType> LeftT(N,N);
      Matrix<RealType> Left(N,N);
      Matrix<RealType> Right(N,N);
      Matrix<RealType> S(N,N);
      RealType H2rescale=optTarget->fillOverlapHamiltonianMatrices(LeftT,Right,S);
      Left=0.0;
//           optTarget->fillOverlapHamiltonianMatrices(Ham2, Ham, Var, S);
//           if (GEVtype=="H2")
//           {
//
//               Left_tmp=Ham;
//               H2rescale=1.0/Ham2(0,0);
//               Right=(1-w_beta)*S + w_beta*H2rescale*Ham2;
//           }
//           else
//           {
//               Right=S;
//               Left_tmp=(1.0-w_beta)*Ham + w_beta*Var;
//           }
      bool apply_inverse(true);
      if(apply_inverse)
      {
        Matrix<RealType> RightT(Right);
        invert_matrix(RightT,false);
        MatrixOperators MO;
        MO.product(RightT,LeftT,Left);
      }
//           //Find largest off-diagonal element compared to diagonal element.
//           //This gives us an idea how well conditioned it is and can be used to stabilize.
//           RealType od_largest(0);
//           for (int i=0; i<N; i++) for (int j=0; j<N; j++)
//             od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
//           app_log()<<"od_largest "<<od_largest<< std::endl;
//           if((nstabilizers>1)and (od_largest>0)) od_largest = std::log(od_largest)/(nstabilizers-1);
//           if (od_largest<=0) od_largest = stabilizerScale;
      app_log()<<"  stabilityBase "<<stabilityBase<< std::endl;
      app_log()<<"  stabilizerScale "<<stabilizerScale<< std::endl;
      RealType safe = Left(0,0);
      for (int stability=0; stability<nstabilizers; stability++)
      {
        for (int i=0; i<N; i++)
          for (int j=0; j<N; j++)
            LeftT(i,j)= Left(j,i);
//             RealType XS(stabilityBase+od_largest*stability);
//             int nms(0);
//             for (int i=0; i<mappedStabilizers.size(); i++) if (mappedStabilizers[i].second==mappedStabilizers[i].second) nms++;
//            if (nms>=3)
//            {
//              RealType estval(0);
//              bool SuccessfulFit(fitMappedStabilizers(mappedStabilizers,XS,estval,100));
//              if (!SuccessfulFit)
//              {
//                if (stability==0)
//                  stabilityBase+=stabilizerScale;
//                else
//                {
//                  RealType maxXS(mappedStabilizers[0].second);
//                  RealType minXS(mappedStabilizers[0].second);
//                  for (int i=1; i<mappedStabilizers.size(); i++)
//                    if (mappedStabilizers[i].second==mappedStabilizers[i].second)
//                    {
//                      maxXS=std::max(maxXS,mappedStabilizers[i].second);
//                      minXS=std::min(minXS,mappedStabilizers[i].second);
//                    }
//                  //resetting XS range
//                  od_largest=(maxXS-minXS)/(nstabilizers-stability+1);
//                  nstabilizers=nstabilizers-stability;
//                  stability=1;
//                  stabilityBase=minXS;
//                  app_log()<<" Resetting XS range new its:"<<stability<<"/"<<nstabilizers;
//                }
//                XS = stabilityBase+od_largest*stability;
//              }
//            }
//             for (int i=1; i<N; i++) LeftT(i,i) += std::exp(XS);
//             app_log()<<"  Using XS:"<<XS<< std::endl;
        RealType XS(stabilityBase+stabilizerScale*stability);
        if (failedTries>0)
        {
          for (int i=1; i<N; i++)
            LeftT(i,i) += std::exp(XS);
          app_log()<<"  Using XS:"<<XS<< std::endl;
        }
        RealType lowestEV(0);
//             myTimers[2]->start();
// //                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
//             lowestEV = getLowestEigenvector(LeftT,currentParameterDirections);
//             myTimers[2]->stop();
        myTimers[2]->start();
        //                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
        lowestEV = getLowestEigenvector(LeftT,currentParameterDirections);
        Lambda = getNonLinearRescale(currentParameterDirections,Right);
        myTimers[2]->stop();
//             Lambda = H2rescale*getNonLinearRescale(currentParameterDirections,S);
        RealType bigVec(0);
        for (int i=0; i<numParams; i++)
          bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
        if (std::abs(Lambda*bigVec)>bigChange)
        {
          app_log()<<"  Failed Step. Magnitude of largest parameter change: "<<std::abs(Lambda*bigVec)<< std::endl;
          if (stability==0)
          {
            failedTries++;
            stability--;
            stabilityBase+=stabilizerScale;
          }
          else
            stability=nstabilizers;
          continue;
//                 mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
        }
        if (MinMethod=="rescale")
        {
          for (int i=0; i<numParams; i++)
            optTarget->Params(i) = currentParameters[i] + Lambda*currentParameterDirections[i+1];
          optTarget->IsValid = true;
        }
        else
        {
          //eigenCG part
          for (int ldi=0; ldi < std::min(eigCG,int(LastDirections.size())); ldi++)
          {
            RealType nrmold(0), ovlpold(0);
            for (int i=1; i<N; i++)
              nrmold += LastDirections[ldi][i]*LastDirections[ldi][i];
            for (int i=1; i<N; i++)
              ovlpold += LastDirections[ldi][i]*currentParameterDirections[i];
            ovlpold*=1.0/nrmold;
            for (int i=1; i<N; i++)
              currentParameterDirections[i] -= ovlpold * LastDirections[ldi][i];
          }
          //if we chose to "freeze" the CSF solutions at their minimum
          //  then we must add them in to the fixed part of the parameter changes
          for (int i=0; i<numParams; i++)
            optparm[i] = currentParameters[i];
          for (int i=0; i<numParams; i++)
            optdir[i] = currentParameterDirections[i+1];
          RealType bigVec(0);
          for (int i=0; i<numParams; i++)
            bigVec = std::max(bigVec,std::abs(optdir[i]));
          TOL = param_tol/bigVec;
          AbsFuncTol=true;
          largeQuarticStep=bigChange/bigVec;
          quadstep = Lambda;
//                  initial guess for line min bracketing
          LambdaMax = quadstep;
          myTimers[3]->start();
          if (MinMethod=="quartic")
          {
            int npts(7);
            quadstep = stepsize*Lambda;
            LambdaMax = Lambda;
            Valid=lineoptimization3(npts,startCost);
          }
          else
            Valid=lineoptimization2();
          myTimers[3]->stop();
          RealType biggestParameterChange = bigVec*std::abs(Lambda);
          if (biggestParameterChange>bigChange)
          {
            app_log()<<"  Failed Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
            failedTries++;
            stability--;
            stabilityBase+=stabilizerScale;
//                   mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
            mappedStabilizers.push_back(make_pair<RealType,RealType>(XS,std::numeric_limits<RealType>::quiet_NaN()));
            continue;
//                     for (int i=0; i<numParams; i++) optTarget->Params(i) = optparm[i];
          }
          else
          {
            for (int i=0; i<numParams; i++)
              optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
            app_log()<<"  Good Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
          }
          //Save this value in here for later
          Lambda = biggestParameterChange;
        }
        //get cost at new minimum
        newCost = optTarget->Cost(false);
        mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(XS,newCost)));
        app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<" Delta Cost:"<<(newCost-lastCost)<< std::endl;
        optTarget->printEstimates();
//                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
        // mmorales
        Valid=optTarget->IsValid;
        if (!ValidCostFunction(Valid))
        {
          app_log()<<"  Good Step, but cost function invalid"<< std::endl;
          failedTries++;
          stability--;
          if(stability<0)
            stabilityBase+=stabilizerScale;
          else
            stability=nstabilizers;
//                   mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
          mappedStabilizers.push_back(make_pair<RealType,RealType>(XS,std::numeric_limits<RealType>::quiet_NaN()));
          continue;
        }
        if (newCost < lastCost)
        {
          //Move was acceptable
          for (int i=0; i<numParams; i++)
            bestParameters[i] = optTarget->Params(i);
          lastCost=newCost;
          BestDirection=currentParameterDirections;
          acceptedOneMove=true;
          deltaPrms= Lambda;
        }
        else
          if ((stability>0) && (newCost>(lastCost-1e-5)))
            stability=nstabilizers;
//             else if ((newCost>lastCost)&&(StabilizerMethod=="fit")&&(mappedStabilizers.size()>2))
//             {
//               stability = std::max(nstabilizers-2,stability);
//               app_log()<<"Small change, moving on to fit."<< std::endl;
//             }
      }
      if (acceptedOneMove)
      {
        for (int i=0; i<numParams; i++)
          optTarget->Params(i) = bestParameters[i];
        currentParameters=bestParameters;
        LastDirections.push_back(BestDirection);
        acceptedOneMove=false;
//               runningStabilityBase = stabilityBase;
//             app_log()<< " Wave Function Parameters updated."<< std::endl;
//             optTarget->reportParameters();
      }
      else
      {
        for (int i=0; i<numParams; i++)
          optTarget->Params(i) = currentParameters[i];
        tries=TotalCGSteps;
      }
    }
  }
  finish();
  return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCFixedSampleLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("yes");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU,"gpu");
  oAttrib.add(vmcMove,"move");
  oAttrib.put(q);
  xmlNodePtr qsave=q;
  xmlNodePtr cur=qsave->children;
  int pid=OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
//         else if (cname == "optimizer")
//         {
//             xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
//             if (att)
//             {
//                 optmethod = (const char*)att;
//             }
//             optNode=cur;
//         }
//         else if (cname == "optimize")
//         {
//             xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
//             if (att)
//             {
//                 optmethod = (const char*)att;
//             }
//         }
    cur=cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers=W.getActiveWalkers();
  //create VMC engine
  if (vmcEngine ==0)
  {
#if defined (QMC_CUDA)
    if (useGPU == "yes")
      vmcEngine = new VMCcuda(W,Psi,H,psiPool);
    else
#endif
      vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool,psiPool);
    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
    vmcEngine->initCommunicator(myComm);
  }
  vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
  vmcEngine->process(qsave);
  bool success=true;
  if (optTarget == 0)
  {
#if defined (QMC_CUDA)
    if (useGPU == "yes")
      optTarget = new QMCCostFunctionCUDA(W,Psi,H,hamPool);
    else
#endif
// #if defined(ENABLE_OPENMP)
//             if (omp_get_max_threads()>1)
//             {
//                 optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
//             }
//             else
// #endif
//                 optTarget = new QMCCostFunctionSingle(W,Psi,H);
      optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
    optTarget->setStream(&app_log());
    success=optTarget->put(q);
  }
  return success;
}

}

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCCorrelatedSamplingLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/QMCCostFunctionBase.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#endif
//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
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

  using MatrixOperators::product;


template<Batching batching>
QMCCorrelatedSamplingLinearOptimize<batching>::QMCCorrelatedSamplingLinearOptimize(MCWalkerConfiguration& w,
										   TrialWaveFunction<batching>& psi,
										   QMCHamiltonian& h,
										   HamiltonianPool<batching>& hpool, WaveFunctionPool& ppool) : QMCLinearOptimize<batching>(w,psi,h,hpool,ppool),
  exp0(-16), nstabilizers(3), stabilizerScale(2.0), bigChange(3), w_beta(0.0), MinMethod("quartic"), GEVtype("mixed")
{
  QDT::IsQMCDriver=false;
  //set the optimization flag
  QDT::QMCDriverMode.set(QDT::QMC_OPTIMIZE,1);
  //read to use vmc output (just in case)
  QDT::RootName = "pot";
  QDT::QMCType ="QMCCorrelatedSamplingLinearOptimize";
  QDT::m_param.add(QLOT::WarmupBlocks,"warmupBlocks","int");
  QDT::m_param.add(stabilizerScale,"stabilizerscale","double");
  QDT::m_param.add(bigChange,"bigchange","double");
  QDT::m_param.add(w_beta,"beta","double");
  QDT::m_param.add(GEVtype,"GEVMethod","string");
  quadstep=-1.0;
  stepsize=0.3;
  QDT::m_param.add(quadstep,"quadstep","double");
  QDT::m_param.add(stepsize,"stepsize","double");
  QDT::m_param.add(exp0,"exp0","double");
  QDT::m_param.add(MinMethod,"MinMethod","string");
  QDT::m_param.add(LambdaMax,"LambdaMax","double");
  QDT::m_param.add(nstabilizers,"nstabilizers","int");
  //Set parameters for line minimization:
}

/** Clean up the vector */
template<Batching batching>
QMCCorrelatedSamplingLinearOptimize<batching>::~QMCCorrelatedSamplingLinearOptimize()
{
}

template<Batching batching>
typename QMCCorrelatedSamplingLinearOptimize<batching>::RealType
QMCCorrelatedSamplingLinearOptimize<batching>::Func(RealType dl)
{
  for (int i=0; i< QLOT::optparm.size(); i++)
    QLOT::optTarget->Params(i) = QLOT::optparm[i] + dl*QLOT::optdir[i];
  QMCT::RealType c = QLOT::optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
//     if (validFuncVal)
  validFuncVal = QLOT::optTarget->IsValid;
  return c;
}

template<Batching batching>
bool QMCCorrelatedSamplingLinearOptimize<batching>::run()
{
  QLOT::start();
//size of matrix
  QLOT::numParams = QLOT::optTarget->NumParams();
  QLOT::N = QLOT::numParams + 1;
//  solve CSFs and other parameters separately then rescale elements accordingly
  int first,last;
  QLOT::getNonLinearRange(first,last);
//     initialize our parameters
  std::vector<RealType> currentParameterDirections(QLOT::N,0);
  std::vector<RealType> currentParameters(QLOT::numParams,0);
  std::vector<RealType> bestParameters(currentParameters);
  QLOT::optdir.resize(QLOT::numParams,0);
  QLOT::optparm.resize(QLOT::numParams,0);
  for (int i=0; i<QLOT::numParams; i++)
    currentParameters[i] = QLOT::optTarget->Params(i);
  //this is the small amount added to the diagonal to stabilize the eigenvalue equation. e^stabilityBase
  RealType stabilityBase(exp0);
  std::vector<std::vector<RealType> > LastDirections;
  RealType deltaPrms(-1.0);
  bool acceptedOneMove(false);
  int tooManyTries(20);
  int failedTries(0);
  Matrix<RealType> Left(QLOT::N,QLOT::N);
  Matrix<RealType> LeftT(QLOT::N,QLOT::N);
  Matrix<RealType> Right(QLOT::N,QLOT::N);
  Right=0;
  LeftT=0;
  Left=0;
  vmcCSEngine->fillOverlapHamiltonianMatrices(LeftT,Right);
  std::vector<std::pair<RealType,RealType> > mappedStabilizers;
  RealType lastCost(0);
  RealType startCost(0);
//    if ((GEVtype!="H2")||(MinMethod!="rescale"))
//    {
  QLOT::myTimers[4]->start();
  startCost = lastCost = QLOT::optTarget->Cost(false);
  QLOT::myTimers[4]->stop();
//    }
  bool apply_inverse(true);
  if(apply_inverse)
  {
    Matrix<RealType> Right_tmp(Right);
    invert_matrix(Right_tmp,false);
    product(Right_tmp,LeftT,Left);
  }
  else
    Left=LeftT;
  //Find largest off-diagonal element compared to diagonal element.
  //This gives us an idea how well conditioned it is and can be used to stabilize.
//     RealType od_largest(0);
//     for (int i=0; i<N; i++) for (int j=0; j<N; j++)
//       od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
//    RealType d_neg(0);
//    for (int i=1; i<N; i++) if (Left(i,i)<d_neg) d_neg=Left(i,i);
//    stabilizerScale = std::max(stabilizerScale*(nstabilizers-1.0),std::log(-d_neg));
//    if(nstabilizers>1)
//      stabilizerScale = stabilizerScale/(nstabilizers-1.0);
  app_log()<<"  stabilityBase "<<stabilityBase<< std::endl;
  app_log()<<"  stabilizerScale "<<stabilizerScale<< std::endl;
  RealType safe = Left(0,0);
  for (int stability=0; stability<nstabilizers; stability++)
  {
    for (int i=0; i<QLOT::N; i++)
      for (int j=0; j<QLOT::N; j++)
      {
        LeftT(i,j) = Left(j,i);
      }
    RealType XS(stabilityBase+stabilizerScale*stability);
    if (failedTries>0)
    {
      for (int i=1; i<QLOT::N; i++)
        LeftT(i,i) += std::exp(XS);
      app_log()<<"  Using XS:"<<XS<< std::endl;
    }
    RealType lowestEV(0);
    RealType bigVec(0);
    QLOT::myTimers[2]->start();
//                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
    if (GEVtype!="sd")
    {
      lowestEV = QLOT::getLowestEigenvector(LeftT,currentParameterDirections);
      Lambda = QLOT::getNonLinearRescale(currentParameterDirections,Right);
    }
    else
    {
      currentParameterDirections[0]=1.0;
      for (int i=0; i<QLOT::numParams; i++)
        bigVec = std::max(bigVec,std::abs(LeftT(i+1,0)));
      //bigVec /= LeftT(0,0);
      RealType rscale(stepsize/bigVec);
      for (int i=0; i<QLOT::numParams; i++)
        currentParameterDirections[i+1]=LeftT(i+1,0)*rscale;
      Lambda=1.0;
    }
    QLOT::myTimers[2]->stop();
    for (int i=0; i<QLOT::numParams; i++)
      bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
    if (std::abs(Lambda*bigVec)>bigChange)
    {
      app_log()<<"  Failed Step. Largest EV parameter change: "<<Lambda*bigVec<< std::endl;
//           if (GEVtype=="H2") continue;
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
      for (int i=0; i<QLOT::numParams; i++)
        bestParameters[i] = QLOT::optTarget->Params(i) = currentParameters[i] + Lambda*currentParameterDirections[i+1];
      if (GEVtype=="H2")
        acceptedOneMove = true;
    }
    else
    {
      for (int i=0; i<QLOT::numParams; i++)
        QLOT::optparm[i] = currentParameters[i];
      for (int i=0; i<QLOT::numParams; i++)
        QLOT::optdir[i] = currentParameterDirections[i+1];
      RealType bigOptVec(0);
      for (int i=0; i<QLOT::numParams; i++)
        bigOptVec = std::max(bigOptVec,std::abs(QLOT::optdir[i]));
      TOL = QLOT::param_tol/bigOptVec;
      AbsFuncTol=true;
      largeQuarticStep=bigChange/bigVec;
      quadstep = stepsize/bigVec;
//                  initial guess for line min bracketing
      LambdaMax = quadstep;
      QLOT::myTimers[3]->start();
      if (MinMethod=="quartic")
      {
        int npts(7);
        quadstep = stepsize*Lambda;
        LambdaMax = Lambda;
        lineoptimization3(npts,startCost);
      }
      else
        lineoptimization2();
      QLOT::myTimers[3]->stop();
      RealType biggestParameterChange = bigOptVec*std::abs(Lambda);
      if (biggestParameterChange>bigChange)
      {
        app_log()<<"  Failed Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
        failedTries++;
        stability--;
        stabilityBase+=stabilizerScale;
//                   mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
        //mappedStabilizers.push_back(make_pair<RealType,RealType>(XS,std::numeric_limits<RealType>::quiet_NaN()));
        mappedStabilizers.push_back(std::pair<RealType,RealType>(XS,std::numeric_limits<RealType>::quiet_NaN()));
        continue;
//                     for (int i=0; i<QLOT::numParams; i++) optTarget->Params(i) = optparm[i];
      }
      else
      {
        for (int i=0; i<QLOT::numParams; i++)
          QLOT::optTarget->Params(i) = QLOT::optparm[i] + Lambda * QLOT::optdir[i];
        app_log()<<"  Largest LM parameter change:"<<biggestParameterChange<< std::endl;
      }
      //Save this value in here for later
      Lambda = biggestParameterChange;
    }
//      if ((GEVtype!="H2")||(MinMethod!="rescale"))
//      {
    //get cost at new minimum
    RealType newCost = QLOT::optTarget->Cost(false);
    mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(XS,newCost)));
    app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<" Delta Cost:"<<(newCost-lastCost)<< std::endl;
    QLOT::optTarget->printEstimates();
    if (newCost < lastCost)
    {
      //Move was acceptable
      for (int i=0; i<QLOT::numParams; i++)
        bestParameters[i] = QLOT::optTarget->Params(i);
      lastCost=newCost;
      acceptedOneMove=true;
      deltaPrms= Lambda;
    }
    else
      if ((stability>0) && (newCost-lastCost>-1e-5))
        stability=nstabilizers;
//    }
  }
  if (acceptedOneMove)
    for (int i=0; i<QLOT::numParams; i++)
      QLOT::optTarget->Params(i) = bestParameters[i];
  else
    for (int i=0; i<QLOT::numParams; i++)
      QLOT::optTarget->Params(i) = currentParameters[i];
  QLOT::finish();
  return (QLOT::optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/

template<Batching batching>
bool QMCCorrelatedSamplingLinearOptimize<batching>::put(xmlNodePtr q)
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
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      QDT::mcwalkerNodePtr.push_back(cur);
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
  if (QDT::W.getActiveWalkers() == 0)
    QDT::addWalkers(omp_get_max_threads());
  QLOT::NumOfVMCWalkers=QDT::W.getActiveWalkers();
  //create VMC engine
  if (QLOT::vmcEngine ==0)
  {
#if defined (QMC_CUDA)
    vmcCSEngine = new VMCcuda(QDT::W,dynamic_cast<TrialWaveFunction<Batching::BATCHED>&>(QDT::Psi),QDT::H,QDT::psiPool);
    vmcCSEngine->setOpt(true);
    QLOT::vmcEngine = vmcCSEngine;
#else
    QLOT::vmcEngine = vmcCSEngine = new VMCLinearOptOMP(QDT::W,QDT::Psi,QDT::H,QLOT::hamPool,QDT::psiPool);
#endif
    QLOT::vmcEngine->setUpdateMode(vmcMove[0] == 'p');
    QLOT::vmcEngine->initCommunicator(QDT::myComm);
  }
  QLOT::vmcEngine->setStatus(QDT::RootName,QDT::h5FileRoot,QDT::AppendRun);
  QLOT::vmcEngine->process(qsave);
  bool success=true;
  if (QLOT::optTarget == 0)
  {
#if defined (QMC_CUDA)
    QLOT::optTarget = new QMCCostFunctionCUDA(QDT::W,QDT::Psi,QDT::H);
#else
    QLOT::optTarget = new QMCCostFunctionOMP(QDT::W,QDT::Psi,QDT::H);
#endif
    QLOT::optTarget->setneedGrads(false);
    QLOT::optTarget->setStream(&app_log());
    success=QLOT::optTarget->put(qsave);
  }
  return success;
}

template class QMCCorrelatedSamplingLinearOptimize<Batching::SINGLE>;
#ifdef QMC_CUDA
template class QMCCorrelatedSamplingLinearOptimize<Batching::BATCHED>;
#endif
  
}

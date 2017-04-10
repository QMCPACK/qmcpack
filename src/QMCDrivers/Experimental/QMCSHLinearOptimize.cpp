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
    
    



#include "QMCDrivers/QMCSHLinearOptimize.h"
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

QMCSHLinearOptimize::QMCSHLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCLinearOptimize(w,psi,h,hpool,ppool),
  bigChange(1), w_beta(0.0), MinMethod("quartic")
{
//     //set the optimization flag
  QMCDriverMode.set(QMC_OPTIMIZE,1);
  //read to use vmc output (just in case)
  RootName = "opt";
  QMCType ="QMCSHLinOpt";
  m_param.add(bigChange,"bigchange","double");
  m_param.add(MinMethod,"MinMethod","double");
  m_param.add(w_beta,"beta","double");
  quadstep=-1.0;
  stepsize=0.75;
  m_param.add(quadstep,"quadstep","double");
  m_param.add(stepsize,"stepsize","double");
  this->add_timers(myTimers);
}

/** Clean up the vector */
QMCSHLinearOptimize::~QMCSHLinearOptimize()
{
}

QMCSHLinearOptimize::RealType QMCSHLinearOptimize::Func(RealType dl)
{
  for (int i=0; i<optparm.size(); i++)
    optTarget->Params(i) = optparm[i] + dl*optdir[i];
  QMCLinearOptimize::RealType c = optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
//     if (validFuncVal)
  validFuncVal = optTarget->IsValid;
  return c;
}

bool QMCSHLinearOptimize::run()
{
  start();
//size of matrix
  numParams = optTarget->NumParams();
  N = numParams + 1;
//     initialize our parameters
  std::vector<RealType> currentOvlp(N,0);
  std::vector<RealType> currentHOvlp(N,0);
  std::vector<RealType> currentParameters(numParams,0);
  RealType E_avg(0);
  std::vector<RealType> bestParameters(currentParameters);
  optdir.resize(numParams,0);
  optparm.resize(numParams,0);
  for (int i=0; i<numParams; i++)
    optparm[i] = currentParameters[i] = optTarget->Params(i);
  bool acceptedOneMove(false);
  int tooManyTries(20);
  int failedTries(0);
  RealType lastCost(0);
  RealType startCost(0);
  startCost = lastCost = optTarget->Cost(false);
  app_log()<<"Starting cost: "<<startCost<< std::endl;
  Matrix<RealType> OM;
  OM.resize(N,N);
  dmcEngine->fillVectors(currentOvlp,currentHOvlp,E_avg,OM);
  for (int i=0; i<numParams; i++)
    optdir[i]=currentOvlp[i];
  std::vector<RealType> dP(N,1);
  for (int i=0; i<numParams; i++)
    dP[i+1]=optdir[i];
  Lambda = getNonLinearRescale(dP,OM);
  app_log()<<"rescaling factor :"<<Lambda<< std::endl;
  RealType bigOptVec(std::abs(optdir[0]));
  for (int i=1; i<numParams; i++)
    bigOptVec =std::max(std::abs(optdir[i]),bigOptVec);
//     app_log()<<"currentOvlp"<< std::endl;
//     for (int i=0; i<numParams; i++)
//     app_log()<<optdir[i]<<" ";
//     app_log()<< std::endl;
//
//     app_log()<<"optparam"<< std::endl;
//     for (int i=0; i<numParams; i++)
//     app_log()<<optparm[i]<<" ";
//     app_log()<< std::endl;
  if (MinMethod=="rescale")
  {
    if (bigOptVec*std::abs(Lambda)>bigChange)
      app_log()<<"  Failed Step. Largest LM parameter change:"<<bigOptVec*std::abs(Lambda)<< std::endl;
    else
    {
      for (int i=0; i<numParams; i++)
        bestParameters[i] = optTarget->Params(i) = currentParameters[i] + Lambda*optdir[i];
      acceptedOneMove = true;
    }
  }
  else
    if (MinMethod=="overlap")
    {
      orthoScale(optdir,OM);
      if (bigOptVec>bigChange)
        app_log()<<"  Failed Step. Largest LM parameter change:"<<bigOptVec<< std::endl;
      else
      {
        for (int i=0; i<numParams; i++)
          bestParameters[i] = optTarget->Params(i) = currentParameters[i] + optdir[i];
        acceptedOneMove = true;
      }
    }
    else
      if (MinMethod=="average")
      {
        for (int i=0; i<numParams; i++)
          bestParameters[i] = optTarget->Params(i) = currentParameters[i] + currentHOvlp[i];
        acceptedOneMove = true;
      }
      else
      {
        TOL = param_tol/bigOptVec;
        AbsFuncTol=true;
        largeQuarticStep=bigChange/bigOptVec;
        quadstep = stepsize*Lambda; //bigOptVec;
//                  initial guess for line min bracketing
        LambdaMax = quadstep;
        myTimers[3]->start();
        if (MinMethod=="quartic")
        {
          int npts(7);
          lineoptimization3(npts,startCost);
        }
        else
          lineoptimization2();
        myTimers[3]->stop();
        RealType biggestParameterChange = bigOptVec*std::abs(Lambda);
        if (biggestParameterChange>bigChange)
        {
          app_log()<<"  Failed Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
          for (int i=0; i<numParams; i++)
            optTarget->Params(i) = bestParameters[i] = currentParameters[i];
        }
        else
        {
          for (int i=0; i<numParams; i++)
            optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
          lastCost = optTarget->Cost(false);
          app_log()<<" Costs: "<<startCost<<" "<<lastCost<< std::endl;
          app_log()<<" Optimal rescaling factor :"<<Lambda<< std::endl;
          if (lastCost<startCost)
            acceptedOneMove = true;
          else
          {
            for (int i=0; i<numParams; i++)
              optTarget->Params(i) = currentParameters[i];
            app_log()<<"  Failed Step. Cost increase "<< std::endl;
          }
        }
      }
//     if (acceptedOneMove)
//         for (int i=0; i<numParams; i++) optTarget->Params(i) = bestParameters[i];
//     else
//         for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i];
//     if (W.getActiveWalkers()>NumOfVMCWalkers)
//     {
//       W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
//       app_log() << "  QMCLinearOptimize::generateSamples removed walkers." << std::endl;
//       app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << std::endl;
//     }
  finish();
  return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCSHLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("no");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
//     oAttrib.add(useGPU,"gpu");
  oAttrib.add(vmcMove,"move");
  oAttrib.put(q);
  optNode=q;
  xmlNodePtr cur=optNode->children;
  int pid=OHMMS::Controller->rank();
//     while (cur != NULL)
//     {
//         std::string cname((const char*)(cur->name));
//         if (cname == "mcwalkerset")
//         {
//             mcwalkerNodePtr.push_back(cur);
//         }
//         cur=cur->next;
//     }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers=W.getActiveWalkers();
  //create VMC engine
  if ((vmcEngine ==0)||(dmcEngine ==0))
  {
#if defined (QMC_CUDA)
    if (useGPU == "yes")
      vmcEngine = new VMCcuda(W,Psi,H,psiPool);
    else
#endif
      vmcEngine = dmcEngine = new DMCOMPOPT(W,Psi,H,hamPool,psiPool);
    dmcEngine->setUpdateMode(vmcMove[0] == 'p');
    dmcEngine->initCommunicator(myComm);
  }
//     app_log()<<RootName<<"   "<<h5FileRoot<< std::endl;
  dmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
  dmcEngine->process(optNode);
//     dmcEngine->setBranchEngine(branchEngine);
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
    optTarget->setneedGrads(false);
    optTarget->setStream(&app_log());
    optTarget->setDMC();
    success=optTarget->put(optNode);
  }
  return success;
}

}

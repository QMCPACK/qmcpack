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


#include <cassert>
#include <iostream>
#include <fstream>
#include "QMCCorrelatedSamplingLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"
#endif
//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionCUDA.h"
#endif

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{
using MatrixOperators::product;


QMCCorrelatedSamplingLinearOptimize::QMCCorrelatedSamplingLinearOptimize(MCWalkerConfiguration& w,
                                                                         TrialWaveFunction& psi,
                                                                         QMCHamiltonian& h,
                                                                         Communicate* comm)
    : QMCLinearOptimize(w, psi, h, comm),
      nstabilizers(3),
      stabilizerScale(2.0),
      bigChange(3),
      exp0(-16),
      MinMethod("quartic"),
      GEVtype("mixed"),
      w_beta(0.0)
{
  IsQMCDriver = false;
  //set the optimization flag
  qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  RootName = "pot";
  m_param.add(stabilizerScale, "stabilizerscale");
  m_param.add(bigChange, "bigchange");
  m_param.add(w_beta, "beta");
  m_param.add(GEVtype, "GEVMethod");
  quadstep = -1.0;
  stepsize = 0.3;
  m_param.add(quadstep, "quadstep");
  m_param.add(stepsize, "stepsize");
  m_param.add(exp0, "exp0");
  m_param.add(MinMethod, "MinMethod");
  m_param.add(LambdaMax, "LambdaMax");
  m_param.add(nstabilizers, "nstabilizers");
  //Set parameters for line minimization:
}

/** Clean up the vector */
QMCCorrelatedSamplingLinearOptimize::~QMCCorrelatedSamplingLinearOptimize() {}

QMCCorrelatedSamplingLinearOptimize::RealType QMCCorrelatedSamplingLinearOptimize::Func(RealType dl)
{
  for (int i = 0; i < optparm.size(); i++)
    optTarget->Params(i) = optparm[i] + dl * optdir[i];
  QMCLinearOptimize::RealType c = optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
  //     if (validFuncVal)
  validFuncVal = optTarget->IsValid;
  return c;
}

bool QMCCorrelatedSamplingLinearOptimize::run()
{
  start();
  //size of matrix
  numParams = optTarget->getNumParams();
  N         = numParams + 1;
  //  solve CSFs and other parameters separately then rescale elements accordingly
  int first, last;
  getNonLinearRange(first, last);
  //     initialize our parameters
  std::vector<RealType> currentParameterDirections(N, 0);
  std::vector<RealType> currentParameters(numParams, 0);
  std::vector<RealType> bestParameters(currentParameters);
  optdir.resize(numParams, 0);
  optparm.resize(numParams, 0);
  for (int i = 0; i < numParams; i++)
    currentParameters[i] = std::real(optTarget->Params(i));
  //this is the small amount added to the diagonal to stabilize the eigenvalue equation. e^stabilityBase
  RealType stabilityBase(exp0);
  std::vector<std::vector<RealType>> LastDirections;
  RealType deltaPrms(-1.0);
  bool acceptedOneMove(false);
  int failedTries(0);
  Matrix<RealType> Left(N, N);
  Matrix<RealType> LeftT(N, N);
  Matrix<RealType> Right(N, N);
  Right = 0;
  LeftT = 0;
  Left  = 0;
  vmcCSEngine->fillOverlapHamiltonianMatrices(LeftT, Right);
  std::vector<std::pair<RealType, RealType>> mappedStabilizers;
  RealType lastCost(0);
  RealType startCost(0);
  //    if ((GEVtype!="H2")||(MinMethod!="rescale"))
  //    {
  cost_function_timer_.start();
  startCost = lastCost = optTarget->Cost(false);
  cost_function_timer_.stop();
  //    }
  bool apply_inverse(true);
  if (apply_inverse)
  {
    Matrix<RealType> Right_tmp(Right);
    invert_matrix(Right_tmp, false);
    product(Right_tmp, LeftT, Left);
  }
  else
    Left = LeftT;
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
  app_log() << "  stabilityBase " << stabilityBase << std::endl;
  app_log() << "  stabilizerScale " << stabilizerScale << std::endl;
  RealType safe = Left(0, 0);
  for (int stability = 0; stability < nstabilizers; stability++)
  {
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
      {
        LeftT(i, j) = Left(j, i);
      }
    RealType XS(stabilityBase + stabilizerScale * stability);
    if (failedTries > 0)
    {
      for (int i = 1; i < N; i++)
        LeftT(i, i) += std::exp(XS);
      app_log() << "  Using XS:" << XS << std::endl;
    }
    RealType lowestEV(0);
    RealType bigVec(0);
    eigenvalue_timer_.start();
    //                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
    if (GEVtype != "sd")
    {
      lowestEV = getLowestEigenvector(LeftT, currentParameterDirections);
      Lambda   = getNonLinearRescale(currentParameterDirections, Right);
    }
    else
    {
      currentParameterDirections[0] = 1.0;
      for (int i = 0; i < numParams; i++)
        bigVec = std::max(bigVec, std::abs(LeftT(i + 1, 0)));
      //bigVec /= LeftT(0,0);
      RealType rscale(stepsize / bigVec);
      for (int i = 0; i < numParams; i++)
        currentParameterDirections[i + 1] = LeftT(i + 1, 0) * rscale;
      Lambda = 1.0;
    }
    eigenvalue_timer_.stop();
    for (int i = 0; i < numParams; i++)
      bigVec = std::max(bigVec, std::abs(currentParameterDirections[i + 1]));
    if (std::abs(Lambda * bigVec) > bigChange)
    {
      app_log() << "  Failed Step. Largest EV parameter change: " << Lambda * bigVec << std::endl;
      //           if (GEVtype=="H2") continue;
      if (stability == 0)
      {
        failedTries++;
        stability--;
        stabilityBase += stabilizerScale;
      }
      else
        stability = nstabilizers;
      continue;
      //                 mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
    }
    if (MinMethod == "rescale")
    {
      for (int i = 0; i < numParams; i++)
      {
        //FIXME This std::real should be removed later when the optimizer starts to work with complex parameters
        optTarget->Params(i) = currentParameters[i] + Lambda * currentParameterDirections[i + 1];
        bestParameters[i]    = std::real(optTarget->Params(i));
      }
      if (GEVtype == "H2")
        acceptedOneMove = true;
    }
    else
    {
      for (int i = 0; i < numParams; i++)
        optparm[i] = currentParameters[i];
      for (int i = 0; i < numParams; i++)
        optdir[i] = currentParameterDirections[i + 1];
      RealType bigOptVec(0);
      for (int i = 0; i < numParams; i++)
        bigOptVec = std::max(bigOptVec, std::abs(optdir[i]));
      TOL              = param_tol / bigOptVec;
      AbsFuncTol       = true;
      largeQuarticStep = bigChange / bigVec;
      quadstep         = stepsize / bigVec;
      //                  initial guess for line min bracketing
      LambdaMax = quadstep;
      line_min_timer_.start();
      if (MinMethod == "quartic")
      {
        int npts(7);
        quadstep  = stepsize * Lambda;
        LambdaMax = Lambda;
        lineoptimization3(npts, startCost);
      }
      else
        lineoptimization2();
      line_min_timer_.stop();
      RealType biggestParameterChange = bigOptVec * std::abs(Lambda);
      if (biggestParameterChange > bigChange)
      {
        app_log() << "  Failed Step. Largest LM parameter change:" << biggestParameterChange << std::endl;
        failedTries++;
        stability--;
        stabilityBase += stabilizerScale;
        //                   mappedStabilizers.push_back(*(new std::pair<RealType,RealType>(std::numeric_limits<RealType>::quiet_NaN(),XS)));
        //mappedStabilizers.push_back(make_pair<RealType,RealType>(XS,std::numeric_limits<RealType>::quiet_NaN()));
        mappedStabilizers.push_back(std::pair<RealType, RealType>(XS, std::numeric_limits<RealType>::quiet_NaN()));
        continue;
        //                     for (int i=0; i<numParams; i++) optTarget->Params(i) = optparm[i];
      }
      else
      {
        for (int i = 0; i < numParams; i++)
          optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
        app_log() << "  Largest LM parameter change:" << biggestParameterChange << std::endl;
      }
      //Save this value in here for later
      Lambda = biggestParameterChange;
    }
    //      if ((GEVtype!="H2")||(MinMethod!="rescale"))
    //      {
    //get cost at new minimum
    RealType newCost = optTarget->Cost(false);
    mappedStabilizers.push_back(*(new std::pair<RealType, RealType>(XS, newCost)));
    app_log() << " OldCost: " << lastCost << " NewCost: " << newCost << " Delta Cost:" << (newCost - lastCost)
              << std::endl;
    optTarget->printEstimates();
    if (newCost < lastCost)
    {
      //Move was acceptable
      for (int i = 0; i < numParams; i++)
        bestParameters[i] = std::real(optTarget->Params(i));
      lastCost        = newCost;
      acceptedOneMove = true;
      deltaPrms       = Lambda;
    }
    else if ((stability > 0) && (newCost - lastCost > -1e-5))
      stability = nstabilizers;
    //    }
  }
  if (acceptedOneMove)
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = bestParameters[i];
  else
    for (int i = 0; i < numParams; i++)
      optTarget->Params(i) = currentParameters[i];
  finish();
  return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool QMCCorrelatedSamplingLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("no");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU, "gpu");
  oAttrib.add(vmcMove, "move");
  oAttrib.put(q);
  xmlNodePtr qsave = q;
  xmlNodePtr cur   = qsave->children;
  int pid          = OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    cur = cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers = W.getActiveWalkers();
  //create VMC engine
  if (vmcEngine == 0)
  {
#if defined(QMC_CUDA)
    vmcEngine   = std::make_unique<VMCcuda>(W, Psi, H, myComm, false);
    vmcCSEngine = dynamic_cast<VMCcuda*>(vmcEngine.get());
    vmcCSEngine->setOpt(true);
#else
    vmcEngine   = std::make_unique<VMCLinearOpt>(W, Psi, H, myComm);
    vmcCSEngine = dynamic_cast<VMCLinearOpt*>(vmcEngine.get());
#endif
    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
  }
  vmcEngine->setStatus(RootName, h5FileRoot, AppendRun);
  vmcEngine->process(qsave);
  bool success = true;
  //allways reset optTarget
#if defined(QMC_CUDA)
  if (useGPU == "yes")
    optTarget = std::make_unique<QMCCostFunctionCUDA>(W, Psi, H, myComm);
  else
#endif
    optTarget = std::make_unique<QMCCostFunction>(W, Psi, H, myComm);
  optTarget->setneedGrads(false);
  optTarget->setStream(&app_log());
  success = optTarget->put(q);

  return success;
}

} // namespace qmcplusplus

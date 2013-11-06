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
#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
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


QMCFixedSampleLinearOptimize::QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCLinearOptimize(w,psi,h,hpool,ppool), Max_iterations(1), exp0(-16), nstabilizers(3),
  stabilizerScale(2.0), bigChange(50), w_beta(0.0),  MinMethod("quartic"), GEVtype("mixed"),
  StabilizerMethod("best"), GEVSplit("no")
{
  IsQMCDriver=false;
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
  m_param.add(MinMethod,"MinMethod","string");
  m_param.add(exp0,"exp0","double");
  stepsize=0.25;
//   stale parameters
//   m_param.add(eigCG,"eigcg","int");
//   m_param.add(TotalCGSteps,"cgsteps","int");
//   m_param.add(w_beta,"beta","double");
//   quadstep=-1.0;
//   m_param.add(quadstep,"quadstep","double");
//   m_param.add(stepsize,"stepsize","double");
//   m_param.add(exp1,"exp1","double");
//   m_param.add(GEVtype,"GEVMethod","string");
//   m_param.add(GEVSplit,"GEVSplit","string");
//   m_param.add(StabilizerMethod,"StabilizerMethod","string");
//   m_param.add(LambdaMax,"LambdaMax","double");
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
//size of matrix
  numParams = optTarget->NumParams();
  N = numParams + 1;
//   where we are and where we are pointing
  vector<RealType> currentParameterDirections(N,0);
  vector<RealType> currentParameters(numParams,0);
  vector<RealType> bestParameters(numParams,0);
  for (int i=0; i<numParams; i++)
    bestParameters[i] = currentParameters[i] = optTarget->Params(i);
//   proposed direction and new parameters
  optdir.resize(numParams,0);
  optparm.resize(numParams,0);

  while (Total_iterations < Max_iterations)
  {
    Total_iterations+=1;
    app_log()<<"Iteration: "<<Total_iterations<<"/"<<Max_iterations<<endl;
    if (!ValidCostFunction(Valid))
      continue;
//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
//     reset params if necessary
    for (int i=0; i<numParams; i++)
      optTarget->Params(i) = currentParameters[i];
    myTimers[4]->start();
    RealType lastCost(optTarget->Cost(true));
    myTimers[4]->start();
//     if cost function is currently invalid continue
    Valid=optTarget->IsValid;
    if (!ValidCostFunction(Valid))
      continue;
    RealType newCost(lastCost);
    RealType startCost(lastCost);
    Matrix<RealType> Left(N,N);
    Matrix<RealType> Right(N,N);
    Matrix<RealType> S(N,N);
//     stick in wrong matrix to reduce the number of matrices we need by 1.( Left is actually stored in Right, & vice-versa)
    optTarget->fillOverlapHamiltonianMatrices(Right,Left,S);
    bool apply_inverse(true);
    if(apply_inverse)
    {
      Matrix<RealType> RightT(Left);
      invert_matrix(RightT,false);
      Left=0;
      product(RightT,Right,Left);
//       Now the left matrix is the Hamiltonian with the inverse of the overlap applied ot it.
    }
    //Find largest off-diagonal element compared to diagonal element.
    //This gives us an idea how well conditioned it is, used to stabilize.
    RealType od_largest(0);
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
    app_log()<<"od_largest "<<od_largest<<endl;
    if(od_largest>0)
      od_largest = std::log(od_largest);
    else
      od_largest = -1e16;
    if (od_largest<stabilityBase)
      stabilityBase=od_largest;
    else
      stabilizerScale = max( 0.2*(od_largest-stabilityBase)/nstabilizers, stabilizerScale);
    app_log()<<"  stabilityBase "<<stabilityBase<<endl;
    app_log()<<"  stabilizerScale "<<stabilizerScale<<endl;
    int failedTries(0);
    bool acceptedOneMove(false);
    for (int stability=0; stability<nstabilizers; stability++)
    {
      bool goodStep(true);
//       store the Hamiltonian matrix in Right
      for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
          Right(i,j)= Left(j,i);
      RealType XS(stabilityBase+stabilizerScale*(failedTries+stability));
      for (int i=1; i<N; i++)
        Right(i,i) += std::exp(XS);
      app_log()<<"  Using XS:"<<XS<<" "<<failedTries<<" "<<stability<<endl;
      RealType lowestEV(0);
      myTimers[2]->start();
      lowestEV = getLowestEigenvector(Right,currentParameterDirections);
      Lambda = getNonLinearRescale(currentParameterDirections,S);
      myTimers[2]->stop();
//       biggest gradient in the parameter direction vector
      RealType bigVec(0);
      for (int i=0; i<numParams; i++)
        bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
//       this can be overwritten during the line minimization
      RealType evaluated_cost(startCost);
      if (MinMethod=="rescale")
      {
        if (std::abs(Lambda*bigVec)>bigChange)
        {
          goodStep=false;
          app_log()<<"  Failed Step. Magnitude of largest parameter change: "<<std::abs(Lambda*bigVec)<<endl;
          if (stability==0)
          {
            failedTries++;
            stability--;
          }
          else
            stability=nstabilizers;
        }
        for (int i=0; i<numParams; i++)
          optTarget->Params(i) = currentParameters[i] + Lambda*currentParameterDirections[i+1];
        optTarget->IsValid = true;
      }
      else
      {
        for (int i=0; i<numParams; i++)
          optparm[i] = currentParameters[i];
        for (int i=0; i<numParams; i++)
          optdir[i] = currentParameterDirections[i+1];
        TOL = param_tol/bigVec;
        AbsFuncTol=true;
        largeQuarticStep=bigChange/bigVec;
        LambdaMax = 0.5*Lambda;
        myTimers[3]->start();
        if (MinMethod=="quartic")
        {
          int npts(7);
          quadstep = stepsize*Lambda;
          largeQuarticStep=bigChange/bigVec;
          Valid=lineoptimization3(npts,evaluated_cost);
        }
        else
          Valid=lineoptimization2();
        myTimers[3]->stop();
        RealType biggestParameterChange = bigVec*std::abs(Lambda);
        if (biggestParameterChange>bigChange)
        {
          goodStep=false;
          failedTries++;
          app_log()<<"  Failed Step. Largest LM parameter change:"<<biggestParameterChange<<endl;
          if (stability==0)
            stability--;
          else
            stability=nstabilizers;
        }
        else
        {
          for (int i=0; i<numParams; i++)
            optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
          app_log()<<"  Good Step. Largest LM parameter change:"<<biggestParameterChange<<endl;
        }
      }
      if (goodStep)
      {
// 	this may have been evaluated allready
// 	newCost=evaluated_cost;
        //get cost at new minimum
        newCost = optTarget->Cost(false);
        app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<" Delta Cost:"<<(newCost-lastCost)<<endl;
        optTarget->printEstimates();
        //                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
        // mmorales
        Valid=optTarget->IsValid;
        if (!ValidCostFunction(Valid))
        {
          goodStep=false;
          app_log()<<"  Good Step, but cost function invalid"<<endl;
          failedTries++;
          if(stability>0)
            stability=nstabilizers;
          else
            stability--;
        }
        if (newCost < lastCost)
        {
          //Move was acceptable
          for (int i=0; i<numParams; i++)
            bestParameters[i] = optTarget->Params(i);
          lastCost=newCost;
          acceptedOneMove=true;
          if(abs(newCost-lastCost)<1e-4)
          {
            failedTries++;
            stability=nstabilizers;
            continue;
          }
        }
        else if (stability>0)
        {
          failedTries++;
          stability=nstabilizers;
          continue;
        }
      }
      app_log().flush();
      app_error().flush();
      if(failedTries>20) break;
        //APP_ABORT("QMCFixedSampleLinearOptimize::run TOO MANY FAILURES");
    }

    if (acceptedOneMove)
    {
      app_log()<<"Setting new Parameters"<<std::endl;
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = bestParameters[i];
    }
    else
    {
      app_log()<<"Revertting to old Parameters"<<std::endl;
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = currentParameters[i];
    }
    app_log().flush();
    app_error().flush();
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
  string useGPU("yes");
  string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU,"gpu");
  oAttrib.add(vmcMove,"move");
  oAttrib.put(q);
  xmlNodePtr qsave=q;
  xmlNodePtr cur=qsave->children;
  int pid=OHMMS::Controller->rank();
  while (cur != NULL)
  {
    string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
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
      optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
    optTarget->setStream(&app_log());
    success=optTarget->put(q);
  }
  return success;
}

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCFixedSampleLinearOptimize.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

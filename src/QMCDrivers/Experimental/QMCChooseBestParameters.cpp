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
#include "QMCDrivers/QMCChooseBestParameters.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/CommOperators.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#endif
#include "QMCDrivers/VMC/VMCSingle.h"
#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "Numerics/Blasf.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif


namespace qmcplusplus
{

QMCChooseBestParameters::QMCChooseBestParameters(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCDriver(w,psi,h,ppool), CloneManager(hpool), optTarget(w,psi,h,hpool), wfNode(NULL)
{
  //set the optimization flag
  QMCDriverMode.set(QMC_OPTIMIZE,1);
  //read to use vmc output (just in case)
  RootName = "pot";
  QMCType ="QMCChooseBestParameters";
}

/** Clean up the vector */
QMCChooseBestParameters::~QMCChooseBestParameters()
{
}


bool QMCChooseBestParameters::run()
{
  optTarget.initCommunicator(myComm);
  optTarget.setRootName(RootName);
  optTarget.setWaveFunctionNode(wfNode);
  optTarget.startOptimization();
  optTarget.recordParametersToPsi(0.0,0.0);
  optTarget.getAvgParameters(naverage);
//       optTarget.resetPsi(true);
  optTarget.reportParameters();
  return true;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCChooseBestParameters::put(xmlNodePtr q)
{
  xmlNodePtr qsave=q;
  xmlNodePtr cur=qsave->children;
  ParameterSet pAttrib;
  pAttrib.add(naverage,"N","scalar");
  pAttrib.put(q);
  optTarget.setStream(&app_log());
  optTarget.put(qsave);
  bool success=true;
  return success;
}

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCChooseBestParameters.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

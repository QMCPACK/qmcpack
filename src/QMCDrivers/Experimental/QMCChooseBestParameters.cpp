//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


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

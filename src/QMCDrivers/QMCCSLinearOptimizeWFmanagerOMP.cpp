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
    
    


#include "QMCDrivers/QMCCSLinearOptimizeWFmanagerOMP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"
//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus
{

QMCCSLinearOptimizeWFmanagerOMP::QMCCSLinearOptimizeWFmanagerOMP(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool):
  QMCCostFunctionBase(w,psi,h), CloneManager(hpool)
{
  app_log()<<" Using QMCCSLinearOptimizeWFmanagerOMP::QMCCSLinearOptimizeWFmanagerOMP"<< std::endl;
}


/** Clean up the vector */
QMCCSLinearOptimizeWFmanagerOMP::~QMCCSLinearOptimizeWFmanagerOMP()
{
}

void QMCCSLinearOptimizeWFmanagerOMP::startOptimization()
{
  ReportCounter=0;
  for (int i=0; i<psiClones.size(); ++i)
    psiClones[i]->startOptimization();
  OptVariablesForPsi.setComputed();
}

void QMCCSLinearOptimizeWFmanagerOMP::resetPsi(bool final_reset)
{
  if (OptVariables.size() < OptVariablesForPsi.size())
  {
    for (int i=0; i<equalVarMap.size(); ++i)
      OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
  }
  else
    for (int i=0; i<OptVariables.size(); ++i)
      OptVariablesForPsi[i]=OptVariables[i];
  if (final_reset)
  {
    for (int i=0; i<psiClones.size(); ++i)
      psiClones[i]->stopOptimization();
  }
  //cout << "######### QMCCSLinearOptimizeWFmanagerOMP::resetPsi " << std::endl;
//     OptVariablesForPsi.print(std::cout);
  //cout << "-------------------------------------- " << std::endl;
  Psi.resetParameters(OptVariablesForPsi);
  for (int i=0; i<psiClones.size(); ++i)
    psiClones[i]->resetParameters(OptVariablesForPsi);
//     for (int i=0; i<psiClones.size(); ++i)
//       psiClones[i]->reportStatus(app_log());
}
}

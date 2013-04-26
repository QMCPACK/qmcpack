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
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  app_log()<<" Using QMCCSLinearOptimizeWFmanagerOMP::QMCCSLinearOptimizeWFmanagerOMP"<<endl;
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
    #pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      for (; it!=it_end; ++it)
        (**it).DataSetForDerivatives.clear();
    }
// is this correct with OMP?
//       MCWalkerConfiguration::iterator it(W.begin());
//       MCWalkerConfiguration::iterator it_end(W.end());
//       for (; it!=it_end; ++it)
//         (**it).DataSetForDerivatives.clear();
  }
  //cout << "######### QMCCSLinearOptimizeWFmanagerOMP::resetPsi " << endl;
//     OptVariablesForPsi.print(cout);
  //cout << "-------------------------------------- " << endl;
  Psi.resetParameters(OptVariablesForPsi);
  for (int i=0; i<psiClones.size(); ++i)
    psiClones[i]->resetParameters(OptVariablesForPsi);
//     for (int i=0; i<psiClones.size(); ++i)
//       psiClones[i]->reportStatus(app_log());
}
}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
* $Id: QMCCSLinearOptimizeWFmanagerOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $
***************************************************************************/

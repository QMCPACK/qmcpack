//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file QMCMain.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCApp/QMCMain.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/Timer.h"
#include "QMCDrivers/QMCDriver.h"
//#include "QMCDrivers/CMC/DummyIonMove.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "qmc_common.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

bool QMCMain::executeDebugSection(xmlNodePtr cur)
{
  app_log() << "QMCMain::executeDebugSection " << endl;
  app_log() << "  Use this to debug new features with <debug/> in the input file " << endl;

  return true;
}

bool QMCMain::executeCMCSection(xmlNodePtr cur)
{
  bool success=true;
  string target("ion0");
  OhmmsAttributeSet a;
  a.add(target,"target");
  a.put(cur);

  MCWalkerConfiguration *ions = ptclPool->getWalkerSet(target);
  TrialWaveFunction* primaryPsi=psiPool->getPrimary();
  QMCHamiltonian* primaryH=hamPool->getPrimary();

  app_log() << "QMCMain::executeCMCSection moving " << target << " by dummy move." << endl;
  //DummyIonMove dummy(*ions,*primaryPsi,*primaryH,*hamPool,*psiPool,qmcDriver);
  //dummy.run();

  int nat = ions->getTotalNum();
  ParticleSet::ParticlePos_t deltaR(nat);

  makeGaussRandomWithEngine(deltaR,Random); //generate random displaement

  //update QMC system
  qmcSystem->update();

  double logpsi1 = primaryPsi->evaluateLog(*qmcSystem);
  cout << "logpsi1 " << logpsi1 << endl;

  double eloc1  = primaryH->evaluate(*qmcSystem);
  cout << "Local Energy " << eloc1 << endl;

  for (int i=0; i<primaryH->sizeOfObservables(); i++)
    app_log() << "  HamTest " << primaryH->getObservableName(i) << " " << primaryH->getObservable(i) << endl;

  for(int iat=0; iat<nat; ++iat)
  {
    ions->R[iat]+=deltaR[iat];

    ions->update(); //update position and distance table of itself 
    primaryH->update_source(*ions);

    qmcSystem->update();
    double logpsi2 = primaryPsi->evaluateLog(*qmcSystem);
    double eloc2  = primaryH->evaluate(*qmcSystem);

    cout << "\nION " << iat << " " << ions->R[iat] << endl;
    cout << "logpsi " << logpsi2 << endl;
    cout << "Local Energy " << eloc2 << endl;
    for (int i=0; i<primaryH->sizeOfObservables(); i++)
      app_log() << "  HamTest " << primaryH->getObservableName(i) << " " << primaryH->getObservable(i) << endl;

    ions->R[iat]-=deltaR[iat];
    ions->update(); //update position and distance table of itself 
    primaryH->update_source(*ions);

    qmcSystem->update();
    double logpsi3 = primaryPsi->evaluateLog(*qmcSystem);
    double eloc3  = primaryH->evaluate(*qmcSystem);

    if(abs(eloc1-eloc3)>1e-12)
    {
      cout << "ERROR Energies are different " << endl;
    }
  }

  //
  //for(int i=0; i<ions->getTotalNum(); ++i)
  //{
  //  ions->R[i]+=deltaR(i);
  //  ions->update();
  //  primaryH->update_source(*ions);
  //  qmcDriver->run();
  //}
  return success;
}

}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author: jnkim $
 * $Revision: 5845 $   $Date: 2013-05-14 16:10:18 -0400 (Tue, 14 May 2013) $
 * $Id: QMCMain.cpp 5845 2013-05-14 20:10:18Z jnkim $
 ***************************************************************************/

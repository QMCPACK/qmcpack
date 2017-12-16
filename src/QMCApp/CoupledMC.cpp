//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
#include "Utilities/Timer.h"
#include "QMCDrivers/QMCDriver.h"
//#include "QMCDrivers/CMC/DummyIonMove.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "qmc_common.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

bool QMCMain::executeDebugSection(xmlNodePtr cur)
{
  app_log() << "QMCMain::executeDebugSection " << std::endl;
  app_log() << "  Use this to debug new features with <debug/> in the input file " << std::endl;

  return true;
}

bool QMCMain::executeCMCSection(xmlNodePtr cur)
{
  bool success=true;
  std::string target("ion0");
  OhmmsAttributeSet a;
  a.add(target,"target");
  a.put(cur);

  MCWalkerConfiguration *ions = ptclPool->getWalkerSet(target);
  TrialWaveFunction* primaryPsi=psiPool->getPrimary();
  QMCHamiltonian* primaryH=hamPool->getPrimary();

  app_log() << "QMCMain::executeCMCSection moving " << target << " by dummy move." << std::endl;
  //DummyIonMove dummy(*ions,*primaryPsi,*primaryH,*hamPool,*psiPool,qmcDriver);
  //dummy.run();

  int nat = ions->getTotalNum();
  ParticleSet::ParticlePos_t deltaR(nat);

  makeGaussRandomWithEngine(deltaR,Random); //generate random displaement

  //update QMC system
  qmcSystem->update();

  double logpsi1 = primaryPsi->evaluateLog(*qmcSystem);
  std::cout << "logpsi1 " << logpsi1 << std::endl;

  double eloc1  = primaryH->evaluate(*qmcSystem);
  std::cout << "Local Energy " << eloc1 << std::endl;

  for (int i=0; i<primaryH->sizeOfObservables(); i++)
    app_log() << "  HamTest " << primaryH->getObservableName(i) << " " << primaryH->getObservable(i) << std::endl;

  for(int iat=0; iat<nat; ++iat)
  {
    ions->R[iat]+=deltaR[iat];

    ions->update(); //update position and distance table of itself 
    primaryH->update_source(*ions);

    qmcSystem->update();
    double logpsi2 = primaryPsi->evaluateLog(*qmcSystem);
    double eloc2  = primaryH->evaluate(*qmcSystem);

    std::cout << "\nION " << iat << " " << ions->R[iat] << std::endl;
    std::cout << "logpsi " << logpsi2 << std::endl;
    std::cout << "Local Energy " << eloc2 << std::endl;
    for (int i=0; i<primaryH->sizeOfObservables(); i++)
      app_log() << "  HamTest " << primaryH->getObservableName(i) << " " << primaryH->getObservable(i) << std::endl;

    ions->R[iat]-=deltaR[iat];
    ions->update(); //update position and distance table of itself 
    primaryH->update_source(*ions);

    qmcSystem->update();
    double logpsi3 = primaryPsi->evaluateLog(*qmcSystem);
    double eloc3  = primaryH->evaluate(*qmcSystem);

    if(std::abs(eloc1-eloc3)>1e-12)
    {
      std::cout << "ERROR Energies are different " << std::endl;
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

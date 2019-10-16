//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCDriverFactory.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCApp/QMCMainState.h"
#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/RMC/RMCFactory.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "QMCDrivers/QMCCorrelatedSamplingLinearOptimize.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include <queue>
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
///initialize the static data member
//ParticleSetPool* QMCMainState::ptclPool = new ParticleSetPool;
QMCMainState::QMCMainState(Communicate* c) : MPIObjectBase(c), curRunType(QMCRunType::DUMMY), qmcSystem(0)
{
  ////create ParticleSetPool
  ptclPool = new ParticleSetPool(myComm);
  //create WaveFunctionPool
  psiPool = new WaveFunctionPool(myComm);
  psiPool->setParticleSetPool(ptclPool);
  //create HamiltonianPool
  hamPool = new HamiltonianPool(myComm);
  hamPool->setParticleSetPool(ptclPool);
  hamPool->setWaveFunctionPool(psiPool);
}

QMCMainState::~QMCMainState()
{
  delete hamPool;
  delete psiPool;
  delete ptclPool;
}

void QMCMainState::putCommunicator(xmlNodePtr cur)
{
  //BROKEN: myComm is ALWAYS initialized by the constructor
  if (myComm)
    return;
  ParameterSet params;
  int nparts = 1;
  params.add(nparts, "groups", "int");
  params.add(nparts, "twistAngles", "int");
  params.put(cur);
  if (nparts > 1)
  {
    app_log() << "  Communiator groups = " << nparts << std::endl;
    myComm = new Communicate(*OHMMS::Controller, nparts);
  }
}

} // namespace qmcplusplus

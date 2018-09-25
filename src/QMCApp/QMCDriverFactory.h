//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#define QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCDrivers/QMCDriver.h"
#include "QMCApp/QMCDriverFactoryInterface.h"
#include "QMCApp/HamiltonianPoolInterface.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include <bitset>

class Communicate;

namespace qmcplusplus
{

extern template class HamiltonianPool<Batching::SINGLE>;
extern template class HamiltonianPool<Batching::BATCHED>;
extern template class QMCOptimize<Batching::SINGLE>;
extern template class QMCOptimize<Batching::BATCHED>;
extern template class WaveFunctionTester<Batching::SINGLE>;
extern template class WaveFunctionTester<Batching::BATCHED>;

  
//forward declaration
class MCWalkerConfiguration;


  
template<Batching batching = Batching::SINGLE>
struct QMCDriverFactory: public QMCDriverFactoryInterface, public MPIObjectBase
{
  ///current QMC mode determined by curQmcModeBits
  unsigned long curQmcMode;

  ///8-bit (ALTERNATE_MODE,SPACEWARP_MODE, MULTIPLE_MODE, UPDATE_MODE)
  std::bitset<QMC_MODE_MAX> curQmcModeBits;

  ///type of qmcdriver
  QMCRunType curRunType;

  ///name of the current QMCriver
  std::string curMethod;

  /** current MCWalkerConfiguration
   */
  MCWalkerConfiguration *qmcSystem;

  /** current QMCDriver
   */
  QMCDriver<batching> *qmcDriver;

  /** ParticleSet Pool
   */
  ParticleSetPool* ptclPool;

  /** TrialWaveFunction Pool
   */
  WaveFunctionPool* psiPool;

  /** QMCHamiltonian Pool
   */
  HamiltonianPool<batching>* hamPool;

  /** default constructor **/
  // QMCDriverFactory(Communicate* c);

  QMCDriverFactory(Communicate* c): MPIObjectBase(c),
				    qmcSystem(0), qmcDriver(0) , curRunType(DUMMY_RUN)
  {
    ////create ParticleSetPool
    ptclPool = new ParticleSetPool(c);
    //create WaveFunctionPool
    psiPool = new WaveFunctionPool(c);
    psiPool->setParticleSetPool(ptclPool);
    //create HamiltonianPool
    hamPool = new HamiltonianPool<batching>(c);
    hamPool->setParticleSetPool(ptclPool);
    hamPool->setWaveFunctionPool(psiPool);
  }

  
  /** set the active qmcDriver */
  void putCommunicator(xmlNodePtr cur);

  /** set the active qmcDriver */
  bool setQMCDriver(int curSeries, xmlNodePtr cur);

  void checkQMCSystem(const std::string& target);
  /** create a new QMCDriver
   */
  void createQMCDriver(xmlNodePtr cur);

  /** virtual destructor **/
  virtual ~QMCDriverFactory();

  void updateQMCSystem() { qmcSystem->update(); };

  //Accessors
  std::string& getMethod() { return curMethod; }
  ParticleSetPool& getParticleSetPool() { return *ptclPool; }
  ParticleSetPool* getParticleSetPoolPtr() { return ptclPool; }
  HamiltonianPoolInterface& getHamiltonianPool() { return *dynamic_cast<HamiltonianPoolInterface*>(hamPool); }
  WaveFunctionPool& getWaveFunctionPool() { return *psiPool; }
  bool driverExists() { return bool(qmcDriver); }
  void setMethod(const std::string& method) { curMethod = method; }
  QMCDriverInterface& getQMCDriver() { return *dynamic_cast<QMCDriverInterface*>(qmcDriver); }
  MCWalkerConfiguration * getQMCSystem() { return qmcSystem; }
};
}
#endif


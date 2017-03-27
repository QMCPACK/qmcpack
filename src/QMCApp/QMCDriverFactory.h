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
#include <bitset>

class Communicate;

namespace qmcplusplus
{

//forward declaration
class MCWalkerConfiguration;
class QMCDriver;
class WaveFunctionPool;
class HamiltonianPool;

struct QMCDriverFactory: public MPIObjectBase
{
  /*! enum for QMC Run Type */
  enum QMCRunType
  {
    DUMMY_RUN, /*!< dummy */
    VMC_RUN, /**< VMC type: vmc, vmc-ptcl, vmc-multiple, vmc-ptcl-multiple */
    CSVMC_RUN,
    DMC_RUN, /**< DMC type: dmc, dmc-ptcl*/
    RMC_RUN, /**< RMC type: rmc, rmc-ptcl */
    OPTIMIZE_RUN,/*!< Optimization */
    VMC_OPT_RUN, /*!< Optimization with vmc blocks */
    LINEAR_OPTIMIZE_RUN,
    CS_LINEAR_OPTIMIZE_RUN,
    WF_TEST_RUN
  };

  /*! enum to set the bit to determine the QMC mode */
  enum QMCModeEnum
  {
    UPDATE_MODE,  /**< bit for move: walker or pbyp */
    MULTIPLE_MODE, /**< bit for multple configuration */
    SPACEWARP_MODE, /**< bit for space-warping */
    ALTERNATE_MODE, /**< bit for performing various analysis and weird qmc methods */
    GPU_MODE,     /**< bit to use GPU driver */
    QMC_MODE_MAX=8
  };

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
  QMCDriver *qmcDriver;

  /** ParticleSet Pool
   */
  ParticleSetPool* ptclPool;

  /** TrialWaveFunction Pool
   */
  WaveFunctionPool* psiPool;

  /** QMCHamiltonian Pool
   */
  HamiltonianPool* hamPool;

  /** default constructor **/
  QMCDriverFactory(Communicate* c);

  /** set the active qmcDriver */
  void putCommunicator(xmlNodePtr cur);

  /** set the active qmcDriver */
  bool setQMCDriver(int curSeries, xmlNodePtr cur);

  /** create a new QMCDriver
   */
  void createQMCDriver(xmlNodePtr cur);

  /** virtual destructor **/
  virtual ~QMCDriverFactory();
};
}
#endif


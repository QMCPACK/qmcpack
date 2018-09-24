//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMCDRIVER_FACTORY_INTERFACE_H
#define QMCPLUSPLUS_QMCDRIVER_FACTORY_INTERFACE_H
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCDrivers/QMCDriver.h"

#include <bitset>

class Communicate;

namespace qmcplusplus
{

//forward declaration
class MCWalkerConfiguration;


class QMCDriverFactoryInterface
{
public:
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

  /** set the active qmcDriver */
  virtual void putCommunicator(xmlNodePtr cur) = 0;

  /** set the active qmcDriver */
  virtual bool setQMCDriver(int curSeries, xmlNodePtr cur) = 0;

  virtual void checkQMCSystem(const std::string&) = 0;

  /** create a new QMCDriver
   */
  virtual void createQMCDriver(xmlNodePtr cur) = 0;
  
  //Accessors
  virtual std::string& getMethod() = 0;
  virtual ParticleSetPool& getParticleSetPool() = 0;
  virtual bool driverExists() = 0;
  // This should be a smart pointer, since some other objects insist on keeping it
  virtual ParticleSetPool* getParticleSetPoolPtr() = 0;
  virtual HamiltonianPoolInterface& getHamiltonianPool() = 0;
  virtual WaveFunctionPool& getWaveFunctionPool() = 0;
  virtual void setMethod(const std::string&) = 0;
  virtual QMCDriverInterface& getQMCDriver() = 0;
  virtual void updateQMCSystem() = 0;
  virtual MCWalkerConfiguration* getQMCSystem() = 0;
};
}
#endif


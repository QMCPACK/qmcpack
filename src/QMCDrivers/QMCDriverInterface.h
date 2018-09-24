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

#ifndef QMCPLUSPLUS_QMCDRIVER_INTERFACE_H
#define QMCPLUSPLUS_QMCDRIVER_INTERFACE_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/PooledData.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerBase.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/BranchIO.h"
#include "QMCWaveFunctions/Batching.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
class Communicate;

namespace qmcplusplus
{
/** @defgroup QMCDrivers QMC Driver group
 * QMC drivers that implement QMC algorithms
 */

/** @defgroup WalkerByWalker QMC Drivers using walker-by-walker update
 * @brief Move all the particles for each walker
 */

/** @defgroup ParticleByParticle QMC Drivers using particle-by-particle update
 * @brief Move particle by particle
 */

/** @defgroup MultiplePsi QMC Drivers for energy differences
 * @brief Umbrella sampling over multiple H/Psi
 *
 * This class of QMC drivers are suitable to evaluate
 * the energy differences of multiple H-Psi pairs.
 */

//forward declarations: Do not include headers if not needed
class MCWalkerConfiguration;
class HDFWalkerOutput;
class TraceManager;


/** @ingroup QMCDrivers
 * @{
 * @brief base class for QMC engines
 */

class QMCDriverInterface
{
public:
  using BranchEngineType = SimpleFixedNodeBranch;
  using QMCT = QMCTraits;
  ///return current step
  virtual int current() const = 0;

  /** set the update mode
   * @param pbyp if true, use particle-by-particle update
   */
  virtual void setUpdateMode(bool pbyp) = 0;

  /** Set the status of the QMCDriver
   * @param aname the root file name
   * @param h5name root name of the master hdf5 file containing previous qmcrun
   * @param append if ture, the run is a continuation of the previous qmc
   *
   * All output files will be of
   * the form "aname.s00X.suffix", where "X" is number
   * of previous QMC runs for the simulation and "suffix"
   * is the suffix for the output file.
   */
  virtual void setStatus(const std::string& aname, const std::string& h5name, bool append) = 0;

  /** add QMCHamiltonian/TrialWaveFunction pair for multiple
   * @param h QMCHamiltonian
   * @param psi TrialWaveFunction
   *
   * *Multiple* drivers use multiple H/Psi pairs to perform correlated sampling
   * for energy difference evaluations.
   */
  virtual void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction<>* psi) = 0;

  /** initialize with xmlNode
   */
  virtual void process(xmlNodePtr cur) = 0;

  /** return a xmlnode with update **/
  virtual xmlNodePtr getQMCNode() = 0;

  virtual void putWalkers(std::vector<xmlNodePtr>& wset) = 0;

  virtual void putTraces(xmlNodePtr txml) = 0;

  virtual bool run() = 0;

  virtual bool put(xmlNodePtr cur) = 0;

  virtual std::string getEngineName() const = 0;

  ///set the BranchEngineType
  virtual void setBranchEngine(BranchEngineType* be) = 0;

  ///return BranchEngineType*
  virtual BranchEngineType* getBranchEngine() = 0;

  virtual int addObservable(const std::string& aname) = 0;

  virtual QMCT::RealType getObservable(int i) = 0;

  virtual void setTau(QMCT::RealType i) = 0;

  ///resetComponents for next run if reusing a driver.
  virtual void resetComponents(xmlNodePtr cur) = 0;

  ///set global offsets of the walkers
  virtual void setWalkerOffsets() = 0;

  //virtual std::vector<RandomGenerator_t*>& getRng() {}

  ///return the random generators
  virtual std::vector<RandomGenerator_t*>& getRng() = 0;

  ///return the i-th random generator
  virtual RandomGenerator_t& getRng(int i) = 0;

  virtual bool putQMCInfo(xmlNodePtr cur) = 0;

  virtual void addWalkers(int nwalkers) = 0;

  /** record the state of the block
   * @param block current block
   *
   * virtual function with a default implementation
   */
  virtual void recordBlock(int block) = 0;

  /** finalize a qmc section
   * @param block current block
   * @param dumpwalkers if true, dump walkers
   *
   * Accumulate energy and weight is written to a hdf5 file.
   * Finialize the estimators
   */
  virtual bool finalize(int block, bool dumpwalkers=true) = 0;

  virtual void adiosCheckpoint(int block) = 0;
  virtual void adiosCheckpointFinal(int block, bool dumpwalkers) = 0;
  virtual std::string getRotationName( std::string RootName) = 0;
  virtual std::string getLastRotationName( std::string RootName) = 0;
  
};
/**@}*/
}

#endif

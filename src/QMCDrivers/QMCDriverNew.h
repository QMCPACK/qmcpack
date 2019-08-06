//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.h
//////////////////////////////////////////////////////////////////////////////////////


/**
 * @file
 * Declaration of QMCDriverNew
 *
 * This will replace QMCDriver once unified drivers are finished
 * the general documentation from QMCDriver.h must be moved before then
 */
#ifndef QMCPLUSPLUS_QMCDRIVERNEW_H
#define QMCPLUSPLUS_QMCDRIVERNEW_H

#include <type_traits>

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/PooledData.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerBase.h"
#include "Particle/MCPopulation.h"
#include "QMCDrivers/Crowd.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/BranchIO.h"
class Communicate;

namespace qmcplusplus
{
//forward declarations: Do not include headers if not needed
class HDFWalkerOutput;
class TraceManager;

/** @ingroup QMCDrivers
 * @{
 * @brief QMCDriverNew Base class for Unified Drivers
 */
class QMCDriverNew : public QMCDriverInterface, public MPIObjectBase
{
public:
  using RealType  = QMCTraits::RealType;
  using IndexType = QMCTraits::IndexType;
  using FullPrecisionRealType = QMCTraits::FullPrecRealType;

  /** @ingroup Type dependent behavior
   * @{
   * @brief use simple metaprogramming in anticipation of single executable
   */
  /** call recompute at the end of each block in the mixed precision case.
   */
  template<typename RT = RealType, typename FPRT = FullPrecisionRealType>
  int defaultBlocksBetweenRecompute() { return 0; }

  template<typename RT = RealType, typename FPRT = FullPrecisionRealType,
	   std::enable_if_t<std::is_same<RT, FPRT>{}>>
  int defaultBlocksBetweenRecompute() { return 1; }
  /** @}
   */
  
  /** separate but similar to QMCModeEnum
   *  
   *  a code smell
   */
  enum
  {
    QMC_UPDATE_MODE,
    QMC_MULTIPLE,
    QMC_OPTIMIZE,
    QMC_WARMUP
  };

  using MCPWalker = MCPopulation::PopulationWalker;
  using Buffer    = MCPWalker::Buffer_t;

  /** bits to classify QMCDriver
   *
   * - qmc_driver_mode[QMC_UPDATE_MODE]? particle-by-particle: walker-by-walker
   * - qmc_driver_mode[QMC_MULTIPLE]? multiple H/Psi : single H/Psi
   * - qmc_driver_mode[QMC_OPTIMIZE]? optimization : vmc/dmc/rmc
   */
  std::bitset<QMC_MODE_MAX> qmc_driver_mode;

  /// whether to allow traces
  bool allow_traces;
  /// traces xml
  xmlNodePtr traces_xml;

  /// Constructor.
  QMCDriverNew(MCPopulation& population,
               TrialWaveFunction& psi,
               QMCHamiltonian& h,
               WaveFunctionPool& ppool,
               Communicate* comm);

  virtual ~QMCDriverNew();

  ///return current step
  inline int current() const { return current_step; }

  /** set the update mode
   * @param pbyp if true, use particle-by-particle update
   */
  inline void setUpdateMode(bool pbyp) { qmc_driver_mode[QMC_UPDATE_MODE] = pbyp; }

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
  void setStatus(const std::string& aname, const std::string& h5name, bool append);

  /** add QMCHamiltonian/TrialWaveFunction pair for multiple
   * @param h QMCHamiltonian
   * @param psi TrialWaveFunction
   *
   * *Multiple* drivers use multiple H/Psi pairs to perform correlated sampling
   * for energy difference evaluations.
   */
  void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi);

  /** initialize with xmlNode
   */
  void process(xmlNodePtr cur);

  /** return a xmlnode with update **/
  xmlNodePtr getQMCNode();

  void putWalkers(std::vector<xmlNodePtr>& wset);

  inline void putTraces(xmlNodePtr txml) { traces_xml = txml; };

  inline void requestTraces(bool traces) { allow_traces = traces; }

  inline std::string getEngineName() const { return QMCType; }

  template<class PDT>
  void setValue(const std::string& aname, PDT x)
  {
    m_param.setValue(aname, x);
  }

  ///set the BranchEngineType
  void setBranchEngine(SimpleFixedNodeBranch* be) { branchEngine = be; }

  ///return BranchEngineType*
  SimpleFixedNodeBranch* getBranchEngine() { return branchEngine; }

  int addObservable(const std::string& aname);

  RealType getObservable(int i);

  void setTau(RealType i) { Tau = i; }

  ///resetComponents for next run if reusing a driver.
  virtual void resetComponents(xmlNodePtr cur)
  {
    qmc_node = cur;
    m_param.put(cur);
  }

  ///set global offsets of the walkers
  void setWalkerOffsets();

  //virtual std::vector<RandomGenerator_t*>& getRng() {}

  ///Observables manager
  EstimatorManagerBase* Estimators;

  ///Traces manager
  TraceManager* Traces;

  ///return the random generators
  inline std::vector<RandomGenerator_t*>& getRng() { return Rng; }

  ///return the i-th random generator
  inline RandomGenerator_t& getRng(int i) { return (*Rng[i]); }

  std::string getEngineName() { return QMCType; }
  unsigned long getDriverMode() { return qmc_driver_mode.to_ulong(); }

protected:
  std::vector<Crowd> crowds_;
  ///branch engine
  SimpleFixedNodeBranch* branchEngine;
  ///drift modifer
  DriftModifierBase* DriftModifier;
  ///randomize it
  bool reset_random;
  ///flag to append or restart the run
  bool append_run;
  ///flag to turn off dumping configurations
  bool dump_config;
  ///true, if the size of population is fixed.
  bool const_population;

  ///the number of blocks to be rolled back
  int RollBackBlocks;
  ///the number to delay updates by
  int k_delay;
  /** period of dumping walker configurations and everything else for restart
   *
   * The unit is a block.
   */
  int check_point_period;
  /** period of dumping walker positions and IDs for Forward Walking
  *
  * The unit is in steps.
  */
  int store_config_period;

  ///Period to recalculate the walker properties from scratch.
  int recalculate_properties_period;

  /** period of recording walker configurations
   *
   * Default is 0 indicating that only the last configuration will be saved.
   */
  int walker_dump_period;

  /** period of recording walker positions and IDs for forward walking afterwards
   *
   */
  int config_dump_period;

  IndexType current_step;
  IndexType max_blocks;
  IndexType max_steps;

  ///number of steps between a step: VMCs do not evaluate energies
  IndexType nSubSteps;

  ///number of warmup steps
  IndexType nWarmupSteps;

  ///counter for number of moves accepted
  IndexType nAccept;

  ///counter for number of moves /rejected
  IndexType nReject;

  /// the number of blocks between recomptePsi
  IndexType nBlocksBetweenRecompute;

  ///the number of walkers
  IndexType nTargetWalkers;
  ///the number of saved samples
  IndexType nTargetSamples;
  ///alternate method of setting QMC run parameters
  IndexType nStepsBetweenSamples;
  ///samples per thread
  RealType nSamplesPerThread;
  ///target population
  RealType nTargetPopulation;

  ///timestep
  RealType Tau;

  ///maximum cpu in secs
  RealType MaxCPUSecs;

  ///Time-step factor \f$ 1/(2\tau)\f$
  RealType m_oneover2tau;
  ///Time-step factor \f$ \sqrt{\tau}\f$
  RealType m_sqrttau;

  ///type of qmc: assigned by subclasses
  std::string QMCType;
  ///the root of h5File
  std::string h5FileRoot;
  ///root of all the output files
  std::string RootName;

  ///store any parameter that has to be read from a file
  ParameterSet m_param;

  ///record engine for walkers
  HDFWalkerOutput* wOut;
  ///the entire (or on node) walker population
  MCPopulation& population_;

  ///trial function
  TrialWaveFunction& Psi;

  WaveFunctionPool& psiPool;

  ///Hamiltonian
  QMCHamiltonian& H;

  ///a list of TrialWaveFunctions for multiple method
  std::vector<TrialWaveFunction*> Psi1;

  ///a list of QMCHamiltonians for multiple method
  std::vector<QMCHamiltonian*> H1;

  ///Random number generators
  std::vector<RandomGenerator_t*> Rng;

  ///a list of mcwalkerset element
  std::vector<xmlNodePtr> mcwalkerNodePtr;

  ///a list of timers
  std::vector<NewTimer*> myTimers;

  ///temporary storage for drift
  ParticleSet::ParticlePos_t drift;

  ///temporary storage for random displacement
  ParticleSet::ParticlePos_t deltaR;


protected:
  // I suspect all of this state is unecessary

  ///pointer to qmc node in xml file
  xmlNodePtr qmc_node;


public:
  ///temporary buffer to accumulate data
  //ostrstream log_buffer;

  //PooledData<RealType> HamPool;

  ///Copy Constructor (disabled).
  QMCDriverNew(const QMCDriverNew&) = delete;
  ///Copy operator (disabled).
  QMCDriverNew& operator=(const QMCDriverNew&) = delete;

  bool putQMCInfo(xmlNodePtr cur);

  void addWalkers(int nwalkers);

  /** record the state of the block
   * @param block current block
   *
   * virtual function with a default implementation
   */
  virtual void recordBlock(int block);

  /** finalize a qmc section
   * @param block current block
   * @param dumpwalkers if true, dump walkers
   *
   * Accumulate energy and weight is written to a hdf5 file.
   * Finialize the estimators
   */
  bool finalize(int block, bool dumpwalkers = true);

  int rotation;
  void adiosCheckpoint(int block);
  void adiosCheckpointFinal(int block, bool dumpwalkers);
  std::string getRotationName(std::string RootName);
  std::string getLastRotationName(std::string RootName);

  NewTimer* checkpointTimer;
};
/**@}*/
} // namespace qmcplusplus

#endif

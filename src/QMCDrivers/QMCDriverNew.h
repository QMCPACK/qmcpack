//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
 *  
 * This driver base class should be generic with respect to precision,
 * value type, device execution, and ...
 * It should contain no typdefs not related to compiler bugs or platform workarounds
 *
 */

#ifndef QMCPLUSPLUS_QMCDRIVERNEW_H
#define QMCPLUSPLUS_QMCDRIVERNEW_H

#include <type_traits>

#include "Configuration.h"
#include "Utilities/PooledData.h"
#include "Utilities/TimerManager.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/EstimatorManagerNew.h"
#include "QMCDrivers/MCPopulation.h"
#include "QMCDrivers/Crowd.h"
#include "QMCDrivers/QMCDriverInterface.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "QMCDrivers/SFNBranch.h"
#include "QMCDrivers/BranchIO.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/ContextForSteps.h"

class Communicate;

namespace qmcplusplus
{
//forward declarations: Do not include headers if not needed
class HDFWalkerOutput;
class TraceManager;
struct SFNBranch;

namespace testing
{
class DMCBatchedTest;
class VMCBatchedTest;
class QMCDriverNewTestWrapper;
} // namespace testing

/** @ingroup QMCDrivers
 * @{
 * @brief QMCDriverNew Base class for Unified Drivers
 *
 * # General Principals
 * * Parameters used unchanged from input object are not copied into class state
 * * The driver state machine should be as minimal as possible.
 * * In non performance critical areas favor clarity over clever optimizations.
 */
class QMCDriverNew : public QMCDriverInterface, public MPIObjectBase
{
public:
  using RealType         = QMCTraits::RealType;
  using IndexType        = QMCTraits::IndexType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
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

  using MCPWalker = MCPopulation::MCPWalker;
  using WFBuffer  = MCPopulation::WFBuffer;

  using SetNonLocalMoveHandler = std::function<void(QMCHamiltonian&)>;
  /** bits to classify QMCDriver
   *
   * - qmc_driver_mode[QMC_UPDATE_MODE]? particle-by-particle: walker-by-walker
   * - qmc_driver_mode[QMC_MULTIPLE]? multiple H/Psi : single H/Psi
   * - qmc_driver_mode[QMC_OPTIMIZE]? optimization : vmc/dmc/rmc
   */
  std::bitset<QMC_MODE_MAX> qmc_driver_mode_;

protected:
  void endBlock();
  /** This is a data structure strictly for QMCDriver and its derived classes
   *
   *  i.e. its nested in scope for a reason
   */
  struct AdjustedWalkerCounts
  {
    IndexType global_walkers;
    std::vector<IndexType> walkers_per_rank;
    std::vector<IndexType> walkers_per_crowd;
    RealType reserve_walkers;
  };

public:
  /// Constructor.
  QMCDriverNew(QMCDriverInput&& input,
               MCPopulation& population,
               TrialWaveFunction& psi,
               QMCHamiltonian& h,
               WaveFunctionPool& ppool,
               const std::string timer_prefix,
               Communicate* comm,
               const std::string& QMC_driver_type,
               SetNonLocalMoveHandler = &QMCDriverNew::defaultSetNonLocalMoveHandler);

  QMCDriverNew(QMCDriverNew&&) = default;

  virtual ~QMCDriverNew();

  ///return current step
  inline IndexType current() const { return current_step_; }

  /** Set the status of the QMCDriver
   * @param aname the root file name
   * @param h5name root name of the master hdf5 file containing previous qmcrun
   * @param append if true, the run is a continuation of the previous qmc
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

  void createRngsStepContexts(int num_crowds);

  void putWalkers(std::vector<xmlNodePtr>& wset);

  /** placate the legacy base class interface
   */
  void setBranchEngine(SimpleFixedNodeBranch* be)
  {
    throw std::runtime_error("You can not use the legacy SimpleFixedNodeBranch class with QMCDriverNew");
  }

  ///set the BranchEngineType
  void setNewBranchEngine(SFNBranch* be) { branch_engine_ = be; }

  /** placate the legacy base class interface
   */
  SimpleFixedNodeBranch* getBranchEngine() { return nullptr; }

  ///return BranchEngineType*
  SFNBranch* getNewBranchEngine() { return branch_engine_; }

  int addObservable(const std::string& aname);

  RealType getObservable(int i);

  ///set global offsets of the walkers
  void setWalkerOffsets();

  std::vector<RandomGenerator_t*> RngCompatibility;

  inline std::vector<RandomGenerator_t*>& getRng() { return RngCompatibility; }

  // ///return the random generators
  //       inline std::vector<std::unique_ptr RandomGenerator_t*>& getRng() { return Rng; }

  ///return the i-th random generator
  inline RandomGenerator_t& getRng(int i) { return (*Rng[i]); }

  std::string getEngineName() { return QMCType; }
  unsigned long getDriverMode() { return qmc_driver_mode_.to_ulong(); }

  IndexType get_living_walkers() const { return population_.get_walkers().size(); }

  /** @ingroup Legacy interface to be dropped
   *  @{
   */
  bool put(xmlNodePtr cur) { return false; };

  /** QMCDriverNew driver second (3rd, 4th...) stage of constructing a valid driver
   *
   *  This is the shared entry point with legacy,
   *  from QMCMain so the API cannot be updated yet
   *
   *  \todo remove cur, the driver and all its child nodes should be completely processed before
   *        this stage of driver initialization is hit.
   */
  virtual void process(xmlNodePtr cur) = 0;

  /** Do common section starting tasks
   *
   *  \todo This should not take xmlNodePtr
   *        It should either take BranchEngineInput and EstimatorInput
   *        And these are the arguments to the branch_engine and estimator_manager
   *        Constructors or these objects should be created elsewhere.
   */
  void startup(xmlNodePtr cur, QMCDriverNew::AdjustedWalkerCounts awc);

  static void initialLogEvaluation(int crowd_id, UPtrVector<Crowd>& crowds, UPtrVector<ContextForSteps>& step_context);


  /** should be set in input don't see a reason to set individually
   * @param pbyp if true, use particle-by-particle update
   */
  inline void setUpdateMode(bool pbyp) { qmc_driver_mode_[QMC_UPDATE_MODE] = pbyp; }

  void putTraces(xmlNodePtr txml) {}
  void requestTraces(bool allow_traces) {}
  /** }@ */

protected:
  /** pure function returning AdjustedWalkerCounts data structure 
   *
   *  The logic is now walker counts is fairly simple.
   *  TotalWalkers trumps all other walker parameters
   *  If TotalWalkers is absent walkers_per_rank is used.
   *  if they are both absent then the default is one walker per crowd,
   *  each rank has crowds walkers.
   *  if crowds aren't specified you get one per main level thread.
   *
   *  You can have crowds or ranks with no walkers.
   *  You cannot have more crowds than threads.
   *
   *  passing num_ranks instead of internally querying comm->size()
   *  makes unit testing much quicker.
   *
   */
  static QMCDriverNew::AdjustedWalkerCounts adjustGlobalWalkerCount(int num_ranks,
                                                                    int rank_id,
                                                                    IndexType desired_count,
                                                                    IndexType walkers_per_rank,
                                                                    RealType reserve_walkers,
                                                                    int num_crowds);

  static void checkNumCrowdsLTNumThreads(const int num_crowds);

  /** The timers for the driver.
   *
   * This cleans up the driver constructor, and a reference to this structure 
   * Takes the timers into thread scope. We assume the timers are threadsafe.
   */
  struct DriverTimers
  {
    NewTimer& checkpoint_timer;
    NewTimer& run_steps_timer;
    NewTimer& init_walkers_timer;
    NewTimer& buffer_timer;
    NewTimer& movepbyp_timer;
    NewTimer& hamiltonian_timer;
    NewTimer& collectables_timer;
    DriverTimers(const std::string& prefix)
        : checkpoint_timer(*timer_manager.createTimer(prefix + "CheckPoint", timer_level_medium)),
          run_steps_timer(*timer_manager.createTimer(prefix + "RunSteps", timer_level_medium)),
          init_walkers_timer(*timer_manager.createTimer(prefix + "InitWalkers", timer_level_medium)),
          buffer_timer(*timer_manager.createTimer(prefix + "Buffer", timer_level_medium)),
          movepbyp_timer(*timer_manager.createTimer(prefix + "MovePbyP", timer_level_medium)),
          hamiltonian_timer(*timer_manager.createTimer(prefix + "Hamiltonian", timer_level_medium)),
          collectables_timer(*timer_manager.createTimer(prefix + "Collectables", timer_level_medium))
    {}
  };

  QMCDriverInput qmcdriver_input_;

  /** @ingroup Driver mutable input values
   *
   *  they should be limited to values that can be changed from input
   *  or are live state.
   *  @{
   */
  RealType max_disp_sq_;
  ///the number of saved samples
  IndexType target_samples_;

  /// the number of blocks between recomptePsi
  IndexType nBlocksBetweenRecompute;

  /**}@*/

  std::vector<std::unique_ptr<Crowd>> crowds_;


  std::string h5_file_root_;

  ///branch engine
  SFNBranch* branch_engine_;
  ///drift modifer
  std::unique_ptr<DriftModifierBase> drift_modifier_;

  ///the number to delay updates by
  int k_delay;

  /** period of recording walker configurations
   *
   * Default is 0 indicating that only the last configuration will be saved.
   */
  int walker_dump_period;


  IndexType current_step_;

  ///counter for number of moves accepted
  IndexType nAccept;

  ///counter for number of moves /rejected
  IndexType nReject;


  ///maximum cpu in secs
  RealType MaxCPUSecs;

  ///Time-step factor \f$ 1/(2\tau)\f$
  RealType m_oneover2tau;
  ///Time-step factor \f$ \sqrt{\tau}\f$
  RealType m_sqrttau;

  ///type of qmc: assigned by subclasses
  const std::string QMCType;
  ///root of all the output files
  std::string root_name_;


  ///the entire (or on node) walker population
  MCPopulation& population_;

  ///trial function
  TrialWaveFunction& Psi;

  ///Hamiltonian
  QMCHamiltonian& H;

  WaveFunctionPool& psiPool;

  /** Observables manager
   *  Has very problematic owner ship and life cycle.
   *  Can be transfered via branch manager one driver to the next indefinitely
   *  TODO:  Modify Branch manager and others to clear this up.
   */
  EstimatorManagerNew* estimator_manager_;

  ///record engine for walkers
  HDFWalkerOutput* wOut;

  /** Per crowd move contexts, this is where the DistanceTables etc. reside
   */
  std::vector<std::unique_ptr<ContextForSteps>> step_contexts_;

  ///a list of TrialWaveFunctions for multiple method
  std::vector<TrialWaveFunction*> Psi1;

  ///a list of QMCHamiltonians for multiple method
  std::vector<QMCHamiltonian*> H1;

  ///Random number generators
  std::vector<std::unique_ptr<RandomGenerator_t>> Rng;

  ///a list of mcwalkerset element
  std::vector<xmlNodePtr> mcwalkerNodePtr;

  ///temporary storage for drift
  ParticleSet::ParticlePos_t drift;

  ///temporary storage for random displacement
  ParticleSet::ParticlePos_t deltaR;

  // ///alternate method of setting QMC run parameters
  // IndexType nStepsBetweenSamples;
  // ///samples per thread
  // IndexType nSamplesPerThread;

  //  TODO: restart
  //  /** period of dumping walker configurations and everything else for restart
  //  *
  //  * The unit is a block.
  //  */
  // int check_point_period_;

  /** }@ */


  DriverTimers timers_;

public:
  ///Copy Constructor (disabled).
  QMCDriverNew(const QMCDriverNew&) = delete;
  ///Copy operator (disabled).
  QMCDriverNew& operator=(const QMCDriverNew&) = delete;

  bool putQMCInfo(xmlNodePtr cur);

  /** Adjust populations local walkers to this number
  * @param nwalkers number of walkers to add
  *
  */
  void makeLocalWalkers(int nwalkers,
                        RealType reserve,
                        const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions);

  DriftModifierBase& get_drift_modifier() const { return *drift_modifier_; }

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
  const std::string& get_root_name() const { return root_name_; }
  std::string getRotationName(std::string RootName);
  std::string getLastRotationName(std::string RootName);

private:
  friend std::ostream& operator<<(std::ostream& o_stream, const QMCDriverNew& qmcd);

  SetNonLocalMoveHandler setNonLocalMoveHandler_;

  static void defaultSetNonLocalMoveHandler(QMCHamiltonian& gold_ham);

  friend class qmcplusplus::testing::VMCBatchedTest;
  friend class qmcplusplus::testing::DMCBatchedTest;
  friend class qmcplusplus::testing::QMCDriverNewTestWrapper;
};
} // namespace qmcplusplus

#endif

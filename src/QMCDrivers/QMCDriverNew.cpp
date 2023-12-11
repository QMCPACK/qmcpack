//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <typeinfo>
#include <cmath>
#include <sstream>
#include <numeric>

#include "QMCDriverNew.h"
#include "Concurrency/ParallelExecutor.hpp"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/FairDivide.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "RandomNumberControl.h"
#include "Estimators/EstimatorManagerNew.h"
#include "hdf/HDFVersion.h"
#include "Utilities/qmc_common.h"
#include "Concurrency/Info.hpp"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/Timer.h"
#include "Message/UniformCommunicateError.h"
#include "EstimatorInputDelegates.h"


namespace qmcplusplus
{
/** Has nasty workaround for RandomNumberControl
 *   
 *  Num crowds must be less than omp_get_max_threads because RandomNumberControl is global c lib function
 *  masquerading as a C++ object.
 */
QMCDriverNew::QMCDriverNew(const ProjectData& project_data,
                           QMCDriverInput&& input,
                           const std::optional<EstimatorManagerInput>& global_emi,
                           WalkerConfigurations& wc,
                           MCPopulation&& population,
                           const std::string timer_prefix,
                           Communicate* comm,
                           const std::string& QMC_driver_type,
                           SetNonLocalMoveHandler snlm_handler)
    : MPIObjectBase(comm),
      qmcdriver_input_(std::move(input)),
      QMCType(QMC_driver_type),
      population_(std::move(population)),
      dispatchers_(!qmcdriver_input_.areWalkersSerialized()),
      estimator_manager_(nullptr),
      timers_(timer_prefix),
      driver_scope_timer_(createGlobalTimer(QMC_driver_type, timer_level_coarse)),
      driver_scope_profiler_(qmcdriver_input_.get_scoped_profiling()),
      project_data_(project_data),
      walker_configs_ref_(wc),
      setNonLocalMoveHandler_(snlm_handler)
{
  // This is done so that the application level input structures reflect the actual input to the code.
  // While the actual simulation objects still take singular input structures at construction.
  auto makeEstimatorManagerInput = [](auto& global_emi, auto& local_emi) -> EstimatorManagerInput {
    if (global_emi.has_value() && local_emi.has_value())
      return {global_emi.value(), local_emi.value()};
    else if (global_emi.has_value())
      return {global_emi.value()};
    else if (local_emi.has_value())
      return {local_emi.value()};
    else
      return {};
  };

  estimator_manager_ =
      std::make_unique<EstimatorManagerNew>(comm,
                                            makeEstimatorManagerInput(global_emi,
                                                                      qmcdriver_input_.get_estimator_manager_input()),
                                            population_.get_golden_hamiltonian(), population.get_golden_electrons(),
                                            population.get_golden_twf());

  drift_modifier_.reset(
      createDriftModifier(qmcdriver_input_.get_drift_modifier(), qmcdriver_input_.get_drift_modifier_unr_a()));

  // This needs to be done here to keep dependency on CrystalLattice out of the QMCDriverInput.
  max_disp_sq_ = input.get_max_disp_sq();
  if (max_disp_sq_ < 0)
  {
    auto& lattice = population.get_golden_electrons().getLattice();
    max_disp_sq_  = lattice.LR_rc * lattice.LR_rc;
  }

  wOut = std::make_unique<HDFWalkerOutput>(population.get_golden_electrons().getTotalNum(), get_root_name(), myComm);
}

// The Rng pointers are transferred from global storage (RandomNumberControl::Children)
// to local storage (Rng) for the duration of QMCDriverNew.
// They are transferred to local storage in createRngsStepContext (called from startup,
// which is usually called from the "process" function in the derived class.)
// The local storage is moved back to the global storage in the destructor.
// In optimization, there are two instances of QMCDriverNew - one for the optimizer and one
// for the vmc engine.   As long as the vmc engine calls process first, it gets valid
// Rng pointers.  The optimizer is called second and gets nullptr, but it doesn't use Rng,
// so it doesn't matter.
// Upon restore, the vmc engine would need to be restored last (otherwise the global storage gets
// the nullptr from the optimizer).  However, the order is fixed by the order the destructors
// are called.
// To work around the issue, check the local pointer for nullptr before restoring to global storage.
QMCDriverNew::~QMCDriverNew()
{
  for (int i = 0; i < Rng.size(); ++i)
    if (Rng[i])
      RandomNumberControl::Children[i].reset(Rng[i].release());
}

void QMCDriverNew::checkNumCrowdsLTNumThreads(const int num_crowds)
{
  int num_threads(Concurrency::maxCapacity<>());
  if (num_crowds > num_threads)
  {
    std::stringstream error_msg;
    error_msg << "Bad Input: num_crowds (" << num_crowds << ") > num_threads (" << num_threads << ")\n";
    throw UniformCommunicateError(error_msg.str());
  }
}

void QMCDriverNew::initializeQMC(const AdjustedWalkerCounts& awc)
{
  ScopedTimer local_timer(timers_.startup_timer);

  app_summary() << QMCType << " Driver running with" << std::endl
                << "             total_walkers     = " << awc.global_walkers << std::endl
                << "             walkers_per_rank  = " << awc.walkers_per_rank << std::endl
                << "             num_crowds        = " << awc.walkers_per_crowd.size() << std::endl
                << "  on rank 0, walkers_per_crowd = " << awc.walkers_per_crowd << std::endl
                << std::endl;

  // set num_global_walkers explicitly and then make local walkers.
  population_.set_num_global_walkers(awc.global_walkers);

  if (qmcdriver_input_.areWalkersSerialized())
  {
    if (estimator_manager_->areThereListeners())
      throw UniformCommunicateError("Serialized walkers ignore multiwalker API's and multiwalker resources and are "
                                    "incompatible with estimators requiring per particle listeners");
  }
  else
  {
    // This needs to happen before walkers are made. i.e. this allows hamiltonian operators to update state
    // the based on the presence of per particle listeners. In the case immediately encountered the operator CoulombPBCAA will
    // call its associated particle set and turnOnPerParticleSK.
    // The design for "initialization" of walker elements is for the golden elements to go through all pre walking state changes
    // and then for the golden elements to be cloned for each walker.
    if (estimator_manager_->areThereListeners())
      population_.get_golden_hamiltonian().informOperatorsOfListener();

    app_debug() << "Creating multi walker shared resources" << std::endl;
    population_.get_golden_electrons().createResource(golden_resource_.pset_res);
    population_.get_golden_twf().createResource(golden_resource_.twf_res);
    population_.get_golden_hamiltonian().createResource(golden_resource_.ham_res);
    app_debug() << "Multi walker shared resources creation completed" << std::endl;
  }

  makeLocalWalkers(awc.walkers_per_rank[myComm->rank()], awc.reserve_walkers);

  crowds_.resize(awc.walkers_per_crowd.size());

  // at this point we can finally construct the Crowd objects.
  for (int i = 0; i < crowds_.size(); ++i)
  {
    crowds_[i] =
        std::make_unique<Crowd>(*estimator_manager_, golden_resource_, population_.get_golden_electrons(),
                                population_.get_golden_twf(), population_.get_golden_hamiltonian(), dispatchers_);
  }

  //now give walkers references to their walkers
  population_.redistributeWalkers(crowds_);

  // Once they are created move contexts can be created.
  createRngsStepContexts(crowds_.size());

  if (qmcdriver_input_.get_measure_imbalance())
    measureImbalance("Startup");
}

/** QMCDriverNew ignores h5name if you want to read and h5 config you have to explicitly
 *  do so.
 */
void QMCDriverNew::setStatus(const std::string& aname, const std::string& h5name, bool append)
{
  app_log() << "\n========================================================="
            << "\n  Start " << QMCType << "\n  File Root " << get_root_name();
  app_log() << "\n=========================================================" << std::endl;

  if (h5name.size())
    h5_file_root_ = h5name;
}

/** Read walker configurations from *.config.h5 files
 * @param wset list of xml elements containing mcwalkerset
 *
 * All this does is look in the walker xml section for the hdf file.
 * It reads that (I think) and if there are active walkers 
 * declares it a restart run.
 *
 * This inferred behavior is asking for trouble.
 * Unified driver will not support until restart feature is
 * re-architected
 */
void QMCDriverNew::putWalkers(std::vector<xmlNodePtr>& wset)
{
  if (wset.empty())
    return;
  const int nfile = wset.size();

  HDFWalkerInputManager W_in(walker_configs_ref_, population_.get_golden_electrons().getTotalNum(), myComm);
  for (int i = 0; i < wset.size(); i++)
    if (W_in.put(wset[i]))
      h5_file_root_ = W_in.getFileRoot();
  //clear the walker set
  wset.clear();
  int nwtot = walker_configs_ref_.getActiveWalkers();
  myComm->bcast(nwtot);
  if (nwtot)
    setWalkerOffsets(walker_configs_ref_, myComm);
}

void QMCDriverNew::recordBlock(int block)
{
  if (qmcdriver_input_.get_dump_config() && block % qmcdriver_input_.get_check_point_period().period == 0)
  {
    ScopedTimer local_timer(timers_.checkpoint_timer);
    population_.saveWalkerConfigurations(walker_configs_ref_);
    setWalkerOffsets(walker_configs_ref_, myComm);
    wOut->dump(walker_configs_ref_, block);
    RandomNumberControl::write(getRngRefs(), get_root_name(), myComm);
  }
}

bool QMCDriverNew::finalize(int block, bool dumpwalkers)
{
  population_.saveWalkerConfigurations(walker_configs_ref_);
  setWalkerOffsets(walker_configs_ref_, myComm);
  app_log() << "  Carry over " << walker_configs_ref_.getGlobalNumWalkers()
            << " walker configurations to the next QMC driver." << std::endl;

  const bool DumpConfig = qmcdriver_input_.get_dump_config();
  if (DumpConfig && dumpwalkers)
    wOut->dump(walker_configs_ref_, block);

  infoSummary.flush();
  infoLog.flush();

  if (DumpConfig)
    RandomNumberControl::write(getRngRefs(), get_root_name(), myComm);

  return true;
}

void QMCDriverNew::makeLocalWalkers(IndexType nwalkers, RealType reserve)
{
  ScopedTimer local_timer(timers_.create_walkers_timer);
  // ensure nwalkers local walkers in population_
  if (population_.get_walkers().size() == 0)
    population_.createWalkers(nwalkers, walker_configs_ref_, reserve);
  else if (population_.get_walkers().size() < nwalkers)
  {
    throw std::runtime_error("Unexpected walker count resulting in dangerous spawning");
    IndexType num_additional_walkers = nwalkers - population_.get_walkers().size();
    for (int i = 0; i < num_additional_walkers; ++i)
      population_.spawnWalker();
  }
  else
  {
    IndexType num_walkers_to_kill = population_.get_walkers().size() - nwalkers;
    for (int i = 0; i < num_walkers_to_kill; ++i)
      population_.killLastWalker();
  }

  // \todo: this could be what is breaking spawned walkers
  for (UPtr<QMCHamiltonian>& ham : population_.get_hamiltonians())
    setNonLocalMoveHandler_(*ham);

  // For the dead ones too. Since this should be on construction but...
  for (UPtr<QMCHamiltonian>& ham : population_.get_dead_hamiltonians())
    setNonLocalMoveHandler_(*ham);
}

/** Creates Random Number generators for crowds and step contexts
 *
 *  This is quite dangerous in that number of crowds can be > omp_get_max_threads()
 *  This is used instead of actually passing number of threads/crowds
 *  controlling threads all over RandomNumberControl.
 */
void QMCDriverNew::createRngsStepContexts(int num_crowds)
{
  step_contexts_.resize(num_crowds);
  Rng.resize(num_crowds);

  if (RandomNumberControl::Children.size() == 0)
  {
    app_warning() << "  Initializing global RandomNumberControl! "
                  << "This message should not be seen in production code but only in unit tests." << std::endl;
    RandomNumberControl::make_seeds();
  }

  for (int i = 0; i < num_crowds; ++i)
  {
    Rng[i].reset(RandomNumberControl::Children[i].release());
    step_contexts_[i] = std::make_unique<ContextForSteps>(*(Rng[i]));
  }
}

void QMCDriverNew::initialLogEvaluation(int crowd_id,
                                        UPtrVector<Crowd>& crowds,
                                        UPtrVector<ContextForSteps>& context_for_steps)
{
  Crowd& crowd = *(crowds[crowd_id]);
  if (crowd.size() == 0)
    return;

  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());
  auto& ps_dispatcher  = crowd.dispatchers_.ps_dispatcher_;
  auto& twf_dispatcher = crowd.dispatchers_.twf_dispatcher_;
  auto& ham_dispatcher = crowd.dispatchers_.ham_dispatcher_;

  const RefVectorWithLeader<ParticleSet> walker_elecs(crowd.get_walker_elecs()[0], crowd.get_walker_elecs());
  const RefVectorWithLeader<TrialWaveFunction> walker_twfs(crowd.get_walker_twfs()[0], crowd.get_walker_twfs());
  const RefVectorWithLeader<QMCHamiltonian> walker_hamiltonians(crowd.get_walker_hamiltonians()[0],
                                                                crowd.get_walker_hamiltonians());

  ResourceCollectionTeamLock<ParticleSet> pset_res_lock(crowd.getSharedResource().pset_res, walker_elecs);
  ResourceCollectionTeamLock<TrialWaveFunction> twfs_res_lock(crowd.getSharedResource().twf_res, walker_twfs);
  ResourceCollectionTeamLock<QMCHamiltonian> hams_res_lock(crowd.getSharedResource().ham_res, walker_hamiltonians);

  auto& walkers = crowd.get_walkers();
  std::vector<bool> recompute_mask(walkers.size(), true);
  ps_dispatcher.flex_loadWalker(walker_elecs, walkers, recompute_mask, true);
  ps_dispatcher.flex_donePbyP(walker_elecs);
  twf_dispatcher.flex_evaluateLog(walker_twfs, walker_elecs);

  // For consistency this should be in ParticleSet as a flex call, but I think its a problem
  // in the algorithm logic and should be removed.
  auto saveElecPosAndGLToWalkers = [](ParticleSet& pset, ParticleSet::Walker_t& walker) { pset.saveWalker(walker); };
  for (int iw = 0; iw < crowd.size(); ++iw)
    saveElecPosAndGLToWalkers(walker_elecs[iw], walkers[iw]);

  std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
      ham_dispatcher.flex_evaluate(walker_hamiltonians, walker_twfs, walker_elecs));

  // \todo rename these are sets not resets.
  auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto local_energy) {
    walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    resetSigNLocalEnergy(walkers[iw], walker_twfs[iw], local_energies[iw]);

  auto evaluateNonPhysicalHamiltonianElements = [](QMCHamiltonian& ham, ParticleSet& pset, MCPWalker& walker) {
    ham.auxHevaluate(pset, walker);
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    evaluateNonPhysicalHamiltonianElements(walker_hamiltonians[iw], walker_elecs[iw], walkers[iw]);

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker) {
    ham.saveProperty(walker.getPropertyBase());
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    savePropertiesIntoWalker(walker_hamiltonians[iw], walkers[iw]);

  auto doesDoinTheseLastMatter = [](MCPWalker& walker) {
    walker.Weight     = 1.;
    walker.wasTouched = false;
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    doesDoinTheseLastMatter(walkers[iw]);
}

std::ostream& operator<<(std::ostream& o_stream, const QMCDriverNew& qmcd)
{
  o_stream << "  time step      = " << qmcd.qmcdriver_input_.get_tau() << '\n';
  o_stream << "  blocks         = " << qmcd.qmcdriver_input_.get_max_blocks() << '\n';
  o_stream << "  steps          = " << qmcd.qmcdriver_input_.get_max_steps() << '\n';
  o_stream << "  substeps       = " << qmcd.qmcdriver_input_.get_sub_steps() << '\n';
  o_stream << "  current        = " << qmcd.current_step_ << '\n';
  o_stream << "  target samples = " << qmcd.target_samples_ << '\n';
  o_stream << "  walkers/mpi    = " << qmcd.population_.get_num_local_walkers() << '\n' << '\n';
  o_stream << "  stepsbetweensamples = " << qmcd.qmcdriver_input_.get_steps_between_samples() << std::endl;
  app_log().flush();

  return o_stream;
}

void QMCDriverNew::defaultSetNonLocalMoveHandler(QMCHamiltonian& ham) {}

QMCDriverNew::AdjustedWalkerCounts QMCDriverNew::adjustGlobalWalkerCount(Communicate& comm,
                                                                         const IndexType current_configs,
                                                                         const IndexType requested_total_walkers,
                                                                         const IndexType requested_walkers_per_rank,
                                                                         const RealType reserve_walkers,
                                                                         int num_crowds)
{
  const int num_ranks = comm.size();
  const int rank_id   = comm.rank();

  // Step 1. set num_crowds by input and Concurrency::maxCapacity<>()
  checkNumCrowdsLTNumThreads(num_crowds);
  if (num_crowds == 0)
    num_crowds = Concurrency::maxCapacity<>();

  AdjustedWalkerCounts awc{0, {}, {}, reserve_walkers};
  awc.walkers_per_rank.resize(num_ranks, 0);

  // Step 2. decide awc.global_walkers and awc.walkers_per_rank based on input values
  if (requested_total_walkers != 0)
  {
    if (requested_total_walkers < num_ranks)
    {
      std::ostringstream error;
      error << "Running on " << num_ranks << " MPI ranks.  The request of " << requested_total_walkers
            << " global walkers cannot be satisfied! Need at least one walker per MPI rank.";
      throw UniformCommunicateError(error.str());
    }
    if (requested_walkers_per_rank != 0 && requested_total_walkers != requested_walkers_per_rank * num_ranks)
    {
      std::ostringstream error;
      error << "Running on " << num_ranks << " MPI ranks, The request of " << requested_total_walkers
            << " global walkers and " << requested_walkers_per_rank << " walkers per rank cannot be satisfied!";
      throw UniformCommunicateError(error.str());
    }
    awc.global_walkers   = requested_total_walkers;
    awc.walkers_per_rank = fairDivide(requested_total_walkers, num_ranks);
  }
  else // requested_total_walkers == 0
  {
    if (requested_walkers_per_rank != 0)
      awc.walkers_per_rank[rank_id] = requested_walkers_per_rank;
    else if (current_configs) // requested_walkers_per_rank == 0 and current_configs > 0
      awc.walkers_per_rank[rank_id] = current_configs;
    else // requested_walkers_per_rank == 0 and current_configs == 0
      awc.walkers_per_rank[rank_id] = num_crowds;
    comm.allreduce(awc.walkers_per_rank);
    awc.global_walkers = std::accumulate(awc.walkers_per_rank.begin(), awc.walkers_per_rank.end(), 0);
  }

  if (awc.global_walkers % num_ranks)
    app_warning() << "TotalWalkers (" << awc.global_walkers << ") not divisible by number of ranks (" << num_ranks
                  << "). This will result in a loss of efficiency.\n";

  // Step 3. decide awc.walkers_per_crowd
  awc.walkers_per_crowd = fairDivide(awc.walkers_per_rank[rank_id], num_crowds);

  if (awc.walkers_per_rank[rank_id] % num_crowds)
    app_warning() << "Walkers per rank (" << awc.walkers_per_rank[rank_id] << ") not divisible by number of crowds ("
                  << num_crowds << "). This will result in a loss of efficiency.\n";

  // \todo some warning if unreasonable number of threads are being used.

  return awc;
}

/** The scalar estimator collection is quite strange
 *
 */
void QMCDriverNew::endBlock()
{
  ScopedTimer local_timer(timers_.endblock_timer);
  RefVector<ScalarEstimatorBase> main_scalar_estimators;

  FullPrecRealType total_block_weight = 0.0;
  // Collect all the ScalarEstimatorsFrom EMCrowds
  unsigned long block_accept = 0;
  unsigned long block_reject = 0;

  std::vector<RefVector<OperatorEstBase>> crowd_operator_estimators;
  // Seems uneeded see EstimatorManagerNew scalar_ests_ documentation.
  std::vector<RefVector<ScalarEstimatorBase>> crowd_scalar_estimators;

  for (const UPtr<Crowd>& crowd : crowds_)
  {
    crowd->stopBlock();
    main_scalar_estimators.push_back(crowd->get_estimator_manager_crowd().get_main_estimator());
    crowd_scalar_estimators.emplace_back(crowd->get_estimator_manager_crowd().get_scalar_estimators());
    total_block_weight += crowd->get_estimator_manager_crowd().get_block_weight();
    block_accept += crowd->get_accept();
    block_reject += crowd->get_reject();

    // This seems altogether easier and more sane.
    crowd_operator_estimators.emplace_back(crowd->get_estimator_manager_crowd().get_operator_estimators());
  }

#ifdef DEBUG_PER_STEP_ACCEPT_REJECT
  app_warning() << "accept: " << block_accept << "   reject: " << block_reject;
  FullPrecRealType total_accept_ratio =
      static_cast<FullPrecRealType>(block_accept) / static_cast<FullPrecRealType>(block_accept + block_reject);
  std::cerr << "   total_accept_ratio: << " << total_accept_ratio << '\n';
#endif
  estimator_manager_->collectMainEstimators(main_scalar_estimators);
  estimator_manager_->collectScalarEstimators(crowd_scalar_estimators);
  estimator_manager_->collectOperatorEstimators(crowd_operator_estimators);

  /// get the average cpu_block time per crowd
  /// cpu_block_time /= crowds_.size();

  estimator_manager_->stopBlock(block_accept, block_reject, total_block_weight);
}

void QMCDriverNew::checkLogAndGL(Crowd& crowd, const std::string_view location)
{
  bool success         = true;
  auto& ps_dispatcher  = crowd.dispatchers_.ps_dispatcher_;
  auto& twf_dispatcher = crowd.dispatchers_.twf_dispatcher_;

  const RefVectorWithLeader<ParticleSet> walker_elecs(crowd.get_walker_elecs()[0], crowd.get_walker_elecs());
  const RefVectorWithLeader<TrialWaveFunction> walker_twfs(crowd.get_walker_twfs()[0], crowd.get_walker_twfs());
  std::vector<TrialWaveFunction::LogValue> log_values(walker_twfs.size());
  std::vector<ParticleSet::ParticleGradient> Gs;
  std::vector<ParticleSet::ParticleLaplacian> Ls;
  Gs.reserve(log_values.size());
  Ls.reserve(log_values.size());

  for (int iw = 0; iw < log_values.size(); iw++)
  {
    log_values[iw] = {walker_twfs[iw].getLogPsi(), walker_twfs[iw].getPhase()};
    Gs.push_back(walker_twfs[iw].G);
    Ls.push_back(walker_twfs[iw].L);
  }

  ps_dispatcher.flex_update(walker_elecs);
  twf_dispatcher.flex_evaluateLog(walker_twfs, walker_elecs);

  RealType threshold;
  // mixed precision can't make this test with cuda direct inversion
  if constexpr (std::is_same<RealType, FullPrecRealType>::value)
    threshold = 100 * std::numeric_limits<float>::epsilon();
  else
    threshold = 500 * std::numeric_limits<float>::epsilon();

  std::ostringstream msg;
  for (int iw = 0; iw < log_values.size(); iw++)
  {
    auto& ref_G = walker_twfs[iw].G;
    auto& ref_L = walker_twfs[iw].L;
    TrialWaveFunction::LogValue ref_log{walker_twfs[iw].getLogPsi(), walker_twfs[iw].getPhase()};
    if (std::abs(std::exp(log_values[iw]) - std::exp(ref_log)) > std::abs(std::exp(ref_log)) * threshold)
    {
      success = false;
      msg << "Logpsi walker[" << iw << "] " << log_values[iw] << " ref " << ref_log << std::endl;
    }

    for (int iel = 0; iel < ref_G.size(); iel++)
    {
      auto grad_diff = ref_G[iel] - Gs[iw][iel];
      if (std::sqrt(std::abs(dot(grad_diff, grad_diff))) > std::sqrt(std::abs(dot(ref_G[iel], ref_G[iel]))) * threshold)
      {
        success = false;
        msg << "walker[" << iw << "] Grad[" << iel << "] ref = " << ref_G[iel] << " wrong = " << Gs[iw][iel]
            << " Delta " << grad_diff << std::endl;
      }

      auto lap_diff = ref_L[iel] - Ls[iw][iel];
      if (std::abs(lap_diff) > std::abs(ref_L[iel]) * threshold)
      {
        // very hard to check mixed precision case, only print, no error out
        if (std::is_same<RealType, FullPrecRealType>::value)
          success = false;
        msg << "walker[" << iw << "] lap[" << iel << "] ref = " << ref_L[iel] << " wrong = " << Ls[iw][iel] << " Delta "
            << lap_diff << std::endl;
      }
    }
  }

  std::cerr << msg.str();
  if (!success)
    throw std::runtime_error(std::string("checkLogAndGL failed at ") + std::string(location) + std::string("\n"));
}

void QMCDriverNew::measureImbalance(const std::string& tag) const
{
  ScopedTimer local_timer(timers_.imbalance_timer);
  Timer only_this_barrier;
  myComm->barrier();
  std::vector<double> my_barrier_time(1, only_this_barrier.elapsed());
  std::vector<double> barrier_time_all_ranks(myComm->size(), 0.0);
  myComm->gather(my_barrier_time, barrier_time_all_ranks, 0);
  if (!myComm->rank())
  {
    auto const count  = static_cast<double>(barrier_time_all_ranks.size());
    const auto max_it = std::max_element(barrier_time_all_ranks.begin(), barrier_time_all_ranks.end());
    const auto min_it = std::min_element(barrier_time_all_ranks.begin(), barrier_time_all_ranks.end());
    app_log() << std::endl
              << tag << " MPI imbalance measured by an additional barrier (slow ranks wait less):" << std::endl
              << "    average wait seconds = "
              << std::accumulate(barrier_time_all_ranks.begin(), barrier_time_all_ranks.end(), 0.0) / count << std::endl
              << "    min wait at rank " << std::distance(barrier_time_all_ranks.begin(), min_it)
              << ", seconds = " << *min_it << std::endl
              << "    max wait at rank " << std::distance(barrier_time_all_ranks.begin(), max_it)
              << ", seconds = " << *max_it << std::endl;
  }
}

void QMCDriverNew::setWalkerOffsets(WalkerConfigurations& walker_configs, Communicate* comm)
{
  std::vector<int> nw(comm->size(), 0);
  std::vector<int> nwoff(comm->size() + 1, 0);
  nw[comm->rank()] = walker_configs.getActiveWalkers();
  comm->allreduce(nw);
  for (int ip = 0; ip < comm->size(); ip++)
    nwoff[ip + 1] = nwoff[ip] + nw[ip];

  walker_configs.setWalkerOffsets(nwoff);
}

} // namespace qmcplusplus

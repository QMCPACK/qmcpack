//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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

#include "QMCDrivers/QMCDriverNew.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/FairDivide.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Estimators/EstimatorManagerNew.h"
#include "HDFVersion.h"
#include "qmc_common.h"
#include "Concurrency/Info.hpp"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"
#include "Utilities/StlPrettyPrint.hpp"

namespace qmcplusplus
{
/** Has nasty workaround for RandomNumberControl
 *   
 *  Num crowds must be less than omp_get_max_threads because RandomNumberControl is global c lib function
 *  masquerading as a C++ object.
 */
QMCDriverNew::QMCDriverNew(QMCDriverInput&& input,
                           MCPopulation& population,
                           TrialWaveFunction& psi,
                           QMCHamiltonian& h,
                           WaveFunctionPool& ppool,
                           const std::string timer_prefix,
                           Communicate* comm,
                           SetNonLocalMoveHandler snlm_handler)
    : MPIObjectBase(comm),
      qmcdriver_input_(input),
      branch_engine_(nullptr),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      estimator_manager_(nullptr),
      wOut(0),
      timers_(timer_prefix),
      setNonLocalMoveHandler_(snlm_handler)
{
  QMCType  = "invalid";
  rotation = 0;

  // This needs to be done here to keep dependency on CrystalLattice out of the QMCDriverInput.
  max_disp_sq_ = input.get_max_disp_sq();
  if (max_disp_sq_ < 0)
  {
    const CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>& lattice = population.get_golden_electrons()->Lattice;
    max_disp_sq_                                              = lattice.LR_rc * lattice.LR_rc;
  }
}

int QMCDriverNew::addObservable(const std::string& aname)
{
  if (estimator_manager_)
    return estimator_manager_->addObservable(aname.c_str());
  else
    return -1;
}

QMCDriverNew::RealType QMCDriverNew::getObservable(int i) { return estimator_manager_->getObservable(i); }


QMCDriverNew::~QMCDriverNew()
{
  for (int i = 0; i < Rng.size(); ++i)
    RandomNumberControl::Children[i] = Rng[i].release();
}

void QMCDriverNew::add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi)
{
  H1.push_back(h);
  Psi1.push_back(psi);
}

void QMCDriverNew::checkNumCrowdsLTNumThreads(const int num_crowds)
{
  int num_threads(Concurrency::maxThreads<>());
  if (num_crowds > num_threads)
  {
    std::stringstream error_msg;
    error_msg << "Bad Input: num_crowds (" << num_crowds << ") > num_threads (" << num_threads << ")\n";
    throw std::runtime_error(error_msg.str());
  }
}

/** process a <qmc/> element
 * @param cur xmlNode with qmc tag
 *
 * This function is called before QMCDriverNew::run and following actions are taken:
 * - Initialize basic data to execute run function.
 * -- distance tables
 * -- resize deltaR and drift with the number of particles
 * -- assign cur to qmcNode
 * - process input file
 *   -- putQMCInfo: <parameter/> s for generic QMC
 *   -- put : extra data by derived classes
 * - initialize branchEngine to accumulate energies
 * - initialize Estimators
 * - initialize Walkers
 */
void QMCDriverNew::startup(xmlNodePtr cur, QMCDriverNew::AdjustedWalkerCounts awc)
{
  app_log() << this->QMCType << " Driver running with target_walkers =" << awc.global_walkers << std::endl
            << "                               walkers_per_rank =" << awc.walkers_per_rank << std::endl
            << "                               num_crowds =" << awc.walkers_per_crowd.size() << std::endl
            << "                    on rank 0, walkers_per_crowd =" << awc.walkers_per_crowd << std::endl
            << std::endl;

  makeLocalWalkers(awc.walkers_per_rank[myComm->rank()], awc.reserve_walkers,
                   ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>(population_.get_num_particles()));

  if (!branch_engine_)
  {
    branch_engine_ = new SFNBranch(qmcdriver_input_.get_tau(), population_.get_num_global_walkers());
  }

  //create and initialize estimator
  estimator_manager_ = branch_engine_->getEstimatorManager();
  if (!estimator_manager_)
  {
    estimator_manager_ = new EstimatorManagerNew(myComm);
    // TODO: remove this when branch engine no longer depends on estimator_mamanger_
    branch_engine_->setEstimatorManager(estimator_manager_);
    // This used to get updated as a side effect of setStatus
    branch_engine_->read(h5_file_root_);
  }
  else
    estimator_manager_->reset();

  if (!drift_modifier_)
    drift_modifier_.reset(createDriftModifier(qmcdriver_input_));

  // I don't think its at all good that the branch engine gets mutated here
  // Carrying the population on is one thing but a branch engine seems like it
  // should be fresh per section.
  branch_engine_->put(cur);
  estimator_manager_->put(H, cur);

  crowds_.resize(awc.walkers_per_crowd.size());

  // at this point we can finally construct the Crowd objects.
  for (int i = 0; i < crowds_.size(); ++i)
  {
    crowds_[i].reset(new Crowd(*estimator_manager_));
  }

  //now give walkers references to their walkers
  population_.distributeWalkers(crowds_);

  // Once they are created move contexts can be created.
  createRngsStepContexts(crowds_.size());

  // if (wOut == 0)
  //   wOut = new HDFWalkerOutput(W, root_name_, myComm);
  branch_engine_->start(root_name_);
  branch_engine_->write(root_name_);

  // PD: not really sure what the point of this is.  Seems to just go to output
  branch_engine_->advanceQMCCounter();
}

/** QMCDriverNew ignores h5name if you want to read and h5 config you have to explicitly
 *  do so.
 */
void QMCDriverNew::setStatus(const std::string& aname, const std::string& h5name, bool append)
{
  root_name_ = aname;
  app_log() << "\n========================================================="
            << "\n  Start " << QMCType << "\n  File Root " << root_name_;
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
  // if (wset.empty())
  //   return;
  // int nfile = wset.size();
  // HDFWalkerInputManager W_in(W, myComm);
  // for (int i = 0; i < wset.size(); i++)
  //   if (W_in.put(wset[i]))
  //     h5FileRoot = W_in.getFileRoot();
  // //clear the walker set
  // wset.clear();
  // int nwtot = W.getActiveWalkers();
  // myComm->bcast(nwtot);
  // if (nwtot)
  // {
  //   int np = myComm->size();
  //   std::vector<int> nw(np, 0), nwoff(np + 1, 0);
  //   nw[myComm->rank()] = W.getActiveWalkers();
  //   myComm->allreduce(nw);
  //   for (int ip = 0; ip < np; ++ip)
  //     nwoff[ip + 1] = nwoff[ip] + nw[ip];
  //   W.setGlobalNumWalkers(nwoff[np]);
  //   W.setWalkerOffsets(nwoff);
  //   qmc_common.is_restart = true;
  // }
  // else
  //   qmc_common.is_restart = false;
}

std::string QMCDriverNew::getRotationName(std::string root_name)
{
  std::string r_RootName;
  if (rotation % 2 == 0)
  {
    r_RootName = root_name;
  }
  else
  {
    r_RootName = root_name + ".bk";
  }
  rotation++;
  return r_RootName;
}

std::string QMCDriverNew::getLastRotationName(std::string root_name)
{
  std::string r_RootName;
  if ((rotation - 1) % 2 == 0)
  {
    r_RootName = root_name;
  }
  else
  {
    r_RootName = root_name + ".bk";
  }
  return r_RootName;
}

void QMCDriverNew::recordBlock(int block)
{
  if (qmcdriver_input_.get_dump_config() && block % qmcdriver_input_.get_check_point_period().period == 0)
  {
    timers_.checkpoint_timer.start();
    branch_engine_->write(root_name_, true); //save energy_history
    RandomNumberControl::write(root_name_, myComm);
    timers_.checkpoint_timer.stop();
  }
}

bool QMCDriverNew::finalize(int block, bool dumpwalkers)
{
  RefVector<MCPWalker> walkers(convertUPtrToRefVector(population_.get_walkers()));
  branch_engine_->finalize(population_.get_num_global_walkers(), walkers);

  if (qmcdriver_input_.get_dump_config())
    RandomNumberControl::write(root_name_, myComm);

  return true;
}

void QMCDriverNew::makeLocalWalkers(IndexType nwalkers,
                                    RealType reserve,
                                    const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  if (population_.get_walkers().size() == 0)
  {
    population_.createWalkers(nwalkers, reserve);
  }
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

  // setWalkerOffsets();
  // ////update the global number of walkers
  // ////int nw=W.getActiveWalkers();
  // ////myComm->allreduce(nw);
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

  TasksOneToOne<> do_per_crowd(num_crowds);

  Rng.resize(num_crowds);

  RngCompatibility.resize(num_crowds);

  if (RandomNumberControl::Children.size() == 0)
  {
    app_warning() << "  Initializing global RandomNumberControl! "
                  << "This message should not be seen in production code but only in unit tests." << std::endl;
    RandomNumberControl::make_seeds();
  }

  for (int i = 0; i < num_crowds; ++i)
  {
    Rng[i].reset(RandomNumberControl::Children[i]);
    // Ye: RandomNumberControl::Children needs to be replaced with unique_ptr and use Rng[i].swap()
    RandomNumberControl::Children[i] = nullptr;
    step_contexts_[i]   = std::make_unique<ContextForSteps>(crowds_[i]->size(), population_.get_num_particles(),
                                                          population_.get_particle_group_indexes(), *(Rng[i]));
    RngCompatibility[i] = Rng[i].get();
  }
}

void QMCDriverNew::initialLogEvaluation(int crowd_id,
                                        UPtrVector<Crowd>& crowds,
                                        UPtrVector<ContextForSteps>& context_for_steps)
{
  Crowd& crowd = *(crowds[crowd_id]);
  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());

  auto& walker_twfs         = crowd.get_walker_twfs();
  auto& mcp_buffers         = crowd.get_mcp_wfbuffers();
  auto& walker_elecs        = crowd.get_walker_elecs();
  auto& walkers             = crowd.get_walkers();
  auto& walker_hamiltonians = crowd.get_walker_hamiltonians();

  crowd.loadWalkers();
  for (ParticleSet& pset : walker_elecs)
    pset.update();

  // Added to match legacy.
  auto cleanDataSet = [](MCPWalker& walker, ParticleSet& pset, TrialWaveFunction& twf) {
    if (walker.DataSet.size())
      walker.DataSet.clear();
    walker.DataSet.rewind();
    walker.registerData();
    twf.registerData(pset, walker.DataSet);
    walker.DataSet.allocate();
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    cleanDataSet(walkers[iw], walker_elecs[iw], walker_twfs[iw]);

  auto copyFrom = [](TrialWaveFunction& twf, ParticleSet& pset, WFBuffer& wfb) { twf.copyFromBuffer(pset, wfb); };
  for (int iw = 0; iw < crowd.size(); ++iw)
    copyFrom(walker_twfs[iw], walker_elecs[iw], mcp_buffers[iw]);

  TrialWaveFunction::flex_evaluateLog(walker_twfs, walker_elecs);

  TrialWaveFunction::flex_updateBuffer(crowd.get_walker_twfs(), crowd.get_walker_elecs(), crowd.get_mcp_wfbuffers());

  // For consistency this should be in ParticleSet as a flex call, but I think its a problem
  // in the algorithm logic and should be removed.
  auto saveElecPosAndGLToWalkers = [](ParticleSet& pset, ParticleSet::Walker_t& walker) { pset.saveWalker(walker); };
  for (int iw = 0; iw < crowd.size(); ++iw)
    saveElecPosAndGLToWalkers(walker_elecs[iw], walkers[iw]);

  std::vector<QMCHamiltonian::FullPrecRealType> local_energies(
      QMCHamiltonian::flex_evaluate(walker_hamiltonians, walker_elecs));
  // This is actually only a partial reset of the walkers properties
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
    walker.ReleasedNodeAge    = 0;
    walker.ReleasedNodeWeight = 0;
    walker.Weight             = 1;
  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    doesDoinTheseLastMatter(walkers[iw]);
}

void QMCDriverNew::setWalkerOffsets()
{
  std::vector<int> nw(myComm->size(), 0), nwoff(myComm->size() + 1, 0);
  //  nw[myComm->rank()] = W.getActiveWalkers();
  myComm->allreduce(nw);
  for (int ip = 0; ip < myComm->size(); ip++)
    nwoff[ip + 1] = nwoff[ip] + nw[ip];
  //  W.setGlobalNumWalkers(nwoff[myComm->size()]);
  //  W.setWalkerOffsets(nwoff);
  long id = nwoff[myComm->rank()];
  for (int iw = 0; iw < nw[myComm->rank()]; ++iw, ++id)
  {
    //    W[iw]->ID       = id;
    //    W[iw]->ParentID = id;
  }
  //  app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples << std::endl;
  //  app_log() << "  Total weight: " << W.EnsembleProperty.Weight << std::endl;
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

QMCDriverNew::AdjustedWalkerCounts QMCDriverNew::adjustGlobalWalkerCount(int num_ranks,
                                                                         int rank_id,
                                                                         IndexType required_total,
                                                                         IndexType walkers_per_rank,
                                                                         RealType reserve_walkers,
                                                                         int num_crowds)
{
  // Step 1. set num_crowds by input and Concurrency::maxThreads<>()
  checkNumCrowdsLTNumThreads(num_crowds);
  if (num_crowds == 0)
    num_crowds = Concurrency::maxThreads<>();

  AdjustedWalkerCounts awc{0, {}, {}, reserve_walkers};

  // Step 2. decide awc.global_walkers and awc.walkers_per_rank based on input values
  if (required_total != 0)
  {
    if (required_total < num_ranks)
    {
      std::ostringstream error;
      error << "Running on " << num_ranks << " MPI ranks and the request of " << required_total
            << " global walkers cannot be satisfied! Need at least one walker per MPI rank.";
      throw std::runtime_error(error.str());
    }
    if (walkers_per_rank != 0 && required_total != walkers_per_rank * num_ranks)
    {
      std::ostringstream error;
      error << "Running on " << num_ranks << " MPI ranks and the request of " << required_total
            << " global walkers and " << walkers_per_rank << " walkers per rank cannot be satisfied!";
      throw std::runtime_error(error.str());
    }
    awc.global_walkers   = required_total;
    awc.walkers_per_rank = fairDivide(required_total, num_ranks);
  }
  else
  {
    if (walkers_per_rank != 0)
      awc.walkers_per_rank = std::vector<IndexType>(num_ranks, walkers_per_rank);
    else
      awc.walkers_per_rank = std::vector<IndexType>(num_ranks, num_crowds);
    awc.global_walkers = awc.walkers_per_rank[0] * num_ranks;
  }

  // Step 3. decide awc.walkers_per_crowd
  awc.walkers_per_crowd = fairDivide(awc.walkers_per_rank[rank_id], num_crowds);

  if (awc.global_walkers % num_ranks)
    app_warning() << "TotalWalkers (" << awc.global_walkers << ") not divisible by number of ranks (" << num_ranks
                  << "). This will result in a loss of efficiency.\n";

  if (awc.walkers_per_rank[rank_id] % num_crowds)
    app_warning() << "Walkers per rank (" << awc.walkers_per_rank[rank_id] << ") not divisible by number of crowds ("
                  << num_crowds << "). This will result in a loss of efficiency.\n";

  // \todo some warning if unreasonable number of threads are being used.

  return awc;
}

void QMCDriverNew::endBlock()
{
  RefVector<ScalarEstimatorBase> all_scalar_estimators;
  FullPrecRealType total_block_weight = 0.0;
  FullPrecRealType total_accept_ratio = 0.0;
  // Collect all the ScalarEstimatorsFrom EMCrowds
  double cpu_block_time = 0.0;
  for (const UPtr<Crowd>& crowd : crowds_)
  {
    crowd->stopBlock();
    auto crowd_sc_est = crowd->get_estimator_manager_crowd().get_scalar_estimators();
    all_scalar_estimators.insert(all_scalar_estimators.end(), std::make_move_iterator(crowd_sc_est.begin()),
                                 std::make_move_iterator(crowd_sc_est.end()));
    total_block_weight += crowd->get_estimator_manager_crowd().get_block_weight();
    total_accept_ratio += crowd->get_accept_ratio() * crowd->get_estimator_manager_crowd().get_block_weight();
    cpu_block_time += crowd->get_estimator_manager_crowd().get_cpu_block_time();
  }
  total_accept_ratio /= total_block_weight;
  estimator_manager_->collectScalarEstimators(all_scalar_estimators);
  cpu_block_time /= crowds_.size();

  estimator_manager_->stopBlockNew(total_accept_ratio, total_block_weight, cpu_block_time);
}

} // namespace qmcplusplus

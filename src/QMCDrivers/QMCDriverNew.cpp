//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: QMCDriver.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <typeinfo>
#include <cmath>

#include "QMCDrivers/QMCDriverNew.h"
#include "Concurrency/TasksOneToOne.hpp"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "HDFVersion.h"
#include "qmc_common.h"
#include "Concurrency/Info.hpp"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"

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
      walkers_per_crowd_(1),
      branch_engine_(nullptr),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      estimator_manager_(nullptr),
      wOut(0),
      timers_(timer_prefix),
      setNonLocalMoveHandler_(snlm_handler)
      // num_crowds_(input.get_num_crowds())
{
  QMCType = "invalid";

  // Avoids segmentation fault when RandomNumberControl::Children is too small, adds surprising behavior
  if(Concurrency::maxThreads() < input.get_num_crowds())
    set_num_crowds(Concurrency::maxThreads(), "RandomNumberControl's maximum children set to omp_get_max_threads()");
  else
    num_crowds_ = input.get_num_crowds();
 
  rotation = 0;

  // This needs to be done here to keep dependency on CrystalLattice out of the QMCDriverInput.
  max_disp_sq_ = input.get_max_disp_sq();
  if(max_disp_sq_ < 0)
  {
    const CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>& lattice = population.get_golden_electrons()->Lattice;
    max_disp_sq_ = lattice.LR_rc * lattice.LR_rc;
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
  for(int i = 0; i < Rng.size(); ++i)
    RandomNumberControl::Children[i] = Rng[i].release();
}

void QMCDriverNew::add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi)
{
  H1.push_back(h);
  Psi1.push_back(psi);
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
void QMCDriverNew::process(xmlNodePtr cur)
{
  // If you really want to persist the MCPopulation it is not the business of QMCDriver to reset it.
  // It could tell it we are starting a new section but shouldn't be pulling internal strings.
  //int numCopies = (H1.empty()) ? 1 : H1.size();
  //W.resetWalkerProperty(numCopies);

  if (!branch_engine_)
  {
    branch_engine_ = new SimpleFixedNodeBranch(qmcdriver_input_.get_tau(), population_.get_num_global_walkers());
  }

  //create and initialize estimator
  estimator_manager_ = branch_engine_->getEstimatorManager();
  if (!estimator_manager_)
  {
    estimator_manager_ = new EstimatorManagerBase(myComm);
    // TODO: remove this when branch engine no longer depends on estimator_mamanger_
    branch_engine_->setEstimatorManager(estimator_manager_);
    // This used to get updated as a side effect of setStatus
    branch_engine_->read(h5_file_root_);
  }

  if (!drift_modifier_)
    drift_modifier_.reset(createDriftModifier(qmcdriver_input_));

  // I don't think its at all good that the branch engine gets mutated here
  // Carrying the population on is one thing but a branch engine seems like it
  // should be fresh per section.
  branch_engine_->put(cur);
  estimator_manager_->put(H, cur);

  crowds_.resize(num_crowds_);
  // at this point we can finally construct the Crowd objects.
  // neglecting first touch for the moment
  // because em cloning is not threadsafe
  for(int i = 0; i < num_crowds_; ++i)
  {
    crowds_[i].reset(new Crowd(*estimator_manager_));
  }//  crowds_.push_back(

  //now give walkers references to their walkers
  auto crowd_start = crowds_.begin();
  auto crowd_end   = crowds_.end();
  std::for_each(crowd_start,
                crowd_end,
                [this](std::unique_ptr<Crowd>& crowd) { crowd->reserve(walkers_per_crowd_); });
  population_.distributeWalkers(crowd_start, crowd_end, walkers_per_crowd_);

  // Once they are created move contexts can be created.
  createRngsStepContexts();

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

void QMCDriverNew::checkNumCrowdsLTNumThreads() const
{
  int num_threads(Concurrency::maxThreads<>());
  if (num_crowds_ > num_threads)
  {
    std::stringstream error_msg;
    error_msg << "Bad Input: num_crowds (" << qmcdriver_input_.get_num_crowds()
       << ") > num_threads (" << num_threads << ")\n";
    throw std::runtime_error(error_msg.str());
  }
}

void QMCDriverNew::set_num_crowds(int num_crowds, const std::string& reason)
{
  num_crowds_ = num_crowds;
  app_warning() << " [INPUT OVERIDDEN] The number of crowds has been set to :  " << num_crowds << '\n';
  app_warning() << " Overiding the input of value of " << qmcdriver_input_.get_num_crowds() << " because " << reason
                << std::endl;
}

void QMCDriverNew::set_walkers_per_rank(int walkers_per_rank, const std::string& reason)
{
  walkers_per_rank_ = walkers_per_rank;
  app_warning() << " [INPUT OVERIDDEN] The number of crowds has been set to :  " << walkers_per_rank << '\n';
  app_warning() << " Overiding the input of value of " << qmcdriver_input_.get_walkers_per_rank() << " because " << reason
                << std::endl;
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

void QMCDriverNew::makeLocalWalkers(int nwalkers, const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  if (population_.get_walkers().size() == 0)
  {
    population_.createWalkers(nwalkers);
  }
  else if(population_.get_walkers().size() < nwalkers)
  {
    throw std::runtime_error("Unexpected walker count resulting in dangerous spawning");
    IndexType num_additional_walkers = nwalkers - population_.get_walkers().size();
    for(int i = 0; i < num_additional_walkers; ++i)
      population_.spawnWalker();
  }
  else
  {
    IndexType num_walkers_to_kill = population_.get_walkers().size() - nwalkers;
    for(int i = 0; i < num_walkers_to_kill; ++i)
      population_.killLastWalker();
  }

  // \todo: this could be what is breaking spawned walkers
  for(UPtr<QMCHamiltonian>& ham : population_.get_hamiltonians())
    setNonLocalMoveHandler_(*ham);

  // For the dead ones too. Since this should be on construction but...
  for(UPtr<QMCHamiltonian>& ham : population_.get_dead_hamiltonians())
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
void QMCDriverNew::createRngsStepContexts()
{
  step_contexts_.resize(num_crowds_);

  TasksOneToOne<> do_per_crowd(num_crowds_);

  Rng.resize(num_crowds_);

  if (RandomNumberControl::Children.size() == 0)
  {
    app_warning() << "  Initializing global RandomNumberControl! "
                  << "This message should not be seen in production code but only in unit tests." << std::endl;
    RandomNumberControl::make_seeds();
  }

  for(int i = 0; i < num_crowds_; ++i)
  {
    Rng[i].reset(RandomNumberControl::Children[i]);
    // Ye: RandomNumberControl::Children needs to be replaced with unique_ptr and use Rng[i].swap()
    RandomNumberControl::Children[i] = nullptr;
    step_contexts_[i] = std::make_unique<ContextForSteps>(crowds_[i]->size(), population_.get_num_particles(),
                                            population_.get_particle_group_indexes(), *(Rng[i]));
  }
}

void QMCDriverNew::initialLogEvaluation(int crowd_id, UPtrVector<Crowd>& crowds, UPtrVector<ContextForSteps>& context_for_steps)
{
  Crowd& crowd = *(crowds[crowd_id]);
  crowd.setRNGForHamiltonian(context_for_steps[crowd_id]->get_random_gen());

  auto& walker_twfs  = crowd.get_walker_twfs();
  auto& mcp_buffers  = crowd.get_mcp_wfbuffers();
  auto& walker_elecs = crowd.get_walker_elecs();
  auto& walkers = crowd.get_walkers();
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

  auto copyFrom = [](TrialWaveFunction& twf, ParticleSet& pset, WFBuffer& wfb){
                      twf.copyFromBuffer(pset,wfb);
                    };
  for (int iw = 0; iw < crowd.size(); ++iw)
    copyFrom(walker_twfs[iw], walker_elecs[iw], mcp_buffers[iw]);

  TrialWaveFunction::flex_evaluateLog(walker_twfs, walker_elecs);

  TrialWaveFunction::flex_updateBuffer(crowd.get_walker_twfs(),
                                       crowd.get_walker_elecs(),
                                       crowd.get_mcp_wfbuffers());

  // For consistency this should be in ParticleSet as a flex call, but I think its a problem
  // in the algorithm logic and should be removed.
  auto saveElecPosAndGLToWalkers = [](ParticleSet& pset, ParticleSet::Walker_t& walker){
                                     pset.saveWalker(walker);};
  for (int iw = 0; iw < crowd.size(); ++iw)
    saveElecPosAndGLToWalkers(walker_elecs[iw], walkers[iw]);

  std::vector<QMCHamiltonian::FullPrecRealType> local_energies(QMCHamiltonian::flex_evaluate(walker_hamiltonians, walker_elecs));
  // This is actually only a partial reset of the walkers properties
  auto resetSigNLocalEnergy = [](MCPWalker& walker, TrialWaveFunction& twf, auto local_energy){
                                walker.resetProperty(twf.getLogPsi(), twf.getPhase(), local_energy);
                                    };
  for (int iw = 0; iw < crowd.size(); ++iw)
    resetSigNLocalEnergy(walkers[iw], walker_twfs[iw], local_energies[iw]);

  auto evaluateNonPhysicalHamiltonianElements = [](QMCHamiltonian& ham, ParticleSet& pset, MCPWalker& walker){
                                                   ham.auxHevaluate(pset, walker);
                                                 };
  for (int iw = 0; iw < crowd.size(); ++iw)
    evaluateNonPhysicalHamiltonianElements(walker_hamiltonians[iw], walker_elecs[iw], walkers[iw]);

  auto savePropertiesIntoWalker = [](QMCHamiltonian& ham, MCPWalker& walker){
                                    ham.saveProperty(walker.getPropertyBase());
                                  };
  for (int iw = 0; iw < crowd.size(); ++iw)
    savePropertiesIntoWalker(walker_hamiltonians[iw], walkers[iw]);

  auto doesDoinTheseLastMatter = [](MCPWalker& walker){
                                   walker.ReleasedNodeAge    = 0;
    walker.ReleasedNodeWeight = 0;
    walker.Weight             = 1;};
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

void QMCDriverNew::defaultSetNonLocalMoveHandler(QMCHamiltonian& ham)
{}


/** static function
 *
 *  
 */
QMCDriverNew::AdjustedWalkerCounts QMCDriverNew::adjustGlobalWalkerCount(Communicate* comm, IndexType desired_count, IndexType walkers_per_rank, int num_crowds)
{
  int ranks = myComm->size();
  AdjustedWalkerCounts awc{0,0,0,0};;
  awc.global_walkers = desired_count;
  awc.walkers_per_rank = walkers_per_rank;
  if ( awc.global_walkers != 0)
  {
    if ( awc.global_walkers % ranks )
    {
      awc.walkers_per_rank = awc.global_walkers / ranks + 1;
      awc.global_walkers = awc.walkers_per_rank * ranks;
      app_warning() << "TargetWalkers not divisible by number of ranks, TargetWalker increased to "
         << awc.global_walkers << '\n';
    }
    else
      awc.walkers_per_rank = awc.global_walkers / ranks;
  }

  // side effect updates Base::walkers_per_crowd_;
  awc = this->calcDefaultLocalWalkers(awc);

  if ( awc.walkers_per_rank * ranks != awc.global_walkers )
  {
    awc.global_walkers = awc.walkers_per_rank * ranks;
    app_warning() << "Walkers per rank number of crowds, TargetWalker increased to "
         << awc.global_walkers << '\n';
  }

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
  // Note: that this is already averaged in crowds and then summed weighted by the walkers 
  // in each crowd (which can be different) and then as a result the average is over
  // local walkers
  total_accept_ratio /= total_block_weight; //population_.get_num_local_walkers();
  estimator_manager_->collectScalarEstimators(all_scalar_estimators);
  cpu_block_time /= crowds_.size();

  // TODO: should be accept rate for block
  estimator_manager_->stopBlockNew(total_accept_ratio, total_block_weight, cpu_block_time);
}

} // namespace qmcplusplus

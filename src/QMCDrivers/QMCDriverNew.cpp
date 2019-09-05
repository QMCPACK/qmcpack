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

/** \file
 * Clean base class for Unified driver
 *  
 * This driver base class should be generic with respect to precision,
 * value type, device execution, and ...
 * It should contain no typdefs not related to compiler bugs or platform workarounds
 *
 * Some integer math is done in non performance critical areas in the clear and
 * not optimized way. Leave these alone.
 */

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
                           MCPopulation&& population,
                           TrialWaveFunction& psi,
                           QMCHamiltonian& h,
                           WaveFunctionPool& ppool,
                           Communicate* comm)
    : MPIObjectBase(comm),
      qmcdriver_input_(input),
      branchEngine(nullptr),
      population_(std::move(population)),
      Psi(psi),
      H(h),
      psiPool(ppool),
      estimator_manager_(nullptr),
      wOut(0),
      walkers_per_crowd_(1)
      // num_crowds_(input.get_num_crowds())
{
  QMCType = "invalid";

  // Avoids segmentation fault when RandomNumberControl::Children is too small, adds surprising behavior
  if(Concurrency::maxThreads() < input.get_num_crowds())
    set_num_crowds(Concurrency::maxThreads(), "RandomNumberControl's maximum children set to omp_get_max_threads()");
  else
    num_crowds_ = input.get_num_crowds();
 
  rotation = 0;

  checkpointTimer = TimerManager.createTimer("checkpoint::recordBlock", timer_level_medium);
}

int QMCDriverNew::addObservable(const std::string& aname)
{
  if (estimator_manager_)
    return estimator_manager_->addObservable(aname.c_str());
  else
    return -1;
}

QMCDriverNew::RealType QMCDriverNew::getObservable(int i) { return estimator_manager_->getObservable(i); }


QMCDriverNew::~QMCDriverNew() {}

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
  // if (qmcdriver_input_.get_reset_random() || RandomNumberControl)
  // {

  // if seeds are not made then neither are the children. So when MoveContexts are created a segfault occurs.
  // For now it is unclear whether get_reset_random should always be true on the first run or what.
  app_log() << "  Regenerate random seeds." << std::endl;
  RandomNumberControl::make_seeds();
  // }
  
  setupWalkers();

  // If you really want to persist the MCPopulation it is not the business of QMCDriver to reset it.
  // It could tell it we are starting a new section but shouldn't be pulling internal strings.
  //int numCopies = (H1.empty()) ? 1 : H1.size();
  //W.resetWalkerProperty(numCopies);

  if (!branchEngine)
  {
    branchEngine = new SimpleFixedNodeBranch(qmcdriver_input_.get_tau(), population_.get_num_global_walkers());
  }

  //create and initialize estimator
  estimator_manager_ = branchEngine->getEstimatorManager();
  if (!estimator_manager_)
  {
    estimator_manager_ = new EstimatorManagerBase(myComm);
    branchEngine->setEstimatorManager(estimator_manager_);
    // This used to get updated as a side effect of setStatus
    branchEngine->read(h5_file_root_);
  }

  if (!drift_modifier_)
    drift_modifier_.reset(createDriftModifier(qmcdriver_input_));

  branchEngine->put(cur);
  estimator_manager_->put(H, cur);

  crowds_.resize(num_crowds_);
  // at this point we can finally construct the Crowd objects.
  // neglecting first touch for the moment
  // because em cloning is not threadsafe
  for(int i = 0; i < num_crowds_; ++i)
  {
    crowds_[i].reset(new Crowd(*estimator_manager_));
  }//  crowds_.push_back(
                        
  // Once they are created move contexts can be created.
  createRngsStepContexts();

  // if (wOut == 0)
  //   wOut = new HDFWalkerOutput(W, root_name_, myComm);
  branchEngine->start(root_name_);
  branchEngine->write(root_name_);

  // PD: not really sure what the point of this is.  Seems to just go to output
  branchEngine->advanceQMCCounter();
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

void QMCDriverNew::set_num_crowds(int num_crowds, const std::string& reason)
{
  num_crowds_ = num_crowds;
  app_warning() << " [INPUT OVERIDDEN] The number of crowds has been set to :  " << num_crowds << '\n';
  app_warning() << " Overiding the input of value of " << qmcdriver_input_.get_num_crowds() << " because " << reason
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
    checkpointTimer->start();
    branchEngine->write(root_name_, true); //save energy_history
    RandomNumberControl::write(root_name_, myComm);
    checkpointTimer->stop();
  }
}

bool QMCDriverNew::finalize(int block, bool dumpwalkers)
{
  //  branchEngine->finalize(W);

  if (qmcdriver_input_.get_dump_config())
    RandomNumberControl::write(root_name_, myComm);

  return true;
}


/** Elements of putQMCInfo that have nothing to do with input
 */
void QMCDriverNew::setupWalkers()
{
  IndexType local_walkers = calc_default_local_walkers();
  addWalkers(local_walkers, ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>(population_.get_num_particles()));
  //now give walkers references to their walkers
  auto crowd_start = crowds_.begin();
  auto crowd_end   = crowds_.end();
  population_.distributeWalkers(crowd_start, crowd_end, walkers_per_crowd_);
}

/** Add walkers to the end of the ensemble of walkers.
 * @param nwalkers number of walkers to add
 *
 */
void QMCDriverNew::addWalkers(int nwalkers, const ParticleAttrib<TinyVector<QMCTraits::RealType, 3>>& positions)
{
  population_.createWalkers(nwalkers, positions);
  // else if (nwalkers < 0)
  // {
  //   W.destroyWalkers(-nwalkers);
  //   app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers()
  //             << std::endl;
  // }
  // else
  // {
  //   app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." << std::endl;
  // }
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

  for(int i = 0; i < num_crowds_; ++i)
  {
    Rng[i].reset(new RandomGenerator_t(*(RandomNumberControl::Children[i])));
    step_contexts_[i].reset(new ContextForSteps(crowds_[i]->size(), population_.get_num_particles(),
                                            population_.get_particle_group_indexes(), *(Rng[i])));
  }
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

} // namespace qmcplusplus

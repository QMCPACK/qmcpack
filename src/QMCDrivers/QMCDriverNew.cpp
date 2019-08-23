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

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"

namespace qmcplusplus
{
QMCDriverNew::QMCDriverNew(QMCDriverInput&& input,
                           MCPopulation& population,
                           TrialWaveFunction& psi,
                           QMCHamiltonian& h,
                           WaveFunctionPool& ppool,
                           Communicate* comm)
    : MPIObjectBase(comm),
      qmcdriver_input_(input),
      branchEngine(nullptr),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      estimator_manager_(nullptr),
      wOut(0),
      walkers_per_crowd_(1),
      num_crowds_(input.get_num_crowds())
{
  QMCType = "invalid";

  ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
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


QMCDriverNew::~QMCDriverNew() { delete_iter(Rng.begin(), Rng.end()); }

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
  if (!qmcdriver_input_.get_append_run())
    current_step_ = 0;
  else
    current_step_ = qmcdriver_input_.get_starting_step();

  setupWalkers();
  // If you really wanter to persist the MCPopulation it is not the business of QMCDriver to reset it.
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
  // if (wOut == 0)
  //   wOut = new HDFWalkerOutput(W, root_name_, myComm);
  branchEngine->start(root_name_);
  branchEngine->write(root_name_);

  if (qmcdriver_input_.get_reset_random())
  {
    app_log() << "  Regenerate random seeds." << std::endl;
    RandomNumberControl::make_seeds();
  }

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
  population_.createWalkers(num_crowds_, walkers_per_crowd_, nwalkers, positions);
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

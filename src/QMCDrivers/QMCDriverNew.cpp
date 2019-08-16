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
QMCDriverNew::QMCDriverNew(QMCDriverInput& input,
                           MCPopulation& population,
                           TrialWaveFunction& psi,
                           QMCHamiltonian& h,
                           WaveFunctionPool& ppool,
                           Communicate* comm)
    : MPIObjectBase(comm),
      qmcdriver_input_(input),
      num_crowds_(input.get_num_crowds()),
      branchEngine(0),
      DriftModifier(0),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      Estimators(0),
      wOut(0)
{
  reset_random = false;

  QMCType = "invalid";

  ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
  rotation = 0;

  checkpointTimer = TimerManager.createTimer("checkpoint::recordBlock", timer_level_medium);
}

int QMCDriverNew::addObservable(const std::string& aname)
{
  if (Estimators)
    return Estimators->addObservable(aname.c_str());
  else
    return -1;
}

QMCDriverNew::RealType QMCDriverNew::getObservable(int i) { return Estimators->getObservable(i); }


QMCDriverNew::~QMCDriverNew()
{
  delete_iter(Rng.begin(), Rng.end());
  if (DriftModifier)
    delete DriftModifier;
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
  //  deltaR.resize(W.getTotalNum());
  //  drift.resize(W.getTotalNum());

  if (!qmcdriver_input_.get_append_run())
    current_step_ = 0;
  else
    current_step_ = qmcdriver_input_.get_starting_step();

  //int numCopies = (H1.empty()) ? 1 : H1.size();
  //W.resetWalkerProperty(numCopies);

  //create branchEngine first
  if (branchEngine == 0)
  {
    branchEngine = new SimpleFixedNodeBranch(qmcdriver_input_.get_tau(), population_.get_num_global_walkers());
  }

  //create and initialize estimator
  Estimators = branchEngine->getEstimatorManager();
  if (Estimators == 0)
  {
    Estimators = new EstimatorManagerBase(myComm);
    branchEngine->setEstimatorManager(Estimators);
    branchEngine->read(h5FileRoot);
  }
  if (DriftModifier == 0)
    DriftModifier = createDriftModifier(cur, myComm);
  DriftModifier->parseXML(cur);
  branchEngine->put(cur);
  // Estimators->put(W, H, cur);
  // if (wOut == 0)
  //   wOut = new HDFWalkerOutput(W, RootName, myComm);
  branchEngine->start(RootName);
  branchEngine->write(RootName);
  //use new random seeds
  if (reset_random)
  {
    app_log() << "  Regenerate random seeds." << std::endl;
    RandomNumberControl::make_seeds();
    reset_random = false;
  }
  //flush the std::ostreams
  infoSummary.flush();
  infoLog.flush();
  //increment QMCCounter of the branch engine
  branchEngine->advanceQMCCounter();
}

void QMCDriverNew::setStatus(const std::string& aname, const std::string& h5name, bool append)
{
  RootName = aname;
  app_log() << "\n========================================================="
            << "\n  Start " << QMCType << "\n  File Root " << RootName;
  app_log() << "\n=========================================================" << std::endl;
  if (h5name.size())
    h5FileRoot = h5name;
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

std::string QMCDriverNew::getRotationName(std::string RootName)
{
  std::string r_RootName;
  if (rotation % 2 == 0)
  {
    r_RootName = RootName;
  }
  else
  {
    r_RootName = RootName + ".bk";
  }
  rotation++;
  return r_RootName;
}

std::string QMCDriverNew::getLastRotationName(std::string RootName)
{
  std::string r_RootName;
  if ((rotation - 1) % 2 == 0)
  {
    r_RootName = RootName;
  }
  else
  {
    r_RootName = RootName + ".bk";
  }
  return r_RootName;
}

void QMCDriverNew::recordBlock(int block)
{
  if (qmcdriver_input_.get_dump_config() && block % qmcdriver_input_.get_check_point_period().period == 0)
  {
    checkpointTimer->start();
    branchEngine->write(RootName, true); //save energy_history
    RandomNumberControl::write(RootName, myComm);
    checkpointTimer->stop();
  }
}

bool QMCDriverNew::finalize(int block, bool dumpwalkers)
{
  //  branchEngine->finalize(W);

  if (qmcdriver_input_.get_dump_config())
    RandomNumberControl::write(RootName, myComm);

  return true;
}


/** Elements of putQMCInfo that have nothing to do with input
 */
void QMCDriverNew::setupWalkers()
{
  //if walkers are initialized via <mcwalkerset/>, use the existing one
  if (qmcdriver_input_.get_qmc_section_count() > 0 || qmc_common.is_restart)
  {
    app_log() << "Using existing walkers " << std::endl;
  }
  else
  { // always reset the walkers
    // Here we do some minimal fixing of walker numbers
    int num_threads(Concurrency::maxThreads<>());
    if (qmcdriver_input_.get_num_crowds() <= 0)
      num_crowds_ = num_threads;
    if (num_crowds_ > num_threads)
    {
      std::stringstream error_msg;
      error_msg << "Bad Input: num_crowds (" << qmcdriver_input_.get_num_crowds() << ") > num_threads (" << num_threads
                << ")\n";
      throw std::runtime_error(error_msg.str());
    }

    // Finding the equal groups that will fit the inputs request

    IndexType local_walkers = calc_default_local_walkers();
    addWalkers(local_walkers);
  }
}

/** Add walkers to the end of the ensemble of walkers.
 * @param nwalkers number of walkers to add
 *
 */
void QMCDriverNew::addWalkers(int nwalkers)
{
  // if (nwalkers > 0)
  // {
  //   //add nwalkers walkers to the end of the ensemble
  //   int nold = W.getActiveWalkers();
  //   app_log() << "  Adding " << nwalkers << " walkers to " << nold << " existing sets" << std::endl;
  //   W.createWalkers(nwalkers);
  //   if (nold)
  //   {
  //     int iw = nold;
  //     for (MCWalkerConfiguration::iterator it = W.begin() + nold; it != W.end(); ++it, ++iw)
  //       (*it)->R = W[iw % nold]->R; //assign existing walker configurations when the number of walkers change
  //   }
  // }
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

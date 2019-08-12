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
#include "Concurrency/TaskBlock.hpp"
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
      branchEngine(0),
      DriftModifier(0),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      Estimators(0),
      Traces(0),
      qmc_node(NULL),
      wOut(0)
{
  reset_random = false;
  append_run   = false;
  dump_config  = false;
  //<parameter name=" "> value </parameter>
  //accept multiple names for the same value
  //recommend using all lower cases for a new parameter
  //m_param.add(storeConfigs,"storeConfigs","int");
  //m_param.add(Period4ConfigDump,"recordConfigs","int");

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
  qmc_node = cur;
  //process common parameters
  putQMCInfo(cur);
  if (!qmcdriver_input_.get_append_run())
    current_step_ = 0;
  else
    current_step_ = qmcdriver_input_.get_starting_step();
  ////set the Tau parameter inside the Hamiltonian
  //H.setTau(Tau);
  //need to initialize properties

  //int numCopies = (H1.empty()) ? 1 : H1.size();
  //W.resetWalkerProperty(numCopies);

  //create branchEngine first
  if (branchEngine == 0)
  {
    branchEngine = new SimpleFixedNodeBranch(Tau, population_.get_num_global_walkers());
  }
  //execute the put function implemented by the derived classes
  put(cur);
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
#if !defined(REMOVE_TRACEMANAGER)
  //create and initialize traces
  if (Traces == 0)
  {
    Traces = new TraceManager(myComm);
  }
  Traces->put(traces_xml, allow_traces, RootName);
#endif
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
  if (append)
    app_log() << " append = yes ";
  else
    app_log() << " append = no ";
  app_log() << "\n=========================================================" << std::endl;
  if (h5name.size())
    h5FileRoot = h5name;
  append_run = append;
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

void QMCDriverNew::adiosCheckpoint(int block) {}

void QMCDriverNew::adiosCheckpointFinal(int block, bool dumpwalkers) {}

void QMCDriverNew::recordBlock(int block)
{
  if (dump_config && block % check_point_period == 0)
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

  if (dump_config)
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
    IndexType rw            = qmcdriver_input_.get_requested_walkers_per_rank();
    walkers_per_crowd_      = (rw % num_crowds_) ? rw / num_crowds_ + 1 : rw / num_crowds_;
    IndexType local_walkers = walkers_per_crowd_ * num_crowds_;

    population_.set_num_local_walkers(local_walkers);
    population_.set_num_global_walkers(local_walkers * num_crowds_ * population_.get_num_ranks());

    int ndiff = 0;

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


xmlNodePtr QMCDriverNew::getQMCNode()
{
  xmlNodePtr newqmc      = xmlCopyNode(qmc_node, 1);
  xmlNodePtr current_ptr = NULL;
  xmlNodePtr cur         = newqmc->children;
  while (cur != NULL && current_ptr == NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "parameter")
    {
      const xmlChar* aptr = xmlGetProp(cur, (const xmlChar*)"name");
      if (aptr)
      {
        if (xmlStrEqual(aptr, (const xmlChar*)"current"))
          current_ptr = cur;
      }
    }
    cur = cur->next;
  }
  if (current_ptr == NULL)
  {
    current_ptr = xmlNewTextChild(newqmc, NULL, (const xmlChar*)"parameter", (const xmlChar*)"0");
    xmlNewProp(current_ptr, (const xmlChar*)"name", (const xmlChar*)"current");
  }
  getContent(current_step_, current_ptr);
  return newqmc;
}

} // namespace qmcplusplus

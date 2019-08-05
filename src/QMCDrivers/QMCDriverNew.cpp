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
 */

#include "QMCDrivers/QMCDriverNew.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "HDFVersion.h"
#include <qmc_common.h>
#include <limits>
#include <typeinfo>

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"

#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

// #ifdef QMC_CUDA
// #include "type_traits/CUDATypes.h"
// #endif

namespace qmcplusplus
{
QMCDriverNew::QMCDriverNew(MCPopulation& population,
                           TrialWaveFunction& psi,
                           QMCHamiltonian& h,
                           WaveFunctionPool& ppool,
                           Communicate* comm)
    : MPIObjectBase(comm),
      branchEngine(0),
      DriftModifier(0),
      population_(population),
      Psi(psi),
      H(h),
      psiPool(ppool),
      Estimators(0),
      Traces(0),
      qmcNode(NULL),
      wOut(0)
{
  ResetRandom     = false;
  AppendRun       = false;
  DumpConfig      = false;
  ConstPopulation = true; //default is a fixed population method
  IsQMCDriver     = true;
  allow_traces    = false;
  MyCounter       = 0;
  //<parameter name=" "> value </parameter>
  //accept multiple names for the same value
  //recommend using all lower cases for a new parameter
  RollBackBlocks = 0;
  m_param.add(RollBackBlocks, "rewind", "int");
  Period4CheckPoint = -1;
  storeConfigs      = 0;
  //m_param.add(storeConfigs,"storeConfigs","int");
  m_param.add(storeConfigs, "storeconfigs", "int");
  m_param.add(storeConfigs, "store_configs", "int");
  Period4CheckProperties = 100;
  m_param.add(Period4CheckProperties, "checkProperties", "int");
  m_param.add(Period4CheckProperties, "checkproperties", "int");
  m_param.add(Period4CheckProperties, "check_properties", "int");
  Period4WalkerDump = 0;
  //m_param.add(Period4WalkerDump,"recordWalkers","int");
  m_param.add(Period4WalkerDump, "record_walkers", "int");
  m_param.add(Period4WalkerDump, "recordwalkers", "int");
  Period4ConfigDump = 0;
  //m_param.add(Period4ConfigDump,"recordConfigs","int");
  m_param.add(Period4ConfigDump, "recordconfigs", "int");
  m_param.add(Period4ConfigDump, "record_configs", "int");
  CurrentStep = 0;
  m_param.add(CurrentStep, "current", "int");
  nBlocks = 1;
  m_param.add(nBlocks, "blocks", "int");
  nSteps = 1;
  m_param.add(nSteps, "steps", "int");
  nSubSteps = 1;
  m_param.add(nSubSteps, "substeps", "int");
  //m_param.add(nSubSteps,"subSteps","int");
  m_param.add(nSubSteps, "sub_steps", "int");
  nWarmupSteps = 0;
  m_param.add(nWarmupSteps, "warmupsteps", "int");
  //m_param.add(nWarmupSteps,"warmupSteps","int");
  m_param.add(nWarmupSteps, "warmup_steps", "int");
  nAccept = 0;
  nReject = 0;
  //nTargetWalkers = W.getActiveWalkers();
  m_param.add(nTargetWalkers, "walkers", "int");
  //sample-related parameters
  //samples will set nTargetPopulation
  nTargetSamples       = 0;
  nStepsBetweenSamples = 1;
  m_param.add(nStepsBetweenSamples, "stepsbetweensamples", "int");
  nSamplesPerThread = 0;
  m_param.add(nSamplesPerThread, "samplesperthread", "real");
  m_param.add(nSamplesPerThread, "dmcwalkersperthread", "real");
  nTargetPopulation = 0;
  m_param.add(nTargetPopulation, "samples", "real");
  Tau = 0.1;
  //m_param.add(Tau,"timeStep","AU");
  m_param.add(Tau, "timestep", "AU");
  m_param.add(Tau, "time_step", "AU");
  //m_param.add(Tau,"Tau","AU");
  m_param.add(Tau, "tau", "AU");
  MaxCPUSecs = 360000; //100 hours
  m_param.add(MaxCPUSecs, "maxcpusecs", "real");
  // by default call recompute at the end of each block in the mixed precision case.
#ifdef QMC_CUDA
  using CTS = CUDAGlobalTypes;
  if (typeid(CTS::RealType) == typeid(float))
  {
    // gpu mixed precision
    nBlocksBetweenRecompute = 1;
  }
  else if (typeid(CTS::RealType) == typeid(double))
  {
    // gpu double precision
    nBlocksBetweenRecompute = 0;
  }
#else
#ifdef MIXED_PRECISION
  // cpu mixed precision
  nBlocksBetweenRecompute = 1;
#else
  // cpu double precision
  nBlocksBetweenRecompute = 0;
#endif
#endif
  m_param.add(nBlocksBetweenRecompute, "blocks_between_recompute", "int");
  QMCType = "invalid";
  ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
  rotation = 0;

  checkpointTimer = TimerManager.createTimer("checkpoint::recordBlock", timer_level_medium);
}

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
  if (ResetRandom)
  {
    app_log() << "  Regenerate random seeds." << std::endl;
    RandomNumberControl::make_seeds();
    ResetRandom = false;
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
  AppendRun = append;
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

void QMCDriverNew::adiosCheckpoint(int block)
{
}

void QMCDriverNew::adiosCheckpointFinal(int block, bool dumpwalkers)
{
}

void QMCDriverNew::recordBlock(int block)
{
  if (DumpConfig && block % Period4CheckPoint == 0)
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

  if (DumpConfig)
    RandomNumberControl::write(RootName, myComm);

  return true;
}

/** Add walkers to the end of the ensemble of walkers.
 * @param nwalkers number of walkers to add
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


/** Parses the xml input file for parameter definitions for a single qmc simulation.
 *
 * Basic parameters are handled here and each driver will perform its own initialization with the input
 * attribute list
 * - checkpoint="-1|0|n" default=-1
 *   -- 1 = do not write anything
 *   -- 0 = dump after the completion of a qmc section
 *   -- n = dump after n blocks
 * - kdelay = "0|1|n" default=0
 */
bool QMCDriverNew::putQMCInfo(xmlNodePtr cur)
{  
    int defaultw = omp_get_max_threads();
    OhmmsAttributeSet aAttrib;
    aAttrib.add(k_delay, "kdelay");
    aAttrib.put(cur);
    if (cur != NULL)
    {
      //initialize the parameter set
      m_param.put(cur);

      xmlNodePtr tcur = cur->children;
      //determine how often to print walkers to hdf5 file
      while (tcur != NULL)
      {
        std::string cname((const char*)(tcur->name));
        if (cname == "record")
        {
          //dump walkers for optimization
          OhmmsAttributeSet rAttrib;
          rAttrib.add(Period4WalkerDump, "stride");
          rAttrib.add(Period4WalkerDump, "period");
          rAttrib.put(tcur);
        }
        else if (cname == "checkpoint")
        {
          OhmmsAttributeSet rAttrib;
          rAttrib.add(Period4CheckPoint, "stride");
          rAttrib.add(Period4CheckPoint, "period");
          rAttrib.put(tcur);
          //DumpConfig=(Period4CheckPoint>0);
        }
        else if (cname == "dumpconfig")
        {
          OhmmsAttributeSet rAttrib;
          rAttrib.add(Period4ConfigDump, "stride");
          rAttrib.add(Period4ConfigDump, "period");
          rAttrib.put(tcur);
        }
        else if (cname == "random")
        {
          ResetRandom = true;
        }
        tcur = tcur->next;
      }
    }
    //set the minimum blocks
    if (nBlocks < 1)
      nBlocks = 1;

    DumpConfig = (Period4CheckPoint >= 0);
    if (Period4CheckPoint < 1)
      Period4CheckPoint = nBlocks;
    //reset CurrentStep to zero if qmc/@continue='no'
    if (!AppendRun)
      CurrentStep = 0;

    //if walkers are initialized via <mcwalkerset/>, use the existing one
    if (qmc_common.qmc_counter || qmc_common.is_restart)
    {
      app_log() << "Using existing walkers " << std::endl;
    }
    else
    { //always reset the walkers
  #ifdef QMC_CUDA
      int nths(1);
  #else
      int nths(omp_get_max_threads());
  #endif
      nTargetWalkers = (std::max(nths, (nTargetWalkers / nths) * nths));
      int nw         = W.getActiveWalkers();
      int ndiff      = 0;
      if (nw)
      {
        // nTargetWalkers == 0, if it is not set by the input file
        ndiff = (nTargetWalkers) ? nTargetWalkers - nw : 0;
      }
      else
      {
        ndiff = (nTargetWalkers) ? nTargetWalkers : defaultw;
      }
      addWalkers(ndiff);
    }

    return (W.getActiveWalkers() > 0);
}

xmlNodePtr QMCDriverNew::getQMCNode()
{
  xmlNodePtr newqmc      = xmlCopyNode(qmcNode, 1);
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
  getContent(CurrentStep, current_ptr);
  return newqmc;
}

} // namespace qmcplusplus

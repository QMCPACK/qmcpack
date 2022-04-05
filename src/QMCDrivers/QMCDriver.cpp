//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDriver.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "RandomNumberControl.h"
#include "hdf/HDFVersion.h"
#include "Utilities/qmc_common.h"
#include <limits>
#include <typeinfo>

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif
#ifdef QMC_CUDA
#include "type_traits/CUDATypes.h"
#endif

namespace qmcplusplus
{
QMCDriver::QMCDriver(MCWalkerConfiguration& w,
                     TrialWaveFunction& psi,
                     QMCHamiltonian& h,
                     Communicate* comm,
                     const std::string& QMC_driver_type,
                     bool enable_profiling)
    : MPIObjectBase(comm),
      DriftModifier(0),
      qmcNode(NULL),
      QMCType(QMC_driver_type),
      W(w),
      Psi(psi),
      H(h),
      driver_scope_timer_(*timer_manager.createTimer(QMC_driver_type, timer_level_coarse)),
      driver_scope_profiler_(enable_profiling)
{
  ResetRandom  = false;
  AppendRun    = false;
  DumpConfig   = false;
  IsQMCDriver  = true;
  allow_traces = false;
  MyCounter    = 0;
  //<parameter name=" "> value </parameter>
  //accept multiple names for the same value
  //recommend using all lower cases for a new parameter
  Period4CheckPoint = -1;
  storeConfigs      = 0;
  //m_param.add(storeConfigs,"storeConfigs");
  m_param.add(storeConfigs, "storeconfigs");
  m_param.add(storeConfigs, "store_configs");
  Period4CheckProperties = 100;
  m_param.add(Period4CheckProperties, "checkProperties");
  m_param.add(Period4CheckProperties, "checkproperties");
  m_param.add(Period4CheckProperties, "check_properties");
  Period4WalkerDump = 0;
  //m_param.add(Period4WalkerDump,"recordWalkers");
  m_param.add(Period4WalkerDump, "record_walkers");
  m_param.add(Period4WalkerDump, "recordwalkers");
  Period4ConfigDump = 0;
  //m_param.add(Period4ConfigDump,"recordConfigs");
  m_param.add(Period4ConfigDump, "recordconfigs");
  m_param.add(Period4ConfigDump, "record_configs");
  CurrentStep = 0;
  m_param.add(CurrentStep, "current");
  nBlocks = 1;
  m_param.add(nBlocks, "blocks");
  nSteps = 1;
  m_param.add(nSteps, "steps");
  nSubSteps = 1;
  m_param.add(nSubSteps, "substeps");
  //m_param.add(nSubSteps,"subSteps");
  m_param.add(nSubSteps, "sub_steps");
  nWarmupSteps = 0;
  m_param.add(nWarmupSteps, "warmupsteps");
  //m_param.add(nWarmupSteps,"warmupSteps");
  m_param.add(nWarmupSteps, "warmup_steps");
  nAccept        = 0;
  nReject        = 0;
  nTargetWalkers = W.getActiveWalkers();
  m_param.add(nTargetWalkers, "walkers");
  //sample-related parameters
  //samples will set nTargetPopulation
  nTargetSamples       = 0;
  nStepsBetweenSamples = 1;
  m_param.add(nStepsBetweenSamples, "stepsbetweensamples");
  nSamplesPerThread = 0;
  m_param.add(nSamplesPerThread, "samplesperthread");
  m_param.add(nSamplesPerThread, "dmcwalkersperthread");

  nTargetPopulation = 0;

  m_param.add(nTargetPopulation, "samples");

  SpinMass = 1.0;
  m_param.add(SpinMass, "SpinMass");

  Tau = 0.1;
  //m_param.add(Tau,"timeStep");
  m_param.add(Tau, "timestep");
  m_param.add(Tau, "time_step");
  //m_param.add(Tau,"Tau");
  m_param.add(Tau, "tau");
  MaxCPUSecs = 360000; //100 hours
  m_param.add(MaxCPUSecs, "maxcpusecs", {}, TagStatus::DEPRECATED);
  m_param.add(MaxCPUSecs, "max_seconds");
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
  m_param.add(nBlocksBetweenRecompute, "blocks_between_recompute");
  ////add each OperatorBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
  rotation = 0;

  checkpointTimer = timer_manager.createTimer("checkpoint::recordBlock", timer_level_medium);
}

QMCDriver::~QMCDriver()
{
  if (DriftModifier)
    delete DriftModifier;
}

void QMCDriver::add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi)
{
  H1.push_back(h);
  Psi1.push_back(psi);
}

/** process a <qmc/> element
 * @param cur xmlNode with qmc tag
 *
 * This function is called before QMCDriver::run and following actions are taken:
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
void QMCDriver::process(xmlNodePtr cur)
{
  deltaR.resize(W.getTotalNum());
  drift.resize(W.getTotalNum());
  qmcNode = cur;
  //process common parameters
  putQMCInfo(cur);
  ////set the Tau parameter inside the Hamiltonian
  //H.setTau(Tau);
  //need to initialize properties
  int numCopies = (H1.empty()) ? 1 : H1.size();
  W.resetWalkerProperty(numCopies);
  //create branchEngine first
  if (!branchEngine)
  {
    branchEngine = std::make_unique<BranchEngineType>(Tau, W.getGlobalNumWalkers());
  }
  //execute the put function implemented by the derived classes
  put(cur);
  //create and initialize estimator
  Estimators = branchEngine->getEstimatorManager();
  if (Estimators == nullptr)
  {
    branchEngine->setEstimatorManager(std::make_unique<EstimatorManagerBase>(myComm));
    Estimators = branchEngine->getEstimatorManager();
    branchEngine->read(h5FileRoot);
  }
  if (DriftModifier == 0)
    DriftModifier = createDriftModifier(cur, myComm);
  DriftModifier->parseXML(cur);
#if !defined(REMOVE_TRACEMANAGER)
  //create and initialize traces
  if (!Traces)
  {
    Traces = std::make_unique<TraceManager>(myComm);
  }
  Traces->put(traces_xml, allow_traces, RootName);
#endif
  branchEngine->put(cur);
  Estimators->put(H, cur);
  if (!wOut)
    wOut = std::make_unique<HDFWalkerOutput>(W.getTotalNum(), RootName, myComm);
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

void QMCDriver::setStatus(const std::string& aname, const std::string& h5name, bool append)
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
 */
void QMCDriver::putWalkers(std::vector<xmlNodePtr>& wset)
{
  if (wset.empty())
    return;
  int nfile = wset.size();
  HDFWalkerInputManager W_in(W, W.getTotalNum(), myComm);
  for (int i = 0; i < wset.size(); i++)
    if (W_in.put(wset[i]))
      h5FileRoot = W_in.getFileRoot();
  //clear the walker set
  wset.clear();
  int nwtot = W.getActiveWalkers();
  myComm->bcast(nwtot);
  if (nwtot)
  {
    int np = myComm->size();
    std::vector<int> nw(np, 0), nwoff(np + 1, 0);
    nw[myComm->rank()] = W.getActiveWalkers();
    myComm->allreduce(nw);
    for (int ip = 0; ip < np; ++ip)
      nwoff[ip + 1] = nwoff[ip] + nw[ip];
    W.setGlobalNumWalkers(nwoff[np]);
    W.setWalkerOffsets(nwoff);
    qmc_common.is_restart = true;
  }
  else
    qmc_common.is_restart = false;
}

std::string QMCDriver::getRotationName(std::string RootName)
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

std::string QMCDriver::getLastRotationName(std::string RootName)
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

void QMCDriver::recordBlock(int block)
{
  if (DumpConfig && block % Period4CheckPoint == 0)
  {
    checkpointTimer->start();
    wOut->dump(W, block);
    branchEngine->write(RootName, true); //save energy_history
    RandomNumberControl::write(RootName, myComm);
    checkpointTimer->stop();
  }
}

bool QMCDriver::finalize(int block, bool dumpwalkers)
{
  if (DumpConfig && dumpwalkers)
    wOut->dump(W, block);
  nTargetWalkers = W.getActiveWalkers();
  MyCounter++;
  infoSummary.flush();
  infoLog.flush();

  branchEngine->finalize(W);

  if (DumpConfig)
    RandomNumberControl::write(RootName, myComm);

  return true;
}

/** Add walkers to the end of the ensemble of walkers.
 * @param nwalkers number of walkers to add
 */
void QMCDriver::addWalkers(int nwalkers)
{
  if (nwalkers > 0)
  {
    //add nwalkers walkers to the end of the ensemble
    int nold = W.getActiveWalkers();
    app_log() << "  Adding " << nwalkers << " walkers to " << nold << " existing sets" << std::endl;
    W.createWalkers(nwalkers);
    if (nold)
    {
      int iw = nold;
      for (MCWalkerConfiguration::iterator it = W.begin() + nold; it != W.end(); ++it, ++iw)
        (*it)->R = W[iw % nold]->R; //assign existing walker configurations when the number of walkers change
    }
  }
  else if (nwalkers < 0)
  {
    W.destroyWalkers(-nwalkers);
    app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers()
              << std::endl;
  }
  else
  {
    app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." << std::endl;
  }
  setWalkerOffsets();

  ////update the global number of walkers
  ////int nw=W.getActiveWalkers();
  ////myComm->allreduce(nw);
}

void QMCDriver::setWalkerOffsets()
{
  std::vector<int> nw(myComm->size(), 0), nwoff(myComm->size() + 1, 0);
  nw[myComm->rank()] = W.getActiveWalkers();
  myComm->allreduce(nw);
  for (int ip = 0; ip < myComm->size(); ip++)
    nwoff[ip + 1] = nwoff[ip] + nw[ip];
  W.setGlobalNumWalkers(nwoff[myComm->size()]);
  W.setWalkerOffsets(nwoff);
  long id = nwoff[myComm->rank()];
  for (int iw = 0; iw < nw[myComm->rank()]; ++iw, ++id)
  {
    W[iw]->ID       = id;
    W[iw]->ParentID = id;
  }
  app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples << std::endl;
  app_log() << "  Total weight: " << W.EnsembleProperty.Weight << std::endl;
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
bool QMCDriver::putQMCInfo(xmlNodePtr cur)
{
  if (!IsQMCDriver)
  {
    app_log() << getName() << "  Skip QMCDriver::putQMCInfo " << std::endl;
    return true;
  }

  ////store the current nSteps and nStepsBetweenSamples
  //int oldStepsBetweenSamples=nStepsBetweenSamples;
  //int oldSteps=nSteps;

  //set the default walker to the number of threads times 10
  Period4CheckPoint = -1;
  // set default for delayed update streak k to zero, meaning use the original Sherman-Morrison rank-1 update
  // if kdelay is set to k (k>1), then the new rank-k scheme is used
#ifdef QMC_CUDA
  kDelay = Psi.getndelay();
#endif
  int defaultw = omp_get_max_threads();
  OhmmsAttributeSet aAttrib;
  aAttrib.add(Period4CheckPoint, "checkpoint");
  aAttrib.add(kDelay, "kdelay");
  aAttrib.put(cur);
#ifdef QMC_CUDA
  W.setkDelay(kDelay);
  kDelay = W.getkDelay(); // in case number is sanitized
#endif
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

  // check input parameters collected by m_param
  //set the minimum blocks
  if (nBlocks < 1)
  {
    app_warning() << "Input parameter \"blocks\" must be positive! Set to 1. User input value " << nBlocks << std::endl;
    nBlocks = 1;
  }

  //set the minimum nSteps
  if (nSteps < 1)
  {
    app_warning() << "Input parameter \"steps\" must be positive! Set to 1. User input value " << nSteps << std::endl;
    nSteps = 1;
  }

  //set the minimum nSubSteps
  if (nSubSteps < 1)
  {
    app_warning() << "Input parameter \"substeps\" must be positive! Set to 1. User input value " << nSubSteps
                  << std::endl;
    nSubSteps = 1;
  }

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
  {
    app_log() << "Resetting walkers" << std::endl;
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

xmlNodePtr QMCDriver::getQMCNode()
{
  xmlNodePtr newqmc      = xmlCopyNode(qmcNode, 1);
  xmlNodePtr current_ptr = NULL;
  xmlNodePtr cur         = newqmc->children;
  while (cur != NULL && current_ptr == NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "parameter")
    {
      const std::string name(getXMLAttributeValue(cur, "name"));
      if (name == "current")
        current_ptr = cur;
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

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCMain.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCMain.h"
#include "DeviceManager.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Platforms/Host/OutputManager.h"
#include "Utilities/Timer.h"
#include "Utilities/TimerManager.h"
#include "Utilities/RunTimeManager.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/InitMolecularSystem.h"
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "Concurrency/OpenMP.h"
#include <queue>
#include <cstring>
#include "hdf/HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/qmc_common.h"
#ifdef BUILD_AFQMC
#include "AFQMC/AFQMCFactory.h"
#endif
#ifdef BUILD_FCIQMC
#include "FCIQMC/App/SQCFactory.h"
#endif

#define STR_VAL(arg) #arg
#define GET_MACRO_VAL(arg) STR_VAL(arg)

namespace qmcplusplus
{
QMCMain::QMCMain(Communicate* c)
    : QMCMainState(c),
      QMCAppBase(),
      FirstQMC(true)
#if !defined(REMOVE_TRACEMANAGER)
      ,
      traces_xml(NULL)
#endif
{
  Communicate node_comm;
  node_comm.initializeAsNodeComm(*OHMMS::Controller);
  // assign accelerators within a node
  DeviceManager::initializeGlobalDeviceManager(node_comm.rank(), node_comm.size());

  app_summary() << "\n=====================================================\n"
                << "                    QMCPACK " << QMCPACK_VERSION_MAJOR << "." << QMCPACK_VERSION_MINOR << "."
                << QMCPACK_VERSION_PATCH << "\n\n"
                << "       (c) Copyright 2003-  QMCPACK developers\n\n"
                << "                    Please cite:\n"
                << " J. Kim et al. J. Phys. Cond. Mat. 30 195901 (2018)\n"
                << "      https://doi.org/10.1088/1361-648X/aab9c3\n";
  qmc_common.print_git_info_if_present(app_summary());
  app_summary() << "=====================================================\n";
  qmc_common.print_options(app_log());
  // clang-format off
  app_summary()
#if !defined(HAVE_MPI)
      << "\n  Built without MPI. Running in serial or with OMP threads." << std::endl
#endif
      << "\n  Total number of MPI ranks = " << OHMMS::Controller->size()
      << "\n  Number of MPI groups      = " << myComm->getNumGroups()
      << "\n  MPI group ID              = " << myComm->getGroupID()
      << "\n  Number of ranks in group  = " << myComm->size()
      << "\n  MPI ranks per node        = " << node_comm.size()
#if defined(ENABLE_OFFLOAD) || defined(ENABLE_CUDA) || defined(ENABLE_ROCM)
      << "\n  Accelerators per node     = " << DeviceManager::getGlobal().getNumDevices()
#endif
      << std::endl;
  // clang-format on

#pragma omp parallel
  {
    const int L1_tid = omp_get_thread_num();
    if (L1_tid == 0)
      app_summary() << "  OMP 1st level threads     = " << omp_get_num_threads() << std::endl;
#pragma omp parallel
    {
      const int L2_tid         = omp_get_thread_num();
      const int L2_num_threads = omp_get_num_threads();
      if (L1_tid == 0 && L2_tid == 0)
      {
        if (L2_num_threads == 1)
          app_summary() << "  OMP nested threading disabled or only 1 thread on the 2nd level" << std::endl;
        else
          app_summary() << "  OMP 2nd level threads     = " << L2_num_threads << std::endl;
      }
    }
  }
  app_summary() << "\n  Precision used in this calculation, see definitions in the manual:"
                << "\n  Base precision      = " << GET_MACRO_VAL(OHMMS_PRECISION)
                << "\n  Full precision      = " << GET_MACRO_VAL(OHMMS_PRECISION_FULL)
#ifdef QMC_CUDA
                << "\n  CUDA base precision = " << GET_MACRO_VAL(CUDA_PRECISION)
                << "\n  CUDA full precision = " << GET_MACRO_VAL(CUDA_PRECISION_FULL)
#endif
                << std::endl;

  // Record features configured in cmake or selected via command-line arguments to the printout
  app_summary() << std::endl;
#if !defined(ENABLE_OFFLOAD) && !defined(ENABLE_CUDA) && !defined(QMC_CUDA) && !defined(ENABLE_ROCM)
  app_summary() << "  CPU only build" << std::endl;
#else
#if defined(ENABLE_OFFLOAD)
  app_summary() << "  OpenMP target offload to accelerators build option is enabled" << std::endl;
#endif
#if defined(ENABLE_CUDA) || defined(QMC_CUDA)
  app_summary() << "  CUDA acceleration build option is enabled" << std::endl;
#endif
#if defined(ENABLE_ROCM)
  app_summary() << "  ROCM acceleration build option is enabled" << std::endl;
#endif
#endif
#ifdef ENABLE_TIMERS
  app_summary() << "  Timer build option is enabled. Current timer level is "
                << timer_manager.get_timer_threshold_string() << std::endl;
#endif
  app_summary() << std::endl;
  app_summary().flush();
}

///destructor
QMCMain::~QMCMain()
{
  // free last_driver before clearing P,Psi,H clones
  last_driver.reset();
  CloneManager::clearClones();
}


bool QMCMain::execute()
{
  Timer t0;
  if (XmlDocStack.empty())
  {
    ERRORMSG("No valid input file exists! Aborting QMCMain::execute")
    return false;
  }

  std::string simulationType = "realspaceQMC";
  { // mmorales: is this necessary??? Don't want to leave xmlNodes lying around unused
    xmlNodePtr cur = XmlDocStack.top()->getRoot();
    OhmmsAttributeSet simType;
    simType.add(simulationType, "type");
    simType.add(simulationType, "name");
    simType.add(simulationType, "method");
    simType.put(cur);
  }

#ifdef BUILD_AFQMC
  if (simulationType == "afqmc")
  {
    NewTimer* t2 = timer_manager.createTimer("Total", timer_level_coarse);
    ScopedTimer t2_scope(*t2);
    app_log() << std::endl
              << "/*************************************************\n"
              << " ********  This is an AFQMC calculation   ********\n"
              << " *************************************************" << std::endl;
    xmlNodePtr cur = XmlDocStack.top()->getRoot();

    xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
    //initialize the random number generator
    xmlNodePtr rptr = myRandomControl.initialize(m_context);

    auto world = boost::mpi3::environment::get_world_instance();
    afqmc::AFQMCFactory afqmc_fac(world);
    if (!afqmc_fac.parse(cur))
    {
      app_log() << " Error in AFQMCFactory::parse() ." << std::endl;
      return false;
    }
    cur = XmlDocStack.top()->getRoot();
    return afqmc_fac.execute(cur);
  }
#else
  if (simulationType == "afqmc")
  {
    app_error() << " Executable not compiled with AFQMC. Recompile with BUILD_AFQMC set to 1." << std::endl;
    return false;
  }
#endif

  NewTimer* t2 = timer_manager.createTimer("Total", timer_level_coarse);
  t2->start();

  NewTimer* t3 = timer_manager.createTimer("Startup", timer_level_coarse);
  t3->start();

  //validate the input file
  bool success = validateXML();
  if (!success)
    myComm->barrier_and_abort("QMCMain::execute. Input document does not contain valid objects");

  //initialize all the instances of distance tables and evaluate them
  ptclPool->reset();
  infoSummary.flush();
  infoLog.flush();
  app_log() << "  Initialization Execution time = " << std::setprecision(4) << t0.elapsed() << " secs" << std::endl;
  //write stuff
  app_log() << "=========================================================\n";
  app_log() << " Summary of QMC systems \n";
  app_log() << "=========================================================\n";
  ptclPool->get(app_log());
  hamPool->get(app_log());
  OHMMS::Controller->barrier();
  if (qmc_common.dryrun)
  {
    app_log() << "  dryrun == 1 Ignore qmc/loop elements " << std::endl;
    APP_ABORT("QMCMain::execute");
  }
  t3->stop();
  Timer t1;
  curMethod              = std::string("invalid");
  qmc_common.qmc_counter = 0;
  for (int qa = 0; qa < m_qmcaction.size(); qa++)
  {
    if (run_time_manager.isStopNeeded())
      break;
    xmlNodePtr cur = m_qmcaction[qa].first;
    std::string cname((const char*)cur->name);
    if (cname == "qmc" || cname == "optimize")
    {
      executeQMCSection(cur);
      qmc_common.qmc_counter++; // increase the counter
    }
    else if (cname == "loop")
    {
      qmc_common.qmc_counter = 0;
      executeLoop(cur);
      qmc_common.qmc_counter = 0;
    }
    else if (cname == "cmc")
    {
      executeCMCSection(cur);
    }
    else if (cname == "debug")
    {
      executeDebugSection(cur);
      app_log() << "  Debug is done. Skip the rest of the input " << std::endl;
      break;
    }
  }
  // free if m_qmcation owns the memory of xmlNodePtr before clearing
  for (auto& qmcactionPair : m_qmcaction)
    if (!qmcactionPair.second)
      xmlFreeNode(qmcactionPair.first);

  m_qmcaction.clear();
  t2->stop();
  app_log() << "  Total Execution time = " << std::setprecision(4) << t1.elapsed() << " secs" << std::endl;
  if (is_manager())
  {
    //generate multiple files
    xmlNodePtr mcptr = NULL;
    if (m_walkerset.size())
      mcptr = m_walkerset[0];
    //remove input mcwalkerset but one
    for (int i = 1; i < m_walkerset.size(); i++)
    {
      xmlUnlinkNode(m_walkerset[i]);
      xmlFreeNode(m_walkerset[i]);
    }
    m_walkerset.clear(); //empty the container
    std::ostringstream np_str, v_str;
    np_str << myComm->size();
    HDFVersion cur_version;
    v_str << cur_version[0] << " " << cur_version[1];
    xmlNodePtr newmcptr = xmlNewNode(NULL, (const xmlChar*)"mcwalkerset");
    xmlNewProp(newmcptr, (const xmlChar*)"fileroot", (const xmlChar*)myProject.CurrentMainRoot().c_str());
    xmlNewProp(newmcptr, (const xmlChar*)"node", (const xmlChar*)"-1");
    xmlNewProp(newmcptr, (const xmlChar*)"nprocs", (const xmlChar*)np_str.str().c_str());
    xmlNewProp(newmcptr, (const xmlChar*)"version", (const xmlChar*)v_str.str().c_str());
    //#if defined(H5_HAVE_PARALLEL)
    xmlNewProp(newmcptr, (const xmlChar*)"collected", (const xmlChar*)"yes");
    //#else
    //      xmlNewProp(newmcptr,(const xmlChar*)"collected",(const xmlChar*)"no");
    //#endif
    if (mcptr == NULL)
    {
      xmlAddNextSibling(lastInputNode, newmcptr);
    }
    else
    {
      xmlReplaceNode(mcptr, newmcptr);
    }
    saveXml();
  }
  return true;
}

void QMCMain::executeLoop(xmlNodePtr cur)
{
  int niter = 1;
  OhmmsAttributeSet a;
  a.add(niter, "max");
  a.put(cur);
  //reset qmc_counter
  qmc_common.qmc_counter = 0;
  app_log() << "Loop execution max-interations = " << niter << std::endl;
  for (int iter = 0; iter < niter; iter++)
  {
    if (run_time_manager.isStopNeeded())
      break;
    xmlNodePtr tcur = cur->children;
    while (tcur != NULL)
    {
      std::string cname((const char*)tcur->name);
      if (cname == "qmc")
      {
        //prevent completed is set
        bool success = executeQMCSection(tcur, true);
        if (!success)
        {
          app_warning() << "  Terminated loop execution. A sub section returns false." << std::endl;
          return;
        }
        qmc_common.qmc_counter++; // increase the counter
      }
      tcur = tcur->next;
    }
  }
  // Destroy the last driver at the end of a loop with no further reuse of a driver needed.
  last_driver.reset(nullptr);
}

bool QMCMain::executeQMCSection(xmlNodePtr cur, bool reuse)
{
  std::string target("e");
  std::string random_test("no");
  OhmmsAttributeSet a;
  a.add(target, "target");
  a.add(random_test, "testrng");
  a.put(cur);
  if (random_test == "yes")
    RandomNumberControl::test();
  if (qmcSystem == 0)
    qmcSystem = ptclPool->getWalkerSet(target);
  bool success = runQMC(cur, reuse);
  FirstQMC     = false;
  return success;
}

/** validate the main document and (read the walker sets !)
 * @return false, if any of the basic objects is not properly created.
 *
 * Current xml schema is changing. Instead validating the input file,
 * we use xpath to create necessary objects. The order follows
 * - project: determine the simulation title and sequence
 * - random: initialize random number generator
 * - particleset: create all the particleset
 * - wavefunction: create wavefunctions
 * - hamiltonian: create hamiltonians
 * Finally, if /simulation/mcwalkerset exists, read the configurations
 * from the external files.
 * TODO: Move this out of what should be a stateless call
 */
bool QMCMain::validateXML()
{
  xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
  OhmmsXPathObject result("//project", m_context);
  myProject.setCommunicator(myComm);
  if (result.empty())
  {
    app_warning() << "Project is not defined" << std::endl;
    myProject.reset();
  }
  else
  {
    try
    {
      myProject.put(result[0]);
    }
    catch (const UniformCommunicateError& ue)
    {
      myComm->barrier_and_abort(ue.what());
    }
  }
  app_summary() << std::endl;
  myProject.get(app_summary());
  app_summary() << std::endl;
  OhmmsXPathObject ham("//hamiltonian", m_context);
  if (ham.empty())
  {
    qmc_common.use_density = true;
  }
  else
  {
    for (int i = 0; i < ham.size(); ++i)
    {
      xmlNodePtr cur = ham[i]->children;
      while (cur != NULL)
      {
        std::string aname = "0";
        OhmmsAttributeSet a;
        a.add(aname, "type");
        a.put(cur);
        if (aname == "mpc" || aname == "MPC")
        {
          qmc_common.use_density = true;
        }
        cur = cur->next;
      }
    }
  }
  if (qmc_common.use_density)
  {
    app_log() << "  hamiltonian has MPC. Will read density if it is found." << std::endl << std::endl;
  }

  //initialize the random number generator
  xmlNodePtr rptr = myRandomControl.initialize(m_context);
  //preserve the input order
  xmlNodePtr cur = XmlDocStack.top()->getRoot()->children;
  lastInputNode  = NULL;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    bool inputnode = true;
    if (cname == "parallel")
    {
      putCommunicator(cur);
    }
    else if (cname == "particleset")
    {
      ptclPool->put(cur);
    }
    else if (cname == "wavefunction")
    {
      psiPool->put(cur);
    }
    else if (cname == "hamiltonian")
    {
      hamPool->put(cur);
    }
    else if (cname == "include")
    {
      //file is provided
      const std::string include_name(getXMLAttributeValue(cur, "href"));
      if (!include_name.empty())
      {
        bool success = pushDocument(include_name);
        if (success)
        {
          inputnode = processPWH(XmlDocStack.top()->getRoot());
          popDocument();
        }
        else
          myComm->barrier_and_abort("Invalid XML document");
      }
      else
        myComm->barrier_and_abort("tag \"include\" must include an \"href\" attribute.");
    }
    else if (cname == "qmcsystem")
    {
      processPWH(cur);
    }
    else if (cname == "init")
    {
      InitMolecularSystem moinit(*ptclPool);
      moinit.put(cur);
    }
#if !defined(REMOVE_TRACEMANAGER)
    else if (cname == "traces")
    {
      traces_xml = cur;
    }
#endif
    else
    {
      //everything else goes to m_qmcaction
      m_qmcaction.push_back(std::pair<xmlNodePtr, bool>(cur, true));
      inputnode = false;
    }
    if (inputnode)
      lastInputNode = cur;
    cur = cur->next;
  }

  if (ptclPool->empty())
    myComm->barrier_and_abort("QMCMain::validateXML. Illegal input. Missing particleset.");

  if (psiPool->empty())
    myComm->barrier_and_abort("QMCMain::validateXML. Illegal input. Missing wavefunction.");

  if (hamPool->empty())
    myComm->barrier_and_abort("QMCMain::validateXML. Illegal input. Missing Hamiltonian.");

  //randomize any particleset with random="yes" && random_source="ion0"
  ptclPool->randomize();

  setMCWalkers(m_context);

  return true;
}


/** grep basic objects and add to Pools
 * @param cur current node
 *
 * Recursive search  all the xml elements with particleset, wavefunction and hamiltonian
 * tags
 */
bool QMCMain::processPWH(xmlNodePtr cur)
{
  //return true and will be ignored
  if (cur == NULL)
    return true;
  bool inputnode = false;
  //save the root to grep @tilematrix
  xmlNodePtr cur_root = cur;
  cur                 = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "simulationcell")
    {
      inputnode = true;
      ptclPool->readSimulationCellXML(cur);
    }
    else if (cname == "particleset")
    {
      inputnode = true;
      ptclPool->put(cur);
    }
    else if (cname == "wavefunction")
    {
      inputnode = true;
      psiPool->put(cur);
    }
    else if (cname == "hamiltonian")
    {
      inputnode = true;
      hamPool->put(cur);
    }
    else
    //add to m_qmcaction
    {
      m_qmcaction.push_back(std::pair<xmlNodePtr, bool>(xmlCopyNode(cur, 1), false));
    }
    cur = cur->next;
  }
  //flush
  app_log().flush();
  return inputnode;
}

/** prepare for a QMC run
 * @param cur qmc element
 * @param reuse if true, the current call is from a loop
 * @return true, if a valid QMCDriver is set.
 */
bool QMCMain::runQMC(xmlNodePtr cur, bool reuse)
{
  std::unique_ptr<QMCDriverInterface> qmc_driver;
  bool append_run = false;

  if (reuse && last_driver)
    qmc_driver = std::move(last_driver);
  else
  {
    QMCDriverFactory driver_factory(myProject);
    try
    {
      QMCDriverFactory::DriverAssemblyState das = driver_factory.readSection(cur);
      qmc_driver = driver_factory.createQMCDriver(cur, das, *qmcSystem, *ptclPool, *psiPool, *hamPool, myComm);
      append_run = das.append_run;
    }
    catch (const UniformCommunicateError& ue)
    {
      myComm->barrier_and_abort(ue.what());
    }
  }

  if (qmc_driver)
  {
    if (last_branch_engine_legacy_driver)
    {
      last_branch_engine_legacy_driver->resetRun(cur);
      qmc_driver->setBranchEngine(std::move(last_branch_engine_legacy_driver));
    }

    //advance the project id
    //if it is NOT the first qmc node and qmc/@append!='yes'
    if (!FirstQMC && !append_run)
      myProject.advance();

    qmc_driver->setStatus(myProject.CurrentMainRoot(), "", append_run);
    // PD:
    // Q: How does m_walkerset_in end up being non empty?
    // A: Anytime that we aren't doing a restart.
    // So put walkers is an exceptional call. This code does not tell a useful
    // story of a QMCDriver's life.
    qmc_driver->putWalkers(m_walkerset_in);
#if !defined(REMOVE_TRACEMANAGER)
    qmc_driver->putTraces(traces_xml);
#endif
    qmc_driver->process(cur);
    infoSummary.flush();
    infoLog.flush();
    Timer qmcTimer;
    qmc_driver->run();
    app_log() << "  QMC Execution time = " << std::setprecision(4) << qmcTimer.elapsed() << " secs" << std::endl;
    // transfer the states of a driver before its destruction
    last_branch_engine_legacy_driver = qmc_driver->getBranchEngine();
    // save the driver in a driver loop
    if (reuse)
      last_driver = std::move(qmc_driver);
    return true;
  }
  else
  {
    // Ye: in which case, the code hits this?
    return false;
  }
}


/** Reads walkers sets from the restart file during XML validation
 *
 *  TODO: Move this, it is not a concern of QMCMain
 */
bool QMCMain::setMCWalkers(xmlXPathContextPtr context_)
{
  OhmmsXPathObject result("/simulation/mcwalkerset", context_);
  for (int iconf = 0; iconf < result.size(); iconf++)
  {
    xmlNodePtr mc_ptr = result[iconf];
    m_walkerset.push_back(mc_ptr);
    m_walkerset_in.push_back(mc_ptr);
  }
  //use the last mcwalkerset to initialize random numbers if possible
  if (result.size())
  {
    std::string fname;
    OhmmsAttributeSet a;
    a.add(fname, "fileroot");
    a.add(fname, "href");
    a.add(fname, "src");
    a.put(result[result.size() - 1]);
    if (fname.size())
      RandomNumberControl::read(fname, myComm);
  }
  return true;
}


} // namespace qmcplusplus

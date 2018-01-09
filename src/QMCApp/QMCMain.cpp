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
#include "QMCApp/QMCMain.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Utilities/OutputManager.h"
#include "Utilities/Timer.h"
#include "Utilities/NewTimer.h"
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Particle/DistanceTable.h"
#include "QMCDrivers/QMCDriver.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/PostProcessor.h"
#endif
#include <queue>
#include <cstring>
#include "HDFVersion.h"
#include "OhmmsData/AttributeSet.h"
#include "qmc_common.h"
#include "qmcpack_version.h"
#ifdef HAVE_ADIOS
#include "ADIOS/ADIOS_config.h"
#include <adios_read.h>
extern "C" {
#include <adios_error.h>
}
#endif
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
  : QMCDriverFactory(c), QMCAppBase(), FirstQMC(true)
#if !defined(REMOVE_TRACEMANAGER)
  , traces_xml(NULL)
#endif
{
  app_summary()
      << "\n=====================================================\n"
      <<  "                    QMCPACK "
      << QMCPACK_VERSION_MAJOR << "." << QMCPACK_VERSION_MINOR << "." << QMCPACK_VERSION_PATCH << " \n"
      << "\n  (c) Copyright 2003-  QMCPACK developers            \n"
#if defined(QMCPACK_GIT_BRANCH)
      << "\n  Git branch: " << QMCPACK_GIT_BRANCH
      << "\n  Last git commit: " << QMCPACK_GIT_HASH
      << "\n  Last commit date: " << QMCPACK_GIT_COMMIT_LAST_CHANGED
#endif
      << "\n=====================================================\n";
  qmc_common.print_options(app_log());
  app_summary()
      << "\n  MPI Nodes            = " << OHMMS::Controller->size()
      << "\n  MPI Nodes per group  = " << myComm->size()
      << "\n  MPI Group ID         = " << myComm->getGroupID()
      << "\n  OMP_NUM_THREADS      = " << omp_get_max_threads()
      << std::endl;
  app_summary()
      << "\n  Precision used in this calculation, see definitions in the manual:"
      << "\n  Base precision      = " << GET_MACRO_VAL(OHMMS_PRECISION)
      << "\n  Full precision      = " << GET_MACRO_VAL(OHMMS_PRECISION_FULL)
#ifdef QMC_CUDA
      << "\n  CUDA base precision = " << GET_MACRO_VAL(CUDA_PRECISION) 
      << "\n  CUDA full precision = " << GET_MACRO_VAL(CUDA_PRECISION_FULL)
#endif
      << std::endl;
  app_summary() << std::endl;
  app_summary().flush();
}

///destructor
QMCMain::~QMCMain()
{
}


bool QMCMain::execute()
{
  Timer t0;
  if(XmlDocStack.empty())
  {
    ERRORMSG("No valid input file exists! Aborting QMCMain::execute")
    return false;
  }

  std::string simulationType = "realspaceQMC";
  {  // mmorales: is this necessary??? Don't want to leave xmlNodes lying around unused 
    xmlNodePtr cur=XmlDocStack.top()->getRoot();
    OhmmsAttributeSet simType;
    simType.add (simulationType, "type");
    simType.add (simulationType, "name");
    simType.add (simulationType, "method");
    simType.put(cur);
  }

#ifdef BUILD_AFQMC
  if(simulationType == "afqmc") {
    NewTimer *t2 = TimerManager.createTimer("Total", timer_level_coarse);
    ScopedTimer t2_scope(t2);
    app_log() << std::endl << "/*************************************************\n"
                      << " ********  This is an AFQMC calculation   ********\n"
                      << " *************************************************" <<std::endl;
    xmlNodePtr cur=XmlDocStack.top()->getRoot(); 

    xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
    //initialize the random number generator
    xmlNodePtr rptr = myRandomControl.initialize(m_context);

    AFQMCFactory afqmc_fac(myComm,myRandomControl);
    if(!afqmc_fac.parse(cur)) {
      app_log()<<" Error in AFQMCFactory::parse() ." <<std::endl;
      return false;
    }
    cur=XmlDocStack.top()->getRoot(); 
    return afqmc_fac.execute(cur);
  } else
#else
  if(simulationType == "afqmc") {
    app_error()<<" Executable not compiled with AFQMC. Recompile with BUILD_AFQMC set to 1." <<std::endl; 
    return false;
  }
#endif

#ifdef BUILD_FCIQMC

  if(simulationType == "fciqmc") {
    app_log() << std::endl << "/*************************************************\n"
                      << " ********  This is a FCIQMC calculation   ********\n"
                      << " *************************************************" <<std::endl;

    xmlNodePtr cur=XmlDocStack.top()->getRoot();

    xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
    //initialize the random number generator
    xmlNodePtr rptr = myRandomControl.initialize(m_context);

    SQCFactory fciqmc_fac(myComm,myRandomControl);
    if(!fciqmc_fac.parse(cur)) {
      app_log()<<" Error in SQCFactory::parse() ." <<std::endl;
      return false;
    }
    cur=XmlDocStack.top()->getRoot();
    return fciqmc_fac.execute(cur);
  }
#else
  if(simulationType == "fciqmc") {
    app_error()<<" Executable not compiled with FCIQMC. Recompile with BUILD_FCIQMC set to 1." <<std::endl; 
    return false;
  }
#endif

  NewTimer *t2 = TimerManager.createTimer("Total", timer_level_coarse);
  t2->start();

  NewTimer *t3 = TimerManager.createTimer("Startup", timer_level_coarse);
  t3->start();

  //validate the input file
  bool success = validateXML();
  if(!success)
  {
    ERRORMSG("Input document does not contain valid objects")
    return false;
  }
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
  if(qmc_common.dryrun)
  {
    app_log() << "  dryrun == 1 Ignore qmc/loop elements " << std::endl;
    APP_ABORT("QMCMain::execute");
  }
  t3->stop();
  Timer t1;
  curMethod = std::string("invalid");
  qmc_common.qmc_counter=0;
  for(int qa=0; qa<m_qmcaction.size(); qa++)
  {
    xmlNodePtr cur=m_qmcaction[qa].first;
    std::string cname((const char*)cur->name);
    if(cname == "postprocess")
    {
#if !defined(REMOVE_TRACEMANAGER)
      postprocess(cur,qa);
#endif
      break;
    }
    else if(cname == "qmc" || cname == "optimize")
    {
      executeQMCSection(cur);
      qmc_common.qmc_counter++; // increase the counter
    }
    else if(cname == "loop")
    {
      qmc_common.qmc_counter=0;
      executeLoop(cur);
      qmc_common.qmc_counter=0;
    }
    else if(cname == "cmc")
    {
      executeCMCSection(cur);
    }
    else if(cname == "debug")
    {
      executeDebugSection(cur);
      app_log() << "  Debug is done. Skip the rest of the input " << std::endl;
      break;
    }
  }
  m_qmcaction.clear();
  t2->stop();
  app_log() << "  Total Execution time = " << std::setprecision(4) << t1.elapsed() << " secs" << std::endl;
  if(is_manager())
  {
    //generate multiple files
    xmlNodePtr mcptr = NULL;
    if(m_walkerset.size())
      mcptr=m_walkerset[0];
    //remove input mcwalkerset but one
    for(int i=1; i<m_walkerset.size(); i++)
    {
      xmlUnlinkNode(m_walkerset[i]);
      xmlFreeNode(m_walkerset[i]);
    }
    m_walkerset.clear();//empty the container
    std::ostringstream np_str, v_str;
    np_str<<myComm->size();
    HDFVersion cur_version;
    v_str << cur_version[0] << " " << cur_version[1];
    xmlNodePtr newmcptr = xmlNewNode(NULL,(const xmlChar*)"mcwalkerset");
    xmlNewProp(newmcptr,(const xmlChar*)"fileroot",(const xmlChar*)myProject.CurrentMainRoot());
    xmlNewProp(newmcptr,(const xmlChar*)"node",(const xmlChar*)"-1");
    xmlNewProp(newmcptr,(const xmlChar*)"nprocs",(const xmlChar*)np_str.str().c_str());
    xmlNewProp(newmcptr,(const xmlChar*)"version",(const xmlChar*)v_str.str().c_str());
//#if defined(H5_HAVE_PARALLEL)
    xmlNewProp(newmcptr,(const xmlChar*)"collected",(const xmlChar*)"yes");
//#else
//      xmlNewProp(newmcptr,(const xmlChar*)"collected",(const xmlChar*)"no");
//#endif
    if(mcptr == NULL)
    {
      xmlAddNextSibling(lastInputNode,newmcptr);
    }
    else
    {
      xmlReplaceNode(mcptr,newmcptr);
    }
    saveXml();
  }
  return true;
}

void QMCMain::executeLoop(xmlNodePtr cur)
{
  int niter=1;
  OhmmsAttributeSet a;
  a.add(niter,"max");
  a.put(cur);
  //reset qmc_counter
  qmc_common.qmc_counter=0;
  app_log() << "Loop execution max-interations = " << niter << std::endl;
  for(int iter=0; iter<niter; iter++)
  {
    xmlNodePtr tcur=cur->children;
    while(tcur != NULL)
    {
      std::string cname((const char*)tcur->name);
      if(cname == "qmc")
      {
        //prevent completed is set
        bool success = executeQMCSection(tcur, false);
        if(!success)
        {
          app_warning() << "  Terminated loop execution. A sub section returns false." << std::endl;
          return;
        }
        qmc_common.qmc_counter++; // increase the counter
      }
      tcur=tcur->next;
    }
  }
}

bool QMCMain::executeQMCSection(xmlNodePtr cur, bool noloop)
{
  std::string target("e");
  std::string random_test("no");
  OhmmsAttributeSet a;
  a.add(target,"target");
  a.add(random_test,"testrng");
  a.put(cur);
  if(random_test=="yes")
    RandomNumberControl::test();
  if(qmcSystem ==0)
    qmcSystem = ptclPool->getWalkerSet(target);
  bool success = runQMC(cur);
  FirstQMC=false;
  return success;
}

/** validate the main document
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
 */
bool QMCMain::validateXML()
{
  xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
#ifdef HAVE_ADIOS
  OhmmsXPathObject ai("//adiosinit",m_context);
  if(ai.empty())
  {
    app_warning()<<"adiosinit is not defined"<< std::endl;
  }
  else
  {
    xmlAttr* curr = ai[0]->properties;
    const char *value = (char *)xmlNodeListGetString(ai[0]->doc, curr->children, 1);
    if(!strncmp((char *)curr->name, "href", 4))
    {
      if (adios_init(value, myComm->getMPI()))
      {
        //fprintf(stderr, "Error: %s %s\n", value, adios_get_last_errmsg());
        APP_ABORT("ADIOS init error. Exiting");
      }
      else
      {
        if (OHMMS::Controller->rank() == 0)
          std::cout << "Adios is initialized" << std::endl;
        ADIOS::set_adios_init(true);
      }
      adios_read_init_method(ADIOS_READ_METHOD_BP, myComm->getMPI(), "verbose=3");
    }
  }
  OhmmsXPathObject io("//checkpoint",m_context);
  if(io.empty())
  {
    app_warning() << "checkpoint IO is not defined, no checkpoint will be written out." << std::endl;
  }
  else{
  
    xmlAttr* curr = io[0]->properties;
    char* value = NULL;
    bool UseADIOS = false;
    bool UseHDF5 = false;
    for(curr; curr; curr = curr->next)
    {
      value = (char *)xmlNodeListGetString(io[0]->doc, curr->children, 1);
      if(!strncmp((char *)curr->name, "adios", 6) && !strncmp(value, "yes", 4))
      {
        UseADIOS = true;
      }
      else if(!strncmp((char *)curr->name, "hdf5", 5) && !strncmp(value, "yes", 4))
      {
        UseHDF5 = true;
      }
      app_log() << "property: " << curr->name << ", value: " << value << std::endl;
    }
    ADIOS::initialize(UseHDF5, UseADIOS);
  }
  OhmmsXPathObject rd("//restart",m_context);
  if(rd.empty())
  {
    app_warning() << "Checkpoint restart read method is not defined. frest start." << std::endl;
  }
  else
  {
    xmlAttr* curr = rd[0]->properties;
    const char *value = (char *)xmlNodeListGetString(rd[0]->doc, curr->children, 1);
    if(!strncmp((char *)curr->name, "method", 6))
    {
      ADIOS::initialize(value);
    }
  }
#endif
  OhmmsXPathObject result("//project",m_context);
  myProject.setCommunicator(myComm);
  if(result.empty())
  {
    app_warning() << "Project is not defined" << std::endl;
    myProject.reset();
  }
  else
  {
    myProject.put(result[0]);
  }
  app_summary() << std::endl;
  myProject.get(app_summary());
  app_summary() << std::endl;
  OhmmsXPathObject ham("//hamiltonian",m_context);
  if(ham.empty())
  {
    qmc_common.use_density=true;
  }
  else
  {
    for(int i=0; i<ham.size(); ++i)
    {
      xmlNodePtr cur=ham[i]->children;
      while(cur != NULL)
      {
        std::string aname="0";
        OhmmsAttributeSet a;
        a.add(aname,"type");
        a.put(cur);
        if(aname == "mpc" || aname == "MPC")
        {
          qmc_common.use_density=true;
        }
        cur = cur->next;
      }
    }
  }
  if(qmc_common.use_density)
  {
    app_log() << "  hamiltonian has MPC. Will read density if it is found." << std::endl;
  }

  //initialize the random number generator
  xmlNodePtr rptr = myRandomControl.initialize(m_context);
  //preserve the input order
  xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;
  lastInputNode = NULL;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    bool inputnode=true;
    if(cname == "parallel")
    {
      putCommunicator(cur);
    }
    else if(cname == "particleset")
    {
      ptclPool->put(cur);
    }
    else if(cname == "wavefunction")
    {
      psiPool->put(cur);
    }
    else if(cname == "hamiltonian")
    {
      hamPool->put(cur);
    }
    else if(cname == "include")
    {
      //file is provided
      const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"href");
      if(a)
      {
        pushDocument((const char*)a);
        inputnode = processPWH(XmlDocStack.top()->getRoot());
        popDocument();
      }
    }
    else if(cname == "qmcsystem")
    {
      processPWH(cur);
    }
    else if(cname == "init")
    {
      InitMolecularSystem moinit(ptclPool);
      moinit.put(cur);
    }
#if !defined(REMOVE_TRACEMANAGER)
    else if(cname == "traces")
    {
      traces_xml = cur;
    }
#endif
    else
    {
      //everything else goes to m_qmcaction
      m_qmcaction.push_back(std::pair<xmlNodePtr,bool>(cur,true));
      inputnode=false;
    }
    if(inputnode)
      lastInputNode=cur;
    cur=cur->next;
  }
  if(ptclPool->empty())
  {
    ERRORMSG("Illegal input. Missing particleset ")
    return false;
  }
  if(psiPool->empty())
  {
    ERRORMSG("Illegal input. Missing wavefunction. ")
    return false;
  }
  if(hamPool->empty())
  {
    ERRORMSG("Illegal input. Missing hamiltonian. ")
    return false;
  }
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
  if(cur == NULL)
    return true;
  bool inputnode=false;
  //save the root to grep @tilematrix
  xmlNodePtr cur_root=cur;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "simulationcell")
    {
      inputnode=true;
      ptclPool->putLattice(cur);
    }
    else if(cname == "particleset")
    {
      inputnode=true;
      ptclPool->putTileMatrix(cur_root);
      ptclPool->put(cur);
    }
    else if(cname == "wavefunction")
    {
      inputnode=true;
      psiPool->put(cur);
    }
    else if(cname == "hamiltonian")
    {
      inputnode=true;
      hamPool->put(cur);
    }
    else
      //add to m_qmcaction
    {
      m_qmcaction.push_back(std::pair<xmlNodePtr,bool>(xmlCopyNode(cur,1),false));
    }
    cur=cur->next;
  }
  //flush
  app_log().flush();
  return inputnode;
}

/** prepare for a QMC run
 * @param cur qmc element
 * @return true, if a valid QMCDriver is set.
 */
bool QMCMain::runQMC(xmlNodePtr cur)
{
  bool append_run = setQMCDriver(myProject.m_series,cur);
  if(qmcDriver)
  {
    //advance the project id
    //if it is NOT the first qmc node and qmc/@append!='yes'
    if(!FirstQMC && !append_run)
      myProject.advance();
    qmcDriver->setStatus(myProject.CurrentMainRoot(),PrevConfigFile, append_run);
    qmcDriver->putWalkers(m_walkerset_in);
#if !defined(REMOVE_TRACEMANAGER)
    qmcDriver->putTraces(traces_xml);
#endif
    qmcDriver->process(cur);
    infoSummary.flush();
    infoLog.flush();
    Timer qmcTimer;
    NewTimer *t1 = TimerManager.createTimer(qmcDriver->getEngineName(), timer_level_coarse);
    t1->start();
    qmcDriver->run();
    t1->stop();
    app_log() << "  QMC Execution time = " << std::setprecision(4) << qmcTimer.elapsed() << " secs " << std::endl;
    //keeps track of the configuration file
    PrevConfigFile = myProject.CurrentMainRoot();
    return true;
  }
  else
  {
    return false;
  }
}

bool QMCMain::setMCWalkers(xmlXPathContextPtr context_)
{
  OhmmsXPathObject result("/simulation/mcwalkerset",context_);
  for(int iconf=0; iconf<result.size(); iconf++)
  {
    xmlNodePtr mc_ptr = result[iconf];
    m_walkerset.push_back(mc_ptr);
    m_walkerset_in.push_back(mc_ptr);
  }
  //use the last mcwalkerset to initialize random numbers if possible
  if(result.size())
  {
    std::string fname;
    OhmmsAttributeSet a;
    a.add(fname,"fileroot");
    a.add(fname,"href");
    a.add(fname,"src");
    a.put(result[result.size()-1]);
    if(fname.size())
      RandomNumberControl::read(fname,myComm);
  }
  return true;
}


void QMCMain::postprocess(xmlNodePtr cur,int qacur)
{
#if !defined(REMOVE_TRACEMANAGER)
  app_log()<<"\nQMCMain::postprocess"<< std::endl;

  int qanext = qacur+1;
  if(qanext>=m_qmcaction.size())
  {
    APP_ABORT("QMCMain::postprocess  no qmc method elements to postprocess");
  }
  else
  {
    app_log()<<"  Assuming same target particleset for all qmc methods"<< std::endl;
    xmlNodePtr next=m_qmcaction[qanext].first;
    std::string target("e");
    OhmmsAttributeSet a;
    a.add(target,"target");
    a.put(next);
    if(qmcSystem ==0)
      qmcSystem = ptclPool->getWalkerSet(target);
  }

  int series_start=myProject.m_series;
  int series_end=series_start;
  for(int qa=qanext;qa<m_qmcaction.size();++qa)
  {
    xmlNodePtr meth=m_qmcaction[qa].first;
    std::string cname((const char*)meth->name);
    if(cname=="qmc"||cname=="optimize"||cname=="cmc")
      series_end++;
    else if(cname=="loop")
    {
      int niter=1;
      OhmmsAttributeSet a;
      a.add(niter,"max");
      a.put(meth);
      series_end+=niter;
    }
  }
  app_log()<<"  Found series";
  for(int s=series_start;s<series_end;++s)
    app_log()<<" "<<s;
  app_log()<< std::endl;

  std::string id = myProject.m_title;

  if(hamPool==0)
    APP_ABORT("QMCMain::postprocess  hamPool is null");
  if(psiPool==0)
    APP_ABORT("QMCMain::postprocess  psiPool is null");
  if(ptclPool==0)
    APP_ABORT("QMCMain::postprocess  ptclPool is null");

  PostProcessor PP(id,series_start,series_end);
  PP.put(cur,*ptclPool,*psiPool,*hamPool);
  PP.postprocess();
#endif
}

}

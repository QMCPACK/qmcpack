//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file QMCMain.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCApp/QMCMain.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/Timer.h"
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Particle/DistanceTable.h"
#include "QMCDrivers/QMCDriver.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include <queue>
#include "HDFVersion.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  QMCMain::QMCMain(Communicate* c): QMCDriverFactory(c), QMCAppBase(), 
  FirstQMC(true) 
  { 

    app_log() 
      << "\n=====================================================\n"
      <<  "                    QMCPACK "
      << QMCPLUSPLUS_VERSION_MAJOR << "." << QMCPLUSPLUS_VERSION_MINOR << "." << QMCPLUSPLUS_VERSION_PATCH << " \n"
      << "\n  (c) Copyright 2003-  QMCPACK developers            \n"
#if defined(QMCPLUSPLUS_BRANCH)
      << "\n  Subversion branch " << QMCPLUSPLUS_BRANCH 
      << "\n  Last modified     " << QMCPLUSPLUS_LAST_CHANGED_DATE
#endif
      << "\n=====================================================\n";

    app_log().flush();
  }

  ///destructor
  QMCMain::~QMCMain() 
  {
  }


  bool QMCMain::execute() {

    if(XmlDocStack.empty()) {
      ERRORMSG("No valid input file exists! Aborting QMCMain::execute")
      return false;
    }

    //validate the input file
    bool success = validateXML();

    if(!success) 
    {
      ERRORMSG("Input document does not contain valid objects")
      return false;
    }

    //initialize all the instances of distance tables and evaluate them 
    ptclPool->reset();

    OhmmsInfo::flush();
    //write stuff
    app_log() << "=========================================================\n";
    app_log() << " Summary of QMC systems \n";
    app_log() << "=========================================================\n";
    ptclPool->get(app_log());
    hamPool->get(app_log());

    OHMMS::Controller->barrier();

    Timer t1;
    curMethod = string("invalid");
    //xmlNodePtr cur=m_root->children;
    for(int qa=0; qa<m_qmcaction.size(); qa++)
    {
      xmlNodePtr cur=m_qmcaction[qa].first;
      string cname((const char*)cur->name);
      if(cname == "qmc" || cname == "optimize")
      {
        executeQMCSection(cur);
      }
      else if(cname == "loop")
      {
        executeLoop(cur);
      }
      //external node, need to free
      //if(m_qmcaction[qa].second) xmlFreeNode(cur); 
    }

    m_qmcaction.clear();

    //xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;
    //while(cur != NULL) 
    //{
    //  string cname((const char*)cur->name);
    //  if(cname == "qmc" || cname == "optimize")
    //  {
    //    if(firstqmc == NULL) firstqmc=cur;
    //    executeQMCSection(cur);
    //  }
    //  else if(cname == "loop")
    //  {
    //    if(firstqmc == NULL) firstqmc=cur;
    //    executeLoop(cur);
    //  }
    //  cur=cur->next;
    //}

    app_log() << "  MPI Nodes            = " << OHMMS::Controller->size() << endl;
    app_log() << "  MPI Nodes per group  = " << myComm->size() << endl;
    app_log() << "  MPI Group ID         = " << myComm->getGroupID() << endl;
    app_log() << "  OMP_NUM_THREADS      = " << omp_get_max_threads() << endl;
    app_log() << "  Total Execution time = " << t1.elapsed() << " secs" << endl;

    //if(OHMMS::Controller->master()) {
    //if(firstqmc != NULL && myComm->master()) { //generate multiple files
    if(is_manager()) 
    { //generate multiple files

      xmlNodePtr mcptr = NULL;
      if(m_walkerset.size()) mcptr=m_walkerset[0];
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

    app_log() << "Loop execution max-interations = " << niter << endl;
    for(int iter=0; iter<niter; iter++)
    {
      xmlNodePtr tcur=cur->children;
      while(tcur != NULL) 
      {
        string cname((const char*)tcur->name);
        if(cname == "qmc")
        {
          //prevent completed is set
          bool success = executeQMCSection(tcur, false);
          if(!success)
          {
            app_warning() << "  Terminated loop execution. A sub section returns false." << endl;
            return;
          }
        }
        tcur=tcur->next;
      }
    }
  }

  bool QMCMain::executeQMCSection(xmlNodePtr cur, bool noloop)
  {
    string target("e");
    string random_test("no");
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
  bool QMCMain::validateXML() {

    xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();

    OhmmsXPathObject result("//project",m_context);

    myProject.setCommunicator(myComm);

    if(result.empty()) 
    {
      app_warning() << "Project is not defined" << endl;
      myProject.reset();
    } else {
      myProject.put(result[0]);
    }

    app_log() << endl;
    myProject.get(app_log());
    app_log() << endl;

    //initialize the random number generator
    xmlNodePtr rptr = myRandomControl.initialize(m_context);

    //preserve the input order
    xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;

    while(cur != NULL) {
      string cname((const char*)cur->name);
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
      {//file is provided
        const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"href");
        if(a) {
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
      else
      { //everything else goes to m_qmcaction
        m_qmcaction.push_back(pair<xmlNodePtr,bool>(cur,true));
        inputnode=false;
      }

      if(inputnode) lastInputNode=cur;
      cur=cur->next;
    }

    if(ptclPool->empty()) {
      ERRORMSG("Illegal input. Missing particleset ")
      return false;
    }

    if(psiPool->empty()) {
      ERRORMSG("Illegal input. Missing wavefunction. ")
      return false;
    }

    if(hamPool->empty()) {
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
  bool QMCMain::processPWH(xmlNodePtr cur) {

    //return true and will be ignored
    if(cur == NULL) return true;

    bool inputnode=true;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "simulationcell") {
        ptclPool->putLattice(cur);
      } else if(cname == "particleset") {
        ptclPool->put(cur);
      } else if(cname == "wavefunction") {
        psiPool->put(cur);
      } else if(cname == "hamiltonian") {
        hamPool->put(cur);
      } else { //add to m_qmcaction
        inputnode=false;
        m_qmcaction.push_back(pair<xmlNodePtr,bool>(xmlCopyNode(cur,1),false));
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
  bool QMCMain::runQMC(xmlNodePtr cur) {

    bool append_run = setQMCDriver(myProject.m_series,cur);

    if(qmcDriver) {

      //advance the project id 
      //if it is NOT the first qmc node and qmc/@append!='yes'
      if(!FirstQMC && !append_run) myProject.advance();

      qmcDriver->setStatus(myProject.CurrentMainRoot(),PrevConfigFile, append_run);
      qmcDriver->putWalkers(m_walkerset_in);
      qmcDriver->process(cur);

      OhmmsInfo::flush();

      Timer qmcTimer;
      qmcDriver->run();
      app_log() << "  QMC Execution time = " << qmcTimer.elapsed() << " secs " << endl;

      //keeps track of the configuration file
      PrevConfigFile = myProject.CurrentMainRoot();
      return true;
    } else {
      return false;
    }
  }

  bool QMCMain::setMCWalkers(xmlXPathContextPtr context_) {

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
      string fname;
      OhmmsAttributeSet a;
      a.add(fname,"fileroot"); a.add(fname,"href"); a.add(fname,"src");
      a.put(result[result.size()-1]);
      if(fname.size()) RandomNumberControl::read(fname,myComm);
    }
    return true;
  }
}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

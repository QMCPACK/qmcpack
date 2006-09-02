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
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Particle/DistanceTable.h"
#include "QMCDrivers/QMCDriver.h"
#include "Message/Communicate.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  QMCMain::QMCMain(int argc, char** argv): QMCAppBase(argc,argv), FirstQMC(true) { 

    app_log() << "\n=========================================================\n"
              <<   "                   qmcplusplus 0.2                       \n"
              << "\n  (c) Copyright 2003-  qmcplusplus developers          \n"
              <<   "=========================================================\n";

    app_log().flush();
  }

  ///destructor
  QMCMain::~QMCMain() {
    DEBUGMSG("QMCMain::~QMCMain")
  }


  bool QMCMain::execute() {

    if(XmlDocStack.empty()) {
      ERRORMSG("No valid input file exists! Aborting QMCMain::execute")
      return false;
    }

    //validate the input file
    bool success = validateXML();

    if(!success) {
      ERRORMSG("Input document does not contain valid objects")
      return false;
    }

    //initialize all the instances of distance tables and evaluate them 
    ptclPool->reset();


    //write stuff
    app_log() << "=========================================================\n";
    app_log() << " Summary of QMC systems \n";
    app_log() << "=========================================================\n";
    ptclPool->get(app_log());
    hamPool->get(app_log());

    curMethod = string("invalid");
    //xmlNodePtr cur=m_root->children;
    xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "qmc" || cname == "vmc" || cname == "dmc" || cname =="rmc" ||
          cname == "optimize") {
        string target("e");
        const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"target");
        if(t) target = (const char*)t;

        bool toRunQMC=true;
        t=xmlGetProp(cur,(const xmlChar*)"completed");
        if(t != NULL) {
          if(xmlStrEqual(t,(const xmlChar*)"yes")) {
            app_log() << "  This " << cname << " section is already executed." << endl;
            toRunQMC=false;
          }
        } else {
          xmlAttrPtr t1=xmlNewProp(cur,(const xmlChar*)"completed", (const xmlChar*)"no");
        }

        if(toRunQMC) {
          qmcSystem = ptclPool->getWalkerSet(target);
          bool good = runQMC(cur);
          if(good) {
            xmlAttrPtr t1=xmlSetProp(cur,(const xmlChar*)"completed", (const xmlChar*)"yes");
            //q.push_back(cur);
            FirstQMC=false;
          } else {
            app_error() << "   " << cname << " failed." << endl;
          }
          t=xmlGetProp(cur,(const xmlChar*)"id");
          if(t == NULL) {
            xmlAttrPtr t1=xmlNewProp(cur,(const xmlChar*)"id", 
                (const xmlChar*)myProject.CurrentMainRoot());
          }
        }
      }
      cur=cur->next;
    }

    if(OHMMS::Controller->master()) {
      int nproc=OHMMS::Controller->ncontexts();
      if(nproc>1) {
        xmlNodePtr t=m_walkerset.back();
        xmlSetProp(t, (const xmlChar*)"node",(const xmlChar*)"0");
        string fname((const char*)xmlGetProp(t,(const xmlChar*)"href"));
        string::size_type ending=fname.find(".p");
        string froot;
        if(ending<fname.size()) 
          froot = string(fname.begin(),fname.begin()+ending);
        else
          froot=fname;
        char pfile[128];
        for(int ip=1; ip<nproc; ip++) {
	  sprintf(pfile,"%s.p%03d",froot.c_str(),ip);
          std::ostringstream ip_str;
          ip_str<<ip;
          t = xmlAddNextSibling(t,xmlCopyNode(t,1));
          xmlSetProp(t, (const xmlChar*)"node",(const xmlChar*)ip_str.str().c_str());
          xmlSetProp(t, (const xmlChar*)"href",(const xmlChar*)pfile);
        }
      }
      saveXml();
    }
    return true;
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

    if(result.empty()) {
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
      if(cname == "particleset") {
        ptclPool->put(cur);
      } else if(cname == "wavefunction") {
        psiPool->put(cur);
      } else if(cname == "hamiltonian") {
        hamPool->put(cur);
      } else if(cname == "include") {//file is provided
        const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"href");
        if(a) {
          pushDocument((const char*)a);
          processPWH(XmlDocStack.top()->getRoot());
          popDocument();
        }
      } else if(cname == "qmcsystem") {
        processPWH(cur);
      } else if(cname == "init") {
        InitMolecularSystem moinit(ptclPool);
        moinit.put(cur);
      }
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

    setMCWalkers(m_context);

    return true;
  }   

  

  /** grep basic objects and add to Pools
   * @param cur current node 
   *
   * Recursive search  all the xml elements with particleset, wavefunction and hamiltonian
   * tags
   */
  void QMCMain::processPWH(xmlNodePtr cur) {

    if(cur == NULL) return;

    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "simulationcell") {
        DistanceTable::createSimulationCell(cur);
      } else if(cname == "particleset") {
        ptclPool->put(cur);
      } else if(cname == "wavefunction") {
        psiPool->put(cur);
      } else if(cname == "hamiltonian") {
        hamPool->put(cur);
      }
      cur=cur->next;
    }
  }

  /** prepare for a QMC run
   * @param cur qmc element
   * @return true, if a valid QMCDriver is set.
   */
  bool QMCMain::runQMC(xmlNodePtr cur) {

    bool append_run = setQMCDriver(myProject.m_series,cur);

    OHMMS::Controller->barrier();

    if(qmcDriver) {
      app_log() << endl;
      myProject.get(app_log());
      app_log() << endl;

      //advance the project id 
      //if it is NOT the first qmc node and qmc/@append!='yes'
      if(!FirstQMC && !append_run) myProject.advance();

      qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
      qmcDriver->putWalkers(m_walkerset_in);
      qmcDriver->process(cur);

      qmcDriver->run();

      //keeps track of the configuration file
      PrevConfigFile = myProject.CurrentRoot();

      //do not change the input href
      //change the content of mcwalkerset/@file attribute
      for(int i=0; i<m_walkerset.size(); i++) {
        xmlSetProp(m_walkerset[i],
            (const xmlChar*)"href", (const xmlChar*)myProject.CurrentRoot());
      }

      //curMethod = what;
      return true;
    } else {
      return false;
    }
  }

  bool QMCMain::setMCWalkers(xmlXPathContextPtr context_) {

    OhmmsXPathObject result("/simulation/mcwalkerset",context_);
    if(result.empty()) {
      if(m_walkerset.empty()) {
        result.put("//qmc",context_);
        xmlNodePtr newnode_ptr = xmlNewNode(NULL,(const xmlChar*)"mcwalkerset");
        if(result.empty()){
          xmlAddChild(XmlDocStack.top()->getRoot(),newnode_ptr);
        } else {
          xmlAddPrevSibling(result[0],newnode_ptr);
        }
        m_walkerset.push_back(newnode_ptr);
      }
    } else {
      for(int iconf=0; iconf<result.size(); iconf++) {
        xmlNodePtr mc_ptr = result[iconf];
        m_walkerset.push_back(mc_ptr);
        m_walkerset_in.push_back(mc_ptr);
      }
    }
    return true;
  }
}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

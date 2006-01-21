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
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCDrivers/DummyQMC.h"
#include "QMCDrivers/VMC.h"
#include "QMCDrivers/VMCParticleByParticle.h"
#include "QMCDrivers/DMCParticleByParticle.h"
#include "QMCDrivers/DMCPbyPOpenMP.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/MolecuDMC.h"
#if !defined(QMCPLUSPLUS_RELEASE)
#include "QMCDrivers/VMCMultiple.h"
#include "QMCDrivers/VMCPbyPMultiple.h"
#include "QMCDrivers/ReptationMC.h"
#include "QMCDrivers/RQMCMultiple.h"
#endif
#include "Utilities/OhmmsInfo.h"
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Message/Communicate.h"
#include "Particle/DistanceTable.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  QMCMain::QMCMain(int argc, char** argv): QMCAppBase(argc,argv), FirstQMC(true),
                                           qmcDriver(0),qmcSystem(0){ 


    //create ParticleSetPool
    ptclPool = new ParticleSetPool;

    //create WaveFunctionPool
    psiPool = new WaveFunctionPool;
    psiPool->setParticleSetPool(ptclPool);

    //create HamiltonianPool
    hamPool = new HamiltonianPool;
    hamPool->setParticleSetPool(ptclPool);
    hamPool->setWaveFunctionPool(psiPool);

    app_log() << "\n=========================================================\n"
              <<   "                   qmcplusplus 0.2                       \n"
              << "\n  (c) Copyright 2003-  qmcplusplus developers          \n"
              <<   "=========================================================\n";

    app_log().flush();
  }

  ///destructor
  QMCMain::~QMCMain() {
    DEBUGMSG("QMCMain::~QMCMain")
    delete hamPool;
    delete psiPool;
    delete ptclPool;
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

    OHMMS::Controller->barrier();

    //write stuff
    app_log() << "=========================================================\n";
    app_log() << " Summary of QMC systems \n";
    app_log() << "=========================================================\n";
    ptclPool->get(app_log());
    hamPool->get(app_log());

    curMethod = string("invalid");
    vector<xmlNodePtr> q;

    //xmlNodePtr cur=m_root->children;
    xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "qmc") {
        string target("e");
        const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"target");
        if(t) target = (const char*)t;
        qmcSystem = ptclPool->getWalkerSet(target);
        bool good = runQMC(cur);
        q.push_back(cur);

        //advance is done in runQMC
        //myProject.advance();
        //set to false
        FirstQMC=false;
      }

      cur=cur->next;
    }

    if(q.size()) { 
      xmlNodePtr newqmc_ptr=qmcDriver->getQMCNode();
      //unlink other qmc node but the last one
      for(int i=0; i<q.size(); i++) {
        xmlUnlinkNode(q[i]);
        xmlFreeNode(q[i]);
      }
      XmlDocStack.top()->addChild(newqmc_ptr);
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

    OHMMS::Controller->barrier();

    string what("invalid");
    string append_tag("no");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(what,"method");
    aAttrib.add(append_tag,"append");
    aAttrib.put(cur);

    bool append_run = (append_tag == "yes");

    if(qmcDriver) {
      if(what != curMethod) {
        delete qmcDriver;
        qmcDriver = 0;
        //if the current qmc method is different from the previous one, append_run is set to false
        append_run = false;
      }
    }

    if(myProject.m_series == 0) append_run = false;

    if(qmcDriver == 0) {
      ///////////////////////////////////////////////
      // get primaryPsi and primaryH
      ///////////////////////////////////////////////
      TrialWaveFunction* primaryPsi= 0;
      QMCHamiltonian* primaryH=0;
      queue<TrialWaveFunction*> targetPsi;//FIFO 
      queue<QMCHamiltonian*> targetH;//FIFO

      xmlNodePtr tcur=cur->children;
      while(tcur != NULL) {
        if(xmlStrEqual(tcur->name,(const xmlChar*)"qmcsystem")) {
          const xmlChar* t= xmlGetProp(tcur,(const xmlChar*)"wavefunction");
          targetPsi.push(psiPool->getWaveFunction((const char*)t));
          t= xmlGetProp(tcur,(const xmlChar*)"hamiltonian");
          targetH.push(hamPool->getHamiltonian((const char*)t));
        }
        tcur=tcur->next;
      }

      if(targetH.empty()) {
        primaryPsi=psiPool->getPrimary();
        primaryH=hamPool->getPrimary();
      } else { //mark the first targetPsi and targetH as the primaries
        primaryPsi=targetPsi.front(); targetPsi.pop();
        primaryH=targetH.front();targetH.pop();
      }
      //set primaryH->Primary
      primaryH->setPrimary(true);
      ///////////////////////////////////////////////
      if (what == "vmc"){
        primaryH->addOperator(new ConservedEnergy,"Flux");
        qmcDriver = new VMC(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = VMC_RUN;
      } else if(what == "vmc-ptcl"){
        primaryH->addOperator(new ConservedEnergy,"Flux");
        qmcDriver = new VMCParticleByParticle(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = VMC_RUN;
      } else if(what == "dmc"){
        MolecuDMC *dmc = new MolecuDMC(*qmcSystem,*primaryPsi,*primaryH);
        //dmc->setBranchInfo(PrevConfigFile);
        qmcDriver=dmc;
        curRunType = DMC_RUN;
      } else if(what == "dmc-ptcl"){
        DMCParticleByParticle *dmc = new DMCParticleByParticle(*qmcSystem,*primaryPsi,*primaryH);
        //dmc->setBranchInfo(PrevConfigFile);
        qmcDriver=dmc;
        curRunType = DMC_RUN;
      } else if(what == "optimize"){
        primaryH->remove("Flux");
        QMCOptimize *opt = new QMCOptimize(*qmcSystem,*primaryPsi,*primaryH);
        //opt->addConfiguration(PrevConfigFile);
        opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("null"));
        qmcDriver=opt;
        curRunType = OPTIMIZE_RUN;
#if !defined(QMCPLUSPLUS_RELEASE)
      } else if(what == "vmc-multi") {
        qmcDriver = new VMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = VMC_RUN;
      } else if(what == "vmc-ptcl-multi") {
        qmcDriver = new VMCPbyPMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = VMC_RUN;
      } else if(what == "rmc") {
        qmcDriver = new ReptationMC(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = RMC_RUN;
      } else if(what == "rmc-multi") {
        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = RMC_RUN;
#endif
      } else if(what == "dmc-omp") {
        DMCPbyPOpenMP *domp = new DMCPbyPOpenMP(*qmcSystem,*primaryPsi,*primaryH);
        domp->makeClones(*hamPool);
        qmcDriver=domp;
        curRunType = DUMMY_RUN; //change this to something else
      } else {
        qmcDriver = new DummyQMC(*qmcSystem,*primaryPsi,*primaryH);
        WARNMSG("Cannot termine what type of qmc to run. Creating DummyQMC for testing")
        curRunType = DUMMY_RUN;
      }
    }

    if(qmcDriver) {

      //advance the project id 
      //if it is NOT the first qmc node and qmc/@append!='yes'
      if(!FirstQMC && !append_run) myProject.advance();

      qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
      qmcDriver->putWalkers(m_walkerset_in);
      qmcDriver->process(cur);

      //set up barrier
      OHMMS::Controller->barrier();

      qmcDriver->run();

      //keeps track of the configuration file
      PrevConfigFile = myProject.CurrentRoot();

      //do not change the input href
      //change the content of mcwalkerset/@file attribute
      for(int i=0; i<m_walkerset.size(); i++) {
        xmlSetProp(m_walkerset[i],
            (const xmlChar*)"href", (const xmlChar*)myProject.CurrentRoot());
      }

      curMethod = what;
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

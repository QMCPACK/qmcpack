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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
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
#include "QMCDrivers/VMCMultiple.h"
#include "QMCDrivers/VMCPbyPMultiple.h"
#include "QMCDrivers/DMCParticleByParticle.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/MolecuDMC.h"
#include "QMCDrivers/ReptationMC.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Message/Communicate.h"
#include "Particle/DistanceTable.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace ohmmsqmc {

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
  }

  ///destructor
  QMCMain::~QMCMain() {
    DEBUGMSG("QMCMain::~QMCMain")
    delete hamPool;
    delete psiPool;
    delete ptclPool;
    //Not a good thing.
    //xmlFreeDoc(m_doc);
  }


  bool QMCMain::execute() {

    //validate the input file
    bool success = validateXML();

    if(!success) {
      ERRORMSG("Input document does not contain valid objects")
      return false;
    }

    //initialize all the instances of distance tables and evaluate them 
    ptclPool->reset();

    curMethod = string("invalid");
    vector<xmlNodePtr> q;

    xmlNodePtr cur=m_root->children;
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

      xmlAddChild(m_root,newqmc_ptr);
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

  /** validate m_doc
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

    xmlXPathContextPtr m_context = xmlXPathNewContext(m_doc);

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"//project",m_context);
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      WARNMSG("Project is not defined")
      myProject.reset();
    } else {
      myProject.put(result->nodesetval->nodeTab[0]);
    }
    xmlXPathFreeObject(result);

    //initialize the random number generator
    xmlNodePtr rptr = myRandomControl.initialize(m_context);

    //preserve the input order
    xmlNodePtr cur=m_root->children;
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
          //root and document are saved
          xmlNodePtr root_save=m_root;
          xmlDocPtr doc_save=m_doc;

          //parse a new file
          m_doc=NULL;
          parse((const char*)a);
          processPWH(xmlDocGetRootElement(m_doc));
          xmlFreeDoc(m_doc);

          //copy the original root and document
          m_root = root_save;
          m_doc = doc_save;
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
    xmlXPathFreeContext(m_context);
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

    string what("invalid");
    string continue_tag("no");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(what,"method");
    aAttrib.add(continue_tag,"continue");
    aAttrib.put(cur);

    bool continue_run = (continue_tag == "yes");

    if(qmcDriver) {
      if(what == curMethod) {
        LOGMSG("Reuse " << what << " driver")
      } else {
        delete qmcDriver;
        qmcDriver = 0;
        //if the current qmc method is different from the previous one, continue_run is set to false
        continue_run = false;
      }
    }

    if(myProject.m_series == 0) continue_run = false;

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
        primaryH->add(new ConservedEnergy,"Flux");
        qmcDriver = new VMC(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = VMC_RUN;
      } else if(what == "vmc-ptcl"){
        primaryH->add(new ConservedEnergy,"Flux");
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
        opt->addConfiguration(PrevConfigFile);
        opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("null"));
        qmcDriver=opt;
        curRunType = OPTIMIZE_RUN;
      } else if(what == "rmc") {
        qmcDriver = new ReptationMC(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = RMC_RUN;
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
      } else if(what == "rmc-multi") {
        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = RMC_RUN;
      } else {
        qmcDriver = new DummyQMC(*qmcSystem,*primaryPsi,*primaryH);
        WARNMSG("Cannot termine what type of qmc to run. Creating DummyQMC for testing")
          curRunType = DUMMY_RUN;
      }
    }

    if(qmcDriver) {

      //advance the project id 
      //if it is NOT the first qmc node and qmc/@continue!='yes'
      if(!FirstQMC && !continue_run) myProject.advance();

      if(continue_run) {
        LOGMSG("Continue a QMC simulation " << what << " File Root = " << myProject.CurrentRoot())
      } else {
        LOGMSG("Starting a QMC simulation " << what << " File Root = " << myProject.CurrentRoot())
      }
      qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, continue_run);
      qmcDriver->putWalkers(m_walkerset_in);
      qmcDriver->process(cur);
      qmcDriver->run();

      //keeps track of the configuration file
      PrevConfigFile = myProject.CurrentRoot();

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

    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"/simulation/mcwalkerset",context_);

    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      if(m_walkerset.empty()) {
        xmlXPathObjectPtr q = xmlXPathEvalExpression((const xmlChar*)"//qmc",context_);
	xmlNodePtr newnode_ptr = xmlNewNode(NULL,(const xmlChar*)"mcwalkerset");
        if(xmlXPathNodeSetIsEmpty(q->nodesetval)) {
	  xmlAddChild(m_root,newnode_ptr);
        } else {
	  xmlAddPrevSibling(q->nodesetval->nodeTab[0],newnode_ptr);
        }
	m_walkerset.push_back(newnode_ptr);
        xmlXPathFreeObject(q);
      }
    } else {
      for(int iconf=0; iconf<result->nodesetval->nodeNr; iconf++) {
	xmlNodePtr mc_ptr = result->nodesetval->nodeTab[iconf];
        m_walkerset.push_back(mc_ptr);
        m_walkerset_in.push_back(mc_ptr);
      }
    }
    xmlXPathFreeObject(result);
    return true;
  }
}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

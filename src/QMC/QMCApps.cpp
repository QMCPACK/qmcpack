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
#include "QMC/QMCApps.h"
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMC/VMC.h"
#include "QMC/VMCParticleByParticle.h"
#include "QMC/DMCParticleByParticle.h"
#include "QMC/VMC_OPT.h"
#include "QMC/MolecuDMC.h"
#include "QMC/ReptationMC.h"
#include "QMC/WaveFunctionTester.h"
#include "QMCHamiltonians/ConservedEnergy.h"

namespace ohmmsqmc {

  QMCApps::QMCApps(int argc, char** argv): 
    m_doc(NULL),m_root(NULL),m_context(NULL),PrevConfigFile("default")
  {
  }

  QMCApps::~QMCApps() {

    xmlXPathFreeContext(m_context);
    xmlFreeDoc(m_doc);

    DEBUGMSG("QMCApps::~QMCApps")
  }

  bool QMCApps::setMCWalkers(xmlNodePtr aroot) {

    bool foundconfig=false;
    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)"/simulation/mcwalkerset",m_context);

    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      if(m_walkerset.empty()) {
	xmlNodePtr newnode_ptr = xmlNewNode(NULL,(const xmlChar*)"mcwalkerset");
	xmlNewProp(newnode_ptr,(const xmlChar*)"file", (const xmlChar*)"nothing");
	xmlAddChild(m_root,newnode_ptr);
	m_walkerset.push_back(newnode_ptr);
      }
    } else {
      int pid=OHMMS::Controller->mycontext(); 
      for(int iconf=0; iconf<result->nodesetval->nodeNr; iconf++) {
	int pid_target=pid;
	xmlNodePtr mc_ptr = result->nodesetval->nodeTab[iconf];
	m_walkerset.push_back(mc_ptr);
	xmlChar* att=xmlGetProp(mc_ptr,(const xmlChar*)"file");
	xmlChar* anode=xmlGetProp(mc_ptr,(const xmlChar*)"node");
	if(anode) {
	  pid_target = atoi((const char*)anode);
	}
	if(pid_target == pid && att) {
	  string cfile((const char*)att);
	  XMLReport("Using previous configuration from " << cfile)

          HDFWalkerInput WO(cfile); 
	  //read only the last ensemble of walkers
	  WO.put(el,-1);

	  PrevConfigFile = cfile;
	  foundconfig=true;
	}
      }
    }
    xmlXPathFreeObject(result);
    return foundconfig;
  }

  /**Run a series of qmc simulations in succession
   *@param aroot the root xmlNode of an input file
   *@return true if successful
   *
   *A XPathContext object is created with the input document.
   *Before running a sequence of qmc runs, the basic objects 
   *are initialized.
   *Expecting mistakes by users, we use xpath to process the
   *basic objects instead of recursive processing as desired.
   *
   *The initialization steps process these elements sequentially.
   *- project : set the name and sequence number of the run
   *- random  : initialization of a random number generator
   *- particleset(s) : initialization of el and ion
   *- wavefunction   : initialization of many-body wavefunctions for el
   *- hamiltonian    : initialization of Hamiltonian
   *- mcwalkerset    : read the walker configuration from a previous run 
   *
   *When the elements above are missing, default values are assigned 
   *if possible.
   *Fatal errors occur when wavefunction is missing or empty because 
   *a requested wavefunction is not available.
   *
   *Once the basic objects are properly initialized, this functon process 
   *the qmc elements sequentially.
   *A run can have a mixture of any number of vmc, opt(itmization) and dmc.
   *After each qmc element, the sequence of the project is incremented and 
   *the file rootname will be set to reflect the change.
   */
  bool QMCApps::run(xmlNodePtr aroot) {
    if(aroot) {
      m_doc=aroot->doc;
      m_root =  aroot;
      m_context = xmlXPathNewContext(m_doc);

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
      if(rptr) {
        xmlAddChild(m_root,rptr);
      }
   
      //call virtual function to initialize the physical entities, such as el.
      if(!init()) return false;
    
      //perform a series of QMC runs
      string methodname("invalid");
      result = xmlXPathEvalExpression((const xmlChar*)"//qmc",m_context);
      if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
	ERRORMSG("QMC Method is not defined." << endl)
        return false;
      } else {
	for(int i=0; i< result->nodesetval->nodeNr; i++) {
	  xmlNodePtr cur=result->nodesetval->nodeTab[i];
	  xmlChar* att=xmlGetProp(cur,(const xmlChar*)"method");
	  if(att) {
	    methodname = (const char*)att;
	  }
	  LOGMSG("Starting a QMC simulation " << methodname)
	  if(methodname == "vmc"){
	    //always check the conserved quantity for vmc
	    XMLReport("Checking the conserved quantity")
	    H.add(new ConservedEnergy,"Flux");
	    VMC vmc(el,Psi,H,cur);
	    vmc.setFileRoot(myProject.CurrentRoot());
	    vmc.run();
	  } else if(methodname == "vmc-ptcl"){
	    //always check the conserved quantity for vmc
	    XMLReport("Checking the conserved quantity")
	    H.add(new ConservedEnergy,"Flux");
	    VMCParticleByParticle vmc(el,Psi,H,cur);
	    vmc.setFileRoot(myProject.CurrentRoot());
	    vmc.run();
	  } else if(methodname == "dmc"){
	    //H.remove("Flux");
	    MolecuDMC dmc(el,Psi,H,cur);
	    dmc.setBranchInfo(PrevConfigFile);
	    dmc.setFileRoot(myProject.CurrentRoot());
	    dmc.run();
	  } else if(methodname == "dmc-ptcl"){
	    DMCParticleByParticle dmc(el,Psi,H,cur);
	    dmc.setBranchInfo(PrevConfigFile);
	    dmc.setFileRoot(myProject.CurrentRoot());
	    dmc.run();
	  } else if(methodname == "optimize"){
	    //Something wrong this! Check this out.
	    //H.remove("Flux");
	    VMC_OPT vmc(el,Psi,H,cur);
	    vmc.addConfiguration(PrevConfigFile);
	    vmc.setFileRoot(myProject.CurrentRoot());
	    vmc.run();
	  } else if(methodname == "rmc") {
	    XMLReport("Running ReptationMC")
	    ReptationMC rmc(el,Psi,H,cur);
	    rmc.setFileRoot(myProject.CurrentRoot());
	    rmc.run();
	  } else if(methodname == "test"){
	    H.remove("Flux");
	    WaveFunctionTester wftest(el,Psi,H,cur);
	    wftest.setFileRoot(myProject.CurrentRoot());
	    wftest.run();
	  } else{
	    ERRORMSG("Options are vmc, dmc or opt")
	  }

          //keeps track of the configuration file
          PrevConfigFile = myProject.CurrentRoot();
	  //change the content of mcwalkerset/@file attribute
	  for(int i=0; i<m_walkerset.size(); i++) {
	    xmlSetProp(m_walkerset[i],(const xmlChar*)"file", 
                       (const xmlChar*)myProject.CurrentRoot());
	  }

	  myProject.advance();
	  //remove this node
	  xmlUnlinkNode(cur);
	  xmlFreeNode(cur);
	}
      }
    } else {
      WARNMSG("No QMC parameters are provided. Peform WaveFunctionTest")
      WaveFunctionTester wftest(el,Psi,H,NULL);
      wftest.setFileRoot(myProject.CurrentRoot());
      wftest.run();
    }

    xmlNodePtr aqmc = xmlNewNode(NULL,(const xmlChar*)"qmc");
    xmlNewProp(aqmc,(const xmlChar*)"method",(const xmlChar*)"dmc");
    xmlNodePtr aparam = xmlNewNode(NULL,(const xmlChar*)"parameter");
    xmlNewProp(aparam,(const xmlChar*)"name",(const xmlChar*)"en_ref");
    char ref_energy[128];
    sprintf(ref_energy,"%15.5e",el.getLocalEnergy());
    xmlNodeSetContent(aparam,(const xmlChar*)ref_energy);
    xmlAddChild(aqmc,aparam);
    xmlAddChild(m_root,aqmc);

    string newxml(myProject.CurrentRoot());
    newxml.append(".cont.xml");
    xmlSaveFormatFile(newxml.c_str(),m_doc,1);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

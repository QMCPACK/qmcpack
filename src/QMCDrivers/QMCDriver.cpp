/////////////////////////////////////////////////////////////////
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
#include "QMCDrivers/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <limits>

namespace qmcplusplus {

  QMCDriver::QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    branchEngine(0), ResetRandom(false), AppendRun(false),
    MyCounter(0), RollBackBlocks(0),
    Period4CheckPoint(1), Period4WalkerDump(0),
    CurrentStep(0), nBlocks(100), nSteps(10), 
    nAccept(0), nReject(0), nTargetWalkers(0),
    Tau(0.01), qmcNode(NULL),
    QMCType("invalid"), 
    qmcComm(0), wOut(0),
    W(w), Psi(psi), H(h), Estimators(0)
  { 

    //use maximum double
    MaxCPUSecs=numeric_limits<RealType>::max();

    m_param.add(nSteps,"steps","int");
    m_param.add(nBlocks,"blocks","int");
    m_param.add(nTargetWalkers,"walkers","int");
    m_param.add(CurrentStep,"current","int");
    m_param.add(Tau,"timeStep","AU");
    m_param.add(Tau,"Tau","AU");
    m_param.add(Tau,"timestep","AU");
    m_param.add(RollBackBlocks,"rewind","int");
    m_param.add(Period4WalkerDump,"recordWalkers","int");
    m_param.add(MaxCPUSecs,"maxcpusecs","real");

    //add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
    H.add2WalkerProperty(W);
  }

  void QMCDriver::setCommunicator(Communicate* c) 
  {
    qmcComm = c ? c : OHMMS::Controller;
  }

  QMCDriver::~QMCDriver() 
  { 
  }

  void QMCDriver::add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi) {
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
  void QMCDriver::process(xmlNodePtr cur) {

    deltaR.resize(W.getTotalNum());
    drift.resize(W.getTotalNum());

    qmcNode=cur;

    putQMCInfo(cur);

    //need to initialize properties
    int numCopies= (H1.empty())?1:H1.size();
    W.resetWalkerProperty(numCopies);

    if(branchEngine==0) 
      branchEngine = new BranchEngineType(Tau,W.getActiveWalkers());

    //execute the put function implemented by the derived classes
    put(cur);

    if(h5FileRoot.size() && RollBackBlocks>1) {
      HDFWalkerInputManager W_in(W);
      W_in.rewind(h5FileRoot,RollBackBlocks);
      RollBackBlocks=0;
    }

    //create and initialize estimator
    Estimators = branchEngine->getEstimatorManager();
    if(Estimators==0)
    {
      Estimators = new EstimatorManager(qmcComm);
      branchEngine->setEstimatorManager(Estimators);
      if(h5FileRoot.size()) branchEngine->read(h5FileRoot);
    }

    branchEngine->put(cur);
    Estimators->put(W,H,cur);

    if(wOut==0) {
      wOut = new HDFWalkerOutput(W,RootName,qmcComm);
      wOut->open();
      branchEngine->write(wOut->getConfigID(),false);
      wOut->close();
      branchEngine->start(RootName,true);
    }
    else
      branchEngine->start(RootName,false);

    //use new random seeds
    if(ResetRandom) {
      app_log() << "  Regenerate random seeds." << endl;
      Random.reset();
      ResetRandom=false;
    }

    //flush the ostreams
    OhmmsInfo::flush();


    //increment QMCCounter of the branch engine
    branchEngine->advanceQMCCounter();
  }

  void QMCDriver::setStatus(const string& aname, const string& h5name, bool append) {
    RootName = aname;
    app_log() << "\n=========================================================" 
              << "\n  Start " << QMCType 
              << "\n  File Root " << RootName;
    if(append) 
      app_log() << " append = yes ";
    else 
      app_log() << " append = no ";
    app_log() << "\n=========================================================" << endl;

    if(h5name.size()) h5FileRoot = h5name;
    AppendRun = append;
  }

  
  /** Read walker configurations from *.config.h5 files
   * @param wset list of xml elements containing mcwalkerset
   */
  void QMCDriver::putWalkers(vector<xmlNodePtr>& wset) {

    if(wset.empty()) return;
    int nfile=wset.size();

    //if(QMCDriverMode[QMC_OPTIMIZE]) 
    //{//for optimization, simply add to ConfigFile
    //  for(int ifile=0; ifile<nfile; ifile++) 
    //    mcwalkerNodePtr.push_back(wset[ifile]);
    //} 
    //else 
    //{
      HDFWalkerInputManager W_in(W);
      if(W_in.put(wset,qmcComm->mycontext())) 
        h5FileRoot = W_in.getLastFile();
    //}

    //clear the walker set
    wset.clear();
  }

  void QMCDriver::recordBlock(int block) {

    if(wOut ==0)
    {
      wOut = new HDFWalkerOutput(W,RootName,qmcComm);
      wOut->open();
      branchEngine->write(wOut->getConfigID(),false);
      wOut->close();
      branchEngine->start(RootName,true);
    }
    //estimator writes
    //Estimators->report(CurrentStep);
    //if Period4WalkerDump>0, record works as the checkpoint
    
    wOut->open();
    if(Period4WalkerDump>0) 
      wOut->append(W,block-1); //block has been incremened already
    else 
      wOut->dump(W);
    const bool overwrite=true;
    branchEngine->write(wOut->getConfigID(),overwrite);
    wOut->close();

    //flush the ostream
    OhmmsInfo::flush();
  }

  bool QMCDriver::finalize(int block) {

    branchEngine->finalize();

    wOut->open();
    const bool overwrite=true;
    branchEngine->write(wOut->getConfigID(),overwrite);
    wOut->close();
    
    delete wOut;
    wOut=0;
    //Estimators->finalize();

    //set the target walkers
    nTargetWalkers = W.getActiveWalkers();

    //increment MyCounter
    MyCounter++;

    //flush the ostream
    OhmmsInfo::flush();

    return true;
  }

  /** Add walkers to the end of the ensemble of walkers.  
   * @param nwalkers number of walkers to add
   */
  void 
  QMCDriver::addWalkers(int nwalkers) {

    if(nwalkers>0) {
      //add nwalkers walkers to the end of the ensemble
      int nold = W.getActiveWalkers();

      app_log() << "  Adding " << nwalkers << " walkers to " << nold << " existing sets" << endl;

      W.createWalkers(nwalkers);
      
      ParticleSet::ParticlePos_t rv(W.getTotalNum());
      MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
      while(it != it_end) {
	(*it)->R=W.R;
	++it;
      }
    } else if(nwalkers<0) {
      W.destroyWalkers(-nwalkers);
      app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers() << endl;
    } else {
      app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." <<  endl;
    }

    //update the global number of walkers
    //int nw=W.getActiveWalkers();
    //qmcComm->allreduce(nw);
    vector<int> nw(qmcComm->ncontexts(),0),nwoff(qmcComm->ncontexts()+1,0);
    nw[qmcComm->mycontext()]=W.getActiveWalkers();
    qmcComm->allreduce(nw);

    for(int ip=0; ip<qmcComm->ncontexts(); ip++) nwoff[ip+1]=nwoff[ip]+nw[ip];
    W.setGlobalNumWalkers(nwoff[qmcComm->ncontexts()]);
    W.setWalkerOffsets(nwoff);

    app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples  <<  endl;
    app_log() << "  Total weight: " << W.EnsembleProperty.Weight  <<  endl;
  }

  
  /** Parses the xml input file for parameter definitions for a single qmc simulation.
   */
  bool 
  QMCDriver::putQMCInfo(xmlNodePtr cur) {
    
    int defaultw = 100;
    int targetw = 0;

    Period4WalkerDump=0;
    Period4CheckPoint=1;
     
    if(cur) {
      //initialize the parameter set
      m_param.put(cur);

      xmlNodePtr tcur=cur->children;
      //determine how often to print walkers to hdf5 file
      while(tcur != NULL) {
	string cname((const char*)(tcur->name));
	if(cname == "record") {
          const xmlChar* aptr=xmlGetProp(tcur,(const xmlChar*)"stride");
          if(aptr) Period4WalkerDump = atoi((const char*)aptr);
          aptr=xmlGetProp(tcur,(const xmlChar*)"period");
          if(aptr) Period4WalkerDump = atoi((const char*)aptr);
	} else if(cname == "checkpoint") {
          const xmlChar* aptr=xmlGetProp(tcur,(const xmlChar*)"period");
          if(aptr) Period4CheckPoint = atoi((const char*)aptr);
        } else if(cname == "random") {
          ResetRandom = true;
        }
	tcur=tcur->next;
      }
    }
    
    //reset CurrentStep to zero if qmc/@continue='no'
    if(!AppendRun) CurrentStep=0;

    app_log() << "  timestep = " << Tau << endl;
    app_log() << "  blocks = " << nBlocks << endl;
    app_log() << "  steps = " << nSteps << endl;
    app_log() << "  current = " << CurrentStep << endl;

    //Need MPI-IO
    app_log() << "  walkers = " << W.getActiveWalkers() << endl;

    /*check to see if the target population is different 
      from the current population.*/ 
    int nw  = W.getActiveWalkers();
    int ndiff = 0;
    if(nw) { // walkers exist
      // nTargetWalkers == 0, if it is not set by the input file
      ndiff = (nTargetWalkers)? nTargetWalkers-nw: 0;
    } else {
      ndiff= (nTargetWalkers)? nTargetWalkers:defaultw;
    }

    addWalkers(ndiff);

    //always true
    return (W.getActiveWalkers()>0);
  }

  xmlNodePtr QMCDriver::getQMCNode() {

    xmlNodePtr newqmc = xmlCopyNode(qmcNode,1);
    //update current
    xmlNodePtr current_ptr=NULL;
    xmlNodePtr cur=newqmc->children;
    while(cur != NULL && current_ptr == NULL) {
      string cname((const char*)(cur->name)); 
      if(cname == "parameter") {
        const xmlChar* aptr= xmlGetProp(cur, (const xmlChar *) "name");
        if(aptr) {
          if(xmlStrEqual(aptr,(const xmlChar*)"current")) current_ptr=cur;
        }
      }
      cur=cur->next;
    }
    if(current_ptr == NULL) {
      current_ptr = xmlNewTextChild(newqmc,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0");
      xmlNewProp(current_ptr,(const xmlChar*)"name",(const xmlChar*)"current");
    } 
    getContent(CurrentStep,current_ptr);
    return newqmc;
  }

}


/***************************************************************************
 * $RCSfile: QMCDriver.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

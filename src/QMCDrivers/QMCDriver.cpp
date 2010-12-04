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
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "HDFVersion.h"
#include <limits>

namespace qmcplusplus {

  QMCDriver::QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool): MPIObjectBase(0),
  branchEngine(0), ResetRandom(false), AppendRun(false), DumpConfig(true),
  MyCounter(0), RollBackBlocks(0),
  Period4CheckPoint(0), Period4WalkerDump(10),Period4ConfigDump(50),
  Period4CheckProperties(100), 
  nBlocks(10), nSteps(0), nTargetWalkers(0), nTargetSamples(0),
  nStepsBetweenSamples(0), nSamplesPerThread(0),
  nAccept(0), nReject(0),  CurrentStep(0), 
  Tau(0.01), qmcNode(NULL),
  QMCType("invalid"), wOut(0), storeConfigs(0),
  W(w), Psi(psi), H(h), psiPool(ppool), Estimators(0)
  { 

    //use maximum double
    MaxCPUSecs=numeric_limits<RealType>::max();

    m_param.add(nSteps,"steps","int");
    m_param.add(nBlocks,"blocks","int");
    m_param.add(nTargetWalkers,"walkers","int");
    m_param.add(CurrentStep,"current","int");
    m_param.add(Tau,"timeStep","AU"); m_param.add(Tau,"timestep","AU"); m_param.add(Tau,"time_step","AU");
    m_param.add(Tau,"Tau","AU"); m_param.add(Tau,"tau","AU");
    m_param.add(RollBackBlocks,"rewind","int");
    m_param.add(MaxCPUSecs,"maxcpusecs","real");
    m_param.add(nTargetSamples,"samples","int");
    m_param.add(Period4WalkerDump,"recordWalkers","int"); m_param.add(Period4WalkerDump,"record_walkers","int"); m_param.add(Period4WalkerDump,"recordwalkers","int");
    m_param.add(Period4ConfigDump,"recordConfigs","int"); m_param.add(Period4ConfigDump,"recordconfigs","int"); m_param.add(Period4ConfigDump,"record_configs","int");
    m_param.add(Period4CheckProperties,"checkProperties","int"); m_param.add(Period4CheckProperties,"checkproperties","int"); m_param.add(Period4CheckProperties,"check_properties","int");
    m_param.add(storeConfigs,"storeConfigs","int"); m_param.add( storeConfigs,"storeconfigs","int"); m_param.add( storeConfigs,"store_configs","int");

    m_param.add(nSamplesPerThread,"samplesperthread","real");
    m_param.add(nStepsBetweenSamples,"stepsbetweensamples","int");
    
    
    ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
    //H.add2WalkerProperty(W);
    if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
  }

  QMCDriver::~QMCDriver() 
  { 
    delete_iter(Rng.begin(),Rng.end());
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

    //process common parameters
    putQMCInfo(cur);
    ////set the Tau parameter inside the Hamiltonian
    //H.setTau(Tau);

    //need to initialize properties
    int numCopies= (H1.empty())?1:H1.size();
    W.resetWalkerProperty(numCopies);

    //create branchEngine first
    if(branchEngine==0) 
      branchEngine = new BranchEngineType(Tau,W.getGlobalNumWalkers());

    //execute the put function implemented by the derived classes
    put(cur);

    //create and initialize estimator
    Estimators = branchEngine->getEstimatorManager();
    if(Estimators==0)
    {
      Estimators = new EstimatorManager(myComm);
      branchEngine->setEstimatorManager(Estimators);
      branchEngine->read(h5FileRoot);
    }

    branchEngine->put(cur);
    Estimators->put(W,H,cur);

    if(wOut==0) wOut = new HDFWalkerOutput(W,RootName,myComm);
    branchEngine->start(RootName);
    branchEngine->write(RootName);

    //use new random seeds
    if(ResetRandom) {
      app_log() << "  Regenerate random seeds." << endl;
      RandomNumberControl::make_seeds();
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

    HDFWalkerInputManager W_in(W,myComm);
    for(int i=0; i<wset.size(); i++)
      if(W_in.put(wset[i])) h5FileRoot = W_in.getFileRoot();

    //clear the walker set
    wset.clear();
  }

  void QMCDriver::recordBlock(int block) {

    ////first dump the data for restart
    if(DumpConfig &&block%Period4CheckPoint == 0)
    {
      wOut->dump(W);
      branchEngine->write(RootName,true); //save energy_history
      if (storeConfigs) wOut->dump( ForwardWalkingHistory);
    }

    //save positions for optimization: this is done within VMC
    //if(QMCDriverMode[QMC_OPTIMIZE]) W.saveEnsemble();
    //if(Period4WalkerDump>0) wOut->append(W);

    //flush the ostream
    //OhmmsInfo::flush();
  }

  bool QMCDriver::finalize(int block) {

    TimerManager.print(myComm);
    TimerManager.reset();

    if(DumpConfig)
    {
      wOut->dump(W);
      branchEngine->finalize(W);
      RandomNumberControl::write(RootName,myComm);
    }

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
      if(nold)
      {
        int iw=nold;
        for(MCWalkerConfiguration::iterator it=W.begin()+nold; it != W.end(); ++it,++iw)
          (*it)->R=W[iw%nold]->R;//assign existing walker configurations when the number of walkers change
      }

    } else if(nwalkers<0) {
      W.destroyWalkers(-nwalkers);
      app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers() << endl;
    } else {
      app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." <<  endl;
    }

    //update the global number of walkers
    //int nw=W.getActiveWalkers();
    //myComm->allreduce(nw);
    vector<int> nw(myComm->size(),0),nwoff(myComm->size()+1,0);
    nw[myComm->rank()]=W.getActiveWalkers();
    myComm->allreduce(nw);

    for(int ip=0; ip<myComm->size(); ip++) nwoff[ip+1]=nwoff[ip]+nw[ip];
    W.setGlobalNumWalkers(nwoff[myComm->size()]);
    W.setWalkerOffsets(nwoff);

    app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples  <<  endl;
    app_log() << "  Total weight: " << W.EnsembleProperty.Weight  <<  endl;
  }

  
  /** Parses the xml input file for parameter definitions for a single qmc simulation.
   */
  bool 
  QMCDriver::putQMCInfo(xmlNodePtr cur) {
    
    //SpeciesSet tspecies(W.getSpeciesSet());
    //RealType mass = tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]));
    //if (mass < 1e-12) {
    //  mass=1.0;
    //  tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]))=1.0;
    //}
    //oneovermass = 1.0/mass;

    //set the default walker to the number of threads times 10
    int defaultw = omp_get_max_threads();
    int targetw = 0;

    //these options are reset for each block
    Period4WalkerDump=10;
    Period4CheckPoint=0;

    OhmmsAttributeSet aAttrib;
    aAttrib.add(Period4CheckPoint,"checkpoint");
    aAttrib.put(cur);
     
    if(cur != NULL) {
      //initialize the parameter set
      m_param.put(cur);
      xmlNodePtr tcur=cur->children;
      
      //determine how often to print walkers to hdf5 file
      while(tcur != NULL) {
	string cname((const char*)(tcur->name));
	if(cname == "record") {
          //dump walkers for optimization
          OhmmsAttributeSet rAttrib;
          rAttrib.add(Period4WalkerDump,"stride");
          rAttrib.add(Period4WalkerDump,"period");
          rAttrib.put(tcur);
	} else if(cname == "checkpoint") {
          OhmmsAttributeSet rAttrib;
          rAttrib.add(Period4CheckPoint,"stride");
          rAttrib.add(Period4CheckPoint,"period");
          rAttrib.put(tcur);
          DumpConfig=(Period4CheckPoint>0);
        }
        else if(cname == "dumpconfig") {
          OhmmsAttributeSet rAttrib; 
          rAttrib.add(Period4ConfigDump,"stride");
          rAttrib.add(Period4ConfigDump,"period");
          rAttrib.put(tcur);
        }
        else if(cname == "random") {
          ResetRandom = true;
        }
	tcur=tcur->next;
      }
    }

    if(Period4CheckPoint==0)  Period4CheckPoint=(nBlocks+1)*nSteps;
    
    int Nthreads = omp_get_max_threads();
    int Nprocs=myComm->size();
    //if target is not give, use whatever it has
     if(nTargetWalkers==0) nTargetWalkers=W.getActiveWalkers();

    //nTargetWalkers is a local quantity.
    nTargetWalkers=std::max(Nthreads,nTargetWalkers);
//     nTargetSamples is set to 
    nTargetSamples=std::max(int(nTargetWalkers*Nprocs*nSamplesPerThread),nTargetSamples);
    
    if(nStepsBetweenSamples)
    {
      int nStepsTotal = (int)std::ceil(RealType(nTargetSamples*nStepsBetweenSamples)/RealType(nTargetWalkers*Nprocs) );
      if (nBlocks<1) nBlocks=1;
      if (nStepsTotal<nBlocks) nBlocks=nStepsTotal;
      nSteps = (int)std::ceil(RealType(nStepsTotal)/RealType(nBlocks));
      nStepsTotal = nSteps*nBlocks;
      nStepsBetweenSamples = (int)std::floor(RealType(nStepsTotal*nTargetWalkers*Nprocs)/RealType(nTargetSamples));
      Period4WalkerDump = nStepsBetweenSamples;
    }
    else if (nTargetSamples)
    {
      int nStepsTotal =  nSteps*nBlocks;
      nStepsBetweenSamples = (int)std::floor(RealType(nStepsTotal*nTargetWalkers*Nprocs)/RealType(nTargetSamples));
      Period4WalkerDump = nStepsBetweenSamples;
    }
    else 
      Period4WalkerDump=(nBlocks+1)*nSteps;


    app_log() << "  Walker Check Points are dumped every " << Period4CheckPoint << " steps." << endl;
    if (Period4WalkerDump>0) app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << endl;
    //reset CurrentStep to zero if qmc/@continue='no'
    if(!AppendRun) CurrentStep=0;

    //target number of walkers is less than the number of threads. Reset it.
    //if(nTargetWalkers && nTargetWalkers<omp_get_max_threads()) 
    //  nTargetWalkers=omp_get_max_threads();

    app_log() << "  timestep = " << Tau << endl;
    app_log() << "  blocks = " << nBlocks << endl;
    app_log() << "  steps = " << nSteps << endl;
//     app_log() << "  mass = " << mass << endl;
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

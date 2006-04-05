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
#include "Message/Communicate.h"

namespace qmcplusplus {

  //initialize the static data
  int QMCDriver::Counter = -1;
  QMCDriver::BranchEngineType* QMCDriver::branchEngine=0;
  
  QMCDriver::QMCDriver(MCWalkerConfiguration& w, 
		       TrialWaveFunction& psi, 
		       QMCHamiltonian& h):
    ResetRandom(false), AppendRun(false),RollBackBlocks(0),
    Period4CheckPoint(1), Period4WalkerDump(0),
    CurrentStep(0), nBlocks(100), nSteps(1000), 
    nAccept(0), nReject(0), nTargetWalkers(0),
    Tau(0.001), qmcNode(NULL),
    QMCType("invalid"), h5FileRoot("invalid"),
    W(w), Psi(psi), H(h), Estimators(0)
  { 
    m_param.add(nSteps,"steps","int");
    m_param.add(nBlocks,"blocks","int");
    m_param.add(nTargetWalkers,"walkers","int");
    m_param.add(CurrentStep,"current","int");
    m_param.add(Tau,"timeStep","AU");
    m_param.add(Tau,"Tau","AU");
    m_param.add(Tau,"timestep","AU");
    m_param.add(RollBackBlocks,"rewind","int");

    //add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
    H.add2WalkerProperty(W);
  }

  QMCDriver::~QMCDriver() { 
    
    if(Estimators) {
      if(Estimators->size()) W.setLocalEnergy(Estimators->average(0));
      delete Estimators;
    }
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
   *   -- putQMCInfo: <parameter/> s for QMC
   *   -- put : extra data by derived classes
   * - initialize branchEngine to accumulate energies
   * - initialize Estimators
   * - initialize Walkers
   * The virtual function put(xmlNodePtr cur) is where QMC algorithm-dependent
   * data are registered and initialized.
   */
  void QMCDriver::process(xmlNodePtr cur) {

    deltaR.resize(W.getTotalNum());
    drift.resize(W.getTotalNum());

    qmcNode=cur;
    putQMCInfo(qmcNode);

    bool firstTime = (branchEngine == 0);
    if(firstTime) {
      branchEngine = new BranchEngineType(Tau,W.getActiveWalkers());
    }

    put(qmcNode);

    if(h5FileRoot.size() && RollBackBlocks>1) {
      HDFWalkerInputManager W_in(W);
      W_in.rewind(h5FileRoot,RollBackBlocks);
      RollBackBlocks=0;
    }

    branchEngine->put(qmcNode);

    if(firstTime && h5FileRoot.size() && h5FileRoot != "invalid") {
      app_log() << "  Initializing BranchEngine with " << h5FileRoot << endl;
      branchEngine->read(h5FileRoot);
    }
    
    //A new run, branchEngine needs to be flushed
    if(!AppendRun) branchEngine->flush(0);

    //create estimator if not allocated
    if(Estimators == 0) Estimators =new ScalarEstimatorManager(H);

    //reset the Properties of all the walkers
    int numCopies= (H1.empty())?1:H1.size();
    W.resetWalkerProperty(numCopies);

    Estimators->put(qmcNode);

    //set the collection mode
    Estimators->setPeriod(nSteps);
    Estimators->setCollectionMode(branchEngine->SwapMode);
    Estimators->resetReportSettings(RootName, AppendRun);

    //initialize the walkers before moving: can be moved to run
    initialize();

    //use new random seeds
    if(ResetRandom) {
      app_log() << "  Regenerate random seeds." << endl;
      Random.reset();
      ResetRandom=false;
    }

    //flush the ostreams
    OhmmsInfo::flush();

    Counter++; 
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

    if(QMCDriverMode[QMC_OPTIMIZE]) {//for optimization, simply add to ConfigFile
      for(int ifile=0; ifile<nfile; ifile++) 
        mcwalkerNodePtr.push_back(wset[ifile]);
    } else {
      HDFWalkerInputManager W_in(W);
      if(W_in.put(wset)) h5FileRoot = W_in.getLastFile();
      //int pid=OHMMS::Controller->mycontext(); 
      //for(int ifile=0; ifile<nfile; ifile++) {
      //  string cfile("invalid"), target("e");
      //  int anode=-1, nwalkers=-1;
      //  OhmmsAttributeSet pAttrib;
      //  pAttrib.add(cfile,"href"); pAttrib.add(cfile,"file"); 
      //  pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 
      //  pAttrib.add(anode,"node");
      //  pAttrib.add(nwalkers,"walkers");
      //  pAttrib.put(wset[ifile]);
      //  int pid_target= (anode<0)?pid:anode;
      //  if(pid_target == pid && cfile != "invalid") {
      //    XMLReport("Using previous configuration of " << target << " from " << cfile)
      //    HDFWalkerInput WO(cfile); 
      //    WO.append(W,nwalkers);
      //    //read random state
      //    WO.getRandomState(true);
      //    h5FileRoot = cfile;
      //  }
      //}
    }

    //clear the walker set
    wset.clear();
  }

  /** Initialize QMCDriver
   *
   * Evaluate the Properties of Walkers when a QMC starts
   */
  void QMCDriver::initialize() {

    //For optimization, do nothing
    if(QMCDriverMode[QMC_OPTIMIZE]) return;

    //For multiple+particle-by-particle, do nothing
    if(QMCDriverMode[QMC_MULTIPLE] && QMCDriverMode[QMC_UPDATE_MODE]) return;

    if(QMCDriverMode[QMC_UPDATE_MODE]) { //using particle-by-particle moves
      bool require_register =  W.createAuxDataSet();
      MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
      if(require_register) {
        while(it != it_end) {
          (*it)->DataSet.rewind();
          W.registerData(**it,(*it)->DataSet);
          ValueType logpsi=Psi.registerData(W,(*it)->DataSet);

          RealType scale=getDriftScale(Tau,W.G);
          (*it)->Drift = scale*W.G;

          RealType ene = H.evaluate(W);
          (*it)->resetProperty(logpsi,Psi.getSign(),ene);
          H.saveProperty((*it)->getPropertyBase());
          ++it;
        } 
      } else {
        updateWalkers(); // simply re-evaluate the values 
      }
    } else { // using walker-by-walker moves

      app_log() << "  Evaluate all the walkers before starting for walker-by-walker update" << endl;

      MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
      while(it != it_end) {
        W.R = (*it)->R;
        //DistanceTable::update(W);
        W.update();

        ValueType logpsi(Psi.evaluateLog(W));
        RealType scale=getDriftScale(Tau,W.G);
        (*it)->Drift = scale*W.G;
        RealType ene = H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getSign(),ene);
        H.saveProperty((*it)->getPropertyBase());
        ++it;
      }
    }

  }

  /** Update walkers
   *
   * Evaluate the properties of all the walkers and update anonyous
   * uffers. Used by particle-by-particle updates.
   */
  void QMCDriver::updateWalkers() {

    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    while(it != it_end) {
      Walker_t& thisWalker(**it);
      Buffer_t& w_buffer(thisWalker.DataSet);
      w_buffer.rewind();
      W.updateBuffer(thisWalker,w_buffer);
      ValueType logpsi=Psi.updateBuffer(W,w_buffer);
      RealType enew= H.evaluate(W);
      thisWalker.resetProperty(logpsi,Psi.getSign(),enew);
      H.saveProperty(thisWalker.getPropertyBase());
      ValueType vsq = Dot(W.G,W.G);
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      thisWalker.Drift = scale*W.G;
      ++it;
    }
  }

  void QMCDriver::recordBlock(int block) {

    //estimator writes
    Estimators->report(CurrentStep);

    //if Period4WalkerDump>0, record works as the checkpoint
    if(Period4WalkerDump>0) {
      if(block%Period4WalkerDump == 0) {
        int now=block/Period4WalkerDump-1;
        bool appendWalker= AppendRun || now>0;
        //HDFWalkerOutput WO(RootName,appendWalker, block-1);
        HDFWalkerOutput WO(RootName,appendWalker, now);
        WO.get(W);
        WO.write(*branchEngine);
      }
    } else {
      if(block%Period4CheckPoint == 0) {
        HDFWalkerOutput WO(RootName,false,0);
        WO.get(W);
        WO.write(*branchEngine);
      }
    }

    //flush the ostream
    OhmmsInfo::flush();
  }

  bool QMCDriver::finalize(int block) {

    branchEngine->update(W.getActiveWalkers(), Estimators->average(0));

    int nconf= (Period4WalkerDump>0) ? block/Period4WalkerDump:1;
    HDFWalkerOutput WOextra(RootName,true,nconf);
    WOextra.write(*branchEngine);

    Estimators->finalize();

    //flush the ostream
    OhmmsInfo::flush();

    return true;
  }

  /** Add walkers to the end of the ensemble of walkers.  
   *@param nwalkers number of walkers to add
   *@return true, if the walker configuration is not empty.
   *
   * Assign positions to any new 
   * walkers \f[ {\bf R}[i] = {\bf R}[i-1] + g{\bf \chi}, \f]
   * where \f$ g \f$ is a constant and \f$ {\bf \chi} \f$
   * is a 3N-dimensional gaussian.
   * As a last step, for each walker calculate 
   * the properties given the new configuration
   <ul>
   <li> Local Energy \f$ E_L({\bf R} \f$
   <li> wavefunction \f$ \Psi({\bf R}) \f$
   <li> wavefunction squared \f$ \Psi^2({\bf R}) \f$
   <li> weight\f$ w({\bf R}) = 1.0\f$  
   <li> drift velocity \f$ {\bf v_{drift}}({\bf R})) \f$
   </ul>
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
    } else {
      app_log() << "  Using the existing " << W.getActiveWalkers() << " walkers" << endl;
    }

  }

  
  /** Parses the xml input file for parameter definitions for a single qmc simulation.
   * \param cur current xmlNode
   *
   Available parameters added to the ParameterSeet
   <ul>
   <li> blocks: number of blocks, default 100
   <li> steps: number of steps per block, default 1000
   <li> walkers: target population of walkers, default 100
   <li> Tau: the timestep, default 0.001
   <li> stride: flag for printing the ensemble of walkers,  default false
   <ul>
   <li> true: print every block
   <li> false: print at the end of the run
   </ul>
   </ul>
   In addition, sets the stride for the scalar estimators
   such that the scalar estimators flush and print to
   file every block and calls the function to initialize
   the walkers.
   *Derived classes can add their parameters.
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
    if(nw) {
      ndiff = nTargetWalkers-nw;
    } else {
      ndiff=(nTargetWalkers)? nTargetWalkers:defaultw;
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
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

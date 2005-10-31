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
    AppendRun(false),pStride(false),CurrentStep(0),
    nBlocks(100), nSteps(1000), 
    nAccept(0), nReject(0), nTargetWalkers(0), AcceptIndex(-1),
    Tau(0.001), qmcNode(NULL),
    QMCType("invalid"), h5FileRoot("invalid"),
    W(w), Psi(psi), H(h), Estimators(0),
    LogOut(0)
  { 
    m_param.add(nSteps,"steps","int");
    m_param.add(nBlocks,"blocks","int");
    m_param.add(nTargetWalkers,"walkers","int");
    m_param.add(CurrentStep,"current","int");
    m_param.add(Tau,"Tau","AU");

    //add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
    H.add2WalkerProperty(W);
  }

  QMCDriver::~QMCDriver() { 
    
    if(Estimators) {
      if(Estimators->size()) W.setLocalEnergy(Estimators->average(0));
      delete Estimators;
    }
    
    if(LogOut) delete LogOut;
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
    put(qmcNode);

    bool firstTime = (branchEngine == 0);
    if(firstTime) {
      branchEngine = new BranchEngineType(Tau,W.getActiveWalkers());
    }

    branchEngine->put(qmcNode,LogOut);

    if(firstTime && h5FileRoot.size() && h5FileRoot != "invalid") {
      LOGMSG("Initializing BranchEngine with " << h5FileRoot)
      HDFWalkerInput wReader(h5FileRoot,0,1);
      wReader.read(*branchEngine);
    }
    
    //A new run, branchEngine needs to be flushed
    if(!AppendRun) branchEngine->flush(0);

    //create estimator if not allocated
    if(Estimators == 0) Estimators =new ScalarEstimatorManager(H);

    //reset the Properties of all the walkers
    int numCopies= (H1.empty())?1:H1.size();
    W.resetWalkerProperty(numCopies);

    Estimators->put(qmcNode);
    //set the stride for the scalar estimators 
    Estimators->setStride(nSteps);

    Estimators->resetReportSettings(RootName, AppendRun);
    AcceptIndex=Estimators->addColumn("AcceptRatio");

    //initialize the walkers before moving: can be moved to run
    initialize();

    Counter++; 
  }

  void QMCDriver::setStatus(const string& aname, const string& h5name, bool append) {

    RootName = aname;
    char logfile[128];
    sprintf(logfile,"%s.%s",RootName.c_str(),QMCType.c_str());
    
    if(LogOut) delete LogOut;
    LogOut = new OhmmsInform(" ",logfile);
    
    LogOut->getStream() << "#Starting a " << QMCType << " run " << endl;

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
      int pid=OHMMS::Controller->mycontext(); 
      for(int ifile=0; ifile<nfile; ifile++) {
        string cfile("invalid"), target("e");
        int anode=-1, nwalkers=-1;
        OhmmsAttributeSet pAttrib;
        pAttrib.add(cfile,"href"); pAttrib.add(cfile,"file"); 
        pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 
        pAttrib.add(anode,"node");
        pAttrib.add(nwalkers,"walkers");
        pAttrib.put(wset[ifile]);
        int pid_target= (anode<0)?pid:anode;
        if(pid_target == pid && cfile != "invalid") {
          XMLReport("Using previous configuration of " << target << " from " << cfile)
          HDFWalkerInput WO(cfile); 
          WO.append(W,nwalkers);
          //read random state
          WO.getRandomState(true);
          h5FileRoot = cfile;
        }
      }
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

          RealType vsq = Dot(W.G,W.G);
          RealType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
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
      LOGMSG("Evaluate all the walkers before starting for walker-by-walker update")

      MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
      while(it != it_end) {
        W.R = (*it)->R;
        //DistanceTable::update(W);
        W.update();

        ValueType logpsi(Psi.evaluateLog(W));
        RealType vsq = Dot(W.G,W.G);
        RealType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
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
      Buffer_t& w_buffer((*it)->DataSet);
      w_buffer.rewind();
      W.updateBuffer(**it,w_buffer);
      ValueType logpsi=Psi.updateBuffer(W,w_buffer);
      RealType enew= H.evaluate(W);
      (*it)->resetProperty(logpsi,Psi.getSign(),enew);
      H.saveProperty((*it)->getPropertyBase());
      ValueType vsq = Dot(W.G,W.G);
      ValueType scale = ((-1.0+sqrt(1.0+2.0*Tau*vsq))/vsq);
      (*it)->Drift = scale*W.G;
      ++it;
    }
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

      LOGMSG("Adding " << nwalkers << " walkers to " << nold << " existing sets")

      W.createWalkers(nwalkers);
      LogOut->getStream() <<"Added " << nwalkers << " walkers" << endl;
      
      ParticleSet::ParticlePos_t rv(W.getTotalNum());
      MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
      while(it != it_end) {
	(*it)->R=W.R;
	++it;
      }
    } else {
      LOGMSG("Using the existing " << W.getActiveWalkers() << " walkers")
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
     
    if(cur) {
      //initialize the parameter set
      m_param.put(cur);

      xmlNodePtr tcur=cur->children;
      //determine how often to print walkers to hdf5 file
      while(tcur != NULL) {
	string cname((const char*)(tcur->name));
	if(cname == "record") {
	  int stemp(-1);
          const xmlChar* aptr=xmlGetProp(tcur,(const xmlChar*)"stride");
          if(aptr) stemp = atoi((const char*)aptr);
	  if(stemp >= 0){
	    pStride = true;
	    LOGMSG(" Append walker ensemble every block.")
	  } else {
	    pStride = false;
	    LOGMSG(" Overwrite walker ensemble every block.")
	  }
	} 
	tcur=tcur->next;
      }
    }
    
    //reset CurrentStep to zero if qmc/@continue='no'
    if(!AppendRun) CurrentStep=0;

    LogOut->getStream() << "#timestep = " << Tau << endl;
    LogOut->getStream() << "#blocks = " << nBlocks << endl;
    LogOut->getStream() << "#steps = " << nSteps << endl;
    LogOut->getStream() << "#walkers = " << W.getActiveWalkers() << endl;
    LogOut->getStream() << "#current = " << CurrentStep << endl;
    m_param.get(cout);

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

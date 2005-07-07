//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jordan Vincent and Jeongnim Kim
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
#include "QMCDrivers/VMC_OPT.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/VMC.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Utilities/Clock.h"
#include "Optimize/Minimize.h"
#include <algorithm>
namespace ohmmsqmc {

  struct VMC_OPT::WalkerData {
    ValueType LogPsi;
    RealType Vloc;
    ParticleSet::ParticleGradient_t G;
    ParticleSet::ParticleLaplacian_t L;
    explicit WalkerData(int n=1):LogPsi(0.0),Vloc(0.0),G(n), L(n){}
  };

  VMC_OPT::VMC_OPT(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), 
    UseWeight(true),
    NumCostCalls(0),
    NumSamples(0),
    cg_tolerance(1.0e-8),
    cg_stepsize(0.01), 
    cg_epsilon(1.0e-6),
    w_en(1.0),
    w_var(1.0)
  { 
    paramList.resize(10); 
    costList.resize(10,0.0); 

    //read to use vmc output (just in case)
    RootName = "vmc";
    QMCType ="opt";

    m_param.add(UseWeight,"useweight","none");
    m_param.add(cg_tolerance,"tolerance","none");
    m_param.add(cg_stepsize,"stepsize","none");
    m_param.add(cg_epsilon,"epsilon","none");

    H_KE.add(H.getHamiltonian("Kinetic"),"Kinetic");
  }

  /** Clean up the vector */
  VMC_OPT::~VMC_OPT() {
    std::vector<WalkerData*>::iterator wit(RefConf.begin());
    std::vector<WalkerData*>::iterator wit_end(RefConf.end());
    while(wit != wit_end) {delete (*wit); ++wit;}
  }

  /** Minimize the cost function with respect to a set of parameters.  
   *
   The choice of the cost function depends on the application.  
   The local energy
   \f[
   E_V() = \frac{ \int d{\bf R}\: \Psi^2({\bf R}) w({\bf R}) E_L({\bf R})}
   { \int d{\bf R}\: \Psi^2({\bf R}) w({\bf R})}
   \f]
   where \f$ E_L \f$ is the local energy and \f$ w \f$ is the weight:
   the subscript of $\Psi$ refers to the \emph{current} parameter set:
   while \f$ \Psi \f$ is a function of the
   \emph{original} parameter set.  In discreet form
   \f[
   E_V() = \frac{\sum_{j=1}^M E_L({\bf R_j}) w({\bf R_j})}
   {\sum_{j=1}^M w({\bf R_j})}
   \f]
   where the summnation is over the \f$ M \f$ configurations 
   \f$ {\bf R} \f$ distributed
   according to the probability distribution \f$ P=|\Psi({\bf R})|^2 \f$.  
   Variance of the local energy
   \f[
   \sigma^2_E() = \frac{ \int d{\bf R}\: \Psi^2({\bf R}) w({\bf R}) 
   \left[ E_L({\bf R})-E_V()\right]^2}{ \int d{\bf R}\:
   \Psi^2({\bf R}) w({\bf R})}
   \f]
   In discreet form
   \f[
   \sigma^2_E() = \frac{\sum_{j=1}^M \left[ E_L({\bf R_j})-E_V()\right]^2
   w({\bf R_j})}{\sum_{j=1}^M w({\bf R_j})}
   \f]
  */

  VMC_OPT::RealType VMC_OPT::evalCost() {

    Estimators->flushreport(1);

    //Estimators::accumulate has been called by correlatedSampling
    RealType eav = Estimators->average(0);
    RealType evar = Estimators->variance(0);
    RealType cost= w_en*eav+w_var*0.5*Tau*evar;


    if(NumCostCalls == 0) {
      log_buffer << NumCostCalls << " ";
      std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<RealType>(log_buffer," "));
      log_buffer << " " << NumSamples;
    }
    
    log_buffer << " " << eav << " " << evar << " " << cost;
    LogOut->getStream() << log_buffer.rdbuf() << endl;
    log_buffer.clear();

    return cost;
  }

  bool
  VMC_OPT::run() {

    Estimators->reportHeader();

    //set the data members to start a new run
    checkConfigurations();

    //estimator has to collect the data over mpi nodes
    Estimators->setCollectionMode(OHMMS::Controller->ncontexts()>1);

    correlatedSampling();

    //   log_buffer << "0 ";
    //   std::copy(OptParams.begin(),OptParams.end(),
    //             ostream_iterator<RealType>(log_buffer," "));
    CostValue = evalCost();

    ConjugateGradient CG;
    CG.Tolerance = cg_tolerance;
    CG.StepSize = cg_stepsize;
    CG.epsilon = cg_epsilon;
    CG.Minimize(*this);

    return true;
  }

  void VMC_OPT::addConfiguration(const string& a) {

    ConfigFile.push_back(a);

  }

  /** Apply constraints on the optimizables. */

  bool VMC_OPT::checkParameters() {

    bool samesign = true;
    //Checking sign change and reject if the sign of the 
    //parameter changes but not seem to do so well.
    //   int i = 0;
    //   while(i<OptParams.size() && samesign) {
    //     if(OptParams[i]*paramList[0][i]<0) samesign = false;
    //     i++;
    //   }
    if(samesign) {  
      paramList.pop_back();
      paramList.push_front(OptParams);
    }

    return samesign;
  }

  VMC_OPT::RealType
  VMC_OPT::Cost() {

    NumCostCalls++;

    if(checkParameters()) {

      resetWaveFunctions();
      ValueType nw_effect = correlatedSampling();

      log_buffer << NumCostCalls << " ";
      std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<RealType>(log_buffer," "));
      log_buffer << " " << nw_effect;

      RealType cost = evalCost();

      costList.pop_back();
      costList.push_front(cost);
      /* stop the optimization of the number of effective walkers
	 is less than 20% the original population */
      if(nw_effect < NumSamples/5) {
      
	ERRORMSG("Number of Effective Walkers is too small " << nw_effect)
	ERRORMSG("Going to stop now.")
        OHMMS::Controller->abort();
	//int index = 0;
	//scalar costmin = costList[index];
	//int i=1;
	//int imin = std::min(NumCostCalls,int(costList.size()));
	//while(i<imin) {
	//  if (costmin > costList[i]){
	//    costmin = costList[i];
	//    index = i;
	//  }
	//  i++;
	//}
	//OptParams = paramList[index];
	//resetWaveFunctions();
	//LogOut->getStream() << "#Run VMC" << endl;
	//run_vmc();
        cost = 1e6;
      } // if nw_effect < NumSamples/5
      CostValue = cost;
    }

    return CostValue;

  }

  /** 
   *@brief Perform the correlated sampling algorthim.
   *
   The correlated sampling algorithim begins by reading in a set of
   configurations (walkers) \f$ {{\bf R}} \f$ that are distributed 
   according to \f$ \Psi(\{{\bf \alpha_0}\})^2 \f$ where 
   \f$ {\{\bf \alpha_0}\} \f$ is the initial
   parameter set. The configurations are created by previous 
   VMC runs and are stored in ".config" 
   files.  For each configuration calculate the local
   energy and weight: 
   \f[ E_L({\bf R}) = \frac{\hat{H}\Psi({\bf R},\{{\bf \alpha_i}\})}
   {\Psi({\bf R},\{{\bf \alpha_0}\})} \f]   
   \f[ w({\bf R}) = \frac{\Psi({\bf R},\{{\bf \alpha_i}\})^2}
   {\Psi({\bf R},\{{\bf \alpha_0}\})^2}, \f]
   where \f$ \{{\bf \alpha_i}\} \f$ is the current parameter set. 
   * \return the number of effective walkers \f$ N_{Eff} = 
   \sum_{i=1}^N w({\bf R_i}) / \sum_{i=1}^{N} w^2({\bf R_i}) \f$
  */

  VMC_OPT::ValueType 
  VMC_OPT::correlatedSampling() {

    // numerator and denominator for number of effective walkers
    //ValueType nw_effect_t = 0.0;
    //ValueType nw_effect_b = 0.0;

    // flush estimators
    Estimators->reset();
    NumSamples = 0;
    //TinyVector<RealType,3> nw_effect(0.0);
    std::vector<RealType> nw_effect(3,0.0);
    std::vector<WalkerData*>::const_iterator wit(RefConf.begin());

    for(int i=0; i<ConfigFile.size(); i++) {

      HDFWalkerInput wReader(ConfigFile[i]);

      while(wReader.put(W)) {

	//accumulate the number of samples      
	NumSamples += W.getActiveWalkers();
	MCWalkerConfiguration::iterator it = W.begin(); 
	MCWalkerConfiguration::iterator it_end = W.end(); 

        while(it != it_end) {

          Walker_t& thisWalker(**it);

          ValueType logpsi0((*wit)->LogPsi);
          ValueType vold((*wit)->Vloc);

	  W.R = thisWalker.R;
	  DistanceTable::update(W);

          ValueType logpsi(Psi.evaluateLog(W,false));
          W.G += (*wit)->G;
          W.L += (*wit)->L;
	  MCWalkerConfiguration::PropertyContainer_t& Properties((*it)->Properties);
          ValueType weight(1.0);
          if(UseWeight) {
            weight = exp(2.0*(logpsi-logpsi0));
	    Properties(LOGPSI) = logpsi;
          } 
	  thisWalker.Properties(LOCALENERGY) = H_KE.evaluate(W)+vold;
	  thisWalker.Weight = weight;
	  nw_effect[0] += weight;
	  nw_effect[1] += weight*weight;
          //move on to next walker
          ++it;++wit;
	}
	Estimators->accumulate(W);
      }
    }

    nw_effect[2] = static_cast<RealType>(NumSamples);
    gsum(nw_effect,0);
    NumSamples = static_cast<IndexType>(nw_effect[2]);
    return (nw_effect[0]*nw_effect[0]/nw_effect[1]);
  }

  void 
  VMC_OPT::checkConfigurations() {
    //getReady();
    RefConf.reserve(ConfigFile.size()*500);
    NumSamples=0;
    for(int i=0; i<ConfigFile.size(); i++) {
      HDFWalkerInput wReader(ConfigFile[i]);
      while(wReader.put(W)) {
	//accumulate the number of samples      
	NumSamples += W.getActiveWalkers();
	MCWalkerConfiguration::iterator it(W.begin()); 
	MCWalkerConfiguration::iterator it_end(W.end()); 
        int nat = W.getTotalNum();
        while(it != it_end) {
	  W.R = (*it)->R;
          WalkerData* awalker=new WalkerData(nat);
	  //update the distance table associated with W
	  DistanceTable::update(W);
	  //evaluate wave function
          //ValueType psi(Psi.evaluate(W));
          ValueType logpsi(Psi.evaluateLog(W,awalker->G, awalker->L));
          RealType e= H.evaluate(W);
          awalker->LogPsi=logpsi;
          awalker->Vloc=H.getLocalPotential();
          RefConf.push_back(awalker);
          ++it;
        }
      }
    }
  }

  /** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}). \f$
   *
   * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
   * set \f$ {{\bf \alpha_{i+1}}}. \f$  
   */
  bool
  VMC_OPT::resetWaveFunctions() {

    int offset = 0;
    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++){
      //locate the id in the variable list for Psi
      int id = Psi.VarList.find(IDtag[i]);
      if(id>= 0) {
	//find the number of variables represented by the id
	int size = Psi.VarList.Sizes[id];
	//loop over the number of variables for the id
	//assign the variable its new value
	RealType* temp = Psi.VarList.Pointers[id];
	for(int j=0; j<size; j++,temp++){
	  *(temp) = OptParams[j+offset];
	}
	//increment offset to account for the fact that id
	//may represent more than one variable
	offset += size;
      } 
    }
    //reset the wavefunction for the new variables
    Psi.reset();

    return true;
  }


  void VMC_OPT::run_vmc() { 

    string vtitle=RootName;
    vtitle.append(".vmc");
  
    VMC vmc(W,Psi,H);
    vmc.setFileRoot(vtitle);
    vmc.process(qmcNode);
    vmc.run();

    ConfigFile.clear();
    ConfigFile.push_back(vtitle);
  }

  /** Parses the xml input file for parameter definitions for the wavefunction optimization.
   * @param q the current xmlNode 
   * @return true if successful
   *
   * Must provide a list of the id's of the variables to be optimized. Other available elements
   <ul>
   <li> mcwalkerset->file: name of configuration file
   <li> optimize->method: method of optimization, right now only 
   congjugate gradient 
   <li> cost->weight: either weight for the energy term, default 1.0, 
   or weight for the variance, default 1.0
   <li> tolerance: for conjugate gradient, default 1e-8
   <li> stepsize: for conjugate gradient, default 0.01
   <li> epsilon: for conjugate gradient, default 1e-6
   </ul>
   */
  bool
  VMC_OPT::put(xmlNodePtr q) {

    xmlNodePtr qsave=q;
    //Estimators.put(q);
    vector<xmlNodePtr> wset,oset,cset;
    xmlNodePtr cur=qsave->children;
    int pid=OHMMS::Controller->mycontext();
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "mcwalkerset") {
        int pid_target=pid;
        xmlChar* anode=xmlGetProp(cur,(const xmlChar*)"node");
        if(anode) {
          pid_target = atoi((const char*)anode);
        }
        if(pid_target == pid) wset.push_back(cur);
      } else if(cname == "optimize") {
	oset.push_back(cur);
      } else if(cname == "cost") {
	cset.push_back(cur);
      }
      cur=cur->next;
    }  

    //erase the existing config files and use new sets
    if(wset.size()) {
      ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
      for(int i=0; i<wset.size(); i++) {
	xmlChar* att= xmlGetProp(wset[i],(const xmlChar*)"file");
	if(att) ConfigFile.push_back((const char*)att);
      }
    }

    if(oset.empty()) {
      ERRORMSG("No Optimization Variables Set")
	return false;
    } else {
      for(int i=0; i<oset.size(); i++) {
	xmlChar* att= xmlGetProp(oset[i],(const xmlChar*)"method");
	if(att) {
	  optmethod = (const char*)att;
	  LogOut->getStream() << "#Optimization: using " << optmethod << " method." << endl;
	} else {
	  optmethod = "cg";
	  LogOut->getStream() << "#Optimization: using " << optmethod << " method." << endl;
	}
	vector<string> idtag;
	putContent(idtag,oset[i]);
	IDtag.reserve(idtag.size());
	for(int id=0; id<idtag.size(); id++) {
	  vector<string>::iterator it 
	    = std::find(IDtag.begin(), IDtag.end(),idtag[id]);
	  if(it == IDtag.end()) IDtag.push_back(idtag[id]);
	}
	std::copy(IDtag.begin(),IDtag.end(), 
		  ostream_iterator<string>(log_buffer," "));
	LogOut->getStream() << "#Optimiziable variables " << log_buffer.rdbuf() << endl;
	log_buffer.clear();
      }
    }
  
    if(cset.empty()){
      LogOut->getStream() 
	<< "#Using Default Cost Function: Cost = <E> + 0.5 Tau <Var>" << endl;
    } else {
      for(int i=0; i<cset.size(); i++){
	string pname;
	RealType wgt=1.0;
	xmlAttrPtr att = cset[i]->properties;
	while(att != NULL) {
	  string aname((const char*)(att->name));
	  const char* vname=(const char*)(att->children->content);
	  if(aname == "name") 
	    pname=vname;
	  else if(aname == "weight") 
	    wgt = atof(vname);
	  att=att->next;
	}
	if(pname == "energy") 
	  w_en = wgt;
	else if(pname == "variance") 
	  w_var = wgt;
	LogOut->getStream() << "#Cost Function: Cost = " << w_en << "*<E> + " 
			    << w_var << "*0.5 Tau <Var>" << endl;
      }
    }  

    //m_param.put(qsave);
    LogOut->getStream() << "#" << optmethod << " values: tolerance = "
			<< cg_tolerance << " stepsize = " << cg_stepsize 
			<< " epsilon = " << cg_epsilon << " Tau = " 
			<< Tau << endl;

    if(!UseWeight){
      LogOut->getStream() << "#All weights set to 1.0" << endl;
      w_en = 0.0;
      w_var = 2.0/Tau;
      LogOut->getStream() << "#Cost Function: Cost = " << w_en << "*<E> + " 
			  << w_var << "*0.5 Tau <Var>" << endl;
    }

    putOptParams();

    LogOut->getStream() 
      << "# list of the configuration files used for optimization" << endl;
    for(int i=0; i<ConfigFile.size(); i++) 
      LogOut->getStream() << "# " << ConfigFile[i] << endl;

    return true;
  }

  bool
  VMC_OPT::putOptParams(){

    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++){
      //locate the id in the variable list for Psi
      int id = Psi.VarList.find(IDtag[i]);
      if(id>= 0) {
	//find the number of variables represented by the id
	int size = Psi.VarList.Sizes[id];
	RealType* temp = Psi.VarList.Pointers[id];
	//loop over the number of variables for the id
	//assign the value to the optimization parameter set
	for(int j=0; j<size; j++){
	  OptParams.push_back(temp[j]);
	}
      } else {
	ERRORMSG("Could not find parameter " << IDtag[i])
	  return false;
	OptParams.push_back(0.0);
      }
    }

    for(int i=0; i<paramList.size(); i++) {
      paramList[i] = OptParams;
    }

    log_buffer << "#Inital Variables ";
    std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<RealType>(log_buffer," "));
    LogOut->getStream() << log_buffer.rdbuf() << endl << endl;
    log_buffer.clear();

    log_buffer << "#Results: index ";
    std::copy(IDtag.begin(),IDtag.end(), ostream_iterator<string>(log_buffer," "));
    log_buffer << " effective-walkers  cost-values ";
    LogOut->getStream() << log_buffer.rdbuf() << endl;
    log_buffer.clear();
    Psi.reset();
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

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
#include "QMCDrivers/QMCOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include <algorithm>
#include <limits>
namespace ohmmsqmc {

  QMCOptimize::QMCOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), 
    UseWeight(false),
    PartID(0),
    NumParts(1),
    NumCostCalls(0),
    NumSamples(0),
    cg_tolerance(1.0e-4),
    cg_stepsize(0.01), 
    cg_epsilon(1.0e-3),
    w_en(0.0),
    w_var(1.0),
    w_abs(1.0)
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

    //H_KE.add(H.getHamiltonian("Kinetic"),"Kinetic");
    for(int i=0; i<H.size(); i++) {
      QMCHamiltonianBase* qH=H.getHamiltonian(i);
      if(qH->UpdateMode[QMCHamiltonianBase::OPTIMIZABLE]) {
        H_KE.add(qH,H.getName(i)); 
      }
    }

    //create an etimator with H_KE
    if(Estimators == 0) Estimators =new ScalarEstimatorManager(H_KE);
    H_KE.add2WalkerProperty(W);

    //check, if H needs to evaluate ratio
    hamiltonianNeedRatio=H.needRatio();
    if(hamiltonianNeedRatio) {
      LOGMSG("QMCHamiltonian requires the ratio evaluations. Use buffer.")
    } else {
      LOGMSG("QMCHamiltonian does not require the ratio evaluations.")
    }
  }

  /** Clean up the vector */
  QMCOptimize::~QMCOptimize() {
  }

  /** Main functio to run a optimization
   */
  bool
  QMCOptimize::run() {


    Estimators->reportHeader();

    dL.resize(W.getTotalNum());

    //set the data members to start a new run
    checkConfigurations();

    //over-write the Etarget by E_T
    if(abs(branchEngine->E_T)>numeric_limits<RealType>::epsilon()) {
      Etarget=branchEngine->E_T;
      w_abs=0.0; w_var=1.0;
    }


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



  void QMCOptimize::addConfiguration(const string& a) {

    ConfigFile.push_back(a);

  }

  /** Apply constraints on the optimizables. */

  bool QMCOptimize::checkParameters() {

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

  QMCOptimize::RealType
  QMCOptimize::Cost() {

    NumCostCalls++;

    if(checkParameters()) {

      //reset the wave function
      resetWaveFunctions();

      //evaluate new local energies
      RealType nw_effect=correlatedSampling();

      log_buffer << setw(5) << NumCostCalls << " ";
      std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<RealType>(log_buffer," "));
      log_buffer << " " << nw_effect;

      //evaluate the cost value
      CostValue = evalCost();
      //if(nw_effect < NumSamples/10) {

      //  ERRORMSG("Number of Effective Walkers is too small " << nw_effect)
      //    ERRORMSG("Going to stop now.")
      //    OHMMS::Controller->abort();
      //  cost = 1e6;
      //} // if nw_effect < NumSamples/5
      //CostValue = cost;
    }

    return CostValue;

  }

  /** Minimize the cost function with respect to a set of parameters.  
  */

  QMCOptimize::RealType QMCOptimize::evalCost() {

    //Estimators::accumulate has been called by correlatedSampling
    RealType eavg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
    RealType evar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-eavg_w*eavg_w;

    RealType wgtinv=1.0/static_cast<RealType>(W.getActiveWalkers());
    RealType eavg = SumValue[SUM_E_BARE]*wgtinv;
    RealType evar = SumValue[SUM_ESQ_BARE]*wgtinv-eavg*eavg; 

    //DIFFERENT COST FUNCTIOn
    //RealType cost= w_en*eavg+w_var*0.5*Tau*evar;
    RealType cost = w_abs*SumValue[SUM_ABSE_BARE]*wgtinv+w_var*evar;

    if(NumCostCalls == 0) {
      log_buffer << setw(5) << NumCostCalls << " ";
      std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<RealType>(log_buffer," "));
      log_buffer <<setw(17) << NumSamples;
    }
    
    log_buffer << setw(20) << eavg_w << setw(20) << eavg << setw(20) 
      << evar_w << setw(20) << evar << setw(20) << cost;
    LogOut->getStream() << log_buffer.rdbuf() << endl;
    log_buffer.clear();

    return cost;
  }

  /**  Perform the correlated sampling algorthim.
   *
   */
  QMCOptimize::RealType QMCOptimize::correlatedSampling() {

    SumValue=0.0;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    RealType eloc_new;
    int iw=0;
    while(it != it_end) {
      Walker_t& thisWalker(**it);
      ValueType* saved = Records[iw];

      //rewind the buffer to get the data from buffer
      thisWalker.DataSet.rewind();
      //W is updated by thisWalker
      W.copyFromBuffer(thisWalker.DataSet);

      ValueType logpsi=0.0;
      if(hamiltonianNeedRatio) {
        logpsi=Psi.evaluateLog(W);
      } else {
        //copy dL from Buffer
        thisWalker.DataSet.get(dL.begin(),dL.end());
        logpsi=Psi.evaluateDeltaLog(W);
        W.G += thisWalker.Drift;
        W.L += dL;
      }

      eloc_new=H_KE.evaluate(W)+saved[ENERGY_FIXED];
      RealType weight = exp(2.0*(logpsi-saved[LOGPSI_FREE]));

      //////////////////////////////////////////
      //THIS WAS TO TEST 
      //////////////////////////////////////////
      //ValueType logpsi2(Psi.evaluateLog(W));
      //RealType et= H.evaluate(W);
      //if(abs(logpsi+saved[LOGPSI_FIXED]-logpsi2)>1e-10 || abs(et-eloc_new)>1e-3)
      //  cout << "Check wfs and energy " << logpsi+saved[LOGPSI_FIXED]-logpsi2 << " " << et-eloc_new << endl;

      saved[ENERGY_NEW]=eloc_new;
      saved[REWEIGHT]=weight;

      RealType delE=abs(eloc_new-Etarget);

      SumValue[SUM_E_BARE] += eloc_new;
      SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
      SumValue[SUM_ABSE_BARE] += delE;
      SumValue[SUM_E_WGT] += eloc_new*weight;
      SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*weight;
      SumValue[SUM_ABSE_WGT] += delE*weight;
      SumValue[SUM_WGT] += weight;
      SumValue[SUM_WGTSQ] += weight*weight;
      ++it;
      ++iw;
    }

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

  /** evaluate everything before optimization */
  void 
  QMCOptimize::checkConfigurations() {

    if(ConfigFile.size()) {
      W.destroyWalkers(W.begin(),W.end());
      for(int i=0; i<ConfigFile.size(); i++) {
        HDFWalkerInput wReader(ConfigFile[i],PartID,NumParts);
        wReader.append(W);
      }

      NumSamples=W.getActiveWalkers();
      Records.resize(NumSamples,6);
      MCWalkerConfiguration::iterator it(W.begin()); 
      MCWalkerConfiguration::iterator it_end(W.end()); 
      int nat = W.getTotalNum();
      int iw=0;
      Etarget=0.0;
      while(it != it_end) {

        Walker_t& thisWalker(**it);

        //clean-up DataSet to save re-used values
        thisWalker.DataSet.clear();
        //rewind the counter
        thisWalker.DataSet.rewind();
        //MCWalkerConfiguraton::registerData add distance-table data
        W.registerData(thisWalker,thisWalker.DataSet);

        ValueType*  saved=Records[iw];

        if(hamiltonianNeedRatio) {
          saved[LOGPSI_FREE]=Psi.evaluateLog(W);
        } else {
          Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], thisWalker.Drift, dL);
          thisWalker.DataSet.add(dL.begin(),dL.end());
        }

        Etarget += saved[ENERGY_TOT] = H.evaluate(W);
        saved[ENERGY_FIXED] = H.getInvariantEnergy();

        ++it;
        ++iw;
      }
      //remove input files
      ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
    }

    Etarget /= static_cast<RealType>(NumSamples);

    LOGMSG("Total number of walkers " << NumSamples << " Eavg = " << Etarget)
  }

  /** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
   *@return true always
   *
   * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
   * set \f$ {{\bf \alpha_{i+1}}}\f$  
   */
  bool
  QMCOptimize::resetWaveFunctions() {

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
  
  void QMCOptimize::WriteStuff() {
    //WRITE A XML FILE
  }

  /** Parses the xml input file for parameter definitions for the wavefunction optimization.
   * @param q current xmlNode 
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
  QMCOptimize::put(xmlNodePtr q) {

    xmlNodePtr qsave=q;
    //Estimators.put(q);
    vector<xmlNodePtr> wset,oset,cset;
    xmlNodePtr cur=qsave->children;
    int pid=OHMMS::Controller->mycontext();
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "mcwalkerset") {
        wset.push_back(cur);
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
      int nfile=wset.size();
      int ng=OHMMS::Controller->ncontexts()/nfile;
      if(ng>=1) {
        NumParts=ng;
        PartID=pid%ng;
        int mygroup=pid/ng;
        string fname("invalid");
        OhmmsAttributeSet pAttrib;
        pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
        pAttrib.put(wset[mygroup]);
        if(fname != "invalid") ConfigFile.push_back(fname);
      } else {//more files than the number of processor
        for(int ifile=0; ifile<nfile; ifile++) {
          int pid_target=pid;
          string fname("invalid");
          OhmmsAttributeSet pAttrib;
          pAttrib.add(pid_target,"node"); 
          pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
          pAttrib.put(wset[ifile]);
          if(pid_target == pid && fname != "invalid") {
            ConfigFile.push_back(fname);
          }
        }
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
      LogOut->getStream() << "#Using Default Cost Function: Cost = <E> + 0.5 Tau <Var>" << endl;
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
        else if(pname == "abs")
          w_abs = wgt;
      }
    }  


    //m_param.put(qsave);
    LogOut->getStream() << "#" << optmethod << " values: tolerance = "
			<< cg_tolerance << " stepsize = " << cg_stepsize 
			<< " epsilon = " << cg_epsilon << " Tau = " 
			<< Tau << endl;

    if(!UseWeight){
      LogOut->getStream() << "#All weights set to 1.0" << endl;
    }

    LogOut->getStream() << "#Cost Function: <E> <Var> |E-E_T| " << w_en << " " << w_var << " " << w_abs << endl;

    putOptParams();

    LogOut->getStream() 
      << "# list of the configuration files used for optimization" << endl;
    for(int i=0; i<ConfigFile.size(); i++) 
      LogOut->getStream() << "# " << ConfigFile[i] << " part " << PartID << "/" << NumParts << endl;

    return true;
  }

  bool
  QMCOptimize::putOptParams(){

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

    log_buffer.setf(ios::scientific, ios::floatfield);
    log_buffer.setf(ios::right,ios::adjustfield);
    log_buffer.precision(10);

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

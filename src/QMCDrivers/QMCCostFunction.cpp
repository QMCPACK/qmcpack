//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#include "QMCDrivers/QMCCostFunction.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/CommOperators.h"
namespace qmcplusplus {

  QMCCostFunction::QMCCostFunction(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    W(w),Psi(psi),H(h),
    UseWeight(false), PowerE(2), NumCostCalls(0), NumSamples(0),
    w_en(0.0), w_var(0.0), w_abs(1.0),
    CorrelationFactor(0.0), m_wfPtr(NULL), m_doc_out(NULL), msg_stream(0)
  { 

    paramList.resize(10); 
    costList.resize(10,0.0); 

    //default: when walkers below 90% stop
    MinNumWalkers = 0.9;

    H_KE.add(H.getHamiltonian("Kinetic"),"Kinetic");
    //move to Optimize class
    ////create an etimator with H_KE
    //if(Estimators == 0) Estimators =new ScalarEstimatorManager(H_KE);
    H_KE.add2WalkerProperty(W);

    IsValid=true;
  }

  /** Clean up the vector */
  QMCCostFunction::~QMCCostFunction() {
    if(m_doc_out != NULL) xmlFreeDoc(m_doc_out);

  }

  void QMCCostFunction::setTargetEnergy(Return_t et) {

    Etarget = et;
    LOGMSG("Etarget (set from previous runs) = " << Etarget)

    //evaluate effective target energy
    EtargetEff=(1.0+CorrelationFactor)*Etarget;

    LOGMSG("Effective Target Energy = " << EtargetEff)
    LOGMSG("Cost Function = " << w_en << "*<E> + " << w_var << "*<Var> + " << w_abs << "*|E-E_T|^" << PowerE)

    if(msg_stream) {
      *msg_stream << "  Total number of walkers          = " << NumSamples << endl;
      *msg_stream << "  Effective Target Energy = " << EtargetEff << endl;
      *msg_stream << "  Cost Function = " << w_en << "*<E> + " << w_var << "*<Var> + " << w_abs << "*|E-E_T|^" << PowerE << endl;
      *msg_stream << "  Optimization report = ";
      *msg_stream << "cost, walkers, eavg/wgt, eavg/walkers, evar/wgt, evar/walkers, evar_abs\n";
      *msg_stream << "  Optimized variables = ";
      std::copy(IDtag.begin(),IDtag.end(),ostream_iterator<string>(*msg_stream,", "));
      *msg_stream << endl;
    }
  }

  QMCCostFunction::Return_t QMCCostFunction::Cost() {

    NumCostCalls++;

    if(checkParameters()) {

      //reset the wave function
      resetWaveFunctions();

      //evaluate new local energies
      NumWalkersEff=correlatedSampling();

      //Estimators::accumulate has been called by correlatedSampling
      //Return_t eavg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
      //Return_t evar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-eavg_w*eavg_w;
      curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
      curVar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-curAvg_w*curAvg_w;

      Return_t wgtinv=1.0/static_cast<Return_t>(W.getActiveWalkers());
      //Return_t eavg = SumValue[SUM_E_BARE]*wgtinv;
      //Return_t evar = SumValue[SUM_ESQ_BARE]*wgtinv-eavg*eavg; 
      curAvg = SumValue[SUM_E_BARE]*wgtinv;
      curVar = SumValue[SUM_ESQ_BARE]*wgtinv-curAvg*curAvg; 

      //DIFFERENT COST FUNCTIOn
      //CostValue = w_en*eavg+w_var*0.5*Tau*evar;
      //CostValue = w_abs*SumValue[SUM_ABSE_BARE]*wgtinv+w_var*evar;
      
      //Return_t evar_abs = SumValue[SUM_ABSE_WGT]/SumValue[SUM_WGT];
      //CostValue = w_abs*evar_abs+w_var*evar;
      curVar_abs = SumValue[SUM_ABSE_WGT]/SumValue[SUM_WGT];
      CostValue = w_abs*curVar_abs+w_var*curVar;

      ////dump to a file
      //log_buffer << setw(5) << NumCostCalls << " ";
      //for(int i=0; i<OptParams.size(); i++) log_buffer << setw(16) << OptParams[i];
      //log_buffer << setw(15) << nw_effect;
      //log_buffer << setw(16) << eavg_w << setw(16) << eavg << setw(16) << evar_w << setw(16) << evar << setw(16) << evar_abs << setw(16) << CostValue;
      //LogOut->getStream() << log_buffer.rdbuf() << endl;
      //log_buffer.clear();

      if(NumWalkersEff < NumSamples*MinNumWalkers) {
        ERRORMSG("Number of Effective Walkers is too small " << NumWalkersEff)
        ERRORMSG("Going to stop now.")
        IsValid=false;
        //OHMMS::Controller->abort();
      } 
    }
    return CostValue;
  }

  /**  Perform the correlated sampling algorthim.
   */
  QMCCostFunction::Return_t QMCCostFunction::correlatedSampling() {

    SumValue=0.0;
    typedef MCWalkerConfiguration::Walker_t Walker_t;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    Return_t eloc_new;
    int iw=0;
    while(it != it_end) {
      Walker_t& thisWalker(**it);
      Return_t* saved = Records[iw];

      //rewind the buffer to get the data from buffer
      thisWalker.DataSet.rewind();
      //W is updated by thisWalker
      W.copyFromBuffer(thisWalker.DataSet);

      Return_t logpsi=0.0;
      //copy dL from Buffer
      thisWalker.DataSet.get(dL.begin(),dL.end());
      logpsi=Psi.evaluateDeltaLog(W);
      W.G += thisWalker.Drift;
      W.L += dL;

      eloc_new=H_KE.evaluate(W)+saved[ENERGY_FIXED];
      Return_t weight = exp(2.0*(logpsi-saved[LOGPSI_FREE]));

      //////////////////////////////////////////
      //THIS WAS TO TEST 
      //////////////////////////////////////////
      //Return_t logpsi2(Psi.evaluateLog(W));
      //Return_t et= H.evaluate(W);
      //if(abs(logpsi+saved[LOGPSI_FIXED]-logpsi2)>1e-10 || abs(et-eloc_new)>1e-3)
      //  cout << "Check wfs and energy " << logpsi+saved[LOGPSI_FIXED]-logpsi2 << " " << et-eloc_new << endl;
      saved[ENERGY_NEW]=eloc_new;
      saved[REWEIGHT]=weight;

      Return_t delE=pow(abs(eloc_new-EtargetEff),PowerE);
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

    //collect everything
    gsum(SumValue,0);

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

  /** evaluate everything before optimization */
  void QMCCostFunction::checkConfigurations(vector<string>& ConfigFile, int partid, int nparts) {

    dL.resize(W.getTotalNum());

    if(ConfigFile.size()) {
      W.destroyWalkers(W.begin(),W.end());
      for(int i=0; i<ConfigFile.size(); i++) {
        HDFWalkerInput wReader(ConfigFile[i],partid,nparts);
        wReader.append(W);
      }

      NumSamples=W.getActiveWalkers();
      Records.resize(NumSamples,6);
      typedef MCWalkerConfiguration::Walker_t Walker_t;
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

        Return_t*  saved=Records[iw];
        Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], thisWalker.Drift, dL);
        thisWalker.DataSet.add(dL.begin(),dL.end());
        Etarget += saved[ENERGY_TOT] = H.evaluate(W);
        saved[ENERGY_FIXED] = H.getInvariantEnergy();

        ++it;
        ++iw;
      }
      //remove input files
      ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
    }

    //Need to sum over the processors
    vector<Return_t> etemp(2);
    etemp[0]=Etarget;
    etemp[1]=static_cast<Return_t>(NumSamples);
    gsum(etemp,0);
    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
    NumSamples = static_cast<int>(etemp[1]);

    setTargetEnergy(Etarget);
    //if(msg_stream) {
    //  *msg_stream << "Etarget (guess from the average) = " << Etarget << endl;
    //}
  }

  /** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
   *@return true always
   *
   * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
   * set \f$ {{\bf \alpha_{i+1}}}\f$  
   */
  bool
  QMCCostFunction::resetWaveFunctions() {

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
	Return_t* temp = Psi.VarList.Pointers[id];
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
  
  void QMCCostFunction::Report() {

    static int writeCounter=0;

    if(OHMMS::Controller->master()) {

      if(m_doc_out == NULL) {
        m_doc_out = xmlNewDoc((const xmlChar*)"1.0");
        xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem"); 
        xmlAddChild(qm_root, m_wfPtr);
        xmlDocSetRootElement(m_doc_out, qm_root);
        m_param_out.resize(IDtag.size(),NULL);
        xmlXPathContextPtr acontext = xmlXPathNewContext(m_doc_out);
        xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//parameter",acontext);
        for(int iparam=0; iparam<result->nodesetval->nodeNr; iparam++) {
          xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
          string aname((const char*)xmlGetProp(cur,(const xmlChar*)"id"));
          vector<string>::iterator it= find(IDtag.begin(), IDtag.end(), aname);
          if(it != IDtag.end()) {
            int item=it-IDtag.begin();
            m_param_out[item]=cur;
          }
        }
        xmlXPathFreeObject(result);
        xmlXPathFreeContext(acontext);
      }

      //update parameters
      for(int item=0; item<m_param_out.size(); item++) {
        if(m_param_out[item] != NULL) {
          getContent(OptParams[item],m_param_out[item]);
        }
      }

      char newxml[128];
      sprintf(newxml,"%s.opt.%d.xml", RootName.c_str(),writeCounter);
      xmlSaveFormatFile(newxml,m_doc_out,1);

      if(msg_stream) {
        *msg_stream << " curCost " 
          << setw(5) << NumCostCalls << setw(16) << CostValue << setw(15) << NumWalkersEff 
          << setw(16) << curAvg_w << setw(16) << curAvg 
          << setw(16) << curVar_w << setw(16) << curVar 
          << setw(16) << curVar_abs << endl;
        *msg_stream << " curVars ";
        for(int i=0; i<OptParams.size(); i++) *msg_stream << setw(16) << OptParams[i];
        *msg_stream << endl;
      }
    }

    writeCounter++;

    OHMMS::Controller->barrier();
  }

  /** Apply constraints on the optimizables. 
   * 
   * Here is where constraints should go
   */ 
  bool QMCCostFunction::checkParameters() {

    bool samesign = true;
    if(samesign) {  
      paramList.pop_back();
      paramList.push_front(OptParams);
    }

    return samesign;
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
  QMCCostFunction::put(xmlNodePtr q) {

    ParameterSet m_param;
    m_param.add(UseWeight,"useweight","none");
    m_param.add(PowerE,"power","int");
    m_param.add(CorrelationFactor,"correlation","scalar");
    m_param.add(MinNumWalkers,"min_walkers","scalar");
    m_param.put(q);

    xmlNodePtr qsave=q;
    //Estimators.put(q);
    vector<xmlNodePtr> oset,cset;
    xmlNodePtr cur=qsave->children;
    int pid=OHMMS::Controller->mycontext();
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "optimize") {
	oset.push_back(cur);
      } else if(cname == "cost") {
	cset.push_back(cur);
      }
      cur=cur->next;
    }  

    //int nfile=mcwalkerNodePtr.size();
    //if(nfile) {
    //  int ng=OHMMS::Controller->ncontexts()/nfile;
    //  if(ng>=1) {
    //    NumParts=ng;
    //    PartID=pid%ng;
    //    int mygroup=pid/ng;
    //    string fname("invalid");
    //    OhmmsAttributeSet pAttrib;
    //    pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
    //    pAttrib.put(mcwalkerNodePtr[mygroup]);
    //    if(fname != "invalid") ConfigFile.push_back(fname);
    //  } else {//more files than the number of processor
    //    for(int ifile=0; ifile<nfile; ifile++) {
    //      int pid_target=pid;
    //      string fname("invalid");
    //      OhmmsAttributeSet pAttrib;
    //      pAttrib.add(pid_target,"node"); 
    //      pAttrib.add(fname,"href"); pAttrib.add(fname,"src"); pAttrib.add(fname,"file");
    //      pAttrib.put(mcwalkerNodePtr[ifile]);
    //      if(pid_target == pid && fname != "invalid") {
    //        ConfigFile.push_back(fname);
    //      }
    //    }
    //  }
    //}

    if(oset.empty()) {
      ERRORMSG("No Optimization Variables Set")
      return false;
    } else {
      for(int i=0; i<oset.size(); i++) {
	vector<string> idtag;
	putContent(idtag,oset[i]);
	IDtag.reserve(idtag.size());
	for(int id=0; id<idtag.size(); id++) {
	  vector<string>::iterator it 
	    = std::find(IDtag.begin(), IDtag.end(),idtag[id]);
	  if(it == IDtag.end()) IDtag.push_back(idtag[id]);
	}
        if(msg_stream) {
          *msg_stream << " Variables to be optimized: ";
          std::copy(IDtag.begin(),IDtag.end(), 
              ostream_iterator<string>(*msg_stream," "));
          *msg_stream << endl;
        }
      }
    }
  
    if(cset.empty()){
      if(msg_stream) *msg_stream << " Using Default Cost Function: Cost = <|E-E_ff|^2>" << endl;
    } else {
      for(int i=0; i<cset.size(); i++){
	string pname;
	Return_t wgt=1.0;
        OhmmsAttributeSet pAttrib;
        pAttrib.add(pname,"name");
        pAttrib.put(cset[i]);
	if(pname == "energy") 
	  putContent(w_en,cset[i]);
	else if(pname == "variance") 
	  putContent(w_var,cset[i]);
        else if(pname == "difference")
	  putContent(w_abs,cset[i]);
      }
    }  
    //    LogOut->getStream() << "#" << optmethod << " values: tolerance = "
    //			<< cg_tolerance << " stepsize = " << cg_stepsize 
    //			<< " epsilon = " << cg_epsilon << " Tau = " 
    //			<< Tau << endl;
    //    if(!UseWeight){
    //      LogOut->getStream() << "#All weights set to 1.0" << endl;
    //    }

    putOptParams();

    //maybe overwritten but will try out
    EtargetEff=(1.0+CorrelationFactor)*Etarget;
    return true;
  }

  bool
  QMCCostFunction::putOptParams(){

    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++){
      //locate the id in the variable list for Psi
      int id = Psi.VarList.find(IDtag[i]);
      if(id>= 0) {
	//find the number of variables represented by the id
	int size = Psi.VarList.Sizes[id];
	Return_t* temp = Psi.VarList.Pointers[id];
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

    if(msg_stream) {
      msg_stream->setf(ios::scientific, ios::floatfield);
      msg_stream->precision(8);
      *msg_stream << " Inital Variables= " ;
      std::copy(OptParams.begin(),OptParams.end(), ostream_iterator<Return_t>(*msg_stream," "));
      *msg_stream << endl;
    }

    Psi.reset();
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

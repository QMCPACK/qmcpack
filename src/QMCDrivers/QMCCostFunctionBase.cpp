//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCCostFunctionBase.h"
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/CommOperators.h"
#include "Optimize/LeastSquaredFit.h"
#include <set>
//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus
{

QMCCostFunctionBase::QMCCostFunctionBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
  MPIObjectBase(0),
  W(w),H(h),Psi(psi),  Write2OneXml(true),
  PowerE(2), NumCostCalls(0), NumSamples(0), MaxWeight(1e6),
  w_en(0.9), w_var(0.1), w_abs(0.0), w_w(0.0), w_beta(0.0), GEVType("mixed"),
  CorrelationFactor(0.0), m_wfPtr(NULL), m_doc_out(NULL), msg_stream(0), debug_stream(0),
  SmallWeight(0), includeNonlocalH("no"),needGrads(true), vmc_or_dmc(2.0),
  targetExcitedStr("no"), targetExcited(false), omega_shift(0.0)
{
  GEVType="mixed";
  //paramList.resize(10);
  //costList.resize(10,0.0);
  //default: don't check fo MinNumWalkers
  MinNumWalkers = 0.3;
  SumValue.resize(SUM_INDEX_SIZE,0.0);
  IsValid=true;
  useNLPPDeriv=false;
#if defined(QMCCOSTFUNCTION_DEBUG)
  char fname[16];
  sprintf(fname,"optdebug.p%d",OHMMS::Controller->mycontext());
  debug_stream = new std::ofstream(fname);
  debug_stream->setf(std::ios::scientific, std::ios::floatfield);
  debug_stream->precision(8);
#endif
}

/** Clean up the vector */
QMCCostFunctionBase::~QMCCostFunctionBase()
{
  delete_iter(dLogPsi.begin(),dLogPsi.end());
  delete_iter(d2LogPsi.begin(),d2LogPsi.end());
//     if (m_doc_out != NULL) xmlFreeDoc(m_doc_out);
  if (debug_stream)
    delete debug_stream;
}

void QMCCostFunctionBase::setRng(std::vector<RandomGenerator_t*>& r)
{
  if(MoverRng.size()<r.size())
  {
    delete_iter(MoverRng.begin(),MoverRng.end());
    delete_iter(RngSaved.begin(),RngSaved.end());
    MoverRng.resize(r.size());
    RngSaved.resize(r.size());
  }
  for(int ip=0; ip<r.size(); ++ip)
    MoverRng[ip]=r[ip];
  for(int ip=0; ip<r.size(); ++ip)
    RngSaved[ip] = new RandomGenerator_t(*MoverRng[ip]);
}

void QMCCostFunctionBase::setTargetEnergy(Return_t et)
{
  //evaluate effective target energy
  EtargetEff= Etarget = et;
//     app_log() << "Effective Target Energy = " << EtargetEff << std::endl;
//     app_log() << "Cost Function = " << w_en << "*<E> + "
//     << w_var << "*<Var> + " << w_w << "*<Var(unreweighted)> " << std::endl;
  //if(UseWeight)
  //app_log() << "Correlated sampling is used." << std::endl;
  //else
  //app_log() << "Weight is set to one." << std::endl;
//     if (msg_stream)
//       {
//         *msg_stream << "  Total number of walkers          = " << NumSamples << std::endl;
//         *msg_stream << "  Effective Target Energy = " << EtargetEff << std::endl;
//         *msg_stream << "  Cost Function = " << w_en << "*<E> + " << w_var << "*<Var> + " << w_w << "*<Var(unreweighted)> " << std::endl;
//         *msg_stream << "  Optimization report = ";
//         *msg_stream << "cost, walkers, eavg/wgt, eavg/walkers, evar/wgt, evar/walkers, evar_abs\n";
//         *msg_stream << "  Optimized variables = " << OptVariables.name(0);
//         for (int i=1; i<OptVariables.size(); ++i) *msg_stream << "," << OptVariables.name(i) ;
//         *msg_stream << std::endl;
//       }
}

QMCCostFunctionBase::Return_t QMCCostFunctionBase::Cost(bool needGrad)
{
  NumCostCalls++;
//reset the wave function
  resetPsi();
//evaluate new local energies
  NumWalkersEff=correlatedSampling(needGrad);
  return computedCost();
}

void QMCCostFunctionBase::printEstimates()
{
  app_log()<<"      Current ene:     " <<curAvg_w<< std::endl;
  app_log()<<"      Current var:     " <<curVar_w<< std::endl;
  app_log()<<"      Current ene_urw:     " <<curAvg<< std::endl;
  app_log()<<"      Current var_urw: " <<curVar<< std::endl;
}

QMCCostFunctionBase::Return_t QMCCostFunctionBase::computedCost()
{
//   Assumes the Sums have been computed all ready
  //Estimators::accumulate has been called by correlatedSampling
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  curVar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-curAvg_w*curAvg_w;
  Return_t wgtinv=1.0/static_cast<Return_t>(NumSamples);
  // app_log() << "wgtinv = " << wgtinv << std::endl;
  curAvg = SumValue[SUM_E_BARE]*wgtinv;
  curVar = SumValue[SUM_ESQ_BARE]*wgtinv-curAvg*curAvg;
  curVar_abs = SumValue[SUM_ABSE_WGT]/SumValue[SUM_WGT];
  // app_log() << "curVar     = " << curVar
  //     << "   curAvg     = " << curAvg
  //     << "   NumWalkersEff     = " << NumWalkersEff << std::endl;
  // app_log() << "SumValue[SUM_WGT] = " << SumValue[SUM_WGT] << std::endl;
  // app_log() << "SumValue[SUM_WGTSQ] = " << SumValue[SUM_WGTSQ] << std::endl;
  // app_log() << "SumValue[SUM_ABSE_WGT] = " << SumValue[SUM_ABSE_WGT] << std::endl;
  // app_log() << "SumValue[SUM_E_WGT] = " << SumValue[SUM_E_WGT] << std::endl;
  // app_log() << "SumValue[SUM_ESQ_WGT] = " << SumValue[SUM_ESQ_WGT] << std::endl;
  Return_t wgt_var = SumValue[SUM_WGTSQ]-SumValue[SUM_WGT]*SumValue[SUM_WGT];
  wgt_var *=wgtinv;
  CostValue = 0.0;
  const Return_t small=1.0e-10;
  if (std::abs(w_abs) > small)
    CostValue += w_abs*curVar_abs;
  if (std::abs(w_var) > small)
    CostValue += w_var*curVar_w;
  if (std::abs(w_en) > small)
    CostValue += w_en*curAvg_w;
  if (std::abs(w_w) > small)
    CostValue += w_w*curVar;
  //CostValue = w_abs*curVar_abs + w_var*curVar_w + w_en*curAvg_w + w_w*curVar;
  // app_log() << "CostValue, NumEffW = " << CostValue <<"  " <<NumWalkersEff << std::endl;
  IsValid=true;
  if(NumWalkersEff < NumSamples*MinNumWalkers)
    //    if (NumWalkersEff < MinNumWalkers)
  {
    ERRORMSG("CostFunction-> Number of Effective Walkers is too small! " << std::endl
          << "  Number of effective walkers (samples) / total number of samples = " << (1.0*NumWalkersEff)/NumSamples << std::endl
          << "  User specified threshold minwalkers = " << MinNumWalkers << std::endl
          << "  If this message appears frequently. You might have to be cautious. " << std::endl
          << "  Find info about parameter \"minwalkers\" in the user manual!");
    // ERRORMSG("Going to stop now.")
    IsValid=false;
  }
  return CostValue;
}

//  /**  Perform the correlated sampling algorthim.
//   */
//  QMCCostFunctionBase::Return_t QMCCostFunctionBase::correlatedSampling() {
//
//    typedef MCWalkerConfiguration::Walker_t Walker_t;
//    MCWalkerConfiguration::iterator it(W.begin());
//    MCWalkerConfiguration::iterator it_end(W.end());
//    Return_t eloc_new;
//    int iw=0;
//    Return_t wgt_tot=0.0;
//
//    while(it != it_end) {
//      Walker_t& thisWalker(**it);
//      Return_t* restrict saved = Records[iw];
//
//      //rewind the buffer to get the data from buffer
//      thisWalker.DataSet.rewind();
//      //W is updated by thisWalker
//      W.copyFromBuffer(thisWalker.DataSet);
//
//      Return_t logpsi=0.0;
//      //copy dL from Buffer
//      thisWalker.DataSet.get(dL.begin(),dL.end());
//      logpsi=Psi.evaluateDeltaLog(W);
//      W.G += thisWalker.Drift;
//      W.L += dL;
//
//      eloc_new=H_KE.evaluate(W)+saved[ENERGY_FIXED];
//      Return_t weight = UseWeight?exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
//
//      saved[ENERGY_NEW]=eloc_new;
//      saved[REWEIGHT]=weight;
//      wgt_tot+=weight;
//      ++it;
//      ++iw;
//    }
//
//    OHMMS::Controller->barrier();
//    //collect the total weight for normalization and apply maximum weight
//    myComm->allreduce(wgt_tot);
//
//    for(int i=0; i<SumValue.size(); i++) SumValue[i]=0.0;
//    wgt_tot=1.0/wgt_tot;
//
//    Return_t wgt_max=MaxWeight*wgt_tot;
//    int nw=W.getActiveWalkers();
//    for(iw=0; iw<nw;iw++) {
//      Return_t* restrict saved = Records[iw];
//      Return_t weight=saved[REWEIGHT]*wgt_tot;
//      Return_t eloc_new=saved[ENERGY_NEW];
//
//      weight = (weight>wgt_max)? wgt_max:weight;
//      Return_t delE=pow(std::abs(eloc_new-EtargetEff),PowerE);
//      SumValue[SUM_E_BARE] += eloc_new;
//      SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
//      SumValue[SUM_ABSE_BARE] += delE;
//      SumValue[SUM_E_WGT] += eloc_new*weight;
//      SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*weight;
//      SumValue[SUM_ABSE_WGT] += delE*weight;
//      SumValue[SUM_WGT] += weight;
//      SumValue[SUM_WGTSQ] += weight*weight;
//    }
//
//    //collect everything
//    myComm->allreduce(SumValue);
//
//    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
//  }

//  void
//  QMCCostFunctionBase::getConfigurations(const std::string& aroot) {
//    if(aroot.size() && aroot != "invalid") {
//      app_log() << "  Reading configurations from the previous qmc block" << std::endl;
//      HDFWalkerInputCollect wReader(aroot);
//      wReader.putSingle(W);
//    }
//
//#if defined(QMCCOSTFUNCTION_DEBUG)
//    if(debug_stream) delete debug_stream;
//    char fname[16];
//    sprintf(fname,"optdebug.p%d",OHMMS::Controller->mycontext());
//    debug_stream = new std::ofstream(fname);
//    debug_stream->setf(std::ios::scientific, std::ios::floatfield);
//    debug_stream->precision(8);
//
//    *debug_stream << "Initial : " << std::endl;
//    for(int i=0; i<OptParams.size(); i++)
//      *debug_stream << " " << IDtag[i] << "=" << OptParams[i];
//    *debug_stream << std::endl;
//#endif
//  }

//  /** evaluate everything before optimization */
//  void
//  QMCCostFunctionBase::getConfigurations(std::vector<std::string>& ConfigFile,
//      int partid, int nparts) {
//    if(ConfigFile.size()) {
//
//      app_log() << "  Reading configurations from mcwalkerset " << std::endl;
//
//      W.destroyWalkers(W.begin(),W.end());
//      for(int i=0; i<ConfigFile.size(); i++) {
//        //JNKIM: needs to change to HDFWalkerInputCollect
//        //HDFWalkerInput0 wReader(ConfigFile[i],partid,nparts);
//        HDFWalkerInputCollect wReader(ConfigFile[i]);
//        wReader.putSingle(W);
//        //wReader.put(W,-1);
//      }
//
//      //remove input files
//      ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
//    }
//  }

//  /** evaluate everything before optimization */
//  void
//  QMCCostFunctionBase::checkConfigurations() {
//
//    dL.resize(W.getTotalNum());
//    int numLocWalkers=W.getActiveWalkers();
//    Records.resize(numLocWalkers,6);
//
//    typedef MCWalkerConfiguration::Walker_t Walker_t;
//    MCWalkerConfiguration::iterator it(W.begin());
//    MCWalkerConfiguration::iterator it_end(W.end());
//    int nat = W.getTotalNum();
//    int iw=0;
//    Etarget=0.0;
//    while(it != it_end) {
//
//      Walker_t& thisWalker(**it);
//
//      //clean-up DataSet to save re-used values
//      thisWalker.DataSet.clear();
//      //rewind the counter
//      thisWalker.DataSet.rewind();
//      //MCWalkerConfiguraton::registerData add distance-table data
//      W.registerData(thisWalker,thisWalker.DataSet);
//
//      Return_t*  saved=Records[iw];
//#if defined(QMC_COMPLEX)
//      app_error() << " Optimization is not working with complex wavefunctions yet" << std::endl;
//      app_error() << "  Needs to fix TrialWaveFunction::evaluateDeltaLog " << std::endl;
//#else
//      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], thisWalker.Drift, dL);
//#endif
//      thisWalker.DataSet.add(dL.first_address(),dL.last_address());
//      Etarget += saved[ENERGY_TOT] = H.evaluate(W);
//      saved[ENERGY_FIXED] = H.getInvariantEnergy();
//
//      ++it;
//      ++iw;
//    }
//
//    //Need to sum over the processors
//    std::vector<Return_t> etemp(2);
//    etemp[0]=Etarget;
//    etemp[1]=static_cast<Return_t>(numLocWalkers);
//    myComm->allreduce(etemp);
//    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
//    NumSamples = static_cast<int>(etemp[1]);
//
//    setTargetEnergy(Etarget);
//
//    ReportCounter=0;
//  }

//  void QMCCostFunctionBase::resetPsi()
//  {
//
//    OptimizableSetType::iterator oit(OptVariables.begin()), oit_end(OptVariables.end());
//    while(oit != oit_end)
//    {
//      Return_t v=(*oit).second;
//      OptVariablesForPsi[(*oit).first]=v;
//      std::map<std::string,std::set<std::string>*>::iterator eit(equalConstraints.find((*oit).first));
//      if(eit != equalConstraints.end())
//      {
//        std::set<std::string>::iterator f((*eit).second->begin()),l((*eit).second->end());
//        while(f != l)
//        {
//          OptVariablesForPsi[(*f)]=v;
//          ++f;
//        }
//      }
//      ++oit;
//    }
//
//    //cout << "QMCCostFunctionBase::resetPsi " << std::endl;
//    //oit=OptVariablesForPsi.begin();
//    //oit_end=OptVariablesForPsi.end();
//    //while(oit != oit_end)
//    //{
//    //  std::cout << (*oit).first << "=" << (*oit).second << " ";
//    //  ++oit;
//    //}
//    //cout << std::endl;
//    Psi.resetParameters(OptVariablesForPsi);
//  }

/** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
 *@return true always
 *
 * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
 * set \f$ {{\bf \alpha_{i+1}}}\f$
 */
//bool
//QMCCostFunctionBase::resetWaveFunctions() {

//  //loop over all the unique id's
//  for(int i=0; i<IDtag.size(); i++)
//  {
//    OptVariables[IDtag[i]]=OptParams[i];
//  }
//  resetPsi();
//  return true;
//}

void QMCCostFunctionBase::Report()
{
  //reset the wavefunction for with the new variables
  resetPsi();
  if (!myComm->rank())
  {
    updateXmlNodes();
    char newxml[128];
    if (Write2OneXml)
      sprintf(newxml,"%s.opt.xml", RootName.c_str());
    else
      sprintf(newxml,"%s.opt.%d.xml", RootName.c_str(),ReportCounter);
    xmlSaveFormatFile(newxml,m_doc_out,1);
    if (msg_stream)
    {
      msg_stream->precision(8);
      *msg_stream << " curCost "
                  << std::setw(5) << ReportCounter << std::setw(16) << CostValue
                  << std::setw(16) << NumWalkersEff
                  << std::setw(16) << curAvg_w << std::setw(16) << curAvg
                  << std::setw(16) << curVar_w << std::setw(16) << curVar
                  << std::setw(16) << curVar_abs << std::endl;
      *msg_stream << " curVars " << std::setw(5) << ReportCounter;
      for (int i=0; i<OptVariables.size(); i++)
        *msg_stream << std::setw(16) << OptVariables[i];
      *msg_stream << std::endl;
    }
    //report the data
    //Psi.reportStatus(*mgs_stream);
  }
#if defined(QMCCOSTFUNCTION_DEBUG)
  *debug_stream << ReportCounter;
  OptVariables.print(*debug_stream);
  *debug_stream << std::endl;
#endif
  ReportCounter++;
  //myComm->barrier();
}

void QMCCostFunctionBase::reportParameters()
{
  //final reset, restoring the OrbitalBase::IsOptimizing to false
  resetPsi(true);
  if (!myComm->rank())
  {
    char newxml[128];
    sprintf(newxml,"%s.opt.xml", RootName.c_str());
    *msg_stream << "  <optVariables href=\"" << newxml << "\">" << std::endl;
    OptVariables.print(*msg_stream);
    *msg_stream << "  </optVariables>" << std::endl;
    updateXmlNodes();
    xmlSaveFormatFile(newxml,m_doc_out,1);
  }
}

/** Apply constraints on the optimizables.
 *
 * Here is where constraints should go
 */
bool QMCCostFunctionBase::checkParameters()
{
  bool samesign = true;
  //if(samesign) {
  //  paramList.pop_back();
  //  paramList.push_front(OptParams);
  //}
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
QMCCostFunctionBase::put(xmlNodePtr q)
{
  std::string writeXmlPerStep("no");
  std::string computeNLPPderiv("no");
  ParameterSet m_param;
  m_param.add(writeXmlPerStep,"dumpXML","string");
  m_param.add(MinNumWalkers,"minwalkers","scalar");
  m_param.add(MaxWeight,"maxWeight","scalar");
  m_param.add(includeNonlocalH,"nonlocalpp","string");
  m_param.add(computeNLPPderiv,"use_nonlocalpp_deriv","string");
  m_param.add(w_beta,"beta","double");
  m_param.add(GEVType,"GEVMethod","string");
  m_param.add(targetExcitedStr,"targetExcited","string");
  m_param.add(omega_shift,"omega","double");
  m_param.put(q);

  tolower(targetExcitedStr);
  targetExcited = ( targetExcitedStr == "yes" );

  if (includeNonlocalH=="yes")
    includeNonlocalH="NonLocalECP";

  if(computeNLPPderiv != "no" && includeNonlocalH != "no")
  {
    app_log() << "   Going to include the derivatives of " << includeNonlocalH << std::endl;
    useNLPPDeriv=true;
  }
  // app_log() << "  QMCCostFunctionBase::put " << std::endl;
  // m_param.get(app_log());
  Write2OneXml = (writeXmlPerStep == "no");
  xmlNodePtr qsave=q;
  //Estimators.put(q);
  std::vector<xmlNodePtr> cset;
  std::vector<std::string> excluded;
  std::map<std::string,std::vector<std::string>*> equalConstraints;
  std::map<std::string,std::vector<std::string>*> negateConstraints;
  std::vector<std::string> idtag;
  xmlNodePtr cur=qsave->children;
  int pid=myComm->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "optimize")
    {
      std::vector<std::string> tmpid;
      putContent(tmpid,cur);
      idtag.insert(idtag.end(),tmpid.begin(),tmpid.end());
    }
    else if (cname == "exclude")
    {
      std::vector<std::string> tmpid;
      putContent(tmpid,cur);
      excluded.insert(excluded.end(),tmpid.begin(),tmpid.end());
    }
    else if (cname == "cost")
    {
      cset.push_back(cur);
    }
    else if (cname == "set")
    {
      std::string ctype("equal");
      std::string s("0");
      OhmmsAttributeSet pAttrib;
      pAttrib.add(ctype,"type");
      pAttrib.add(s,"name");
      pAttrib.put(cur);
      if (ctype == "equal" || ctype == "=")
      {
        std::map<std::string,std::vector<std::string>*>::iterator eit(equalConstraints.find(s));
        std::vector<std::string>* eqSet=0;
        if (eit == equalConstraints.end())
        {
          eqSet = new std::vector<std::string>;
          equalConstraints[s]=eqSet;
        }
        else
          eqSet = (*eit).second;
        std::vector<std::string> econt;
        putContent(econt,cur);
        eqSet->insert(eqSet->end(),econt.begin(),econt.end());
      }
    }
    cur=cur->next;
  }
  //build optimizables from the wavefunction
  OptVariablesForPsi.clear();
  Psi.checkInVariables(OptVariablesForPsi);
  OptVariablesForPsi.resetIndex();
  //synchronize OptVariables and OptVariablesForPsi
  OptVariables=OptVariablesForPsi;
  //first disable <exclude>.... </exclude> from the big list used by a TrialWaveFunction
  OptVariablesForPsi.disable(excluded.begin(),excluded.end(),false);
  //now, set up the real variable list for optimization
  //check <equal>
  int nc=0;
  if (equalConstraints.size())
  {
    std::map<std::string,std::vector<std::string>*>::iterator eit(equalConstraints.begin());
    while (eit != equalConstraints.end())
    {
      nc+=(*eit).second->size();
      //actiave the active variable even though it is probably unnecessary
      OptVariablesForPsi.activate((*eit).second->begin(),(*eit).second->end(),false);
      excluded.insert(excluded.end(),(*eit).second->begin(),(*eit).second->end());
      ++eit;
    }
  }
  //build OptVariables which is equal to or identical to OptVariablesForPsi
  //disable the variables that are equal to a variable
  OptVariables.disable(excluded.begin(),excluded.end(),false);
  //set up OptVariables and OptVariablesForPsi
  OptVariables.activate(idtag.begin(),idtag.end(),true);
  OptVariablesForPsi.activate(idtag.begin(),idtag.end(),true);
  //found constraints build equalVarMap
  if (nc>0)
  {
    equalVarMap.reserve(nc+OptVariables.size());
    //map the basic lists from the active list
    for (int i=0; i<OptVariables.size(); ++i)
    {
      int bigloc=OptVariablesForPsi.getIndex(OptVariables.name(i));
      if (bigloc<0)
        continue;
      equalVarMap.push_back(TinyVector<int,2>(bigloc,i));
    }
    //add <equal/>
    std::map<std::string,std::vector<std::string>*>::iterator eit(equalConstraints.begin());
    while (eit != equalConstraints.end())
    {
      int loc=OptVariables.getIndex((*eit).first);
      if (loc>=0)
      {
        const std::vector<std::string>& elist(*((*eit).second));
        for (int i=0; i<elist.size(); ++i)
        {
          int bigloc=OptVariablesForPsi.getIndex(elist[i]);
          if (bigloc<0)
            continue;
          equalVarMap.push_back(TinyVector<int,2>(bigloc,loc));
        }
      }
      //remove std::vector<std::string>
      delete(*eit).second;
      ++eit;
    }
  }
  //get the indices
  Psi.checkOutVariables(OptVariablesForPsi);
  NumOptimizables=OptVariables.size();
  if (NumOptimizables == 0)
  {
    APP_ABORT("QMCCostFunctionBase::put No valid optimizable variables are found.");
  }
//     app_log() << "<active-optimizables> " << std::endl;
//     OptVariables.print(app_log());
//     app_log() << "</active-optimizables>" << std::endl;
  if (msg_stream)
    msg_stream->setf(std::ios::scientific, std::ios::floatfield);
  resetPsi();
//     Psi.reportStatus(app_log());
  //set the cost function
  if (cset.empty())
  {
    if (msg_stream)
      *msg_stream << " Using Default Cost Function: Cost = <|E-E_ff|^2>" << std::endl;
  }
  else
    resetCostFunction(cset);
  //maybe overwritten but will try out
  EtargetEff=Etarget;
  return true;
}

void QMCCostFunctionBase::resetCostFunction(std::vector<xmlNodePtr>& cset)
{
  for (int i=0; i<cset.size(); i++)
  {
    std::string pname;
    Return_t wgt=1.0;
    OhmmsAttributeSet pAttrib;
    pAttrib.add(pname,"name");
    pAttrib.put(cset[i]);
    if (pname == "energy")
      putContent(w_en,cset[i]);
    else if ((pname == "variance") || (pname == "unreweightedvariance"))
      putContent(w_w,cset[i]);
    else if (pname == "difference")
      putContent(w_abs,cset[i]);
    else if ((pname == "reweightedVariance") || (pname == "reweightedvariance"))
      putContent(w_var,cset[i]);
  }
}


void QMCCostFunctionBase::updateXmlNodes()
{
  if (m_doc_out == NULL) //first time, create a document tree and get parameters and attributes to be updated
  {
    m_doc_out = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    xmlAddChild(qm_root, m_wfPtr);
    xmlDocSetRootElement(m_doc_out, qm_root);
    xmlXPathContextPtr acontext = xmlXPathNewContext(m_doc_out);
    //check var
    xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//var",acontext);
    for (int iparam=0; iparam<result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
      const xmlChar* iptr=xmlGetProp(cur,(const xmlChar*)"id");
      if (iptr == NULL)
        continue;
      std::string aname((const char*)iptr);
      opt_variables_type::iterator oit(OptVariablesForPsi.find(aname));
      if (oit != OptVariablesForPsi.end())
      {
        paramNodes[aname]=cur;
      }
    }
    xmlXPathFreeObject(result);
    //check radfunc
    result = xmlXPathEvalExpression((const xmlChar*)"//radfunc",acontext);
    for (int iparam=0; iparam<result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
      const xmlChar* iptr=xmlGetProp(cur,(const xmlChar*)"id");
      if (iptr == NULL)
        continue;
      std::string aname((const char*)iptr);
      std::string expID=aname+"_E";
      xmlAttrPtr aptr=xmlHasProp(cur,(const xmlChar*)"exponent");
      opt_variables_type::iterator oit(OptVariablesForPsi.find(expID));
      if (aptr != NULL && oit != OptVariablesForPsi.end())
      {
        attribNodes[expID]= std::pair<xmlNodePtr,std::string>(cur,"exponent");
      }
      std::string cID=aname+"_C";
      aptr=xmlHasProp(cur,(const xmlChar*)"contraction");
      oit=OptVariablesForPsi.find(cID);
      if (aptr != NULL && oit != OptVariablesForPsi.end())
      {
        attribNodes[cID]= std::pair<xmlNodePtr,std::string>(cur,"contraction");
      }
    }
    xmlXPathFreeObject(result);
    //check ci
    result = xmlXPathEvalExpression((const xmlChar*)"//ci",acontext);
    for (int iparam=0; iparam<result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
      const xmlChar* iptr=xmlGetProp(cur,(const xmlChar*)"id");
      if (iptr == NULL)
        continue;
      std::string aname((const char*)iptr);
      xmlAttrPtr aptr=xmlHasProp(cur,(const xmlChar*)"coeff");
      opt_variables_type::iterator oit(OptVariablesForPsi.find(aname));
      if (aptr != NULL && oit != OptVariablesForPsi.end())
      {
        attribNodes[aname]= std::pair<xmlNodePtr,std::string>(cur,"coeff");
      }
    }
    xmlXPathFreeObject(result);
    //check csf
    result = xmlXPathEvalExpression((const xmlChar*)"//csf",acontext);
    for (int iparam=0; iparam<result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
      const xmlChar* iptr=xmlGetProp(cur,(const xmlChar*)"id");
      if (iptr == NULL)
        continue;
      std::string aname((const char*)iptr);
      xmlAttrPtr aptr=xmlHasProp(cur,(const xmlChar*)"coeff");
      opt_variables_type::iterator oit(OptVariablesForPsi.find(aname));
      if (aptr != NULL && oit != OptVariablesForPsi.end())
      {
        attribNodes[aname]= std::pair<xmlNodePtr,std::string>(cur,"coeff");
      }
    }
    xmlXPathFreeObject(result);
    addCoefficients(acontext, "//coefficient");
    addCoefficients(acontext, "//coefficients");
    xmlXPathFreeContext(acontext);
  }
  //     Psi.reportStatus(app_log());
  std::map<std::string,xmlNodePtr>::iterator pit(paramNodes.begin()), pit_end(paramNodes.end());
  while (pit != pit_end)
  {
    Return_t v=OptVariablesForPsi[(*pit).first];
    getContent(v,(*pit).second);
    //         vout <<(*pit).second<< std::endl;
    ++pit;
  }
  std::map<std::string,std::pair<xmlNodePtr,std::string> >::iterator ait(attribNodes.begin()), ait_end(attribNodes.end());
  while (ait != ait_end)
  {
    std::ostringstream vout;
    vout.setf(std::ios::scientific, std::ios::floatfield);
    vout.precision(16);
    vout << OptVariablesForPsi[(*ait).first];
    xmlSetProp((*ait).second.first, (const xmlChar*)(*ait).second.second.c_str(),(const xmlChar*)vout.str().c_str());
    ++ait;
  }
  std::map<std::string,xmlNodePtr>::iterator cit(coeffNodes.begin()), cit_end(coeffNodes.end());
  while (cit != cit_end)
  {
    std::string rname((*cit).first);
    OhmmsAttributeSet cAttrib;
    std::string datatype("none");
    std::string aname("0");
    cAttrib.add(datatype,"type");
    cAttrib.add(aname,"id");
    cAttrib.put((*cit).second);
    if (datatype == "Array")
    {
      //
      aname.append("_");
      opt_variables_type::iterator vit(OptVariablesForPsi.begin());
      std::vector<Return_t> c;
      while (vit != OptVariablesForPsi.end())
      {
        if ((*vit).first.find(aname) == 0)
        {
          c.push_back((*vit).second);
        }
        ++vit;
      }
      xmlNodePtr contentPtr = cit->second;
      if (xmlNodeIsText(contentPtr->children))
        contentPtr = contentPtr->children;
      getContent(c,contentPtr);
    }
    else
    {
      xmlNodePtr cur=(*cit).second->children;
      while (cur != NULL)
      {
        std::string cname((const char*)(cur->name));
        if (cname == "lambda")
        {
          int i=0;
          int j=-1;
          OhmmsAttributeSet pAttrib;
          pAttrib.add(i,"i");
          pAttrib.add(j,"j");
          pAttrib.put(cur);
          char lambda_id[32];
          if (j<0)
            sprintf(lambda_id,"%s_%d",rname.c_str(),i);
          else
            sprintf(lambda_id,"%s_%d_%d",rname.c_str(),i,j);
          opt_variables_type::iterator vTarget(OptVariablesForPsi.find(lambda_id));
          if (vTarget != OptVariablesForPsi.end())
          {
            std::ostringstream vout;
            vout.setf(std::ios::scientific, std::ios::floatfield);
            vout.precision(16);
            vout << (*vTarget).second;
            xmlSetProp(cur, (const xmlChar*)"c",(const xmlChar*)vout.str().c_str());
          }
        }
        cur=cur->next;
      }
    }
    ++cit;
  }
}

/** add coefficient or coefficients
 * @param acontext context from which xpath cname is searched
 * @param cname xpath
 */
void QMCCostFunctionBase::addCoefficients(xmlXPathContextPtr acontext, const char* cname)
{
  xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)cname,acontext);
  for (int iparam=0; iparam<result->nodesetval->nodeNr; iparam++)
  {
    xmlNodePtr cur= result->nodesetval->nodeTab[iparam];
    OhmmsAttributeSet cAttrib;
    std::string aname("0");
    std::string optimize("yes");
    std::string datatype("none");
    cAttrib.add(aname,"id");
    cAttrib.add(aname,"name");
    cAttrib.add(datatype,"type");
    cAttrib.add(optimize, "optimize");
    cAttrib.put(cur);
    if (aname[0] == '0')
      continue;
    if (datatype == "Array")
    {
      if(optimize == "yes")
        coeffNodes[aname]=cur;
    }
    else
    {
      //check if any optimizables contains the id of coefficients
      bool notlisted=true;
      opt_variables_type::iterator oit(OptVariablesForPsi.begin()),oit_end(OptVariablesForPsi.end());
      while (notlisted && oit != oit_end)
      {
        const std::string& oname((*oit).first);
        notlisted=(oname.find(aname)>=oname.size());
        ++oit;
      }
      if (!notlisted)
      {
        coeffNodes[aname]=cur;
      }
    }
  }
  xmlXPathFreeObject(result);
}

bool
QMCCostFunctionBase::lineoptimization(const std::vector<Return_t>& x0, const std::vector<Return_t>& gr, Return_t val0,
                                      Return_t&  dl, Return_t& val_proj, Return_t& lambda_max)
{
  return false;
  const int maxclones=3;
  const int max_poly=3;
  //Matrix<Return_t> js(maxclones+1,x0.size());
  Vector<Return_t> y(maxclones+1);
  Vector<Return_t> sigma(maxclones+1);
  Matrix<Return_t> A(maxclones+1,max_poly);
  sigma=1e-6; //a small value
  Return_t gr_norm=0.0;
  for (int j=0; j<x0.size(); ++j)
  {
    //js(0,j)=x0[j];
    gr_norm+=gr[j]*gr[j];
  }
  Return_t nw=1.0/static_cast<QMCTraits::RealType>(NumSamples);
  //Return_t MaxDispl=0.04;
  gr_norm=std::sqrt(gr_norm);
  Return_t dx=lambda_max/gr_norm;
  dx=std::min((QMCTraits::RealType)0.25,dx);
  if (val0<1e12)
  {
    y[0]=val0;
    sigma[0]=std::sqrt(val0)*nw;
  }
  else
  {
    for (int j=0; j<x0.size(); ++j)
      Params(j)=x0[j];
    Return_t costval=Cost();
    y[0]=costval;
    sigma[0]=std::sqrt(costval)*nw;
  }
  app_log() << "  lineoptimization (" << 0.0 << "," << y[0] << ")";
  for (int k=0; k<max_poly; ++k)
    A(0,k)=0.0;
  Return_t dxmax=0.0;
  for (int i=1; i<=maxclones; ++i)
  {
    dxmax+=dx;
    //set OptParams to vary
    for (int j=0; j<x0.size(); ++j)
    {
      //js(i,j)=OptParams[j]=x0[j]+dx*gr[j];
      Params(j)=x0[j]+dxmax*gr[j];
    }
    Return_t costval=Cost();
    y[i]=costval;
    sigma[i]=std::sqrt(costval)*nw;
    for (int k=0; k<max_poly; ++k)
      A(i,k)=std::pow(dxmax,k);
    app_log() << " (" << dxmax << "," << y[i] << ")";
  }
  app_log() << std::endl;
  Vector<QMCTraits::RealType> polys(max_poly);
  Vector<QMCTraits::RealType> errors(max_poly);
  LeastSquaredFitLU(y,sigma,A,polys,errors);
  dl=-polys[1]/polys[2]*0.5;
  val_proj=polys[0]+dl*(polys[1]+dl*polys[2]);
  if (dl<dx*0.25)
    lambda_max *=0.5; // narrow the bracket
  if (dl>dxmax*5.0)
    lambda_max *=2.0; //widen the bracket
  app_log() << "    minimum at " << dl << " estimated=" << val_proj << " LambdaMax " << lambda_max;
  for (int j=0; j<x0.size(); ++j)
    Params(j)=x0[j]+dl*gr[j];
  val_proj=Cost();
  app_log() << "  cost = " << val_proj << std::endl;
  //Psi.reportStatus(app_log());
  //return false;
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  If the LMYEngine is available, returns the cost function calculated by the engine.
///         Otherwise, returns the usual cost function.
///
/// \param[in]      needDeriv             whether derivative vectors should be computed
///
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LMY_ENGINE
QMCCostFunctionBase::Return_t QMCCostFunctionBase::LMYEngineCost(const bool needDeriv, cqmc::engine::LMYEngine * EngineObj) {

  // prepare local energies, weights, and possibly derivative vectors, and compute standard cost
  const Return_t standardCost = this->Cost(needDeriv);

  // if we are using the LMYEngine, compute and return it's cost function value
  return this->LMYEngineCost_detail(EngineObj);

  // otherwise return the standard cost function
  return standardCost;

}
#endif

}

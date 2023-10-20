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


#include "QMCCostFunctionBase.h"
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/XMLParsingString.h"
#include "Message/CommOperators.h"
#include "Message/UniformCommunicateError.h"

#include <array>

//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus
{
QMCCostFunctionBase::QMCCostFunctionBase(ParticleSet& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm)
    : MPIObjectBase(comm),
      reportH5(false),
      CI_Opt(false),
      W(w),
      Psi(psi),
      H(h),
      Write2OneXml(true),
      PowerE(2),
      NumCostCalls(0),
      NumSamples(0),
      w_en(0.9),
      w_var(0.1),
      w_abs(0.0),
      w_w(0.0),
      MaxWeight(1e6),
      w_beta(0.0),
      GEVType("mixed"),
      vmc_or_dmc(2.0),
      needGrads(true),
      targetExcitedStr("no"),
      targetExcited(false),
      omega_shift(0.0),
      msg_stream(nullptr),
      m_wfPtr(nullptr),
      m_doc_out(nullptr),
      do_override_output(true)
{
  //default: don't check fo MinNumWalkers
  MinNumWalkers = 0.3;
  SumValue.resize(SUM_INDEX_SIZE, 0.0);
  IsValid = true;
#if defined(QMCCOSTFUNCTION_DEBUG)
  std::array<char, 16> fname;
  int length = std::snprintf(fname.data(), fname.size(), "optdebug.p%d", OHMMS::Controller->rank());
  if (length < 0)
    throw std::runtime_error("Error generating filename");
  debug_stream = std::make_unique<std::ofstream>(fname.data());
  debug_stream->setf(std::ios::scientific, std::ios::floatfield);
  debug_stream->precision(8);
#endif
}

/** Clean up the vector */
QMCCostFunctionBase::~QMCCostFunctionBase()
{
  delete_iter(dLogPsi.begin(), dLogPsi.end());
  delete_iter(d2LogPsi.begin(), d2LogPsi.end());
  if (m_doc_out != nullptr)
    xmlFreeDoc(m_doc_out);
  debug_stream.reset();
}

void QMCCostFunctionBase::setRng(RefVector<RandomBase<FullPrecRealType>> r)
{
  if (MoverRng.size() < r.size())
  {
    MoverRng.resize(r.size());
    RngSaved.resize(r.size());
  }
  for (int ip = 0; ip < r.size(); ++ip)
    MoverRng[ip] = &r[ip].get();
  for (int ip = 0; ip < r.size(); ++ip)
    RngSaved[ip] = r[ip].get().makeClone();
}

void QMCCostFunctionBase::setTargetEnergy(Return_rt et)
{
  //evaluate effective target energy
  EtargetEff = Etarget = et;
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

QMCCostFunctionBase::Return_rt QMCCostFunctionBase::Cost(bool needGrad)
{
  NumCostCalls++;
  //reset the wave function
  resetPsi();
  //evaluate new local energies
  EffectiveWeight effective_weight = correlatedSampling(needGrad);
  IsValid                          = isEffectiveWeightValid(effective_weight);
  return computedCost();
}

void QMCCostFunctionBase::printEstimates()
{
  app_log() << "      Current ene:     " << curAvg_w << std::endl;
  app_log() << "      Current var:     " << curVar_w << std::endl;
  app_log() << "      Current ene_urw:     " << curAvg << std::endl;
  app_log() << "      Current var_urw: " << curVar << std::endl;
}

QMCCostFunctionBase::Return_rt QMCCostFunctionBase::computedCost()
{
  //   Assumes the Sums have been computed all ready
  //Estimators::accumulate has been called by correlatedSampling
  curAvg_w         = SumValue[SUM_E_WGT] / SumValue[SUM_WGT];
  curVar_w         = SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT] - curAvg_w * curAvg_w;
  Return_rt wgtinv = 1.0 / static_cast<Return_rt>(NumSamples);
  // app_log() << "wgtinv = " << wgtinv << std::endl;
  curAvg     = SumValue[SUM_E_BARE] * wgtinv;
  curVar     = SumValue[SUM_ESQ_BARE] * wgtinv - curAvg * curAvg;
  curVar_abs = SumValue[SUM_ABSE_WGT] / SumValue[SUM_WGT];
  // app_log() << "curVar     = " << curVar
  //     << "   curAvg     = " << curAvg << std::endl;
  // app_log() << "SumValue[SUM_WGT] = " << SumValue[SUM_WGT] << std::endl;
  // app_log() << "SumValue[SUM_WGTSQ] = " << SumValue[SUM_WGTSQ] << std::endl;
  // app_log() << "SumValue[SUM_ABSE_WGT] = " << SumValue[SUM_ABSE_WGT] << std::endl;
  // app_log() << "SumValue[SUM_E_WGT] = " << SumValue[SUM_E_WGT] << std::endl;
  // app_log() << "SumValue[SUM_ESQ_WGT] = " << SumValue[SUM_ESQ_WGT] << std::endl;
  CostValue             = 0.0;
  const Return_rt small = 1.0e-10;
  if (std::abs(w_abs) > small)
    CostValue += w_abs * curVar_abs;
  if (std::abs(w_var) > small)
    CostValue += w_var * curVar_w;
  if (std::abs(w_en) > small)
    CostValue += w_en * curAvg_w;
  if (std::abs(w_w) > small)
    CostValue += w_w * curVar;
  return CostValue;
}


void QMCCostFunctionBase::Report()
{
  //reset the wavefunction for with the new variables
  resetPsi();
  if (!myComm->rank())
  {
    updateXmlNodes();
    std::array<char, 128> newxml;
    int length{0};
    if (Write2OneXml)
      length = std::snprintf(newxml.data(), newxml.size(), "%s.opt.xml", RootName.c_str());
    else
      length = std::snprintf(newxml.data(), newxml.size(), "%s.opt.%d.xml", RootName.c_str(), ReportCounter);
    if (length < 0)
      throw std::runtime_error("Error generating fileroot");
    xmlSaveFormatFile(newxml.data(), m_doc_out, 1);
    if (msg_stream)
    {
      msg_stream->precision(8);
      *msg_stream << " curCost " << std::setw(5) << ReportCounter << std::setw(16) << CostValue << std::setw(16)
                  << curAvg_w << std::setw(16) << curAvg << std::setw(16) << curVar_w << std::setw(16) << curVar
                  << std::setw(16) << curVar_abs << std::endl;
      *msg_stream << " curVars " << std::setw(5) << ReportCounter;
      for (int i = 0; i < OptVariables.size(); i++)
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
  //final reset
  resetPsi(true);
  if (!myComm->rank())
  {
    std::ostringstream vp_filename;
    vp_filename << RootName << ".vp.h5";
    hdf_archive hout;
    OptVariables.writeToHDF(vp_filename.str(), hout);

    UniqueOptObjRefs opt_obj_refs = Psi.extractOptimizableObjectRefs();
    for (auto opt_obj : opt_obj_refs)
      opt_obj.get().writeVariationalParameters(hout);

    std::string newxml = RootName + ".opt.xml";
    *msg_stream << "  <optVariables href=\"" << newxml << "\">" << std::endl;
    OptVariables.print(*msg_stream);
    *msg_stream << "  </optVariables>" << std::endl;
    updateXmlNodes();
    xmlSaveFormatFile(newxml.c_str(), m_doc_out, 1);
  }
}
/** This function stores optimized CI coefficients in HDF5 
  * Other parameters (Jastors, orbitals), can also be stored to H5 from this function
  * This function should be called before reportParameters()
  * Since it is needed to call updateXmlNodes() and xmlSaveFormatFile()
  * While it is possible to call updateXmlNodes() from QMCLinearOptimize.cpp 
  * It is not clean to call xmlSaveFormatFile() from QMCLinearOptimize.cpp 
  *
  * @param OptVariables.size() OptVariables.name(i) OptVariables[i]
  * 
  * OptVariables.size(): contains the total number of Optimized variables (not Only CI coeff)
  * OptVariables.name(i): the tag of the optimized variable. To store we use the name of the variable 
  * OptVariables[i]: The new value of the optimized variable 
*/
void QMCCostFunctionBase::reportParametersH5()
{
  if (!myComm->rank())
  {
    int ci_size = 0;
    std::vector<opt_variables_type::value_type> CIcoeff;
    for (int i = 0; i < OptVariables.size(); i++)
    {
      std::array<char, 128> Coeff;
      if (std::snprintf(Coeff.data(), Coeff.size(), "CIcoeff_%d", ci_size + 1) < 0)
        throw std::runtime_error("Error generating fileroot");
      if (OptVariables.name(i) != Coeff.data())
      {
        if (ci_size > 0)
          break;
        else
          continue;
      }

      CIcoeff.push_back(OptVariables[i]);
      ci_size++;
    }
    if (ci_size > 0)
    {
      CI_Opt = true;
      newh5 = RootName + ".opt.h5";
      *msg_stream << "  <Ci Coeffs saved in opt_coeffs=\"" << newh5 << "\">" << std::endl;
      hdf_archive hout;
      hout.create(newh5, H5F_ACC_TRUNC);
      hout.push("MultiDet", true);
      hout.write(ci_size, "NbDet");
      hout.write(CIcoeff, "Coeff");
      hout.close();
    }
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
bool QMCCostFunctionBase::put(xmlNodePtr q)
{
  std::string includeNonlocalH;
  std::string writeXmlPerStep("no");
  std::string computeNLPPderiv;
  astring variational_subset_str;
  ParameterSet m_param;
  m_param.add(writeXmlPerStep, "dumpXML");
  m_param.add(MinNumWalkers, "minwalkers");
  m_param.add(MaxWeight, "maxWeight");
  m_param.add(includeNonlocalH, "nonlocalpp", {}, TagStatus::DEPRECATED);
  m_param.add(computeNLPPderiv, "use_nonlocalpp_deriv", {}, TagStatus::DEPRECATED);
  m_param.add(w_beta, "beta");
  m_param.add(GEVType, "GEVMethod");
  m_param.add(targetExcitedStr, "targetExcited");
  m_param.add(omega_shift, "omega");
  m_param.add(do_override_output, "output_vp_override", {true});
  m_param.add(variational_subset_str, "variational_subset");
  m_param.put(q);

  if (!includeNonlocalH.empty())
    app_warning() << "'nonlocalpp' no more affects any part of the execution. Please remove it from your input file."
                  << std::endl;
  if (!computeNLPPderiv.empty())
    app_warning()
        << "'use_nonlocalpp_deriv' no more affects any part of the execution. Please remove it from your input file."
        << std::endl;

  targetExcitedStr = lowerCase(targetExcitedStr);
  targetExcited    = (targetExcitedStr == "yes");

  variational_subset_names = convertStrToVec<std::string>(variational_subset_str.s);

  // app_log() << "  QMCCostFunctionBase::put " << std::endl;
  // m_param.get(app_log());
  Write2OneXml     = (writeXmlPerStep == "no");
  xmlNodePtr qsave = q;
  //Estimators.put(q);
  std::vector<xmlNodePtr> cset;
  std::vector<std::string> excluded;
  std::map<std::string, std::vector<std::string>*> equalConstraints;
  std::map<std::string, std::vector<std::string>*> negateConstraints;
  std::vector<std::string> idtag;
  xmlNodePtr cur = qsave->children;
  int pid        = myComm->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "optimize")
    {
      std::vector<std::string> tmpid;
      putContent(tmpid, cur);
      idtag.insert(idtag.end(), tmpid.begin(), tmpid.end());
    }
    else if (cname == "exclude")
    {
      std::vector<std::string> tmpid;
      putContent(tmpid, cur);
      excluded.insert(excluded.end(), tmpid.begin(), tmpid.end());
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
      pAttrib.add(ctype, "type");
      pAttrib.add(s, "name");
      pAttrib.put(cur);
      if (ctype == "equal" || ctype == "=")
      {
        std::map<std::string, std::vector<std::string>*>::iterator eit(equalConstraints.find(s));
        std::vector<std::string>* eqSet = 0;
        if (eit == equalConstraints.end())
        {
          eqSet               = new std::vector<std::string>;
          equalConstraints[s] = eqSet;
        }
        else
          eqSet = (*eit).second;
        std::vector<std::string> econt;
        putContent(econt, cur);
        eqSet->insert(eqSet->end(), econt.begin(), econt.end());
      }
    }
    cur = cur->next;
  }

  UniqueOptObjRefs opt_obj_refs = extractOptimizableObjects(Psi);
  app_log() << " TrialWaveFunction \"" << Psi.getName() << "\" has " << opt_obj_refs.size()
            << " optimizable objects:" << std::endl;
  for (OptimizableObject& obj : opt_obj_refs)
    app_log() << "   '" << obj.getName() << "'" << (obj.isOptimized() ? " optimized" : " fixed") << std::endl;

  //build optimizables from the wavefunction
  OptVariablesForPsi.clear();
  for (OptimizableObject& obj : opt_obj_refs)
    if (obj.isOptimized())
      obj.checkInVariablesExclusive(OptVariablesForPsi);
  OptVariablesForPsi.resetIndex();
  app_log() << " Variational subset selects " << OptVariablesForPsi.size() << " parameters." << std::endl;

  //synchronize OptVariables and OptVariablesForPsi
  OptVariables  = OptVariablesForPsi;
  InitVariables = OptVariablesForPsi;
  //first disable <exclude>.... </exclude> from the big list used by a TrialWaveFunction
  OptVariablesForPsi.disable(excluded.begin(), excluded.end(), false);
  //now, set up the real variable list for optimization
  //check <equal>
  int nc = 0;
  if (equalConstraints.size())
  {
    std::map<std::string, std::vector<std::string>*>::iterator eit(equalConstraints.begin());
    while (eit != equalConstraints.end())
    {
      nc += (*eit).second->size();
      //actiave the active variable even though it is probably unnecessary
      OptVariablesForPsi.activate((*eit).second->begin(), (*eit).second->end(), false);
      excluded.insert(excluded.end(), (*eit).second->begin(), (*eit).second->end());
      ++eit;
    }
  }
  //build OptVariables which is equal to or identical to OptVariablesForPsi
  //disable the variables that are equal to a variable
  OptVariables.disable(excluded.begin(), excluded.end(), false);
  //set up OptVariables and OptVariablesForPsi
  OptVariables.activate(idtag.begin(), idtag.end(), true);
  OptVariablesForPsi.activate(idtag.begin(), idtag.end(), true);
  //found constraints build equalVarMap
  if (nc > 0)
  {
    equalVarMap.reserve(nc + OptVariables.size());
    //map the basic lists from the active list
    for (int i = 0; i < OptVariables.size(); ++i)
    {
      int bigloc = OptVariablesForPsi.getIndex(OptVariables.name(i));
      if (bigloc < 0)
        continue;
      equalVarMap.push_back(TinyVector<int, 2>(bigloc, i));
    }
    //add <equal/>
    std::map<std::string, std::vector<std::string>*>::iterator eit(equalConstraints.begin());
    while (eit != equalConstraints.end())
    {
      int loc = OptVariables.getIndex((*eit).first);
      if (loc >= 0)
      {
        const std::vector<std::string>& elist(*((*eit).second));
        for (int i = 0; i < elist.size(); ++i)
        {
          int bigloc = OptVariablesForPsi.getIndex(elist[i]);
          if (bigloc < 0)
            continue;
          equalVarMap.push_back(TinyVector<int, 2>(bigloc, loc));
        }
      }
      //remove std::vector<std::string>
      delete (*eit).second;
      ++eit;
    }
  }
  //get the indices
  Psi.checkOutVariables(OptVariablesForPsi);
  NumOptimizables = OptVariables.size();
  if (NumOptimizables == 0)
  {
    APP_ABORT("QMCCostFunctionBase::put No valid optimizable variables are found.");
  }
  else
    app_log() << " In total " << NumOptimizables << " parameters being optimized after applying constraints."
              << std::endl;
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
  EtargetEff = Etarget;
  return true;
}

void QMCCostFunctionBase::resetCostFunction(std::vector<xmlNodePtr>& cset)
{
  for (int i = 0; i < cset.size(); i++)
  {
    std::string pname;
    OhmmsAttributeSet pAttrib;
    pAttrib.add(pname, "name");
    pAttrib.put(cset[i]);
    if (pname == "energy")
      putContent(w_en, cset[i]);
    else if ((pname == "variance") || (pname == "unreweightedvariance"))
      putContent(w_w, cset[i]);
    else if (pname == "difference")
      putContent(w_abs, cset[i]);
    else if ((pname == "reweightedVariance") || (pname == "reweightedvariance"))
      putContent(w_var, cset[i]);
  }
}


void QMCCostFunctionBase::updateXmlNodes()
{
  if (m_doc_out == NULL) //first time, create a document tree and get parameters and attributes to be updated
  {
    m_doc_out          = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    xmlNodePtr wf_root = xmlAddChild(qm_root, xmlCopyNode(m_wfPtr, 1));
    xmlDocSetRootElement(m_doc_out, qm_root);
    xmlXPathContextPtr acontext = xmlXPathNewContext(m_doc_out);

    if (do_override_output)
    {
      std::ostringstream vp_filename;
      vp_filename << RootName << ".vp.h5";

      OhmmsXPathObject vp_file_nodes("//override_variational_parameters", acontext);
      if (vp_file_nodes.empty())
      {
        // Element is not present. Create a new one.
        xmlNodePtr vp_file_node = xmlNewNode(NULL, BAD_CAST "override_variational_parameters");
        xmlSetProp(vp_file_node, BAD_CAST "href", BAD_CAST vp_filename.str().c_str());
        xmlAddChild(wf_root, vp_file_node);
      }
      else
      {
        // Element is present. Rewrite the href attribute.
        for (int iparam = 0; iparam < vp_file_nodes.size(); iparam++)
        {
          xmlSetProp(vp_file_nodes[iparam], BAD_CAST "href", BAD_CAST vp_filename.str().c_str());
        }
      }
    }

    //check var
    xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//var", acontext);
    for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
      std::string aname(getXMLAttributeValue(cur, "id"));
      if (aname.empty())
        continue;
      if (auto oit = OptVariablesForPsi.find(aname); oit != OptVariablesForPsi.end())
        paramNodes[aname] = cur;
    }
    xmlXPathFreeObject(result);
    //check radfunc
    result = xmlXPathEvalExpression((const xmlChar*)"//radfunc", acontext);
    for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
      std::string aname(getXMLAttributeValue(cur, "id"));
      if (aname.empty())
        continue;
      if (xmlAttrPtr aptr = xmlHasProp(cur, (const xmlChar*)"exponent"); aptr != nullptr)
      {
        std::string expID = aname + "_E";
        if (auto oit = OptVariablesForPsi.find(expID); oit != OptVariablesForPsi.end())
          attribNodes[expID] = std::pair<xmlNodePtr, std::string>(cur, "exponent");
      }
      std::string cID = aname + "_C";
      if (xmlAttrPtr aptr = xmlHasProp(cur, (const xmlChar*)"contraction"); aptr != nullptr)
        if (auto oit = OptVariablesForPsi.find(cID); oit != OptVariablesForPsi.end())
          attribNodes[cID] = std::pair<xmlNodePtr, std::string>(cur, "contraction");
    }
    xmlXPathFreeObject(result);
    //check ci
    result = xmlXPathEvalExpression((const xmlChar*)"//ci", acontext);
    for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
      std::string aname(getXMLAttributeValue(cur, "id"));
      if (aname.empty())
        continue;
      xmlAttrPtr aptr = xmlHasProp(cur, (const xmlChar*)"coeff");
      opt_variables_type::iterator oit(OptVariablesForPsi.find(aname));
      if (xmlAttrPtr aptr = xmlHasProp(cur, (const xmlChar*)"coeff"); aptr != NULL)
        if (auto oit = OptVariablesForPsi.find(aname); oit != OptVariablesForPsi.end())
          attribNodes[aname] = std::pair<xmlNodePtr, std::string>(cur, "coeff");
    }
    xmlXPathFreeObject(result);
    //check csf
    result = xmlXPathEvalExpression((const xmlChar*)"//csf", acontext);
    for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
    {
      xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
      std::string aname(getXMLAttributeValue(cur, "id"));
      if (aname.empty())
        continue;
      if (xmlAttrPtr aptr = xmlHasProp(cur, (const xmlChar*)"coeff"); aptr != nullptr)
        if (auto oit = OptVariablesForPsi.find(aname); oit != OptVariablesForPsi.end())
          attribNodes[aname] = std::pair<xmlNodePtr, std::string>(cur, "coeff");
    }
    xmlXPathFreeObject(result);
    if (CI_Opt)
    {
      //check multidet
      result = xmlXPathEvalExpression((const xmlChar*)"//detlist", acontext);
      for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
      {
        xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
        xmlSetProp(cur, (const xmlChar*)"opt_coeffs", (const xmlChar*)newh5.c_str());
      }
      xmlXPathFreeObject(result);
    }

    addCoefficients(acontext, "//coefficient");
    addCoefficients(acontext, "//coefficients");
    addCJParams(acontext, "//jastrow");
    xmlXPathFreeContext(acontext);
  }
  for (const auto& [pname, pptr] : paramNodes)
  {
    //FIXME real value is forced here to makde sure that the code builds
    Return_rt v = std::real(OptVariablesForPsi[pname]);
    getContent(v, pptr);
  }
  for (const auto& [aname, attrib] : attribNodes)
  {
    std::ostringstream vout;
    vout.setf(std::ios::scientific, std::ios::floatfield);
    vout.precision(16);
    vout << OptVariablesForPsi[aname];
    xmlSetProp(attrib.first, (const xmlChar*)attrib.second.c_str(), (const xmlChar*)vout.str().c_str());
  }
  for (const auto& [cname, cptr] : coeffNodes)
  {
    std::string rname(cname);
    OhmmsAttributeSet cAttrib;
    std::string datatype("none");
    std::string aname("0");
    cAttrib.add(datatype, "type");
    cAttrib.add(aname, "id");
    cAttrib.put(cptr);
    if (datatype == "Array")
    {
      //
      aname.append("_");
      std::vector<Return_rt> c;
      for (const auto& [name, value] : OptVariablesForPsi)
      {
        if (name.find(aname) == 0)
        {
          //FIXME real value is forced here to makde sure that the code builds
          c.push_back(std::real(value));
        }
      }
      xmlNodePtr contentPtr = cptr;
      if (xmlNodeIsText(contentPtr->children))
        contentPtr = contentPtr->children;
      getContent(c, contentPtr);
    }
    // counting jastrow variables
    else if (rname.find("cj_") == 0)
    {
      printCJParams(cptr, rname);
    }
    else
    {
      xmlNodePtr cur = cptr->children;
      while (cur != NULL)
      {
        std::string childName((const char*)(cur->name));
        if (childName == "lambda")
        {
          int i = 0;
          int j = -1;
          OhmmsAttributeSet pAttrib;
          pAttrib.add(i, "i");
          pAttrib.add(j, "j");
          pAttrib.put(cur);
          std::array<char, 32> lambda_id;
          int length{0};
          if (j < 0)
            length = std::snprintf(lambda_id.data(), lambda_id.size(), "%s_%d", rname.c_str(), i);
          else
            length = std::snprintf(lambda_id.data(), lambda_id.size(), "%s_%d_%d", rname.c_str(), i, j);
          if (length < 0)
            throw std::runtime_error("Error generating lambda_id");
          opt_variables_type::iterator vTarget(OptVariablesForPsi.find(lambda_id.data()));
          if (vTarget != OptVariablesForPsi.end())
          {
            std::ostringstream vout;
            vout.setf(std::ios::scientific, std::ios::floatfield);
            vout.precision(16);
            vout << (*vTarget).second;
            xmlSetProp(cur, (const xmlChar*)"c", (const xmlChar*)vout.str().c_str());
          }
        }
        cur = cur->next;
      }
    }
  }
}

/** add coefficient or coefficients
 * @param acontext context from which xpath cname is searched
 * @param cname xpath
 */
void QMCCostFunctionBase::addCoefficients(xmlXPathContextPtr acontext, const char* cname)
{
  xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)cname, acontext);
  for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
  {
    xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
    OhmmsAttributeSet cAttrib;
    std::string aname("0");
    std::string optimize("yes");
    std::string datatype("none");
    cAttrib.add(aname, "id");
    cAttrib.add(aname, "name");
    cAttrib.add(datatype, "type");
    cAttrib.add(optimize, "optimize");
    cAttrib.put(cur);
    if (aname[0] == '0')
      continue;
    if (datatype == "Array")
    {
      if (optimize == "yes")
        coeffNodes[aname] = cur;
    }
    else
    {
      //check if any optimizables contains the id of coefficients
      bool notlisted = true;
      opt_variables_type::iterator oit(OptVariablesForPsi.begin()), oit_end(OptVariablesForPsi.end());
      while (notlisted && oit != oit_end)
      {
        const std::string& oname((*oit).first);
        notlisted = (oname.find(aname) >= oname.size());
        ++oit;
      }
      if (!notlisted)
      {
        coeffNodes[aname] = cur;
      }
    }
  }
  xmlXPathFreeObject(result);
}

void QMCCostFunctionBase::addCJParams(xmlXPathContextPtr acontext, const char* cname)
{
  xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)cname, acontext);
  for (int iparam = 0; iparam < result->nodesetval->nodeNr; iparam++)
  {
    // check that we're in a counting jastrow tag space
    xmlNodePtr cur = result->nodesetval->nodeTab[iparam];
    OhmmsAttributeSet cAttrib;
    std::string aname("none");
    std::string atype("none");
    cAttrib.add(aname, "name");
    cAttrib.add(atype, "type");
    cAttrib.put(cur);
    if (atype == "Counting")
    {
      // iterate through children
      xmlNodePtr cur2 = cur->xmlChildrenNode;
      int Fnum        = 0;
      bool Ftag_found = false;
      // find the <var name="F" /> tag, save the location
      while (cur2 != NULL)
      {
        std::string tname2((const char*)cur2->name);
        OhmmsAttributeSet cAttrib2;
        std::string aname2("none");
        std::string aopt2("yes");
        std::string ref_id("g0");
        cAttrib2.add(aname2, "name");
        cAttrib2.add(aopt2, "opt");
        cAttrib2.add(ref_id, "reference_id");
        cAttrib2.put(cur2);
        // find the <var name="F" /> tags and register the tag location first
        if (tname2 == "var" && aname2 == "F" && (aopt2 == "yes" || aopt2 == "true"))
        {
          coeffNodes["cj_F"] = cur2;
          Ftag_found         = true;
        }
        cur2 = cur2->next;
      }

      // count the total number of registered F matrix variables
      opt_variables_type::iterator oit(OptVariables.begin()), oit_end(OptVariables.end());
      for (; oit != oit_end; ++oit)
      {
        const std::string& oname((*oit).first);
        if (oname.find("F_") == 0)
          ++Fnum;
      }
      ++Fnum;

      // F variable tag isn't found; build the tag.
      if (!Ftag_found)
      {
        xmlNodeAddContent(cur, (const xmlChar*)"\n      ");
        xmlNodePtr F_tag = xmlNewChild(cur, NULL, (const xmlChar*)"var", NULL);
        xmlSetProp(F_tag, (const xmlChar*)"name", (const xmlChar*)"F");
        xmlSetProp(F_tag, (const xmlChar*)"opt", (const xmlChar*)"true");
        xmlNodeSetContent(F_tag, (const xmlChar*)" ");
        coeffNodes["cj_F"] = F_tag;
        xmlNodeAddContent(cur, (const xmlChar*)"\n    ");
      }

      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL)
      {
        std::string tname2((const char*)cur2->name);
        OhmmsAttributeSet cAttrib2;
        std::string aname2("none");
        std::string atype2("normalized_gaussian");
        std::string aopt2("yes");
        std::string ref_id("g0");
        cAttrib2.add(aname2, "name");
        cAttrib2.add(atype2, "type");
        cAttrib2.add(aopt2, "opt");
        cAttrib2.add(ref_id, "reference_id");
        cAttrib2.put(cur2);
        // find <region /> tag
        if (tname2 == "region" && (aopt2 == "true" || aopt2 == "yes") && Fnum != 0)
        {
          // if type is voronoi, reconstruct the tag structure
          if (atype2 == "voronoi")
          {
            // set property to normalized_gaussian, add reference, remove source
            xmlUnsetProp(cur, (const xmlChar*)"source");
            xmlSetProp(cur2, (const xmlChar*)"type", (const xmlChar*)"normalized_gaussian");
            xmlSetProp(cur2, (const xmlChar*)"reference_id", (const xmlChar*)ref_id.c_str());
            // remove existing children
            xmlNodePtr child = cur2->xmlChildrenNode;
            while (child != NULL)
            {
              xmlUnlinkNode(child);
              child = cur2->xmlChildrenNode;
            }
            // get Fdim
            int Fdim = (std::sqrt(1 + 8 * Fnum) - 1) / 2;
            for (int i = 0; i < Fdim; ++i)
            {
              std::ostringstream os;
              os << "g" << i;
              std::string gid = os.str();
              std::string opt = (gid != ref_id) ? "true" : "false";

              // create function tag, set id
              xmlNodeAddContent(cur2, (const xmlChar*)"\n        ");
              xmlNodePtr function_tag = xmlNewChild(cur2, NULL, (const xmlChar*)"function", NULL);
              xmlNewProp(function_tag, (const xmlChar*)"id", (const xmlChar*)gid.c_str());
              xmlNodeAddContent(function_tag, (const xmlChar*)"\n          ");
              // create A tag
              xmlNodePtr A_tag = xmlNewChild(function_tag, NULL, (const xmlChar*)"var", NULL);
              xmlNodeAddContent(function_tag, (const xmlChar*)"\n          ");
              xmlNewProp(A_tag, (const xmlChar*)"name", (const xmlChar*)"A");
              xmlNewProp(A_tag, (const xmlChar*)"opt", (const xmlChar*)opt.c_str());
              os.str("");
              if (gid == ref_id)
              {
                os << std::setprecision(6) << std::scientific;
                os << InitVariables.find(gid + "_A_xx")->second << " " << InitVariables.find(gid + "_A_xy")->second
                   << " " << InitVariables.find(gid + "_A_xz")->second << " "
                   << InitVariables.find(gid + "_A_yy")->second << " " << InitVariables.find(gid + "_A_yz")->second
                   << " " << InitVariables.find(gid + "_A_zz")->second;
              }
              else
                os << " ";
              xmlNodeSetContent(A_tag, (const xmlChar*)(os.str().c_str()));

              // create B tag
              xmlNodePtr B_tag = xmlNewChild(function_tag, NULL, (const xmlChar*)"var", NULL);
              xmlNodeAddContent(function_tag, (const xmlChar*)"\n          ");
              xmlNewProp(B_tag, (const xmlChar*)"name", (const xmlChar*)"B");
              xmlNewProp(B_tag, (const xmlChar*)"opt", (const xmlChar*)opt.c_str());
              os.str("");
              if (gid == ref_id)
              {
                os << std::setprecision(6) << std::scientific;
                os << InitVariables.find(gid + "_B_x")->second << " " << InitVariables.find(gid + "_B_y")->second << " "
                   << InitVariables.find(gid + "_B_z")->second;
              }
              else
                os << " ";
              xmlNodeSetContent(B_tag, (const xmlChar*)(os.str().c_str()));

              // create C tag
              xmlNodePtr C_tag = xmlNewChild(function_tag, NULL, (const xmlChar*)"var", NULL);
              xmlNodeAddContent(function_tag, (const xmlChar*)"\n        ");
              xmlNewProp(C_tag, (const xmlChar*)"name", (const xmlChar*)"C");
              xmlNewProp(C_tag, (const xmlChar*)"opt", (const xmlChar*)opt.c_str());
              os.str("");
              if (gid == ref_id)
              {
                os << std::setprecision(6) << std::scientific;
                os << InitVariables.find(gid + "_C")->second;
              }
              else
                os << " ";
              xmlNodeSetContent(C_tag, (const xmlChar*)(os.str().c_str()));
            }
            xmlNodeAddContent(cur2, (const xmlChar*)"\n      ");
          }

          // register the optimizable variables for each function
          // find <function /> tags
          xmlNodePtr cur3 = cur2->xmlChildrenNode;
          while (cur3 != NULL)
          {
            std::string tname3((const char*)cur3->name);
            OhmmsAttributeSet cAttrib3;
            std::string fid("none");
            cAttrib3.add(fid, "id");
            cAttrib3.put(cur3);
            // for each function tag, register coeffNodes[gid_A/B/C] as each variable location
            if (tname3 == "function")
            {
              // find <var name="A/B/C" /> tags
              xmlNodePtr cur4 = cur3->xmlChildrenNode;
              while (cur4 != NULL)
              {
                std::string tname4((const char*)cur4->name);
                OhmmsAttributeSet cAttrib4;
                std::string aname4("none");
                std::string aopt4("false");
                cAttrib4.add(aname4, "name");
                cAttrib4.add(aopt4, "opt");
                cAttrib4.put(cur4);
                if (tname4 == "var" && (aopt4 == "true" || aopt4 == "yes"))
                {
                  std::string varname = "cj_" + fid + "_" + aname4;
                  bool notlisted      = true;
                  opt_variables_type::iterator oit(OptVariables.begin()), oit_end(OptVariables.end());
                  while (notlisted && oit != oit_end)
                  {
                    const std::string& oname((*oit).first);
                    notlisted = (oname.find(varname.substr(3)) >= oname.size());
                    ++oit;
                  }
                  if (!notlisted)
                  {
                    coeffNodes[varname] = cur4;
                  }
                }
                cur4 = cur4->next;
              }
            }
            cur3 = cur3->next;
          }
        }
        cur2 = cur2->next;
      }
    }
  }
  xmlXPathFreeObject(result);
}


void QMCCostFunctionBase::printCJParams(xmlNodePtr cur, std::string& rname)
{
  opt_variables_type::iterator vit(OptVariables.begin());
  // F matrix variables
  if (rname.find("cj_F") < rname.size())
  {
    // get a vector of pairs with f matrix names, values
    std::vector<Return_rt> f_vals;
    for (auto vit = OptVariables.begin(); vit != OptVariables.end(); ++vit)
    {
      if ((*vit).first.find("F_") == 0)
      {
        Return_rt fval = std::real(vit->second);
        f_vals.push_back(fval);
      }
    }
    // manually add final element = 0
    f_vals.push_back(0);
    // get the dimensions of the f matrix
    int Fdim = (std::sqrt(1 + 8 * f_vals.size()) - 1) / 2;
    std::ostringstream os;
    std::string pad_str = std::string(14, ' ');
    os << std::endl;
    os << std::setprecision(6) << std::scientific;
    // print out the values with the proper indentation
    auto fit = f_vals.begin();
    for (int i = 0; i < Fdim; ++i)
    {
      os << pad_str;
      for (int j = 0; j < Fdim; ++j)
      {
        if (j < i)
          os << pad_str; //print blank
        else             // print value
        {
          os << "  " << (*fit);
          ++fit;
        }
      }
      // line break
      if (i < Fdim - 1)
        os << std::endl;
    }
    // assign to tag
    xmlNodePtr cur2 = cur->children;
    xmlNodeSetContent(cur2, (const xmlChar*)(os.str().c_str()));
    xmlNodeAddContent(cur2, (const xmlChar*)"\n      ");
  }
  // gaussian function parameters A, B, C
  else
  {
    std::string var_prefix = rname.substr(3);
    std::vector<Return_rt> vals;

    for (auto vit = OptVariablesForPsi.begin(); vit != OptVariablesForPsi.end(); ++vit)
    {
      if (vit->first.find(var_prefix) == 0)
      {
        Return_rt val = std::real(vit->second);
        vals.push_back(val);
      }
    }
    std::ostringstream os;
    os << std::setprecision(6) << std::scientific;
    for (int i = 0; i < vals.size(); ++i)
    {
      os << vals[i];
      if (i < vals.size() - 1)
        os << " ";
    }
    xmlNodePtr cur2 = cur->children;
    xmlNodeSetContent(cur2, (const xmlChar*)(os.str().c_str()));
  }
}

bool QMCCostFunctionBase::isEffectiveWeightValid(EffectiveWeight effective_weight) const
{
  app_log() << "Effective weight of all the samples measured by correlated sampling is " << effective_weight
            << std::endl;
  if (effective_weight < MinNumWalkers)
  {
    WARNMSG("    Smaller than the user specified threshold \"minwalkers\" = "
            << MinNumWalkers << std::endl
            << "  If this message appears frequently. You might have to be cautious. " << std::endl
            << "  Find info about parameter \"minwalkers\" in the user manual!");
    return false;
  }

  return true;
}

UniqueOptObjRefs QMCCostFunctionBase::extractOptimizableObjects(TrialWaveFunction& psi) const
{
  const auto& names(variational_subset_names);
  // survey all the optimizable objects
  const auto opt_obj_refs = psi.extractOptimizableObjectRefs();
  // check if input names are valid
  for (auto& name : names)
    if (std::find_if(opt_obj_refs.begin(), opt_obj_refs.end(),
                     [&name](const OptimizableObject& obj) { return name == obj.getName(); }) == opt_obj_refs.end())
    {
      std::ostringstream msg;
      msg << "Variational subset entry '" << name << "' doesn't exist in the trial wavefunction which contains";
      for (OptimizableObject& obj : opt_obj_refs)
        msg << " '" << obj.getName() << "'";
      msg << "." << std::endl;
      throw UniformCommunicateError(msg.str());
    }

  for (OptimizableObject& obj : opt_obj_refs)
    obj.setOptimization(names.empty() || std::find_if(names.begin(), names.end(), [&obj](const std::string& name) {
                                           return name == obj.getName();
                                         }) != names.end());
  return opt_obj_refs;
}

void QMCCostFunctionBase::resetOptimizableObjects(TrialWaveFunction& psi, const opt_variables_type& opt_variables) const
{
  const auto opt_obj_refs = extractOptimizableObjects(psi);
  for (OptimizableObject& obj : opt_obj_refs)
    if (obj.isOptimized())
      obj.resetParametersExclusive(opt_variables);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  If the LMYEngine is available, returns the cost function calculated by the engine.
///         Otherwise, returns the usual cost function.
///
/// \param[in]      needDeriv             whether derivative vectors should be computed
///
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LMY_ENGINE
QMCCostFunctionBase::Return_rt QMCCostFunctionBase::LMYEngineCost(const bool needDeriv,
                                                                  cqmc::engine::LMYEngine<Return_t>* EngineObj)
{
  // prepare local energies, weights, and possibly derivative vectors, and compute standard cost
  const Return_rt standardCost = this->Cost(needDeriv);

  // since we are using the LMYEngine, compute and return it's cost function value
  return this->LMYEngineCost_detail(EngineObj);
}
#endif

} // namespace qmcplusplus

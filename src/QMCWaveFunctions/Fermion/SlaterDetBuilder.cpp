//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/SPOSetScanner.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"

#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminant.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#if defined(QMC_CUDA)
#include "QMCWaveFunctions/Fermion/DiracDeterminantCUDA.h"
#endif
#include "QMCWaveFunctions/Fermion/BackflowBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDetWithBackflow.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include <vector>
//#include "QMCWaveFunctions/Fermion/ci_node.h"
#include "QMCWaveFunctions/Fermion/ci_configuration.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "QMCWaveFunctions/Fermion/SPOSetProxy.h"
#include "QMCWaveFunctions/Fermion/SPOSetProxyForMSD.h"

#include <bitset>
#include <unordered_map>

namespace qmcplusplus
{
SlaterDetBuilder::SlaterDetBuilder(Communicate* comm, ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets)
    : WaveFunctionComponentBuilder(comm, els),
      targetPsi(psi),
      ptclPool(psets),
      slaterdet_0(nullptr),
      multislaterdet_0(nullptr),
      multislaterdetfast_0(nullptr)
{
  ClassName   = "SlaterDetBuilder";
  BFTrans     = 0;
  UseBackflow = false;
}

/** process <determinantset>
 *
 * determinantset
 * - basiset 0..1
 * - sposet 0..*
 * - slaterdeterminant
 *   - determinant 0..*
 * - ci
 */
WaveFunctionComponent* SlaterDetBuilder::buildComponent(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "put(xmlNodePtr)");
  ///save the current node
  xmlNodePtr curRoot = cur;
  xmlNodePtr BFnode = nullptr;
  bool success = true, FastMSD = true;
  std::string cname, tname;
  std::map<std::string, SPOSetPtr> spomap;
  bool multiDet = false;

  if (!mySPOSetBuilderFactory)
  { //always create one, using singleton and just to access the member functions
    mySPOSetBuilderFactory = std::make_unique<SPOSetBuilderFactory>(myComm, targetPtcl, ptclPool);
    mySPOSetBuilderFactory->createSPOSetBuilder(curRoot);
  }

  //check the basis set
  cur = curRoot->children;
  while (cur != NULL) //check the basis set
  {
    getNodeName(cname, cur);
    if (cname == basisset_tag)
    {
      mySPOSetBuilderFactory->loadBasisSetFromXML(cur);
    }
    else if (cname == sposet_tag)
    {
      app_log() << "Creating SPOSet in SlaterDetBuilder::put(xmlNodePtr cur).\n";
      std::string spo_name;
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add(spo_name, "name");
      spoAttrib.put(cur);
      app_log() << "spo_name = " << spo_name << std::endl;
      SPOSetPtr spo = mySPOSetBuilderFactory->createSPOSet(cur);
      //spo->put(cur, spomap);
      if (spomap.find(spo_name) != spomap.end())
      {
        app_error() << "SPOSet name \"" << spo_name << "\" is already in use.\n";
        abort();
      }
      spomap[spo_name] = spo;
      spo->objectName  = spo_name;
      assert(spomap.find(spo_name) != spomap.end());
      //	slaterdet_0->add(spo,spo_name);
    }
    else if (cname == backflow_tag)
    {
      app_log() << "Creating Backflow transformation in SlaterDetBuilder::put(xmlNodePtr cur).\n";
      // to simplify the logic inside DiracDeterminantWithBackflow,
      // I'm requiring that only a single <backflow> block appears
      // in the xml file
      if (BFTrans != 0)
      {
        APP_ABORT("Only a single backflow block is allowed in the xml. Please collect all transformations into a "
                  "single block. \n");
      }
      UseBackflow = true;
      // creating later due to problems with ParticleSets
      //BFTrans = new BackflowTransformation(targetPtcl,ptclPool);
      BFTrans = nullptr;
      BFnode  = cur;
      // read xml later, in case some ParticleSets are read from hdf5 file.
      //BFTrans->put(cur);
    }
    cur = cur->next;
  }

  //sposet_builder is defined outside <determinantset/>
  if (spomap.empty())
  {
    cur = curRoot->children;
    while (cur != NULL) //check the basis set
    {
      getNodeName(cname, cur);
      if (cname == sd_tag)
      {
        xmlNodePtr cur1 = cur->children;
        while (cur1 != NULL)
        {
          getNodeName(tname, cur1);
          if (tname == det_tag)
          {
            std::string aspo, did;
            OhmmsAttributeSet a;
            a.add(did, "id");
            a.add(aspo, "sposet");
            a.put(cur1);
            if (aspo.empty())
              aspo = did;
            SPOSet* aset = get_sposet(aspo);
            if (aset)
              spomap[aspo] = aset;
            else
            {
              mySPOSetBuilderFactory->createSPOSetBuilder(curRoot);
              aset = mySPOSetBuilderFactory->createSPOSet(cur1);
              if (aset)
                spomap[aspo] = aset;
            }
          }
          cur1 = cur1->next;
        }
      }
      else if (cname == multisd_tag)
      {
        std::string spo_alpha;
        std::string spo_beta;
        OhmmsAttributeSet a;
        a.add(spo_alpha, "spo_up");
        a.add(spo_beta, "spo_dn");
        a.put(cur);
        SPOSet* alpha = get_sposet(spo_alpha);
        SPOSet* beta  = get_sposet(spo_beta);
        if (alpha && beta)
        {
          spomap[spo_alpha] = alpha;
          spomap[spo_beta]  = beta;
        }
      }
      cur = cur->next;
    }
  }

  if (spomap.empty())
  {
    APP_ABORT_TRACE(__FILE__, __LINE__, " No sposet is found to build slaterdeterminant or multideterminant");
  }

  cur = curRoot->children;
  while (cur != NULL) //check the basis set
  {
    getNodeName(cname, cur);
    if (cname == sposcanner_tag)
    {
      if (myComm->rank() == 0)
      {
        SPOSetScanner ascanner(spomap, targetPtcl, ptclPool);
        ascanner.put(cur);
      }
    }
    else if (cname == sd_tag)
    {
      multiDet = false;
      if (slaterdet_0)
      {
        APP_ABORT("slaterdet is already instantiated.");
      }
      if (UseBackflow)
        slaterdet_0 = new SlaterDetWithBackflow(targetPtcl, BFTrans);
      else
        slaterdet_0 = new SlaterDeterminant_t(targetPtcl);

      // Copy any entries in sposetmap into slaterdet_0
      std::map<std::string, SPOSetPtr>::iterator iter;
      for (iter = spomap.begin(); iter != spomap.end(); iter++)
        slaterdet_0->add(iter->second, iter->first);
      size_t spin_group = 0;
      xmlNodePtr tcur   = cur->children;
      while (tcur != NULL)
      {
        getNodeName(tname, tcur);
        if (tname == det_tag || tname == rn_tag)
        {
          if (putDeterminant(tcur, spin_group))
            spin_group++;
        }
        tcur = tcur->next;
      }
    }
    else if (cname == multisd_tag)
    {
      multiDet = true;
      if (slaterdet_0)
      {
        APP_ABORT("can't combine slaterdet with multideterminant.");
      }
      if (multislaterdet_0 || multislaterdetfast_0)
      {
        APP_ABORT("multideterminant is already instantiated.");
      }
      std::string spo_alpha;
      std::string spo_beta;
      std::string fastAlg("yes");
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add(spo_alpha, "spo_up");
      spoAttrib.add(spo_beta, "spo_dn");
      spoAttrib.add(fastAlg, "Fast");
      spoAttrib.put(cur);
      if (spo_alpha == spo_beta)
      {
        app_error()
            << "In SlaterDetBuilder: In MultiSlaterDeterminant construction, SPO sets must be different. spo_up: "
            << spo_alpha << "  spo_dn: " << spo_beta << "\n";
        abort();
      }
      if (spomap.find(spo_alpha) == spomap.end())
      {
        app_error() << "In SlaterDetBuilder: SPOSet \"" << spo_alpha
                    << "\" is not found. Expected for MultiSlaterDeterminant.\n";
        abort();
      }
      if (spomap.find(spo_beta) == spomap.end())
      {
        app_error() << "In SlaterDetBuilder: SPOSet \"" << spo_beta
                    << "\" is not found. Expected for MultiSlaterDeterminant.\n";
        abort();
      }
      FastMSD = (fastAlg == "yes");
      if (FastMSD)
      {
        if (UseBackflow)
        {
          APP_ABORT("Backflow is not implemented with multi determinants.");
        }
        app_log() << "Using Bryan's algorithm for MultiSlaterDeterminant expansion. \n";
        MultiDiracDeterminant* up_det = 0;
        MultiDiracDeterminant* dn_det = 0;
        app_log() << "Creating base determinant (up) for MSD expansion. \n";
        up_det = new MultiDiracDeterminant((SPOSetPtr)spomap.find(spo_alpha)->second, 0);
        app_log() << "Creating base determinant (down) for MSD expansion. \n";
        dn_det = new MultiDiracDeterminant((SPOSetPtr)spomap.find(spo_beta)->second, 1);

        multislaterdetfast_0 = new MultiSlaterDeterminantFast(targetPtcl, up_det, dn_det);
        success              = createMSDFast(multislaterdetfast_0, cur);
      }
      else
      {
        SPOSetProxyForMSD* spo_up;
        SPOSetProxyForMSD* spo_dn;
        spo_up = new SPOSetProxyForMSD(spomap.find(spo_alpha)->second, targetPtcl.first(0), targetPtcl.last(0));
        spo_dn = new SPOSetProxyForMSD(spomap.find(spo_beta)->second, targetPtcl.first(1), targetPtcl.last(1));
        if (UseBackflow)
        {
          app_log() << "Multi-Slater Determinant expansion with Backflow. \n";
          app_log() << "Using a list of dirac determinants for MultiSlaterDeterminant expansion. \n";
          multislaterdet_0 = new MultiSlaterDeterminantWithBackflow(targetPtcl, spo_up, spo_dn, BFTrans);
          success          = createMSD(multislaterdet_0, cur);
        }
        else
        {
          app_log() << "Using a list of dirac determinants for MultiSlaterDeterminant expansion. \n";
          multislaterdet_0 = new MultiSlaterDeterminant(targetPtcl, spo_up, spo_dn);
          success          = createMSD(multislaterdet_0, cur);
        }
      }
    }
    cur = cur->next;
  }


  if (!multiDet && !slaterdet_0)
  {
    //fatal
    PRE.error("Failed to create a SlaterDeterminant.", true);
    return nullptr;
  }
  if (multiDet && (!multislaterdet_0 && !multislaterdetfast_0))
  {
    //fatal
    PRE.error("Failed to create a MultiSlaterDeterminant.", true);
    return nullptr;
  }
  // change DistanceTables if using backflow
  if (UseBackflow)
  {
    BackflowBuilder bfbuilder(targetPtcl, ptclPool);
    BFTrans = bfbuilder.buildBackflowTransformation(BFnode);
    if (multiDet)
    {
      if (FastMSD)
      {
        multislaterdetfast_0->setBF(BFTrans);
        multislaterdetfast_0->resetTargetParticleSet(BFTrans->QP);
        //           if(BFTrans->isOptimizable()) multislaterdetfast_0->Optimizable = true;
      }
      else
      {
        multislaterdet_0->setBF(BFTrans);
        multislaterdet_0->resetTargetParticleSet(BFTrans->QP);
        //           if(BFTrans->isOptimizable()) multislaterdet_0->Optimizable = true;
      }
    }
    else
    {
      slaterdet_0->setBF(BFTrans);
      slaterdet_0->resetTargetParticleSet(targetPtcl);
      if (BFTrans->isOptimizable())
        slaterdet_0->Optimizable = true;
    }
  }
  //only single slater determinant
  if (multiDet)
  {
    if (FastMSD)
      return multislaterdetfast_0;
    else
      return multislaterdet_0;
  }
  else
    return slaterdet_0;
}


/** process determiment element
 *
 * determinant has
 * - id unique name
 * - sposet reference to the pre-defined sposet; when missing, use id
 * - group electron species name, u or d
magnetic system
 * Extra attributes to handled the original released-node case
 */
bool SlaterDetBuilder::putDeterminant(xmlNodePtr cur, int spin_group)
{
  ReportEngine PRE(ClassName, "putDeterminant(xmlNodePtr,int)");

  SpeciesSet& myspecies = targetPtcl.mySpecies;

  std::string spin_name = myspecies.speciesName[spin_group];
  std::string sposet;
  std::string basisName("invalid");
  std::string detname("0"), refname("0");
  std::string s_detSize("0");

  OhmmsAttributeSet aAttrib;
  aAttrib.add(basisName, basisset_tag);
  aAttrib.add(detname, "id");
  aAttrib.add(sposet, "sposet");
  aAttrib.add(refname, "ref");
  aAttrib.add(s_detSize, "DetSize");

  std::string s_cutoff("0.0");
  std::string s_radius("0.0");
  int s_smallnumber(-999999);
  int rntype(0);
  aAttrib.add(s_cutoff, "Cutoff");
  aAttrib.add(s_radius, "Radius");
  aAttrib.add(s_smallnumber, "smallnumber");
  aAttrib.add(s_smallnumber, "eps");
  aAttrib.add(rntype, "primary");
  aAttrib.add(spin_name, "group");
  aAttrib.put(cur);

  // whether to use an optimizable slater determinant
  std::string optimize("no");
#if defined(ENABLE_CUDA)
  std::string useGPU("yes");
#else
  std::string useGPU("no");
#endif
  int delay_rank(0);
  OhmmsAttributeSet sdAttrib;
  sdAttrib.add(delay_rank, "delay_rank");
  sdAttrib.add(optimize, "optimize");
  sdAttrib.add(useGPU, "gpu");
  sdAttrib.put(cur->parent);

  { //check determinant@group
    int spin_group_in = spin_group;
    if (isdigit(spin_name[0]))
      spin_group_in = atoi(spin_name.c_str());
    else
      spin_group_in = myspecies.findSpecies(spin_name);
    if (spin_group_in < myspecies.size() && spin_group_in != spin_group)
    {
      spin_group = spin_group_in;
      app_log() << "  Overwrite group = " << spin_group << std::endl;
    }
  }

  //old input does not have sposet
  if (sposet.empty())
    sposet = detname;

  app_log() << "  Creating a determinant " << detname << " group=" << spin_group << " sposet=" << sposet << std::endl;

  std::map<std::string, SPOSetPtr>& spo_ref(slaterdet_0->mySPOSet);
  std::map<std::string, SPOSetPtr>::iterator lit(spo_ref.find(sposet));
  SPOSetPtr psi = 0;
  if (lit == spo_ref.end())
  {
    psi = get_sposet(sposet); //check if the named sposet exists
    if (psi == 0)
    {
      //SPOSet[detname]=psi;
      app_log() << "  Create a new SPO set " << sposet << std::endl;
      psi = mySPOSetBuilderFactory->createSPOSet(cur);
    }
    //psi->put(cur);
    psi->checkObject();
    slaterdet_0->add(psi, detname);
  }
  else
  {
    app_log() << "  Reusing a SPO set " << sposet << std::endl;
    psi = (*lit).second;
  }

  int firstIndex = targetPtcl.first(spin_group);
  int lastIndex  = targetPtcl.last(spin_group);
  if (firstIndex == lastIndex)
    return true;
  std::string dname;
  getNodeName(dname, cur);
  DiracDeterminantBase* adet = 0;
  {
#if defined(QMC_CUDA)
    adet = new DiracDeterminantCUDA(psi, firstIndex);
#else
    if (UseBackflow)
      adet = new DiracDeterminantWithBackflow(targetPtcl, psi, BFTrans, firstIndex);
#if defined(ENABLE_CUDA)
    else if (useGPU == "yes")
    {
      app_log() << "  Using DiracDeterminant with DelayedUpdateCUDA engine" << std::endl;
      adet = new DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>(psi,firstIndex);
    }
#endif
    else
    {
      app_log() << "  Using DiracDeterminant with DelayedUpdate engine" << std::endl;
      adet = new DiracDeterminant<>(psi, firstIndex);
    }
#endif
  }

  if (delay_rank < 0 || delay_rank > lastIndex - firstIndex)
  {
    std::ostringstream err_msg;
    err_msg << "SlaterDetBuilder::putDeterminant delay_rank must be positive "
            << "and no larger than the electron count within a determinant!\n"
            << "Acceptable value [1," << lastIndex - firstIndex << "], "
            << "user input " + std::to_string(delay_rank);
    APP_ABORT(err_msg.str());
  }
  else if (delay_rank == 0)
  {
    app_log() << "  Setting delay_rank by default!" << std::endl;
    if (lastIndex - firstIndex >= 192)
      delay_rank = 32;
    else
      delay_rank = 1;
  }

  if (delay_rank > 1)
    app_log() << "  Using rank-" << delay_rank << " delayed update" << std::endl;
  else
    app_log() << "  Using rank-1 Sherman-Morrison Fahy update" << std::endl;
  adet->set(firstIndex, lastIndex - firstIndex, delay_rank);
#ifdef QMC_CUDA
  targetPsi.setndelay(delay_rank);
#endif
  slaterdet_0->add(adet, spin_group);
  if (psi->isOptimizable())
    slaterdet_0->Optimizable = true;

  app_log() << std::endl;
  app_log().flush();
  return true;
}

bool SlaterDetBuilder::createMSDFast(MultiSlaterDeterminantFast* multiSD, xmlNodePtr cur)
{
  bool success = true;
  std::vector<ci_configuration> uniqueConfg_up, uniqueConfg_dn;
  std::vector<std::string> CItags;
  std::string HDF5Path(""), cname;

  bool optimizeCI;
  int nels_up = multiSD->nels_up;
  int nels_dn = multiSD->nels_dn;
  multiSD->initialize();
  //Check id multideterminants are in HDF5

  xmlNodePtr curTemp = cur, DetListNode = nullptr;
  curTemp            = curTemp->children;
  OhmmsAttributeSet ciAttrib;
  while (curTemp != NULL) //check the basis set
  {
    getNodeName(cname, curTemp);
    if (cname == "detlist")
    {
      DetListNode = curTemp;
    }
    curTemp = curTemp->next;
  }
  ciAttrib.add(HDF5Path, "href");
  ciAttrib.put(DetListNode);
  if (HDF5Path != "")
  {
    app_log() << "Found Multideterminants in H5 File" << std::endl;
    #ifndef QMC_COMPLEX
    success = readDetListH5(cur,
                            uniqueConfg_up,
                            uniqueConfg_dn,
                            *(multiSD->C2node_up),
                            *(multiSD->C2node_dn),
                            CItags,
                            *(multiSD->C),
                            optimizeCI,
                            nels_up,
                            nels_dn);
     #else
      APP_ABORT("Build Multideterminant in H5 File is not implemented for complex CI coefficients");
     #endif
  }
  else
    success = readDetList(cur,
                          uniqueConfg_up,
                          uniqueConfg_dn,
                          *(multiSD->C2node_up),
                          *(multiSD->C2node_dn),
                          CItags,
                          *(multiSD->C),
                          optimizeCI,
                          nels_up,
                          nels_dn,
                          *(multiSD->CSFcoeff),
                          *(multiSD->DetsPerCSF),
                          *(multiSD->CSFexpansion),
                          multiSD->usingCSF);
  if (!success)
    return false;
  // you should choose the det with highest weight for reference
  multiSD->Dets[0]->ReferenceDeterminant  = 0; // for now
  multiSD->Dets[0]->NumDets               = uniqueConfg_up.size();
  std::vector<ci_configuration2>& list_up = *(multiSD->Dets[0]->ciConfigList);
  list_up.resize(uniqueConfg_up.size());
  for (int i = 0; i < list_up.size(); i++)
  {
    list_up[i].occup.resize(nels_up);
    int cnt = 0;
    for (int k = 0; k < uniqueConfg_up[i].occup.size(); k++)
      if (uniqueConfg_up[i].occup[k])
        list_up[i].occup[cnt++] = k;
    if (cnt != nels_up)
    {
      APP_ABORT("Error in SlaterDetBuilder::createMSDFast, problems with ci configuration list. \n");
    }
  }
  multiSD->Dets[0]->set(multiSD->FirstIndex_up, nels_up, multiSD->Dets[0]->Phi->getOrbitalSetSize());
  multiSD->Dets[1]->ReferenceDeterminant  = 0; // for now
  multiSD->Dets[1]->NumDets               = uniqueConfg_dn.size();
  std::vector<ci_configuration2>& list_dn = *(multiSD->Dets[1]->ciConfigList);
  list_dn.resize(uniqueConfg_dn.size());
  for (int i = 0; i < list_dn.size(); i++)
  {
    list_dn[i].occup.resize(nels_dn);
    int cnt = 0;
    for (int k = 0; k < uniqueConfg_dn[i].occup.size(); k++)
      if (uniqueConfg_dn[i].occup[k])
        list_dn[i].occup[cnt++] = k;
    if (cnt != nels_dn)
    {
      APP_ABORT("Error in SlaterDetBuilder::createMSDFast, problems with ci configuration list. \n");
    }
  }
  multiSD->Dets[1]->set(multiSD->FirstIndex_dn, nels_dn, multiSD->Dets[1]->Phi->getOrbitalSetSize());
  if (multiSD->CSFcoeff->size() == 1)
    optimizeCI = false;
  if (optimizeCI)
  {
    app_log() << "CI coefficients are optimizable. \n";
    std::string resetCI("no");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add(resetCI, "reset_coeff");
    spoAttrib.put(cur);
    if (resetCI == "yes")
    {
      if (multiSD->usingCSF)
        for (int i = 1; i < multiSD->CSFcoeff->size(); i++)
          (*(multiSD->CSFcoeff))[i] = 0;
      else
        for (int i = 1; i < multiSD->C->size(); i++)
          (*(multiSD->C))[i] = 0;
      app_log() << "CI coefficients are reset. \n";
    }
    multiSD->Optimizable = multiSD->CI_Optimizable = true;
    if (multiSD->usingCSF)
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->CSFcoeff[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD->CSFcoeff->size(); i++)
      {
        //std::stringstream sstr;
        //sstr << "CIcoeff" << "_" << i;
        multiSD->myVars->insert(CItags[i], (*(multiSD->CSFcoeff))[i], true, optimize::LINEAR_P);
      }
    }
    else
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->C[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD->C->size(); i++)
      {
        multiSD->myVars->insert(CItags[i], (*(multiSD->C))[i], true, optimize::LINEAR_P);         
      }
    }
  }
  else
  {
    app_log() << "CI coefficients are not optimizable. \n";
    multiSD->CI_Optimizable = false;
  }

  if (multiSD->Dets[0]->Optimizable == true || multiSD->Dets[1]->Optimizable == true)
  {
    //safety checks for orbital optimization
    if (multiSD->Dets[0]->Optimizable != multiSD->Dets[1]->Optimizable)
      APP_ABORT("Optimizing the SPOSet of only one spin is not supported!\n");
    if (multiSD->usingCSF)
      APP_ABORT("Currently, Using CSF is not available with MSJ Orbital Optimization!\n");

    //checks that the hartree fock determinant is the first in the multislater expansion
    for (int i = 0; i < nels_up; i++)
    {
      if ((uniqueConfg_up[0].occup[i] != true) || (uniqueConfg_dn[0].occup[i] != true))
        APP_ABORT(
            "The Hartree Fock Reference Determinant must be first in the Multi-Slater expansion for the input!\n");
    }

    app_log() << "WARNING: Unrestricted Orbital Optimization will be performed. Spin symmetry is not guaranteed to be "
                 "preserved!\n";

    // mark the overall optimization flag
    multiSD->Optimizable = true;

    // The primary purpose of this function is to create all the optimizable orbital rotation parameters.
    // But if orbital rotation parameters were supplied by the user it will also apply a unitary transformation
    // and then remove the orbital rotation parameters
    multiSD->buildOptVariables();
  }

  return success;
}
/*
  bool SlaterDetBuilder::createMSD(MultiSlaterDeterminant* multiSD, xmlNodePtr cur)
  {
     bool success=true;

     std::vector<ci_configuration> uniqueConfg_up, uniqueConfg_dn;
     std::vector<std::string> CItags;
     bool optimizeCI;

     success = readDetList(cur,uniqueConfg_up,uniqueConfg_dn,multiSD->C2node_up, multiSD->C2node_dn,CItags,multiSD->C,optimizeCI,multiSD->nels_up,multiSD->nels_dn);
     if(!success) return false;

     multiSD->resize(uniqueConfg_up.size(),uniqueConfg_dn.size());
     SPOSetProxyForMSD* spo = multiSD->spo_up;
     spo->occup.resize(uniqueConfg_up.size(),multiSD->nels_up);
     for(int i=0; i<uniqueConfg_up.size(); i++)
     {
       int nq=0;
       ci_configuration& ci = uniqueConfg_up[i];
       for(int k=0; k<ci.occup.size(); k++) {
         if(ci.occup[k]) {
           spo->occup(i,nq++) = k;
         }
       }
       DiracDeterminant* adet = new DiracDeterminant((SPOSetPtr) spo,0);
       adet->set(multiSD->FirstIndex_up,multiSD->nels_up);
       multiSD->dets_up.push_back(adet);
     }
     spo = multiSD->spo_dn;
     spo->occup.resize(uniqueConfg_dn.size(),multiSD->nels_dn);
     for(int i=0; i<uniqueConfg_dn.size(); i++)
     {
       int nq=0;
       ci_configuration& ci = uniqueConfg_dn[i];
       for(int k=0; k<ci.occup.size(); k++) {
         if(ci.occup[k]) {
           spo->occup(i,nq++) = k;
         }
       }
       DiracDeterminant* adet = new DiracDeterminant((SPOSetPtr) spo,0);
       adet->set(multiSD->FirstIndex_dn,multiSD->nels_dn);
       multiSD->dets_dn.push_back(adet);
     }

     if(optimizeCI) {
       app_log() <<"CI coefficients are optimizable. ";
       multiSD->Optimizable=true;
       multiSD->myVars.insert(CItags[0],multiSD->C[0],false,optimize::LINEAR_P);
       for(int i=1; i<multiSD->C.size(); i++) {
         //std::stringstream sstr;
         //sstr << "CIcoeff" << "_" << i;
         multiSD->myVars.insert(CItags[i],multiSD->C[i],true,optimize::LINEAR_P);
       }
     }
     return success;
  }
*/
bool SlaterDetBuilder::createMSD(MultiSlaterDeterminant* multiSD, xmlNodePtr cur)
{
  bool success = true;
  std::vector<ci_configuration> uniqueConfg_up, uniqueConfg_dn;
  std::vector<std::string> CItags;
  bool optimizeCI;
  int nels_up = multiSD->nels_up;
  int nels_dn = multiSD->nels_dn;
  success     = readDetList(cur,
                        uniqueConfg_up,
                        uniqueConfg_dn,
                        multiSD->C2node_up,
                        multiSD->C2node_dn,
                        CItags,
                        multiSD->C,
                        optimizeCI,
                        nels_up,
                        nels_dn,
                        multiSD->CSFcoeff,
                        multiSD->DetsPerCSF,
                        multiSD->CSFexpansion,
                        multiSD->usingCSF);
  if (!success)
    return false;
  multiSD->resize(uniqueConfg_up.size(), uniqueConfg_dn.size());
  SPOSetProxyForMSD* spo = multiSD->spo_up;
  spo->occup.resize(uniqueConfg_up.size(), multiSD->nels_up);
  for (int i = 0; i < uniqueConfg_up.size(); i++)
  {
    int nq               = 0;
    ci_configuration& ci = uniqueConfg_up[i];
    for (int k = 0; k < ci.occup.size(); k++)
    {
      if (ci.occup[k])
      {
        spo->occup(i, nq++) = k;
      }
    }
    DiracDeterminantBase* adet;
    if (UseBackflow)
    {
      adet = new DiracDeterminantWithBackflow(targetPtcl, (SPOSetPtr)spo, 0, 0);
    }
    else
    {
      adet = new DiracDeterminant<>((SPOSetPtr)spo, 0);
    }
    adet->set(multiSD->FirstIndex_up, multiSD->nels_up);
    multiSD->dets_up.push_back(adet);
  }
  spo = multiSD->spo_dn;
  spo->occup.resize(uniqueConfg_dn.size(), multiSD->nels_dn);
  for (int i = 0; i < uniqueConfg_dn.size(); i++)
  {
    int nq               = 0;
    ci_configuration& ci = uniqueConfg_dn[i];
    for (int k = 0; k < ci.occup.size(); k++)
    {
      if (ci.occup[k])
      {
        spo->occup(i, nq++) = k;
      }
    }
    DiracDeterminantBase* adet;
    if (UseBackflow)
    {
      adet = new DiracDeterminantWithBackflow(targetPtcl, (SPOSetPtr)spo, 0, 0);
    }
    else
    {
      adet = new DiracDeterminant<>((SPOSetPtr)spo, 0);
    }
    adet->set(multiSD->FirstIndex_dn, multiSD->nels_dn);
    multiSD->dets_dn.push_back(adet);
  }
  if (multiSD->CSFcoeff.size() == 1 || multiSD->C.size() == 1)
    optimizeCI = false;
  if (optimizeCI)
  {
    app_log() << "CI coefficients are optimizable. \n";
    std::string resetCI("no");
    OhmmsAttributeSet spoAttrib;
    spoAttrib.add(resetCI, "reset_coeff");
    spoAttrib.put(cur);
    if (resetCI == "yes")
    {
      if (multiSD->usingCSF)
        for (int i = 1; i < multiSD->CSFcoeff.size(); i++)
          multiSD->CSFcoeff[i] = 0;
      else
        for (int i = 1; i < multiSD->C.size(); i++)
          multiSD->C[i] = 0;
      app_log() << "CI coefficients are reset. \n";
    }
    multiSD->Optimizable = true;
    if (multiSD->usingCSF)
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->CSFcoeff[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD->CSFcoeff.size(); i++)
      {
        //std::stringstream sstr;
        //sstr << "CIcoeff" << "_" << i;
        multiSD->myVars.insert(CItags[i], multiSD->CSFcoeff[i], true, optimize::LINEAR_P);
      }
    }
    else
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->C[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD->C.size(); i++)
      {
        //std::stringstream sstr;
        //sstr << "CIcoeff" << "_" << i;
        multiSD->myVars.insert(CItags[i], multiSD->C[i], true, optimize::LINEAR_P);
      }
    }
  }
  else
  {
    app_log() << "CI coefficients are not optimizable. \n";
    multiSD->Optimizable = false;
  }
  return success;
}


bool SlaterDetBuilder::readDetList(xmlNodePtr cur,
                                   std::vector<ci_configuration>& uniqueConfg_up,
                                   std::vector<ci_configuration>& uniqueConfg_dn,
                                   std::vector<size_t>& C2node_up,
                                   std::vector<size_t>& C2node_dn,
                                   std::vector<std::string>& CItags,
                                   std::vector<ValueType>& coeff,
                                   bool& optimizeCI,
                                   int nels_up,
                                   int nels_dn,
                                   std::vector<ValueType>& CSFcoeff,
                                   std::vector<size_t>& DetsPerCSF,
                                   std::vector<RealType>& CSFexpansion,
                                   bool& usingCSF)
{
  bool success = true;
  uniqueConfg_up.clear();
  uniqueConfg_dn.clear();
  C2node_up.clear();
  C2node_dn.clear();
  CItags.clear();
  coeff.clear();
  CSFcoeff.clear();
  DetsPerCSF.clear();
  CSFexpansion.clear();
  std::vector<ci_configuration> confgList_up;
  std::vector<ci_configuration> confgList_dn;
  std::string optCI    = "no";
  RealType cutoff      = 0.0;
  RealType zero_cutoff = 0.0;
  OhmmsAttributeSet ciAttrib;
  ciAttrib.add(optCI, "optimize");
  ciAttrib.add(optCI, "Optimize");
  ciAttrib.put(cur);
  optimizeCI         = (optCI == "yes");
  xmlNodePtr curRoot = cur, DetListNode = nullptr;
  std::string cname, cname0;
  cur = curRoot->children;
  while (cur != NULL) //check the basis set
  {
    getNodeName(cname, cur);
    if (cname == "detlist")
    {
      DetListNode = cur;
      app_log() << "Found determinant list. \n";
    }
    cur = cur->next;
  }
  size_t NCA, NCB, NEA, NEB, nstates, ndets = 0, count = 0, cnt0 = 0;
  std::string Dettype   = "DETS";
  std::string CSFChoice = "qchem_coeff";
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(NCA, "nca");
  spoAttrib.add(NCB, "ncb");
  spoAttrib.add(NEA, "nea");
  spoAttrib.add(NEB, "neb");
  spoAttrib.add(ndets, "size");
  spoAttrib.add(Dettype, "type");
  spoAttrib.add(nstates, "nstates");
  spoAttrib.add(cutoff, "cutoff");
  spoAttrib.add(zero_cutoff, "zero_cutoff");
  spoAttrib.add(zero_cutoff, "zerocutoff");
  spoAttrib.add(CSFChoice, "sortby");
  spoAttrib.put(DetListNode);
  if (ndets == 0)
  {
    APP_ABORT("size==0 in detlist is not allowed. Use slaterdeterminant in this case.\n");
  }
  if (Dettype == "DETS" || Dettype == "Determinants")
    usingCSF = false;
  else if (Dettype == "CSF")
    usingCSF = true;
  else
  {
    APP_ABORT("Only allowed type in detlist is DETS or CSF.\n");
  }
  if (zero_cutoff > 0)
    app_log() << "  Initializing CI coeffs less than " << zero_cutoff << " to zero." << std::endl;
  // cheating until I fix the converter
  NCA = nels_up - NEA;
  NCB = nels_dn - NEB;
  // mmorales: a little messy for now, clean up later
  cur = DetListNode->children;
  ci_configuration dummyC_alpha;
  ci_configuration dummyC_beta;
  dummyC_alpha.occup.resize(NCA + nstates, false);
  for (size_t i = 0; i < NCA + NEA; i++)
    dummyC_alpha.occup[i] = true;
  dummyC_beta.occup.resize(NCB + nstates, false);
  for (size_t i = 0; i < NCB + NEB; i++)
    dummyC_beta.occup[i] = true;
  RealType sumsq_qc = 0.0;
  RealType sumsq    = 0.0;
  //app_log() <<"alpha reference: \n" <<dummyC_alpha;
  //app_log() <<"beta reference: \n" <<dummyC_beta;
  if (usingCSF)
  {
    app_log() << "Reading CSFs." << std::endl;
    while (cur != NULL) //check the basis set
    {
      getNodeName(cname, cur);
      if (cname == "csf")
      {
        RealType exctLvl, qc_ci = 0.0;
        OhmmsAttributeSet confAttrib;
        std::string tag, OccString;
        #ifdef QMC_COMPLEX
          RealType ci_real = 0.0, ci_imag = 0.0;
          confAttrib.add(ci_real,"coeff_real");
          confAttrib.add(ci_imag,"coeff_imag");
        #else
          RealType ci = 0.0;
          confAttrib.add(ci, "coeff");
        #endif
        confAttrib.add(qc_ci, "qchem_coeff");
        confAttrib.add(tag, "id");
        confAttrib.add(OccString, "occ");
        confAttrib.add(exctLvl, "exctLvl");
        confAttrib.put(cur);
        if (qc_ci == 0.0)
          #ifdef QMC_COMPLEX
            qc_ci = ci_real;
          #else
            qc_ci = ci;
          #endif

          #ifdef QMC_COMPLEX
            ValueType ci(ci_real, ci_imag);
          #endif
        //Can discriminate based on any of 3 criterion
        if (((std::abs(qc_ci) < cutoff) && (CSFChoice == "qchem_coeff")) ||
            ((CSFChoice == "exctLvl") && (exctLvl > cutoff)) || ((CSFChoice == "coeff") && (std::abs(ci) < cutoff)))
        {
          cur = cur->next;
          cnt0++;
          continue;
        }
        cnt0++;
        if (std::abs(qc_ci) < zero_cutoff)
          #ifdef QMC_COMPLEX
            ci_real = 0.0;
          #else
            ci = 0.0;
          #endif
        CSFcoeff.push_back(ci);
        sumsq_qc += qc_ci * qc_ci;
        DetsPerCSF.push_back(0);
        CItags.push_back(tag);
        count++;
        xmlNodePtr csf = cur->children;
        while (csf != NULL)
        {
          getNodeName(cname0, csf);
          if (cname0 == "det")
          {
            std::string alpha, beta, tag0;
            RealType coef = 0.0;
            OhmmsAttributeSet detAttrib;
            detAttrib.add(tag0, "id");
            detAttrib.add(coef, "coeff");
            detAttrib.add(beta, "beta");
            detAttrib.add(alpha, "alpha");
            detAttrib.put(csf);
            size_t nq = 0;
            if (alpha.size() < nstates)
            {
              std::cerr << "alpha: " << alpha << std::endl;
              APP_ABORT("Found incorrect alpha determinant label. size < nca+nstates");
            }
            for (size_t i = 0; i < nstates; i++)
            {
              if (alpha[i] != '0' && alpha[i] != '1')
              {
                std::cerr << alpha << std::endl;
                APP_ABORT("Found incorrect determinant label.");
              }
              if (alpha[i] == '1')
                nq++;
            }
            if (nq != NEA)
            {
              std::cerr << "alpha: " << alpha << std::endl;
              APP_ABORT("Found incorrect alpha determinant label. noccup != nca+nea");
            }
            nq = 0;
            if (beta.size() < nstates)
            {
              std::cerr << "beta: " << beta << std::endl;
              APP_ABORT("Found incorrect beta determinant label. size < ncb+nstates");
            }
            for (size_t i = 0; i < nstates; i++)
            {
              if (beta[i] != '0' && beta[i] != '1')
              {
                std::cerr << beta << std::endl;
                APP_ABORT("Found incorrect determinant label.");
              }
              if (beta[i] == '1')
                nq++;
            }
            if (nq != NEB)
            {
              std::cerr << "beta: " << beta << std::endl;
              APP_ABORT("Found incorrect beta determinant label. noccup != ncb+neb");
            }
            //app_log() <<" <ci id=\"coeff_" <<ntot++ <<"\" coeff=\"" <<ci*coef <<"\" alpha=\"" <<alpha <<"\" beta=\"" <<beta <<"\" />" << std::endl;
            DetsPerCSF.back()++;
            CSFexpansion.push_back(coef);
            coeff.push_back(coef * ci);
            confgList_up.push_back(dummyC_alpha);
            for (size_t i = 0; i < NCA; i++)
              confgList_up.back().occup[i] = true;
            for (size_t i = NCA; i < NCA + nstates; i++)
              confgList_up.back().occup[i] = (alpha[i - NCA] == '1');
            confgList_dn.push_back(dummyC_beta);
            for (size_t i = 0; i < NCB; i++)
              confgList_dn.back().occup[i] = true;
            for (size_t i = NCB; i < NCB + nstates; i++)
              confgList_dn.back().occup[i] = (beta[i - NCB] == '1');
          } // if(name=="det")
          csf = csf->next;
        } // csf loop
        if (DetsPerCSF.back() == 0)
        {
          APP_ABORT("Found empty CSF (no det blocks).");
        }
      } // if (name == "csf")
      cur = cur->next;
    }
    if (cnt0 != ndets)
    {
      std::cerr << "count, ndets: " << cnt0 << "  " << ndets << std::endl;
      APP_ABORT("Problems reading determinant ci_configurations. Found a number of determinants inconsistent with xml "
                "file size parameter.\n");
    }
    //if(!usingCSF)
    //  if(confgList_up.size() != ndets || confgList_dn.size() != ndets || coeff.size() != ndets) {
    //    APP_ABORT("Problems reading determinant ci_configurations.");
    //  }
    C2node_up.resize(coeff.size());
    C2node_dn.resize(coeff.size());
    app_log() << "Found " << coeff.size() << " terms in the MSD expansion.\n";
    RealType sumsq = 0.0;
    for (size_t i = 0; i < coeff.size(); i++)
      sumsq += std::abs(coeff[i] * coeff[i]);
    app_log() << "Norm of ci vector (sum of ci^2): " << sumsq << std::endl;
    app_log() << "Norm of qchem ci vector (sum of qchem_ci^2): " << sumsq_qc << std::endl;
    for (size_t i = 0; i < confgList_up.size(); i++)
    {
      bool found = false;
      size_t k   = -1;
      for (size_t j = 0; j < uniqueConfg_up.size(); j++)
      {
        if (confgList_up[i] == uniqueConfg_up[j])
        {
          found = true;
          k     = j;
          break;
        }
      }
      if (found)
      {
        C2node_up[i] = k;
      }
      else
      {
        uniqueConfg_up.push_back(confgList_up[i]);
        C2node_up[i] = uniqueConfg_up.size() - 1;
      }
    }
    for (size_t i = 0; i < confgList_dn.size(); i++)
    {
      bool found = false;
      size_t k   = -1;
      for (size_t j = 0; j < uniqueConfg_dn.size(); j++)
      {
        if (confgList_dn[i] == uniqueConfg_dn[j])
        {
          found = true;
          k     = j;
          break;
        }
      }
      if (found)
      {
        C2node_dn[i] = k;
      }
      else
      {
        uniqueConfg_dn.push_back(confgList_dn[i]);
        C2node_dn[i] = uniqueConfg_dn.size() - 1;
      }
    }
  }
  else
  {
    app_log() << "Reading CI expansion." << std::endl;

    int cntup = 0;
    int cntdn = 0;
    std::unordered_map<std::string, int> MyMapUp;
    std::unordered_map<std::string, int> MyMapDn;
    while (cur != NULL) //check the basis set
    {
      getNodeName(cname, cur);
      if (cname == "configuration" || cname == "ci")
      {
        RealType qc_ci = 0.0;
        std::string alpha, beta, tag;
        OhmmsAttributeSet confAttrib;
        #ifdef QMC_COMPLEX
        RealType ci_real = 0.0, ci_imag = 0.0;
        confAttrib.add(ci_real, "coeff_real");
        confAttrib.add(ci_imag, "coeff_imag");
        ValueType ci(ci_real, ci_imag);
        #else
        RealType ci = 0.0;
        confAttrib.add(ci, "coeff");
        #endif
        confAttrib.add(qc_ci, "qchem_coeff");
        confAttrib.add(alpha, "alpha");
        confAttrib.add(beta, "beta");
        confAttrib.add(tag, "id");
        confAttrib.put(cur);


        //Will always loop through the whole determinant set as no assumption on the order of the determinant is made
        if (std::abs(ci) < cutoff)
        {
          cur = cur->next;
          continue;
        }

        for (size_t i = 0; i < nstates; i++)
        {
          if (alpha[i] != '0' && alpha[i] != '1')
          {
            std::cerr << alpha << std::endl;
            APP_ABORT("Found incorrect determinant label.");
          }
          if (beta[i] != '0' && beta[i] != '1')
          {
            std::cerr << beta << std::endl;
            APP_ABORT("Found incorrect determinant label.");
          }
        }

        if (alpha.size() < nstates)
        {
          std::cerr << "alpha: " << alpha << std::endl;
          APP_ABORT("Found incorrect alpha determinant label. size < nca+nstates");
        }
        if (beta.size() < nstates)
        {
          std::cerr << "beta: " << beta << std::endl;
          APP_ABORT("Found incorrect beta determinant label. size < ncb+nstates");
        }

        coeff.push_back(ci);
        CItags.push_back(tag);

        std::unordered_map<std::string, int>::const_iterator gotup = MyMapUp.find(alpha);


        if (gotup == MyMapUp.end())
        {
          uniqueConfg_up.push_back(dummyC_alpha);
          uniqueConfg_up.back().add_occupation(alpha);
          C2node_up.push_back(cntup);
          MyMapUp.insert(std::pair<std::string, int>(alpha, cntup));
          cntup++;
        }
        else
        {
          C2node_up.push_back(gotup->second);
        }

        std::unordered_map<std::string, int>::const_iterator gotdn = MyMapDn.find(beta);

        if (gotdn == MyMapDn.end())
        {
          uniqueConfg_dn.push_back(dummyC_beta);
          uniqueConfg_dn.back().add_occupation(beta);
          C2node_dn.push_back(cntdn);
          MyMapDn.insert(std::pair<std::string, int>(beta, cntdn));
          cntdn++;
        }
        else
        {
          C2node_dn.push_back(gotdn->second);
        }

        //if (qc_ci == 0.0)
        //  qc_ci = ci;

        cnt0++;
        sumsq_qc += qc_ci * qc_ci;
        sumsq += std::abs(ci * ci);
      }
      cur = cur->next;
    }

    app_log() << "Found " << coeff.size() << " terms in the MSD expansion.\n";
    app_log() << "Norm of ci vector (sum of ci^2): " << sumsq << std::endl;
    app_log() << "Norm of qchem ci vector (sum of qchem_ci^2): " << sumsq_qc << std::endl;

  } //usingCSF

  app_log() << "Found " << uniqueConfg_up.size() << " unique up determinants.\n";
  app_log() << "Found " << uniqueConfg_dn.size() << " unique down determinants.\n";

  return success;
}


bool SlaterDetBuilder::readDetListH5(xmlNodePtr cur,
                                     std::vector<ci_configuration>& uniqueConfg_up,
                                     std::vector<ci_configuration>& uniqueConfg_dn,
                                     std::vector<size_t>& C2node_up,
                                     std::vector<size_t>& C2node_dn,
                                     std::vector<std::string>& CItags,
                                     std::vector<RealType>& coeff,
                                     bool& optimizeCI,
                                     int nels_up,
                                     int nels_dn)
{
  bool success = true;
  uniqueConfg_up.clear();
  uniqueConfg_dn.clear();
  C2node_up.clear();
  C2node_dn.clear();
  CItags.clear();
  coeff.clear();
  std::string CICoeffH5path("");
  std::vector<ci_configuration> confgList_up;
  std::vector<ci_configuration> confgList_dn;
  std::vector<RealType> CIcoeff;
  std::vector<std::string> ConfigTag;
  std::string optCI = "no";
  RealType cutoff   = 0.0;
  OhmmsAttributeSet ciAttrib;
  ciAttrib.add(optCI, "optimize");
  ciAttrib.put(cur);
  optimizeCI         = (optCI == "yes");
  xmlNodePtr curRoot = cur, DetListNode = nullptr;
  std::string cname, cname0, multidetH5path;
  cur = curRoot->children;
  while (cur != NULL) //check the basis set
  {
    getNodeName(cname, cur);
    if (cname == "detlist")
    {
      DetListNode = cur;
      app_log() << "Found determinant list. \n";
    }
    cur = cur->next;
  }
  size_t NCA, NCB, NEA, NEB, nstates, ndets = 0;
  size_t H5_ndets, H5_nstates;
  int N_int;
  const int bit_kind  = 64;
  std::string Dettype = "DETS";
  RealType sumsq      = 0.0;
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(NCA, "nca");
  spoAttrib.add(NCB, "ncb");
  spoAttrib.add(NEA, "nea");
  spoAttrib.add(NEB, "neb");
  spoAttrib.add(ndets, "size");
  spoAttrib.add(Dettype, "type");
  spoAttrib.add(nstates, "nstates");
  spoAttrib.add(cutoff, "cutoff");
  spoAttrib.add(multidetH5path, "href");
  spoAttrib.add(CICoeffH5path, "opt_coeffs");
  spoAttrib.put(DetListNode);
  if (ndets == 0)
  {
    APP_ABORT("size==0 in detlist is not allowed. Use slaterdeterminant in this case.\n");
  }
  if (Dettype != "DETS" && Dettype != "Determinants")
    APP_ABORT("Reading from HDF5 is only enabled for CI DETS. Must be accessed through (type=\"DETS\") or "
              "(type=\"Determinants\") .\n");
  app_log() << "Reading CI expansion from HDF5:" << multidetH5path << std::endl;

  hdf_archive hin;
  if (!hin.open(multidetH5path.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }

  if (!hin.push("MultiDet"))
  {
    std::cerr << "Could not open Multidet Group in H5 file" << std::endl;
    abort();
  }

  hin.read(H5_ndets, "NbDet");
  if (ndets != H5_ndets)
  {
    std::cerr << "Number of determinants in H5 file (" << H5_ndets << ") different from number of dets in XML ("
              << ndets << ")" << std::endl;
    abort();
  }

  hin.read(H5_nstates, "nstate");
  if (nstates != H5_nstates)
  {
    std::cerr << "Number of states/orbitals in H5 file (" << H5_nstates
              << ") different from number of states/orbitals in XML (" << nstates << ")" << std::endl;
    abort();
  }

  hin.read(N_int, "Nbits");
  CIcoeff.resize(ndets);
  ConfigTag.resize(ndets);

  hin.read(CIcoeff, "Coeff");

  ///IF OPTIMIZED COEFFICIENTS ARE PRESENT IN opt_coeffs Path
  ///THEY ARE READ FROM DIFFERENT HDF5 the replace the previous coeff
  ///It is important to still read all old coeffs and only replace the optimized ones
  ///in order to keep coherence with the cutoff on the number of determinants 
  ///REMEMBER!! FIRST COEFF IS FIXED. THEREFORE WE DO NOT REPLACE IT!!!
  if(CICoeffH5path!="")
  {

    int OptCiSize=0;
    std::vector<RealType> CIcoeffopt;
    hdf_archive coeffin;
    if (!coeffin.open(CICoeffH5path.c_str(), H5F_ACC_RDONLY))
    {
      std::cerr << "Could not open H5 file containing Optimized Coefficients" << std::endl;
      abort();
    }

    if (!coeffin.push("MultiDet"))
    {
      std::cerr << "Could not open Multidet Group in H5 file" << std::endl;
      abort();
    }
    coeffin.read(OptCiSize, "NbDet");
    CIcoeffopt.resize(OptCiSize);
    
    coeffin.read(CIcoeffopt, "Coeff");
    coeffin.close();

    for (int i=0; i<OptCiSize;i++)
          CIcoeff[i+1]=CIcoeffopt[i];
  
    app_log()<<"The first "<<OptCiSize<<" Optimized coefficients were substituted to the original set of coefficients." <<std::endl;
  }

  Matrix<long int> tempAlpha(ndets, N_int);
  hin.read(tempAlpha, "CI_Alpha");

  Matrix<long int> tempBeta(ndets, N_int);
  hin.read(tempBeta, "CI_Beta");

  std::string MyCIAlpha, MyCIBeta;
  MyCIAlpha.resize(nstates);
  MyCIBeta.resize(nstates);

  ci_configuration dummyC_alpha;
  ci_configuration dummyC_beta;

  dummyC_alpha.occup.resize(nstates, false);
  dummyC_beta.occup.resize(nstates, false);

  for (size_t i = 0; i < nstates; i++)
    dummyC_alpha.occup[i] = true;

  for (size_t i = 0; i < nstates; i++)
    dummyC_beta.occup[i] = true;

  hin.close();
  app_log() << " Done reading CIs from H5!!" << std::endl;

  int cntup = 0;
  int cntdn = 0;
  std::unordered_map<std::string, int> MyMapUp;
  std::unordered_map<std::string, int> MyMapDn;

  app_log() << " Sorting unique CIs" << std::endl;
  for (int ni = 0; ni < ndets; ni++)
  {
    if (std::abs(CIcoeff[ni]) < cutoff)
      continue;

    int j = 0;
    for (int k = 0; k < N_int; k++)
    {
      long int a               = tempAlpha[ni][k];
      std::bitset<bit_kind> a2 = a;

      auto b  = tempBeta[ni][k];
      auto b2 = std::bitset<bit_kind>(b);

      for (int i = 0; i < bit_kind; i++)
      {
        if (j < nstates)
        {
          MyCIAlpha[j] = a2[i] ? '1' : '0';
          MyCIBeta[j]  = b2[i] ? '1' : '0';
          j++;
        }
      }
    }
    coeff.push_back(CIcoeff[ni]);
    std::ostringstream h5tag;
    h5tag << "CIcoeff_" << ni;
    CItags.push_back(h5tag.str());

    std::unordered_map<std::string, int>::const_iterator gotup = MyMapUp.find(MyCIAlpha);

    if (gotup == MyMapUp.end())
    {
      uniqueConfg_up.push_back(dummyC_alpha);
      uniqueConfg_up.back().add_occupation(MyCIAlpha);
      C2node_up.push_back(cntup);
      MyMapUp.insert(std::pair<std::string, int>(MyCIAlpha, cntup));
      cntup++;
    }
    else
    {
      C2node_up.push_back(gotup->second);
    }

    std::unordered_map<std::string, int>::const_iterator gotdn = MyMapDn.find(MyCIBeta);

    if (gotdn == MyMapDn.end())
    {
      uniqueConfg_dn.push_back(dummyC_beta);
      uniqueConfg_dn.back().add_occupation(MyCIBeta);
      C2node_dn.push_back(cntdn);
      MyMapDn.insert(std::pair<std::string, int>(MyCIBeta, cntdn));
      cntdn++;
    }
    else
    {
      C2node_dn.push_back(gotdn->second);
    }

    sumsq += CIcoeff[ni] * CIcoeff[ni];
  }

  app_log() << " Done Sorting unique CIs" << std::endl;
  app_log() << "Found " << coeff.size() << " terms in the MSD expansion.\n";
  app_log() << "Norm of ci vector (sum of ci^2): " << sumsq << std::endl;

  app_log() << "Found " << uniqueConfg_up.size() << " unique up determinants.\n";
  app_log() << "Found " << uniqueConfg_dn.size() << " unique down determinants.\n";


  return success;
}
} // namespace qmcplusplus

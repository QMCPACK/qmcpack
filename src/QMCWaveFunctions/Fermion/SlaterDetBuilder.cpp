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

#include <type_traits>
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "SlaterDetBuilder.h"
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
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
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
SlaterDetBuilder::SlaterDetBuilder(Communicate* comm,
                                   SPOSetBuilderFactory& factory,
                                   ParticleSet& els,
                                   TrialWaveFunction& psi,
                                   PtclPoolType& psets)
    : WaveFunctionComponentBuilder(comm, els), sposet_builder_factory_(factory), targetPsi(psi), ptclPool(psets)
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
std::unique_ptr<WaveFunctionComponent> SlaterDetBuilder::buildComponent(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "put(xmlNodePtr)");
  ///save the current node
  xmlNodePtr curRoot = cur;
  xmlNodePtr BFnode  = nullptr;
  std::string cname, tname;
  std::map<std::string, SPOSetPtr> spomap;
  bool multiDet = false;
  std::string msd_algorithm;

  if (sposet_builder_factory_.empty())
  { //always create one, using singleton and just to access the member functions
    app_warning() << "!!!!!!! Deprecated input style: creating SPO set inside determinantset. Support for this usage "
                     "will soon be removed. SPO sets should be built outside."
                  << std::endl;
    sposet_builder_factory_.createSPOSetBuilder(curRoot);
  }

  //check the basis set
  cur = curRoot->children;
  while (cur != NULL) //check the basis set
  {
    getNodeName(cname, cur);
    if (cname == sposet_tag)
    {
      app_warning() << "!!!!!!! Deprecated input style: creating SPO set inside determinantset. Support for this usage "
                       "will soon be removed. SPO sets should be built outside."
                    << std::endl;
      app_log() << "Creating SPOSet in SlaterDetBuilder::put(xmlNodePtr cur).\n";
      std::string spo_name;
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add(spo_name, "name");
      spoAttrib.put(cur);
      app_log() << "spo_name = " << spo_name << std::endl;
      SPOSetPtr spo = sposet_builder_factory_.createSPOSet(cur);
      if (spomap.find(spo_name) != spomap.end())
      {
        app_error() << "SPOSet name \"" << spo_name << "\" is already in use.\n";
        abort();
      }
      spomap[spo_name] = spo;
      spo->setName(spo_name);
      assert(spomap.find(spo_name) != spomap.end());
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

  cur = curRoot->children;
  while (cur != NULL)
  {
    getNodeName(cname, cur);
    if (cname == sd_tag)
    {
      app_summary() << std::endl;
      app_summary() << "   Single Slater determinant" << std::endl;
      app_summary() << "   -------------------------" << std::endl;
      app_summary() << std::endl;

      multiDet = false;
      if (slaterdet_0)
      {
        APP_ABORT("slaterdet is already instantiated.");
      }
      if (UseBackflow)
      {
        app_summary() << "    Using backflow transformation." << std::endl;
        slaterdet_0 = std::make_unique<SlaterDetWithBackflow>(targetPtcl, BFTrans);
      }
      else
        slaterdet_0 = std::make_unique<SlaterDeterminant_t>(targetPtcl);

      size_t spin_group = 0;
      xmlNodePtr tcur   = cur->children;
      while (tcur != NULL)
      {
        getNodeName(tname, tcur);
        if (tname == det_tag || tname == rn_tag)
        {
          if (spin_group >= targetPtcl.groups())
          {
            std::ostringstream err_msg;
            err_msg << "Need only " << targetPtcl.groups() << " determinant input elements. Found more." << std::endl;
            throw std::runtime_error(err_msg.str());
          }
          if (putDeterminant(tcur, spin_group))
            spin_group++;
        }
        tcur = tcur->next;
      }
      if (spin_group < targetPtcl.groups())
      {
        std::ostringstream err_msg;
        err_msg << "Not enough determinant input elements. "
                << "Need " << targetPtcl.groups() << " but only found " << spin_group << "." << std::endl;
        throw std::runtime_error(err_msg.str());
      }
    }
    else if (cname == multisd_tag)
    {
      app_summary() << std::endl;
      app_summary() << "   Multi Slater determinants" << std::endl;
      app_summary() << "   -------------------------" << std::endl;
      app_summary() << std::endl;

      multiDet = true;

      if (slaterdet_0)
        APP_ABORT("can't combine slaterdet with multideterminant.");

      if (multislaterdet_0 || multislaterdetfast_0)
        APP_ABORT("multideterminant is already instantiated.");

      const int nGroups = targetPtcl.groups();
      std::vector<std::string> spoNames(nGroups);

      std::string fastAlg;
      OhmmsAttributeSet spoAttrib;
      for (int grp = 0; grp < nGroups; grp++)
        spoAttrib.add(spoNames[grp], "spo_" + std::to_string(grp));
      if (nGroups == 2)
      { //for legacy
        spoAttrib.add(spoNames[0], "spo_up");
        spoAttrib.add(spoNames[1], "spo_dn");
      }
      spoAttrib.add(fastAlg, "Fast", {"", "yes", "no"}, TagStatus::DELETED);
      spoAttrib.add(msd_algorithm, "algorithm", {"precomputed_table_method", "table_method", "all_determinants"});
      spoAttrib.put(cur);

      //new format
      std::vector<std::unique_ptr<SPOSet>> spo_clones;

      for (int grp = 0; grp < nGroups; grp++)
      {
        SPOSetPtr spo_tmp = sposet_builder_factory_.getSPOSet(spoNames[grp]);
        if (spo_tmp == nullptr)
        {
          std::stringstream err_msg;
          err_msg << "In SlaterDetBuilder: SPOSet \"" << spoNames[grp]
                  << "\" is not found. Expected for MultiSlaterDeterminant." << std::endl;
          myComm->barrier_and_abort(err_msg.str());
        }
        spo_clones.emplace_back(spo_tmp->makeClone());
      }


      if (msd_algorithm == "precomputed_table_method" || msd_algorithm == "table_method")
      {
        app_summary() << "    Using Bryan's table method." << std::endl;
        if (UseBackflow)
        {
          APP_ABORT("Backflow is not implemented with the table method.");
        }

        bool spinor = targetPtcl.is_spinor_;
        std::vector<std::unique_ptr<MultiDiracDeterminant>> dets;
        for (int grp = 0; grp < nGroups; grp++)
        {
          app_log() << "      Creating base determinant (" << grp << ") for MSD expansion. \n";
          dets.emplace_back(std::make_unique<MultiDiracDeterminant>(std::move(spo_clones[grp]), spinor));
        }

        if (msd_algorithm == "precomputed_table_method")
        {
          app_summary() << "    Using the table method with precomputing. Faster" << std::endl;
          multislaterdetfast_0 = std::make_unique<MultiSlaterDeterminantFast>(targetPtcl, std::move(dets), true);
        }
        else
        {
          app_summary() << "    Using the table method without precomputing. Slower." << std::endl;
          multislaterdetfast_0 = std::make_unique<MultiSlaterDeterminantFast>(targetPtcl, std::move(dets), false);
        }

        multislaterdetfast_0->initialize();
        createMSDFast(multislaterdetfast_0->Dets, *multislaterdetfast_0->C2node, *multislaterdetfast_0->C,
                      *multislaterdetfast_0->CSFcoeff, *multislaterdetfast_0->DetsPerCSF,
                      *multislaterdetfast_0->CSFexpansion, multislaterdetfast_0->usingCSF,
                      *multislaterdetfast_0->myVars, multislaterdetfast_0->Optimizable,
                      multislaterdetfast_0->CI_Optimizable, cur);

        // The primary purpose of this function is to create all the optimizable orbital rotation parameters.
        // But if orbital rotation parameters were supplied by the user it will also apply a unitary transformation
        // and then remove the orbital rotation parameters
        multislaterdetfast_0->buildOptVariables();
      }
      else
      {
        if (nGroups != 2)
        {
          PRE.error("MSD using all_determinants algorithm requires two particle species.");
          return nullptr;
        }
        app_summary() << "    Using a list of determinants for multi-deterimant expansion." << std::endl;
        auto spo_up =
            std::make_unique<SPOSetProxyForMSD>(std::move(spo_clones[0]), targetPtcl.first(0), targetPtcl.last(0));
        auto spo_dn =
            std::make_unique<SPOSetProxyForMSD>(std::move(spo_clones[1]), targetPtcl.first(1), targetPtcl.last(1));
        if (UseBackflow)
        {
          app_summary() << "    Using backflow transformation." << std::endl;
          multislaterdet_0 = std::make_unique<MultiSlaterDeterminantWithBackflow>(targetPtcl, std::move(spo_up),
                                                                                  std::move(spo_dn), BFTrans);
          createMSD(*multislaterdet_0, cur);
        }
        else
        {
          multislaterdet_0 = std::make_unique<MultiSlaterDeterminant>(targetPtcl, std::move(spo_up), std::move(spo_dn));
          createMSD(*multislaterdet_0, cur);
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
      if (msd_algorithm == "all_determinants")
        multislaterdet_0->setBF(BFTrans);
      else
        myComm->barrier_and_abort("Backflow is not supported by Multi-Slater determinants using the table method!");
    }
    else
    {
      slaterdet_0->setBF(BFTrans);
      if (BFTrans->isOptimizable())
        slaterdet_0->Optimizable = true;
    }
  }
  //only single slater determinant
  if (multiDet)
  {
    if (msd_algorithm == "all_determinants")
      return std::move(multislaterdet_0);
    else
      return std::move(multislaterdetfast_0);
  }
  else
    return std::move(slaterdet_0);
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
  std::string sposet_name;
  std::string basisName("invalid");
  std::string detname("0"), refname("0");
  std::string s_detSize("0");

  OhmmsAttributeSet aAttrib;
  aAttrib.add(basisName, "basisset");
  aAttrib.add(detname, "id");
  aAttrib.add(sposet_name, "sposet");
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
#if defined(ENABLE_OFFLOAD)
  std::string use_batch("yes");
#else
  std::string use_batch("no");
#endif
  std::string useGPU;
  int delay_rank(0);

  OhmmsAttributeSet sdAttrib;
  sdAttrib.add(delay_rank, "delay_rank");
  sdAttrib.add(optimize, "optimize");
  sdAttrib.add(use_batch, "batch");
#if defined(ENABLE_CUDA) || defined(ENABLE_OFFLOAD)
  sdAttrib.add(useGPU, "gpu", {"yes", "no"});
#endif
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
  if (sposet_name.empty())
    sposet_name = detname;

  app_summary() << std::endl;
  app_summary() << "     Determinant" << std::endl;
  app_summary() << "     -----------" << std::endl;
  app_summary() << "      Name: " << detname << "   Spin group: " << spin_group << "   SPO name: " << sposet_name
                << std::endl;
  app_summary() << std::endl;

  SPOSetPtr psi = sposet_builder_factory_.getSPOSet(sposet_name);
  //check if the named sposet exists
  if (psi == 0)
  {
    app_warning() << "!!!!!!! Deprecated input style: creating SPO set inside determinantset. Support for this usage "
                     "will soon be removed. SPO sets should be built outside."
                  << std::endl;
    app_log() << "      Create a new SPO set " << sposet_name << std::endl;
    psi = sposet_builder_factory_.createSPOSet(cur);
  }
  psi->checkObject();
  std::unique_ptr<SPOSet> psi_clone(psi->makeClone());

  int firstIndex = targetPtcl.first(spin_group);
  int lastIndex  = targetPtcl.last(spin_group);
  if (firstIndex == lastIndex)
    return true;

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
    if (lastIndex - firstIndex >= 192)
      delay_rank = 32;
    else
      delay_rank = 1;
    app_summary() << "      Setting delay_rank to default value " << delay_rank << std::endl;
  }

  if (delay_rank > 1)
    app_summary() << "      Using rank-" << delay_rank << " delayed update" << std::endl;
  else
    app_summary() << "      Using rank-1 Sherman-Morrison Fahy update (SM1)" << std::endl;

  DiracDeterminantBase* adet = 0;

  //TODO: the switch logic should be improved as we refine the input tags.
#if defined(QMC_CUDA)
  app_summary() << "      Using legacy CUDA acceleration." << std::endl;
  adet = new DiracDeterminantCUDA(std::move(psi_clone), firstIndex);
#else
  if (UseBackflow)
  {
    app_summary() << "      Using backflow transformation." << std::endl;
    adet = new DiracDeterminantWithBackflow(targetPtcl, std::move(psi_clone), BFTrans, firstIndex);
  }
  else
  {
    if (use_batch == "yes")
    {
      app_summary() << "      Using walker batching." << std::endl;
#if defined(ENABLE_CUDA) && defined(ENABLE_OFFLOAD)
      if (useGPU == "yes")
      {
        app_summary() << "      Running on an NVIDIA GPU via CUDA acceleration and OpenMP offload." << std::endl;
        adet = new DiracDeterminantBatched<
            MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>(std::move(psi_clone),
                                                                                         firstIndex);
      }
      else
#endif
      {
        app_summary() << "      Running on an accelerator via OpenMP offload. Only SM1 update is supported. "
                         "delay_rank is ignored."
                      << std::endl;
        adet = new DiracDeterminantBatched<>(std::move(psi_clone), firstIndex);
      }
    }
    else
    {
#if defined(ENABLE_CUDA)
      if (useGPU == "yes")
      {
        app_summary() << "      Running on an NVIDIA GPU via CUDA acceleration." << std::endl;
        adet = new DiracDeterminant<DelayedUpdateCUDA<ValueType, QMCTraits::QTFull::ValueType>>(std::move(psi_clone),
                                                                                                firstIndex);
      }
      else
#endif
      {
        app_summary() << "      Running on CPU." << std::endl;
        adet = new DiracDeterminant<>(std::move(psi_clone), firstIndex);
      }
    }
  }
#endif

  adet->set(firstIndex, lastIndex - firstIndex, delay_rank);
#ifdef QMC_CUDA
  targetPsi.setndelay(delay_rank);
#endif
  slaterdet_0->add(adet, spin_group);

  app_log() << std::endl;
  app_log().flush();
  return true;
}

bool SlaterDetBuilder::createMSDFast(std::vector<std::unique_ptr<MultiDiracDeterminant>>& Dets,
                                     std::vector<std::vector<size_t>>& C2nodes,
                                     std::vector<ValueType>& C,
                                     std::vector<ValueType>& CSFcoeff,
                                     std::vector<size_t>& DetsPerCSF,
                                     std::vector<RealType>& CSFexpansion,
                                     bool& usingCSF,
                                     opt_variables_type& myVars,
                                     bool& Optimizable,
                                     bool& CI_Optimizable,
                                     xmlNodePtr cur)
{
  bool success = true;


  bool optimizeCI;

  const int nGroups = targetPtcl.groups();
  assert(nGroups == Dets.size());
  std::vector<int> nptcls(nGroups);
  for (int grp = 0; grp < nGroups; grp++)
    nptcls[grp] = targetPtcl.groupsize(grp);

  std::vector<std::vector<ci_configuration>> uniqueConfgs(nGroups);
  std::vector<std::string> CItags;

  //Check id multideterminants are in HDF5

  xmlNodePtr curTemp = cur, DetListNode = nullptr;
  curTemp = curTemp->children;
  while (curTemp != NULL) //check the basis set
  {
    std::string cname;
    getNodeName(cname, curTemp);
    if (cname == "detlist")
    {
      DetListNode = curTemp;
    }
    curTemp = curTemp->next;
  }
  XMLAttrString HDF5Path(DetListNode, "href");
  if (HDF5Path != "")
  {
    app_log() << "Found Multideterminants in H5 File" << std::endl;
    success = readDetListH5(cur, uniqueConfgs, C2nodes, CItags, C, optimizeCI, nptcls);
  }
  else
  {
    success = readDetList(cur, uniqueConfgs, C2nodes, CItags, C, optimizeCI, nptcls, CSFcoeff, DetsPerCSF, CSFexpansion,
                          usingCSF);
  }
  if (!success)
    return false;

  for (int grp = 0; grp < nGroups; grp++)
  {
    std::vector<ci_configuration2>& list = Dets[grp]->getCIConfigList();
    list.resize(uniqueConfgs[grp].size());
    for (int i = 0; i < list.size(); i++)
    {
      list[i].occup.resize(nptcls[grp]);
      int cnt = 0;
      for (int k = 0; k < uniqueConfgs[grp][i].occup.size(); k++)
        if (uniqueConfgs[grp][i].occup[k])
          list[i].occup[cnt++] = k;
      if (cnt != nptcls[grp])
      {
        APP_ABORT("Error in SlaterDetBuilder::createMSDFast for ptcl group "
                  << grp << ", problems with ci configuration list. \n");
      }
    }
    // you should choose the det with highest weight for reference. for now choosing 0
    Dets[grp]->set(targetPtcl.first(grp), nptcls[grp], 0);
  }

  if (CSFcoeff.size() == 1)
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
      if (usingCSF)
        for (int i = 1; i < CSFcoeff.size(); i++)
          CSFcoeff[i] = 0;
      else
        for (int i = 1; i < C.size(); i++)
          C[i] = 0;
      app_log() << "CI coefficients are reset. \n";
    }
    Optimizable = CI_Optimizable = true;
    if (usingCSF)
      for (int i = 1; i < CSFcoeff.size(); i++)
        myVars.insert(CItags[i], CSFcoeff[i], true, optimize::LINEAR_P);
    else
      for (int i = 1; i < C.size(); i++)
        myVars.insert(CItags[i], C[i], true, optimize::LINEAR_P);
  }
  else
  {
    app_log() << "CI coefficients are not optimizable. \n";
    CI_Optimizable = false;
  }

  bool any_optimizable = false;
  for (int grp = 0; grp < nGroups; grp++)
  {
    if (Dets[grp]->Optimizable == true)
    {
      any_optimizable = true;
      break;
    }
  }
  if (any_optimizable)
  {
    for (int grp = 0; grp < nGroups; grp++)
    {
      if (Dets[grp]->Optimizable != true)
        APP_ABORT("Optimizing the SPOSet of only only species is not supported!\n");
    }
    if (usingCSF)
      APP_ABORT("Currently, Using CSF is not available with MSJ Orbital Optimization!\n");

    for (int grp = 0; grp < nGroups; grp++)
    {
      for (int i = 0; i < nptcls[grp]; i++)
      {
        if (uniqueConfgs[grp][0].occup[i] != true)
          APP_ABORT(
              "The Hartee Fock Reference Determinant must be the first in the Multi-Slater expansion for the input!\n");
      }
    }
    app_warning() << "Unrestricted Orbital Optimization will be performed. Spin symmetry is not guaranteed to be "
                     "preserved!\n";

    Optimizable = true;
  }

  return success;
}

bool SlaterDetBuilder::createMSD(MultiSlaterDeterminant& multiSD, xmlNodePtr cur)
{
  bool success = true;
  std::vector<std::vector<ci_configuration>> uniqueConfgs(2);
  std::vector<std::string> CItags;
  bool optimizeCI;
  std::vector<int> nels(2);
  nels[0] = multiSD.nels_up;
  nels[1] = multiSD.nels_dn;
  std::vector<std::vector<size_t>> C2nodes(2);

  //Check id multideterminants are in HDF5

  xmlNodePtr curTemp = cur, DetListNode = nullptr;
  curTemp = curTemp->children;
  while (curTemp != NULL) //check the basis set
  {
    std::string cname;
    getNodeName(cname, curTemp);
    if (cname == "detlist")
    {
      DetListNode = curTemp;
    }
    curTemp = curTemp->next;
  }
  XMLAttrString HDF5Path(DetListNode, "href");
  if (HDF5Path != "")
  {
    app_log() << "Found Multideterminants in H5 File" << std::endl;
    success = readDetListH5(cur, uniqueConfgs, C2nodes, CItags, multiSD.C, optimizeCI, nels);
  }
  else
    success = readDetList(cur, uniqueConfgs, C2nodes, CItags, multiSD.C, optimizeCI, nels, multiSD.CSFcoeff,
                          multiSD.DetsPerCSF, multiSD.CSFexpansion, multiSD.usingCSF);
  if (!success)
    return false;

  multiSD.C2node_up = C2nodes[0];
  multiSD.C2node_dn = C2nodes[1];
  multiSD.resize(uniqueConfgs[0].size(), uniqueConfgs[1].size());
  // alpha dets
  {
    auto& spo = multiSD.spo_up;
    spo->occup.resize(uniqueConfgs[0].size(), multiSD.nels_up);
    multiSD.dets_up.reserve(uniqueConfgs[0].size());
    for (int i = 0; i < uniqueConfgs[0].size(); i++)
    {
      int nq               = 0;
      ci_configuration& ci = uniqueConfgs[0][i];
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
        adet = new DiracDeterminantWithBackflow(targetPtcl, std::static_pointer_cast<SPOSet>(spo), 0, 0);
      }
      else
      {
        adet = new DiracDeterminant<>(std::static_pointer_cast<SPOSet>(spo), 0);
      }
      adet->set(multiSD.FirstIndex_up, multiSD.nels_up);
      multiSD.dets_up.emplace_back(adet);
    }
  }
  // beta dets
  {
    auto& spo = multiSD.spo_dn;
    spo->occup.resize(uniqueConfgs[1].size(), multiSD.nels_dn);
    multiSD.dets_dn.reserve(uniqueConfgs[1].size());
    for (int i = 0; i < uniqueConfgs[1].size(); i++)
    {
      int nq               = 0;
      ci_configuration& ci = uniqueConfgs[1][i];
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
        adet = new DiracDeterminantWithBackflow(targetPtcl, std::static_pointer_cast<SPOSet>(spo), 0, 0);
      }
      else
      {
        adet = new DiracDeterminant<>(std::static_pointer_cast<SPOSet>(spo), 0);
      }
      adet->set(multiSD.FirstIndex_dn, multiSD.nels_dn);
      multiSD.dets_dn.emplace_back(adet);
    }
  }
  if (multiSD.CSFcoeff.size() == 1 || multiSD.C.size() == 1)
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
      if (multiSD.usingCSF)
        for (int i = 1; i < multiSD.CSFcoeff.size(); i++)
          multiSD.CSFcoeff[i] = 0;
      else
        for (int i = 1; i < multiSD.C.size(); i++)
          multiSD.C[i] = 0;
      app_log() << "CI coefficients are reset. \n";
    }
    multiSD.Optimizable = true;
    if (multiSD.usingCSF)
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->CSFcoeff[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD.CSFcoeff.size(); i++)
      {
        //std::stringstream sstr;
        //sstr << "CIcoeff" << "_" << i;
        multiSD.myVars.insert(CItags[i], multiSD.CSFcoeff[i], true, optimize::LINEAR_P);
      }
    }
    else
    {
      //          multiSD->myVars.insert(CItags[0],multiSD->C[0],false,optimize::LINEAR_P);
      for (int i = 1; i < multiSD.C.size(); i++)
      {
        //std::stringstream sstr;
        //sstr << "CIcoeff" << "_" << i;
        multiSD.myVars.insert(CItags[i], multiSD.C[i], true, optimize::LINEAR_P);
      }
    }
  }
  else
  {
    app_log() << "CI coefficients are not optimizable. \n";
    multiSD.Optimizable = false;
  }
  return success;
}

bool SlaterDetBuilder::readDetList(xmlNodePtr cur,
                                   std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                                   std::vector<std::vector<size_t>>& C2nodes,
                                   std::vector<std::string>& CItags,
                                   std::vector<ValueType>& coeff,
                                   bool& optimizeCI,
                                   std::vector<int>& nptcls,
                                   std::vector<ValueType>& CSFcoeff,
                                   std::vector<size_t>& DetsPerCSF,
                                   std::vector<RealType>& CSFexpansion,
                                   bool& usingCSF)
{
  bool success = true;

  const int nGroups = uniqueConfgs.size();
  for (int grp = 0; grp < nGroups; grp++)
  {
    uniqueConfgs[grp].clear();
    C2nodes[grp].clear();
  }
  CItags.clear();
  coeff.clear();
  CSFcoeff.clear();
  DetsPerCSF.clear();
  CSFexpansion.clear();
  std::vector<std::vector<ci_configuration>> confgLists(nGroups);
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

  std::vector<size_t> NCs(nGroups);
  std::vector<size_t> NEs(nGroups);
  size_t nstates        = 0;
  size_t ndets          = 0;
  size_t count          = 0;
  size_t cnt0           = 0;
  std::string Dettype   = "DETS";
  std::string CSFChoice = "qchem_coeff";
  OhmmsAttributeSet spoAttrib;
  for (int grp = 0; grp < nGroups; grp++)
  {
    spoAttrib.add(NCs[grp], "nc" + std::to_string(grp));
    spoAttrib.add(NEs[grp], "ne" + std::to_string(grp));
  }
  if (nGroups == 2)
  { //for legacy
    spoAttrib.add(NCs[0], "nca");
    spoAttrib.add(NCs[1], "ncb");
    spoAttrib.add(NEs[0], "nea");
    spoAttrib.add(NEs[1], "neb");
  }
  spoAttrib.add(ndets, "size");
  spoAttrib.add(nstates, "nstates");
  spoAttrib.add(Dettype, "type");
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

  for (int grp = 0; grp < nGroups; grp++)
  {
    if (NEs[grp] == 0)
      NEs[grp] = nptcls[grp] - NCs[grp];
    else if (NEs[grp] != nptcls[grp] - NCs[grp])
      throw std::runtime_error("ne is not equal to n - nc for group " + std::to_string(grp));
  }

  cur = DetListNode->children;
  std::vector<ci_configuration> dummyCs(nGroups);
  for (int grp = 0; grp < nGroups; grp++)
  {
    dummyCs[grp].occup.resize(NCs[grp] + nstates, false);
    for (size_t i = 0; i < NCs[grp] + NEs[grp]; i++)
      dummyCs[grp].occup[i] = true;
  }
  RealType sumsq_qc = 0.0;
  RealType sumsq    = 0.0;
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
        RealType ci_real = 0.0, ci_imag = 0.0;
        confAttrib.add(ci_real, "coeff");
        confAttrib.add(ci_real, "coeff_real");
        confAttrib.add(ci_imag, "coeff_imag");
        confAttrib.add(qc_ci, "qchem_coeff");
        confAttrib.add(tag, "id");
        confAttrib.add(OccString, "occ");
        confAttrib.add(exctLvl, "exctLvl");
        confAttrib.put(cur);

        if (qc_ci == 0.0)
          qc_ci = ci_real;

        ValueType ci;
#ifdef QMC_COMPLEX
        ci = ValueType(ci_real, ci_imag);
#else
        if (ci_imag != RealType(0))
          myComm->barrier_and_abort(
              "SlaterDetBuilder::readDetList. Build with QMC_COMPLEX if using complex CI expansion coefficients.");
        ci = ci_real;
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
          ci = 0.0;
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
            std::vector<std::string> occs(nGroups);
            std::string tag0;
            RealType coef = 0.0;
            OhmmsAttributeSet detAttrib;
            detAttrib.add(tag0, "id");
            detAttrib.add(coef, "coeff");
            for (int grp = 0; grp < nGroups; grp++)
              detAttrib.add(occs[grp], "occ" + std::to_string(grp));
            if (nGroups == 2)
            { //for legacy
              detAttrib.add(occs[0], "alpha");
              detAttrib.add(occs[1], "beta");
            }
            detAttrib.put(csf);
            for (int grp = 0; grp < nGroups; grp++)
            {
              size_t nq = 0;
              if (occs[grp].size() < nstates)
              {
                std::cerr << "occ" << grp << ": " << occs[grp] << std::endl;
                APP_ABORT("Found incorrect group" + std::to_string(grp) + " determinant label. size < nc+nstates");
              }
              for (size_t i = 0; i < nstates; i++)
              {
                if (occs[grp][i] != '0' && occs[grp][i] != '1')
                {
                  std::cerr << occs[grp] << std::endl;
                  APP_ABORT("Found incorrect determinant label.");
                }
                if (occs[grp][i] == '1')
                  nq++;
              }
              if (nq != NEs[grp])
              {
                std::cerr << "occ" << grp << ": " << occs[grp] << std::endl;
                APP_ABORT("Found incorrect group" + std::to_string(grp) + " determinant label. nocc != nc+ne");
              }
            }
            DetsPerCSF.back()++;
            CSFexpansion.push_back(coef);
            coeff.push_back(coef * ci);
            for (int grp = 0; grp < nGroups; grp++)
            {
              confgLists[grp].push_back(dummyCs[grp]);
              for (size_t i = 0; i < NCs[grp]; i++)
                confgLists[grp].back().occup[i] = true;
              for (size_t i = NCs[grp]; i < NCs[grp] + nstates; i++)
                confgLists[grp].back().occup[i] = (occs[grp][i - NCs[grp]] == '1');
            }
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

    for (int grp = 0; grp < nGroups; grp++)
      C2nodes[grp].resize(coeff.size());
    app_log() << "Found " << coeff.size() << " terms in the MSD expansion.\n";
    RealType sumsq = 0.0;
    for (size_t i = 0; i < coeff.size(); i++)
      sumsq += std::abs(coeff[i] * coeff[i]);
    app_log() << "Norm of ci vector (sum of ci^2): " << sumsq << std::endl;
    app_log() << "Norm of qchem ci vector (sum of qchem_ci^2): " << sumsq_qc << std::endl;
    for (int grp = 0; grp < nGroups; grp++)
    {
      for (size_t i = 0; i < confgLists[grp].size(); i++)
      {
        bool found = false;
        size_t k   = -1;
        for (size_t j = 0; j < uniqueConfgs[grp].size(); j++)
        {
          if (confgLists[grp][i] == uniqueConfgs[grp][j])
          {
            found = true;
            k     = j;
            break;
          }
        }
        if (found)
        {
          C2nodes[grp][i] = k;
        }
        else
        {
          uniqueConfgs[grp].push_back(confgLists[grp][i]);
          C2nodes[grp][i] = uniqueConfgs[grp].size() - 1;
        }
      }
    }
  }
  else
  {
    app_log() << "Reading CI expansion." << std::endl;

    std::vector<int> cnts(nGroups, 0);
    std::vector<std::unordered_map<std::string, int>> MyMaps(nGroups);
    while (cur != NULL) //check the basis set
    {
      getNodeName(cname, cur);
      if (cname == "configuration" || cname == "ci")
      {
        RealType qc_ci = 0.0;
        std::vector<std::string> occs(nGroups);
        std::string tag;
        OhmmsAttributeSet confAttrib;
        RealType ci_real = 0.0, ci_imag = 0.0;
        confAttrib.add(ci_real, "coeff");
        confAttrib.add(ci_real, "coeff_real");
        confAttrib.add(ci_imag, "coeff_imag");
        confAttrib.add(qc_ci, "qchem_coeff");
        for (int grp = 0; grp < nGroups; grp++)
          confAttrib.add(occs[grp], "occ" + std::to_string(grp));
        if (nGroups == 2)
        { //legacy
          confAttrib.add(occs[0], "alpha");
          confAttrib.add(occs[1], "beta");
        }
        confAttrib.add(tag, "id");
        confAttrib.put(cur);

        ValueType ci;
#ifdef QMC_COMPLEX
        ci = ValueType(ci_real, ci_imag);
#else
        if (ci_imag != RealType(0))
          myComm->barrier_and_abort(
              "SlaterDetBuilder::readDetList. Build with QMC_COMPLEX if using complex CI expansion coefficients.");
        ci = ci_real;
#endif

        //Will always loop through the whole determinant set as no assumption on the order of the determinant is made
        if (std::abs(ci) < cutoff)
        {
          cur = cur->next;
          continue;
        }

        for (int grp = 0; grp < nGroups; grp++)
        {
          for (size_t i = 0; i < nstates; i++)
          {
            if (occs[grp][i] != '0' && occs[grp][i] != '1')
            {
              std::cerr << occs[grp] << std::endl;
              APP_ABORT("Found incorrect determinant label for group " + std::to_string(grp));
            }
          }
          if (occs[grp].size() < nstates)
          {
            std::cerr << "occ" << grp << ": " << occs[grp] << std::endl;
            APP_ABORT("Found incorrect group" + std::to_string(grp) + " determinant label. size < nc+nstates");
          }
        }

        coeff.push_back(ci);
        CItags.push_back(tag);

        for (int grp = 0; grp < nGroups; grp++)
        {
          std::unordered_map<std::string, int>::const_iterator got = MyMaps[grp].find(occs[grp]);
          if (got == MyMaps[grp].end())
          {
            uniqueConfgs[grp].push_back(dummyCs[grp]);
            uniqueConfgs[grp].back().add_occupation(occs[grp]);
            C2nodes[grp].push_back(cnts[grp]);
            MyMaps[grp].insert(std::pair<std::string, int>(occs[grp], cnts[grp]));
            cnts[grp]++;
          }
          else
          {
            C2nodes[grp].push_back(got->second);
          }
        }

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

  for (int grp = 0; grp < nGroups; grp++)
    app_log() << "Found " << uniqueConfgs[grp].size() << " unique group" << grp << " determinants.\n";

  return success;
}

bool SlaterDetBuilder::readDetListH5(xmlNodePtr cur,
                                     std::vector<std::vector<ci_configuration>>& uniqueConfgs,
                                     std::vector<std::vector<size_t>>& C2nodes,
                                     std::vector<std::string>& CItags,
                                     std::vector<ValueType>& coeff,
                                     bool& optimizeCI,
                                     std::vector<int>& nptcls)
{
  bool success = true;
  int extlevel(0);
  const int nGroups = uniqueConfgs.size();
  for (int grp = 0; grp < nGroups; grp++)
  {
    uniqueConfgs[grp].clear();
    C2nodes[grp].clear();
  }
  CItags.clear();
  coeff.clear();
  std::string CICoeffH5path("");
  std::vector<std::vector<ci_configuration>> confgLists(nGroups);
  std::vector<ValueType> CIcoeff;
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

  app_log() << "  H5 code path implicitly assumes NC0 = NC1 = ...  = 0" << std::endl;
  for (int grp = 0; grp < nGroups; grp++)
    app_log() << "NE" << grp << " = " << nptcls[grp] << ", ";
  app_log() << std::endl;
  size_t nstates = 0;
  size_t ndets   = 0;
  size_t H5_ndets, H5_nstates;
  /// 64 bit fixed width integer
  const unsigned bit_kind = 64;
  static_assert(bit_kind == sizeof(int64_t) * 8, "Must be 64 bit fixed width integer");
  /// the number of 64 bit integers which represent the binary string for occupation
  int N_int;
  std::string Dettype = "DETS";
  ValueType sumsq     = 0.0;
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(ndets, "size");
  spoAttrib.add(Dettype, "type");
  spoAttrib.add(nstates, "nstates");
  spoAttrib.add(extlevel, "ext_level");
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
  if (nstates == 0)
    nstates = H5_nstates;
  else if (nstates != H5_nstates)
  {
    std::cerr << "Number of states/orbitals in H5 file (" << H5_nstates
              << ") different from number of states/orbitals in XML (" << nstates << ")" << std::endl;
    abort();
  }

  hin.read(N_int, "Nbits");
  CIcoeff.resize(ndets);
  ConfigTag.resize(ndets);

  readCoeffs(hin, CIcoeff, ndets, extlevel);

  ///IF OPTIMIZED COEFFICIENTS ARE PRESENT IN opt_coeffs Path
  ///THEY ARE READ FROM DIFFERENT HDF5 the replace the previous coeff
  ///It is important to still read all old coeffs and only replace the optimized ones
  ///in order to keep coherence with the cutoff on the number of determinants
  ///REMEMBER!! FIRST COEFF IS FIXED. THEREFORE WE DO NOT REPLACE IT!!!
  if (CICoeffH5path != "")
  {
    int OptCiSize = 0;
    std::vector<ValueType> CIcoeffopt;
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

    readCoeffs(coeffin, CIcoeffopt, ndets, extlevel);

    coeffin.close();

    for (int i = 0; i < OptCiSize; i++)
      CIcoeff[i + 1] = CIcoeffopt[i];

    app_log() << "The first " << OptCiSize
              << " Optimized coefficients were substituted to the original set of coefficients." << std::endl;
  }

  std::vector<Matrix<int64_t>> temps;
  for (int grp = 0; grp < nGroups; grp++)
  {
    Matrix<int64_t> tmp(ndets, N_int);
    temps.push_back(tmp);
    std::string ci_str;

    if (!hin.readEntry(temps[grp], "CI_" + std::to_string(grp)))
    {
      //for backwards compatibility
      if (grp == 0)
        hin.read(temps[grp], "CI_Alpha");
      else if (grp == 1)
        hin.read(temps[grp], "CI_Beta");
      else
        APP_ABORT("Unknown HDF5 CI format");
    }
  }

  std::vector<std::string> MyCIs(nGroups);
  std::vector<ci_configuration> dummyCs(nGroups);
  for (int grp = 0; grp < nGroups; grp++)
  {
    MyCIs[grp].resize(nstates);
    dummyCs[grp].occup.resize(nstates, false);
    for (size_t i = 0; i < nstates; i++)
      dummyCs[grp].occup[i] = true;
  }

  hin.close();
  app_log() << " Done reading " << ndets << " CIs from H5!" << std::endl;

  std::vector<int> cnts(nGroups, 0);
  std::vector<std::unordered_map<std::string, int>> MyMaps(nGroups);

  app_log() << " Sorting unique CIs" << std::endl;
  ///This loop will find all unique Determinants in and store them "unsorted" in a new container uniqueConfg_up
  ///and uniqueConfg_dn. The sorting is not done here
  for (int ni = 0; ni < ndets; ni++)
  {
    if (std::abs(CIcoeff[ni]) < cutoff)
      continue;

    int j = 0;
    for (int k = 0; k < N_int; k++)
    {
      std::vector<std::bitset<bit_kind>> a2s(nGroups);
      for (int grp = 0; grp < nGroups; grp++)
      {
        int64_t a = temps[grp][ni][k];
        a2s[grp]  = a;
      }

      for (int i = 0; i < bit_kind; i++)
      {
        if (j < nstates)
        {
          for (int grp = 0; grp < nGroups; grp++)
            MyCIs[grp][j] = a2s[grp][i] ? '1' : '0';
          j++;
        }
      }
    }
    coeff.push_back(CIcoeff[ni]);
    std::ostringstream h5tag;
    h5tag << "CIcoeff_" << ni;
    CItags.push_back(h5tag.str());

    for (int grp = 0; grp < nGroups; grp++)
    {
      std::unordered_map<std::string, int>::const_iterator got = MyMaps[grp].find(MyCIs[grp]);

      if (got == MyMaps[grp].end())
      {
        uniqueConfgs[grp].push_back(dummyCs[grp]);
        uniqueConfgs[grp].back().add_occupation(MyCIs[grp]);
        C2nodes[grp].push_back(cnts[grp]);
        MyMaps[grp].insert(std::pair<std::string, int>(MyCIs[grp], cnts[grp]));
        cnts[grp]++;
      }
      else
      {
        C2nodes[grp].push_back(got->second);
      }
    }
    sumsq += CIcoeff[ni] * CIcoeff[ni];
  }

  app_log() << " Done Sorting unique CIs" << std::endl;
  app_log() << "Found " << coeff.size() << " terms in the MSD expansion.\n";
  app_log() << "Norm of ci vector (sum of ci^2): " << sumsq << std::endl;

  for (int grp = 0; grp < nGroups; grp++)
    app_log() << "Found " << uniqueConfgs[grp].size() << " unique group" << grp << " determinants.\n";

  return success;
}
} // namespace qmcplusplus

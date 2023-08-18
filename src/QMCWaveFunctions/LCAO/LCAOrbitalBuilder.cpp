//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include "LCAOrbitalBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "MultiQuinticSpline1D.h"
#include "Numerics/SoaCartesianTensor.h"
#include "Numerics/SoaSphericalTensor.h"
#include "SoaAtomicBasisSet.h"
#include "SoaLocalizedBasisSet.h"
#include "LCAOrbitalSet.h"
#include "AOBasisBuilder.h"
#include "MultiFunctorAdapter.h"
#if !defined(QMC_COMPLEX)
#include "LCAOrbitalSetWithCorrection.h"
#include "CuspCorrectionConstruction.h"
#endif
#include "hdf/hdf_archive.h"
#include "Message/CommOperators.h"
#include "Utilities/ProgressReportEngine.h"
#include "CPU/math.hpp"

#include <array>

namespace qmcplusplus
{
/** traits for a localized basis set; used by createBasisSet
   *
   * T radial function value type
   * ORBT orbital value type, can be complex
   * ROT {0=numuerica;, 1=gto; 2=sto}
   * SH {0=cartesian, 1=spherical}
   * If too confusing, inroduce enumeration.
   */
template<typename T, typename ORBT, int ROT, int SH>
struct ao_traits
{};

/** specialization for numerical-cartesian AO */
template<typename T, typename ORBT>
struct ao_traits<T, ORBT, 0, 0>
{
  using radial_type  = MultiQuinticSpline1D<T>;
  using angular_type = SoaCartesianTensor<T>;
  using ao_type      = SoaAtomicBasisSet<radial_type, angular_type>;
  using basis_type   = SoaLocalizedBasisSet<ao_type, ORBT>;
};

/** specialization for numerical-spherical AO */
template<typename T, typename ORBT>
struct ao_traits<T, ORBT, 0, 1>
{
  using radial_type  = MultiQuinticSpline1D<T>;
  using angular_type = SoaSphericalTensor<T>;
  using ao_type      = SoaAtomicBasisSet<radial_type, angular_type>;
  using basis_type   = SoaLocalizedBasisSet<ao_type, ORBT>;
};

/** specialization for GTO-cartesian AO */
template<typename T, typename ORBT>
struct ao_traits<T, ORBT, 1, 0>
{
  using radial_type  = MultiFunctorAdapter<GaussianCombo<T>>;
  using angular_type = SoaCartesianTensor<T>;
  using ao_type      = SoaAtomicBasisSet<radial_type, angular_type>;
  using basis_type   = SoaLocalizedBasisSet<ao_type, ORBT>;
};

/** specialization for GTO-cartesian AO */
template<typename T, typename ORBT>
struct ao_traits<T, ORBT, 1, 1>
{
  using radial_type  = MultiFunctorAdapter<GaussianCombo<T>>;
  using angular_type = SoaSphericalTensor<T>;
  using ao_type      = SoaAtomicBasisSet<radial_type, angular_type>;
  using basis_type   = SoaLocalizedBasisSet<ao_type, ORBT>;
};

/** specialization for STO-spherical AO */
template<typename T, typename ORBT>
struct ao_traits<T, ORBT, 2, 1>
{
  using radial_type  = MultiFunctorAdapter<SlaterCombo<T>>;
  using angular_type = SoaSphericalTensor<T>;
  using ao_type      = SoaAtomicBasisSet<radial_type, angular_type>;
  using basis_type   = SoaLocalizedBasisSet<ao_type, ORBT>;
};


inline bool is_same(const xmlChar* a, const char* b) { return !strcmp((const char*)a, b); }

using BasisSet_t = LCAOrbitalSet::basis_type;

LCAOrbitalBuilder::LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("LCAO", comm),
      targetPtcl(els),
      sourcePtcl(ions),
      h5_path(""),
      SuperTwist(0.0),
      doCuspCorrection(false)
{
  ClassName = "LCAOrbitalBuilder";
  ReportEngine PRE(ClassName, "createBasisSet");

  std::string cuspC("no"); // cusp correction
  OhmmsAttributeSet aAttrib;
  aAttrib.add(cuspC, "cuspCorrection");
  aAttrib.add(h5_path, "href");
  aAttrib.add(PBCImages, "PBCimages");
  aAttrib.add(SuperTwist, "twist");
  aAttrib.put(cur);

  if (cuspC == "yes")
    doCuspCorrection = true;
  //Evaluate the Phase factor. Equals 1 for OBC.
  EvalPeriodicImagePhaseFactors(SuperTwist, PeriodicImagePhaseFactors);

  // no need to wait but load the basis set
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "basisset")
    {
      std::string basisset_name_input(getXMLAttributeValue(element, "name"));
      std::string basisset_name(basisset_name_input.empty() ? "LCAOBSet" : basisset_name_input);
      if (basisset_map_.find(basisset_name) != basisset_map_.end())
      {
        std::ostringstream err_msg;
        err_msg << "Cannot create basisset " << basisset_name << " which already exists." << std::endl;
        throw std::runtime_error(err_msg.str());
      }
      if (h5_path != "")
        basisset_map_[basisset_name] = loadBasisSetFromH5(element);
      else
        basisset_map_[basisset_name] = loadBasisSetFromXML(element, cur);
    }
  });

  // deprecated h5 basis set handling when basisset element is missing
  if (basisset_map_.size() == 0 && h5_path != "")
  {
    app_warning() << "!!!!!!! Deprecated input style: missing basisset element. "
                  << "LCAO needs an explicit basisset XML element. "
                  << "Fallback on loading an implicit one." << std::endl;
    basisset_map_["LCAOBSet"] = loadBasisSetFromH5(cur);
  }

  if (basisset_map_.size() == 0)
    throw std::runtime_error("No basisset found in the XML input!");
}

LCAOrbitalBuilder::~LCAOrbitalBuilder()
{
  //properly cleanup
}

int LCAOrbitalBuilder::determineRadialOrbType(xmlNodePtr cur) const
{
  std::string keyOpt;
  std::string transformOpt;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(keyOpt, "keyword");
  aAttrib.add(keyOpt, "key");
  aAttrib.add(transformOpt, "transform");
  aAttrib.put(cur);

  int radialOrbType = -1;
  if (transformOpt == "yes" || keyOpt == "NMO")
    radialOrbType = 0;
  else
  {
    if (keyOpt == "GTO")
      radialOrbType = 1;
    if (keyOpt == "STO")
      radialOrbType = 2;
  }
  return radialOrbType;
}

std::unique_ptr<BasisSet_t> LCAOrbitalBuilder::loadBasisSetFromXML(xmlNodePtr cur, xmlNodePtr parent)
{
  ReportEngine PRE(ClassName, "loadBasisSetFromXML(xmlNodePtr)");
  int ylm = -1;
  {
    xmlNodePtr cur1 = cur->xmlChildrenNode;
    while (cur1 != NULL && ylm < 0)
    {
      if (is_same(cur1->name, "atomicBasisSet"))
      {
        std::string sph;
        OhmmsAttributeSet att;
        att.add(sph, "angular");
        att.put(cur1);
        ylm = (sph == "cartesian") ? 0 : 1;
      }
      cur1 = cur1->next;
    }
  }

  if (ylm < 0)
    PRE.error("Missing angular attribute of atomicBasisSet.", true);

  int radialOrbType = determineRadialOrbType(cur);
  if (radialOrbType < 0)
  {
    app_warning() << "Radial orbital type cannot be determined based on the attributes of basisset line. "
                  << "Trying the parent element." << std::endl;
    radialOrbType = determineRadialOrbType(parent);
  }

  if (radialOrbType < 0)
    PRE.error("Unknown radial function for LCAO orbitals. Specify keyword=\"NMO/GTO/STO\" .", true);

  BasisSet_t* myBasisSet = nullptr;
  /** process atomicBasisSet per ion species */
  switch (radialOrbType)
  {
  case (0): //numerical
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiQuintic," << ylm << ">" << std::endl;
    if (ylm)
      myBasisSet = createBasisSet<0, 1>(cur);
    else
      myBasisSet = createBasisSet<0, 0>(cur);
    break;
  case (1): //gto
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiGTO," << ylm << ">" << std::endl;
    if (ylm)
      myBasisSet = createBasisSet<1, 1>(cur);
    else
      myBasisSet = createBasisSet<1, 0>(cur);
    break;
  case (2): //sto
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiSTO," << ylm << ">" << std::endl;
    myBasisSet = createBasisSet<2, 1>(cur);
    break;
  default:
    PRE.error("Cannot construct SoaAtomicBasisSet<ROT,YLM>.", true);
    break;
  }

  return std::unique_ptr<BasisSet_t>(myBasisSet);
}

std::unique_ptr<BasisSet_t> LCAOrbitalBuilder::loadBasisSetFromH5(xmlNodePtr parent)
{
  ReportEngine PRE(ClassName, "loadBasisSetFromH5()");

  hdf_archive hin(myComm);
  int ylm = -1;
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      PRE.error("Could not open H5 file", true);

    hin.push("basisset", false);

    std::string sph;
    std::string ElemID0 = "atomicBasisSet0";

    hin.push(ElemID0.c_str(), false);

    if (!hin.readEntry(sph, "angular"))
      PRE.error("Could not find name of  basisset group in H5; Probably Corrupt H5 file", true);
    ylm = (sph == "cartesian") ? 0 : 1;
    hin.close();
  }

  myComm->bcast(ylm);
  if (ylm < 0)
    PRE.error("Missing angular attribute of atomicBasisSet.", true);

  int radialOrbType = determineRadialOrbType(parent);
  if (radialOrbType < 0)
    PRE.error("Unknown radial function for LCAO orbitals. Specify keyword=\"NMO/GTO/STO\" .", true);

  BasisSet_t* myBasisSet = nullptr;
  /** process atomicBasisSet per ion species */
  switch (radialOrbType)
  {
  case (0): //numerical
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiQuintic," << ylm << ">" << std::endl;
    if (ylm)
      myBasisSet = createBasisSetH5<0, 1>();
    else
      myBasisSet = createBasisSetH5<0, 0>();
    break;
  case (1): //gto
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiGTO," << ylm << ">" << std::endl;
    if (ylm)
      myBasisSet = createBasisSetH5<1, 1>();
    else
      myBasisSet = createBasisSetH5<1, 0>();
    break;
  case (2): //sto
    app_log() << "  LCAO: SoaAtomicBasisSet<MultiSTO," << ylm << ">" << std::endl;
    myBasisSet = createBasisSetH5<2, 1>();
    break;
  default:
    PRE.error("Cannot construct SoaAtomicBasisSet<ROT,YLM>.", true);
    break;
  }
  return std::unique_ptr<BasisSet_t>(myBasisSet);
}


template<int I, int J>
LCAOrbitalBuilder::BasisSet_t* LCAOrbitalBuilder::createBasisSet(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createBasisSet(xmlNodePtr)");

  using ao_type    = typename ao_traits<RealType, ValueType, I, J>::ao_type;
  using basis_type = typename ao_traits<RealType, ValueType, I, J>::basis_type;

  basis_type* mBasisSet = new basis_type(sourcePtcl, targetPtcl);

  //list of built centers
  std::vector<std::string> ao_built_centers;

  /** process atomicBasisSet per ion species */
  cur = cur->xmlChildrenNode;
  while (cur != NULL) //loop over unique ioons
  {
    std::string cname((const char*)(cur->name));

    if (cname == "atomicBasisSet")
    {
      std::string elementType;
      std::string sph;
      OhmmsAttributeSet att;
      att.add(elementType, "elementType");
      att.put(cur);

      if (elementType.empty())
        PRE.error("Missing elementType attribute of atomicBasisSet.", true);

      auto it = std::find(ao_built_centers.begin(), ao_built_centers.end(), elementType);
      if (it == ao_built_centers.end())
      {
        AOBasisBuilder<ao_type> any(elementType, myComm);
        any.put(cur);
        auto aoBasis = any.createAOSet(cur);
        if (aoBasis)
        {
          //add the new atomic basis to the basis set
          int activeCenter = sourcePtcl.getSpeciesSet().findSpecies(elementType);
          mBasisSet->add(activeCenter, std::move(aoBasis));
        }
        ao_built_centers.push_back(elementType);
      }
    }
    cur = cur->next;
  } // done with basis set
  mBasisSet->setBasisSetSize(-1);
  mBasisSet->setPBCParams(PBCImages, SuperTwist, PeriodicImagePhaseFactors);
  return mBasisSet;
}


template<int I, int J>
LCAOrbitalBuilder::BasisSet_t* LCAOrbitalBuilder::createBasisSetH5()
{
  ReportEngine PRE(ClassName, "createBasisSetH5(xmlNodePtr)");

  using ao_type    = typename ao_traits<RealType, ValueType, I, J>::ao_type;
  using basis_type = typename ao_traits<RealType, ValueType, I, J>::basis_type;

  basis_type* mBasisSet = new basis_type(sourcePtcl, targetPtcl);

  //list of built centers
  std::vector<std::string> ao_built_centers;

  int Nb_Elements(0);
  std::string basiset_name;

  /** process atomicBasisSet per ion species */
  app_log() << "Reading BasisSet from HDF5 file:" << h5_path << std::endl;

  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      PRE.error("Could not open H5 file", true);

    hin.push("basisset", false);

    hin.read(Nb_Elements, "NbElements");
  }

  myComm->bcast(Nb_Elements);
  if (Nb_Elements < 1)
    PRE.error("Missing elementType attribute of atomicBasisSet.", true);

  for (int i = 0; i < Nb_Elements; i++)
  {
    std::string elementType, dataset;
    std::stringstream tempElem;
    std::string ElemID0 = "atomicBasisSet", ElemType;
    tempElem << ElemID0 << i;
    ElemType = tempElem.str();

    if (myComm->rank() == 0)
    {
      hin.push(ElemType.c_str(), false);

      if (!hin.readEntry(basiset_name, "name"))
        PRE.error("Could not find name of  basisset group in H5; Probably Corrupt H5 file", true);
      if (!hin.readEntry(elementType, "elementType"))
        PRE.error("Could not read elementType in H5; Probably Corrupt H5 file", true);
    }
    myComm->bcast(basiset_name);
    myComm->bcast(elementType);

    auto it = std::find(ao_built_centers.begin(), ao_built_centers.end(), elementType);
    if (it == ao_built_centers.end())
    {
      AOBasisBuilder<ao_type> any(elementType, myComm);
      any.putH5(hin);
      auto aoBasis = any.createAOSetH5(hin);
      if (aoBasis)
      {
        //add the new atomic basis to the basis set
        int activeCenter = sourcePtcl.getSpeciesSet().findSpecies(elementType);
        mBasisSet->add(activeCenter, std::move(aoBasis));
      }
      ao_built_centers.push_back(elementType);
    }

    if (myComm->rank() == 0)
      hin.pop();
  }

  if (myComm->rank() == 0)
  {
    hin.pop();
    hin.close();
  }
  mBasisSet->setBasisSetSize(-1);
  mBasisSet->setPBCParams(PBCImages, SuperTwist, PeriodicImagePhaseFactors);
  return mBasisSet;
}


std::unique_ptr<SPOSet> LCAOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createSPO(xmlNodePtr)");
  std::string spo_name(""), cusp_file(""), optimize("no");
  std::string basisset_name("LCAOBSet");
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(spo_name, "name");
  spoAttrib.add(spo_name, "id");
  spoAttrib.add(cusp_file, "cuspInfo");
  spoAttrib.add(basisset_name, "basisset");
  spoAttrib.put(cur);

  std::unique_ptr<BasisSet_t> myBasisSet;
  if (basisset_map_.find(basisset_name) == basisset_map_.end())
    myComm->barrier_and_abort("basisset \"" + basisset_name + "\" cannot be found\n");
  else
    myBasisSet.reset(basisset_map_[basisset_name]->makeClone());

  std::unique_ptr<SPOSet> sposet;
  if (doCuspCorrection)
  {
#if defined(QMC_COMPLEX)
    myComm->barrier_and_abort(
        "LCAOrbitalBuilder::createSPOSetFromXML cusp correction is not supported on complex LCAO.");
#else
    app_summary() << "        Using cusp correction." << std::endl;
    auto lcwc = std::make_unique<LCAOrbitalSetWithCorrection>(spo_name, sourcePtcl, targetPtcl, std::move(myBasisSet));
    loadMO(lcwc->lcao, cur);
    lcwc->setOrbitalSetSize(lcwc->lcao.getOrbitalSetSize());
    sposet = std::move(lcwc);
#endif
  }
  else
  {
    auto lcos = std::make_unique<LCAOrbitalSet>(spo_name, std::move(myBasisSet));
    loadMO(*lcos, cur);
    sposet = std::move(lcos);
  }

#if !defined(QMC_COMPLEX)
  if (doCuspCorrection)
  {
    // Create a temporary particle set to use for cusp initialization.
    // The particle coordinates left at the end are unsuitable for further computations.
    // The coordinates get set to nuclear positions, which leads to zero e-N distance,
    // which causes a NaN in SoaAtomicBasisSet.h
    // This problem only appears when the electron positions are specified in the input.
    // The random particle placement step executes after this part of the code, overwriting
    // the leftover positions from the cusp initialization.
    ParticleSet tmp_targetPtcl(targetPtcl);

    const int num_centers = sourcePtcl.getTotalNum();
    auto& lcwc            = dynamic_cast<LCAOrbitalSetWithCorrection&>(*sposet);

    const int orbital_set_size = lcwc.getOrbitalSetSize();
    Matrix<CuspCorrectionParameters> info(num_centers, orbital_set_size);

    // set a default file name if not given
    if (cusp_file.empty())
      cusp_file = spo_name + ".cuspInfo.xml";

    bool file_exists(myComm->rank() == 0 && std::ifstream(cusp_file).good());
    myComm->bcast(file_exists);
    app_log() << "  Cusp correction file " << cusp_file << (file_exists ? " exits." : " doesn't exist.") << std::endl;

    // validate file if it exists
    if (file_exists)
    {
      bool valid = 0;
      if (myComm->rank() == 0)
        valid = readCuspInfo(cusp_file, spo_name, orbital_set_size, info);
      myComm->bcast(valid);
      if (!valid)
        myComm->barrier_and_abort("Invalid cusp correction file " + cusp_file);
#ifdef HAVE_MPI
      for (int orb_idx = 0; orb_idx < orbital_set_size; orb_idx++)
        for (int center_idx = 0; center_idx < num_centers; center_idx++)
          broadcastCuspInfo(info(center_idx, orb_idx), *myComm, 0);
#endif
    }
    else
    {
      generateCuspInfo(info, tmp_targetPtcl, sourcePtcl, lcwc.lcao, spo_name, *myComm);
      if (myComm->rank() == 0)
        saveCusp(cusp_file, info, spo_name);
    }

    applyCuspCorrection(info, tmp_targetPtcl, sourcePtcl, lcwc.lcao, lcwc.cusp, spo_name);
  }
#endif

  return sposet;
}


/** Parse the xml file for information on the Dirac determinants.
   *@param cur the current xmlNode
   */
bool LCAOrbitalBuilder::loadMO(LCAOrbitalSet& spo, xmlNodePtr cur)
{
#undef FunctionName
#define FunctionName                                      \
  printf("Calling FunctionName from %s\n", __FUNCTION__); \
  FunctionNameReal
  //Check if HDF5 present
  ReportEngine PRE("LCAOrbitalBuilder", "put(xmlNodePtr)");

  //initialize the number of orbital by the basis set size
  int norb = spo.getBasisSetSize();
  std::string debugc("no");
  double orbital_mix_magnitude = 0.0;
  bool PBC                     = false;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb, "orbitals");
  aAttrib.add(norb, "size");
  aAttrib.add(debugc, "debug");
  aAttrib.add(orbital_mix_magnitude, "orbital_mix_magnitude");
  aAttrib.put(cur);
  xmlNodePtr occ_ptr   = NULL;
  xmlNodePtr coeff_ptr = NULL;
  cur                  = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "occupation")
    {
      occ_ptr = cur;
    }
    else if (cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var")
    {
      coeff_ptr = cur;
    }
    cur = cur->next;
  }
  if (coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << std::endl;
    return true;
  }
  spo.setOrbitalSetSize(norb);
  bool success = putOccupation(spo, occ_ptr);
  if (h5_path == "")
    success = putFromXML(spo, coeff_ptr);
  else
  {
    hdf_archive hin(myComm);

    if (myComm->rank() == 0)
    {
      if (!hin.open(h5_path, H5F_ACC_RDONLY))
        APP_ABORT("LCAOrbitalBuilder::putFromH5 missing or incorrect path to H5 file.");

      try
      {
        hin.push("PBC", false);
        PBC = true;
      }
      catch (const std::exception& e)
      {
        app_debug() << e.what() << std::endl;
        PBC = false;
      }

      if (PBC)
        hin.read(PBC, "PBC");

      hin.close();
    }
    myComm->bcast(PBC);
    if (PBC)
      success = putPBCFromH5(spo, coeff_ptr);
    else
      success = putFromH5(spo, coeff_ptr);
  }

  // Ye: used to construct cusp correction
  //bool success2 = transformSPOSet();
  if (debugc == "yes")
  {
    app_log() << "   Single-particle orbital coefficients dims=" << spo.C->rows() << " x " << spo.C->cols()
              << std::endl;
    app_log() << *spo.C << std::endl;
  }

  return success;
}

bool LCAOrbitalBuilder::putFromXML(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr)
{
  int norbs = 0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norbs, "size");
  aAttrib.add(norbs, "orbitals");
  aAttrib.put(coeff_ptr);
  if (norbs < spo.getOrbitalSetSize())
  {
    return false;
    APP_ABORT("LCAOrbitalBuilder::putFromXML missing or incorrect size");
  }
  if (norbs)
  {
    std::vector<ValueType> Ctemp;
    int BasisSetSize = spo.getBasisSetSize();
    Ctemp.resize(norbs * BasisSetSize);
    putContent(Ctemp, coeff_ptr);
    int n = 0, i = 0;
    std::vector<ValueType>::iterator cit(Ctemp.begin());
    while (i < spo.getOrbitalSetSize())
    {
      if (Occ[n] > std::numeric_limits<RealType>::epsilon())
      {
        std::copy(cit, cit + BasisSetSize, (*spo.C)[i]);
        i++;
      }
      n++;
      cit += BasisSetSize;
    }
  }
  return true;
}

/** read data from a hdf5 file
   * @param norb number of orbitals to be initialized
   * @param coeff_ptr xmlnode for coefficients
   */
bool LCAOrbitalBuilder::putFromH5(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr)
{
  int neigs  = spo.getBasisSetSize();
  int setVal = -1;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal, "spindataset");
  aAttrib.add(neigs, "size");
  aAttrib.add(neigs, "orbitals");
  aAttrib.put(coeff_ptr);
  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      APP_ABORT("LCAOrbitalBuilder::putFromH5 missing or incorrect path to H5 file.");

    Matrix<RealType> Ctemp;
    std::array<char, 72> name;


    //This is to make sure of Backward compatibility with previous tags.
    int name_len = std::snprintf(name.data(), name.size(), "%s%d", "/Super_Twist/eigenset_", setVal);
    if (name_len < 0)
      throw std::runtime_error("Error generating name");
    std::string setname(name.data(), name_len);
    if (!hin.readEntry(Ctemp, setname))
    {
      name_len = std::snprintf(name.data(), name.size(), "%s%d", "/KPTS_0/eigenset_", setVal);
      if (name_len < 0)
        throw std::runtime_error("Error generating name");
      setname = std::string(name.data(), name_len);
      hin.read(Ctemp, setname);
    }
    hin.close();

    if (Ctemp.cols() != spo.getBasisSetSize())
    {
      std::ostringstream err_msg;
      err_msg << "Basis set size " << spo.getBasisSetSize() << " mismatched the number of MO coefficients columns "
              << Ctemp.cols() << " from h5." << std::endl;
      myComm->barrier_and_abort(err_msg.str());
    }

    int norbs = spo.getOrbitalSetSize();
    if (Ctemp.rows() < norbs)
    {
      std::ostringstream err_msg;
      err_msg << "Need " << norbs << " orbitals. Insufficient rows of MO coefficients " << Ctemp.rows() << " from h5."
              << std::endl;
      myComm->barrier_and_abort(err_msg.str());
    }

    int n = 0, i = 0;
    while (i < norbs)
    {
      if (Occ[n] > 0.0)
      {
        std::copy(Ctemp[n], Ctemp[n + 1], (*spo.C)[i]);
        i++;
      }
      n++;
    }
  }
  myComm->bcast(spo.C->data(), spo.C->size());
  return true;
}


/** read data from a hdf5 file
   * @param norb number of orbitals to be initialized
   * @param coeff_ptr xmlnode for coefficients
   */
bool LCAOrbitalBuilder::putPBCFromH5(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr)
{
  ReportEngine PRE("LCAOrbitalBuilder", "LCAOrbitalBuilder::putPBCFromH5");
  int norbs      = spo.getOrbitalSetSize();
  int neigs      = spo.getBasisSetSize();
  int setVal     = -1;
  bool IsComplex = false;
  bool MultiDet  = false;
  PosType SuperTwist(0.0);
  PosType SuperTwistH5(0.0);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal, "spindataset");
  aAttrib.add(neigs, "size");
  aAttrib.add(neigs, "orbitals");
  aAttrib.put(coeff_ptr);
  hdf_archive hin(myComm);

  xmlNodePtr curtemp = coeff_ptr;

  std::string xmlTag("determinantset");
  std::string MSDTag("sposet");
  std::string SDTag("determinant");
  std::string EndTag("qmcsystem");
  std::string curname;

  do
  {
    std::stringstream ss;
    curtemp = curtemp->parent;
    ss << curtemp->name;
    ss >> curname;
    if (curname == MSDTag)
      MultiDet = true; ///Used to know if running an MSD calculation - needed for order of Orbitals.
    if (curname == SDTag)
      MultiDet = false;

  } while ((xmlTag != curname) && (curname != EndTag));
  if (curname == EndTag)
  {
    APP_ABORT(
        "Could not find in wf file the \"sposet\" or \"determinant\" tags. Please verify input or contact developers");
  }

  aAttrib.add(SuperTwist, "twist");
  aAttrib.put(curtemp);

  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      APP_ABORT("LCAOrbitalBuilder::putFromH5 missing or incorrect path to H5 file.");
    hin.push("parameters");
    hin.read(IsComplex, "IsComplex");
    hin.pop();

    std::string setname("/Super_Twist/Coord");
    hin.read(SuperTwistH5, setname);
    if (std::abs(SuperTwistH5[0] - SuperTwist[0]) >= 1e-6 || std::abs(SuperTwistH5[1] - SuperTwist[1]) >= 1e-6 ||
        std::abs(SuperTwistH5[2] - SuperTwist[2]) >= 1e-6)
    {
      app_log() << "Super Twist in XML : " << SuperTwist[0] << "    In H5:" << SuperTwistH5[0] << std::endl;
      app_log() << "                     " << SuperTwist[1] << "          " << SuperTwistH5[1] << std::endl;
      app_log() << "                     " << SuperTwist[2] << "          " << SuperTwistH5[2] << std::endl;
      app_log() << "Diff in Coord     x :" << std::abs(SuperTwistH5[0] - SuperTwist[0]) << std::endl;
      app_log() << "                  y :" << std::abs(SuperTwistH5[1] - SuperTwist[1]) << std::endl;
      app_log() << "                  z :" << std::abs(SuperTwistH5[2] - SuperTwist[2]) << std::endl;
      APP_ABORT("Requested Super Twist in XML and Super Twist in HDF5 do not Match!!! Aborting.");
    }
    //SuperTwist=SuperTwistH5;
    Matrix<ValueType> Ctemp;
    LoadFullCoefsFromH5(hin, setVal, SuperTwist, Ctemp, MultiDet);

    int n = 0, i = 0;
    while (i < norbs)
    {
      if (Occ[n] > 0.0)
      {
        std::copy(Ctemp[n], Ctemp[n + 1], (*spo.C)[i]);
        i++;
      }
      n++;
    }

    hin.close();
  }
#ifdef HAVE_MPI
  myComm->comm.broadcast_n(spo.C->data(), spo.C->size());
#endif
  return true;
}


bool LCAOrbitalBuilder::putOccupation(LCAOrbitalSet& spo, xmlNodePtr occ_ptr)
{
  //die??
  if (spo.getBasisSetSize() == 0)
  {
    APP_ABORT("LCAOrbitalBuilder::putOccupation detected ZERO BasisSetSize");
    return false;
  }
  Occ.resize(std::max(spo.getBasisSetSize(), spo.getOrbitalSetSize()));
  Occ = 0.0;
  for (int i = 0; i < spo.getOrbitalSetSize(); i++)
    Occ[i] = 1.0;
  std::vector<int> occ_in;
  std::string occ_mode("table");
  if (occ_ptr == NULL)
  {
    occ_mode = "ground";
  }
  else
  {
    const std::string o(getXMLAttributeValue(occ_ptr, "mode"));
    if (!o.empty())
      occ_mode = o;
  }
  //Do nothing if mode == ground
  if (occ_mode == "excited")
  {
    putContent(occ_in, occ_ptr);
    for (int k = 0; k < occ_in.size(); k++)
    {
      if (occ_in[k] < 0) //remove this, -1 is to adjust the base
        Occ[-occ_in[k] - 1] = 0.0;
      else
        Occ[occ_in[k] - 1] = 1.0;
    }
  }
  else if (occ_mode == "table")
  {
    putContent(Occ, occ_ptr);
  }
  return true;
}

void LCAOrbitalBuilder::readRealMatrixFromH5(hdf_archive& hin,
                                             const std::string& setname,
                                             Matrix<LCAOrbitalBuilder::RealType>& Creal) const
{
  hin.read(Creal, setname);
}

void LCAOrbitalBuilder::LoadFullCoefsFromH5(hdf_archive& hin,
                                            int setVal,
                                            PosType& SuperTwist,
                                            Matrix<std::complex<RealType>>& Ctemp,
                                            bool MultiDet)
{
  Matrix<RealType> Creal;
  Matrix<RealType> Ccmplx;

  std::array<char, 72> name;
  int name_len{0};
  ///When running Single Determinant calculations, MO coeff loaded based on occupation and lowest eingenvalue.
  ///However, for solids with multideterminants, orbitals are order by kpoints; first all MOs for kpoint 1, then 2 etc
  /// The multideterminants occupation is specified in the input/HDF5 and theefore as long as there is consistency between
  /// the order in which we read the orbitals and the occupation, we are safe. In the case of Multideterminants generated
  /// by pyscf and Quantum Package, They are stored in the same order as generated for quantum package and one should use
  /// the orbitals labelled eigenset_unsorted.

  if (MultiDet == false)
    name_len = std::snprintf(name.data(), name.size(), "%s%d", "/Super_Twist/eigenset_", setVal);
  else
    name_len = std::snprintf(name.data(), name.size(), "%s%d", "/Super_Twist/eigenset_unsorted_", setVal);
  if (name_len < 0)
    throw std::runtime_error("Error generating name");

  std::string setname(name.data(), name_len);
  readRealMatrixFromH5(hin, setname, Creal);

  bool IsComplex = true;
  hin.read(IsComplex, "/parameters/IsComplex");
  if (IsComplex == false)
  {
    Ccmplx.resize(Creal.rows(), Creal.cols());
    Ccmplx = 0.0;
  }
  else
  {
    setname += "_imag";
    readRealMatrixFromH5(hin, setname, Ccmplx);
  }

  Ctemp.resize(Creal.rows(), Creal.cols());
  for (int i = 0; i < Ctemp.rows(); i++)
    for (int j = 0; j < Ctemp.cols(); j++)
      Ctemp[i][j] = std::complex<RealType>(Creal[i][j], Ccmplx[i][j]);
}

void LCAOrbitalBuilder::LoadFullCoefsFromH5(hdf_archive& hin,
                                            int setVal,
                                            PosType& SuperTwist,
                                            Matrix<RealType>& Creal,
                                            bool MultiDet)
{
  bool IsComplex = false;
  hin.read(IsComplex, "/parameters/IsComplex");
  if (IsComplex &&
      (std::abs(SuperTwist[0]) >= 1e-6 || std::abs(SuperTwist[1]) >= 1e-6 || std::abs(SuperTwist[2]) >= 1e-6))
  {
    std::string setname("This Wavefunction is Complex and you are using the real version of QMCPACK. "
                        "Please re-run this job with the Complex build of QMCPACK.");
    APP_ABORT(setname.c_str());
  }

  std::array<char, 72> name;
  int name_len{0};
  bool PBC = false;
  hin.read(PBC, "/PBC/PBC");
  if (MultiDet && PBC)
    name_len = std::snprintf(name.data(), name.size(), "%s%d", "/Super_Twist/eigenset_unsorted_", setVal);
  else
    name_len = std::snprintf(name.data(), name.size(), "%s%d", "/Super_Twist/eigenset_", setVal);
  if (name_len < 0)
    throw std::runtime_error("Error generating name");

  readRealMatrixFromH5(hin, std::string(name.data(), name_len), Creal);
}

/// Periodic Image Phase Factors computation to be determined
void LCAOrbitalBuilder::EvalPeriodicImagePhaseFactors(PosType SuperTwist,
                                                      std::vector<RealType>& LocPeriodicImagePhaseFactors)
{
  const int NbImages = (PBCImages[0] + 1) * (PBCImages[1] + 1) * (PBCImages[2] + 1);
  LocPeriodicImagePhaseFactors.resize(NbImages);
  for (size_t i = 0; i < NbImages; i++)
    LocPeriodicImagePhaseFactors[i] = 1.0;
}

void LCAOrbitalBuilder::EvalPeriodicImagePhaseFactors(PosType SuperTwist,
                                                      std::vector<std::complex<RealType>>& LocPeriodicImagePhaseFactors)
{
  // Allow computation to continue with no HDF file if the system has open boundary conditions.
  // The complex build is usually only used with open BC for testing.
  bool usesOpenBC = PBCImages[0] == 0 && PBCImages[1] == 0 && PBCImages[2] == 0;

  ///Exp(ik.g) where i is imaginary, k is the supertwist and g is the translation vector PBCImage.
  if (h5_path != "" && !usesOpenBC)
  {
    hdf_archive hin(myComm);
    if (myComm->rank() == 0)
    {
      if (!hin.open(h5_path, H5F_ACC_RDONLY))
        APP_ABORT("Could not open H5 file");

      hin.push("Cell", false);

      hin.read(Lattice, "LatticeVectors");
      hin.close();
    }
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        myComm->bcast(Lattice(i, j));
  }
  else if (!usesOpenBC)
  {
    APP_ABORT("Attempting to run PBC LCAO with no HDF5 support. Behaviour is unknown. Safer to exit");
  }

  int phase_idx = 0;
  int TransX, TransY, TransZ;
  RealType phase;

  for (int i = 0; i <= PBCImages[0]; i++) //loop Translation over X
  {
    TransX = ((i % 2) * 2 - 1) * ((i + 1) / 2);
    for (int j = 0; j <= PBCImages[1]; j++) //loop Translation over Y
    {
      TransY = ((j % 2) * 2 - 1) * ((j + 1) / 2);
      for (int k = 0; k <= PBCImages[2]; k++) //loop Translation over Z
      {
        TransZ = ((k % 2) * 2 - 1) * ((k + 1) / 2);
        RealType s, c;
        PosType Val;
        Val[0] = TransX * Lattice(0, 0) + TransY * Lattice(1, 0) + TransZ * Lattice(2, 0);
        Val[1] = TransX * Lattice(0, 1) + TransY * Lattice(1, 1) + TransZ * Lattice(2, 1);
        Val[2] = TransX * Lattice(0, 2) + TransY * Lattice(1, 2) + TransZ * Lattice(2, 2);

        phase = dot(SuperTwist, Val);
        qmcplusplus::sincos(phase, &s, &c);

        LocPeriodicImagePhaseFactors.emplace_back(c, s);
      }
    }
  }
}
} // namespace qmcplusplus

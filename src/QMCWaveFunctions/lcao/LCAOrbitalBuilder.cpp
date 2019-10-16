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


#include "OhmmsData/AttributeSet.h"
#include <QMCWaveFunctions/SPOSet.h>
#include <QMCWaveFunctions/lcao/NGFunctor.h>
#include <QMCWaveFunctions/lcao/MultiQuinticSpline1D.h>
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"
#include "QMCWaveFunctions/lcao/SoaSphericalTensor.h"
#include "QMCWaveFunctions/lcao/SoaAtomicBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
//#include "QMCWaveFunctions/lcao/RadialOrbitalSetBuilder.h"
#include "QMCWaveFunctions/lcao/AOBasisBuilder.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#include "QMCWaveFunctions/lcao/MultiFunctorAdapter.h"
#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/lcao/LCAOrbitalSetWithCorrection.h"
#include "QMCWaveFunctions/lcao/CuspCorrectionConstruction.h"
#endif
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"
#include "Utilities/ProgressReportEngine.h"

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


LCAOrbitalBuilder::LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder(comm), targetPtcl(els), sourcePtcl(ions), myBasisSet(nullptr), h5_path(""), doCuspCorrection(false)
{
  ClassName = "LCAOrbitalBuilder";
  ReportEngine PRE(ClassName, "createBasisSet");

  std::string keyOpt("NMO");       // Numerical Molecular Orbital
  std::string transformOpt("yes"); // Numerical Molecular Orbital
  std::string cuspC("no");         // cusp correction
  PosType SuperTwist(0.0);         // Supertwist coordinates
  cuspInfo = "";                   // file with precalculated cusp correction info
  OhmmsAttributeSet aAttrib;
  aAttrib.add(keyOpt, "keyword");
  aAttrib.add(keyOpt, "key");
  aAttrib.add(transformOpt, "transform");
  aAttrib.add(cuspC, "cuspCorrection");
  aAttrib.add(cuspInfo, "cuspInfo");
  aAttrib.add(h5_path, "href");
  aAttrib.add(PBCImages, "PBCimages");
  aAttrib.add(SuperTwist, "twist");
  aAttrib.put(cur);

  if (cur != NULL)
    aAttrib.put(cur);

  if (std::abs(SuperTwist[0] - 0.0) >= 1e-6 || std::abs(SuperTwist[1] - 0.0) >= 1e-6 ||
      std::abs(SuperTwist[2] - 0.0) >= 1e-6)
  {
    std::string error_msg("You are attempting to use a Super Twist other than Gamma. "
                          "This feature is being implemented but not supported yet. "
                          "Please contact developers for more details !!! Aborting.");
    APP_ABORT(error_msg.c_str());
  }

  radialOrbType = -1;
  if (transformOpt == "yes")
    radialOrbType = 0;
  else
  {
    if (keyOpt == "GTO")
      radialOrbType = 1;
    if (keyOpt == "STO")
      radialOrbType = 2;
  }

  if (radialOrbType < 0)
    PRE.error("Unknown radial function for LCAO orbitals. Specify keyword=\"NMO/GTO/STO\" .", true);

  if (cuspC == "yes")
    doCuspCorrection = true;
  //Evaluate the Phase factor. Equals 1 for OBC.
  EvalPeriodicImagePhaseFactors(SuperTwist);
  // no need to wait but load the basis set
  if (h5_path != "")
    loadBasisSetFromH5();
}

LCAOrbitalBuilder::~LCAOrbitalBuilder()
{
  //properly cleanup
}

void LCAOrbitalBuilder::loadBasisSetFromXML(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "loadBasisSetFromXML(xmlNodePtr)");
  if (myBasisSet)
  {
    app_log() << "Reusing previously loaded BasisSet." << std::endl;
    return;
  }

  if (!is_same(cur->name, "basisset"))
  { //heck to handle things like <sposet_builder>
    xmlNodePtr cur1 = cur->xmlChildrenNode;
    while (cur1 != NULL)
    {
      if (is_same(cur1->name, "basisset"))
        cur = cur1;
      cur1 = cur1->next;
    }
  }

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
}

void LCAOrbitalBuilder::loadBasisSetFromH5()
{
  ReportEngine PRE(ClassName, "loadBasisSetFromH5()");
  if (myBasisSet)
  {
    app_log() << "Reusing previously loaded BasisSet." << std::endl;
    return;
  }

  hdf_archive hin(myComm);
  int ylm = -1;
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      PRE.error("Could not open H5 file", true);
    if (!hin.push("basisset"))
      PRE.error("Could not open basisset group in H5; Probably Corrupt H5 file", true);

    std::string sph;
    std::string ElemID0 = "atomicBasisSet0";
    if (!hin.push(ElemID0.c_str()))
      PRE.error("Could not open  group Containing atomic Basis set in H5; Probably Corrupt H5 file", true);
    if (!hin.readEntry(sph, "angular"))
      PRE.error("Could not find name of  basisset group in H5; Probably Corrupt H5 file", true);
    ylm = (sph == "cartesian") ? 0 : 1;
    hin.close();
  }

  myComm->bcast(ylm);
  if (ylm < 0)
    PRE.error("Missing angular attribute of atomicBasisSet.", true);

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
        ao_type* aoBasis = any.createAOSet(cur);
        if (aoBasis)
        {
          //add the new atomic basis to the basis set
          int activeCenter = sourcePtcl.getSpeciesSet().findSpecies(elementType);
          mBasisSet->add(activeCenter, aoBasis);
        }
        ao_built_centers.push_back(elementType);
      }
    }
    cur = cur->next;
  } // done with basis set
  mBasisSet->setBasisSetSize(-1);
  mBasisSet->setPBCParams(PBCImages, PeriodicImagePhaseFactors);
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
    if (!hin.push("basisset"))
      PRE.error("Could not open basisset group in H5; Probably Corrupt H5 file", true);
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
      if (!hin.push(ElemType.c_str()))
        PRE.error("Could not open  group Containing atomic Basis set in H5; Probably Corrupt H5 file", true);
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
      ao_type* aoBasis = any.createAOSetH5(hin);
      if (aoBasis)
      {
        //add the new atomic basis to the basis set
        int activeCenter = sourcePtcl.getSpeciesSet().findSpecies(elementType);
        mBasisSet->add(activeCenter, aoBasis);
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
  mBasisSet->setPBCParams(PBCImages, PeriodicImagePhaseFactors);
  return mBasisSet;
}


SPOSet* LCAOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createSPO(xmlNodePtr)");
  std::string spo_name(""), id, cusp_file(""), optimize("no");
  OhmmsAttributeSet spoAttrib;
  spoAttrib.add(spo_name, "name");
  spoAttrib.add(id, "id");
  spoAttrib.add(cusp_file, "cuspInfo");
  spoAttrib.add(optimize, "optimize");
  spoAttrib.put(cur);

  if (myBasisSet == nullptr)
    PRE.error("Missing basisset.", true);

  if (optimize == "yes")
    app_log() << "  SPOSet " << spo_name << " is optimizable\n";

  LCAOrbitalSet* lcos = nullptr;
#if !defined(QMC_COMPLEX)
  LCAOrbitalSetWithCorrection* lcwc = nullptr;
  if (doCuspCorrection)
    lcos = lcwc = new LCAOrbitalSetWithCorrection(sourcePtcl, targetPtcl, myBasisSet, optimize == "yes");
  else
    lcos = new LCAOrbitalSet(myBasisSet, optimize == "yes");
#else
  lcos = new LCAOrbitalSet(myBasisSet, optimize == "yes");
#endif
  loadMO(*lcos, cur);

#if !defined(QMC_COMPLEX)
  if (doCuspCorrection)
  {
    int num_centers = sourcePtcl.getTotalNum();

    // Sometimes sposet attribute is 'name' and sometimes it is 'id'
    if (id == "")
      id = spo_name;

    int orbital_set_size = lcos->getOrbitalSetSize();
    Matrix<CuspCorrectionParameters> info(num_centers, orbital_set_size);

    bool valid = false;
    if (myComm->rank() == 0)
    {
      valid = readCuspInfo(cusp_file, id, orbital_set_size, info);
    }
#ifdef HAVE_MPI
    myComm->comm.broadcast_value(valid);
    if (valid)
    {
      for (int orb_idx = 0; orb_idx < orbital_set_size; orb_idx++)
      {
        for (int center_idx = 0; center_idx < num_centers; center_idx++)
        {
          broadcastCuspInfo(info(center_idx, orb_idx), *myComm, 0);
        }
      }
    }
#endif
    if (!valid)
    {
      generateCuspInfo(orbital_set_size, num_centers, info, targetPtcl, sourcePtcl, *lcwc, id, *myComm);
    }

    applyCuspCorrection(info, num_centers, orbital_set_size, targetPtcl, sourcePtcl, *lcwc, id);
  }
#endif

  return lcos;
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
  spo.setOrbitalSetSize(norb);
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
    else if (cname == "opt_vars")
    {
      spo.params_supplied = true;
      putContent(spo.params, cur);
    }
    cur = cur->next;
  }
  if (coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << std::endl;
    return spo.setIdentity(true);
  }
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
      //TO REVIEWERS:: IDEAL BEHAVIOUR SHOULD BE:
      /*
         if(!hin.push("PBC")
             PBC=false;
         else
            if (!hin.readEntry(PBC,"PBC"))
                APP_ABORT("Could not read PBC dataset in H5 file. Probably corrupt file!!!.");
        // However, it always succeeds to enter the if condition even if the group does not exists...
        */
      hin.push("PBC");
      PBC = false;
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

  //init_LCOrbitalSetOpt(orbital_mix_magnitude);

  return success;
}

bool LCAOrbitalBuilder::putFromXML(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr)
{
  spo.Identity = true;
  int norbs    = 0;
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
    spo.setIdentity(false);
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
#if defined(HAVE_LIBHDF5)
  int norbs  = spo.getOrbitalSetSize();
  int neigs  = spo.getBasisSetSize();
  int setVal = -1;
  std::string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal, "spindataset");
  aAttrib.add(neigs, "size");
  aAttrib.add(neigs, "orbitals");
  aAttrib.put(coeff_ptr);
  spo.setIdentity(false);
  hdf_archive hin(myComm);
  if (myComm->rank() == 0)
  {
    if (!hin.open(h5_path, H5F_ACC_RDONLY))
      APP_ABORT("LCAOrbitalBuilder::putFromH5 missing or incorrect path to H5 file.");

    Matrix<RealType> Ctemp(neigs, spo.getBasisSetSize());
    char name[72];

    //This is to make sure of Backward compatibility with previous tags.
    sprintf(name, "%s%d", "/Super_Twist/eigenset_", setVal);
    setname = name;
    if (!hin.readEntry(Ctemp, setname))
    {
      sprintf(name, "%s%d", "/KPTS_0/eigenset_", setVal);
      setname = name;
      if (!hin.readEntry(Ctemp, setname))
      {
        setname = "LCAOrbitalBuilder::putFromH5 Missing " + setname + " from HDF5 File.";
        APP_ABORT(setname.c_str());
      }
    }
    hin.close();

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
#else
  APP_ABORT("LCAOrbitalBuilder::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


/** read data from a hdf5 file
   * @param norb number of orbitals to be initialized
   * @param coeff_ptr xmlnode for coefficients
   */
bool LCAOrbitalBuilder::putPBCFromH5(LCAOrbitalSet& spo, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  ReportEngine PRE("LCAOrbitalBuilder", "LCAOrbitalBuilder::putPBCFromH5");
  int norbs      = spo.getOrbitalSetSize();
  int neigs      = spo.getBasisSetSize();
  int setVal     = -1;
  bool IsComplex = false;
  PosType SuperTwist(0.0);
  PosType SuperTwistH5(0.0);
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setVal, "spindataset");
  aAttrib.add(neigs, "size");
  aAttrib.add(neigs, "orbitals");
  aAttrib.put(coeff_ptr);
  spo.setIdentity(false);
  hdf_archive hin(myComm);

  xmlNodePtr curtemp = coeff_ptr->parent->parent->parent;
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

    Matrix<ValueType> Ctemp(neigs, spo.getBasisSetSize());
    LoadFullCoefsFromH5(hin, setVal, SuperTwist, Ctemp);

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

#else
  APP_ABORT("LCAOrbitalBuilder::putFromH5 HDF5 is disabled.")
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
    const XMLAttrString o(occ_ptr, "mode");
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

void readRealMatrixFromH5(hdf_archive& hin, const std::string& setname, Matrix<LCAOrbitalBuilder::RealType>& Creal)
{
  if (!hin.readEntry(Creal, setname))
  {
    std::string error_msg = "LCAOrbitalBuilder::readRealMatrixFromH5 Missing " + setname + " from HDF5 File.";
    APP_ABORT(error_msg.c_str());
  }
}

void LCAOrbitalBuilder::LoadFullCoefsFromH5(hdf_archive& hin,
                                            int setVal,
                                            PosType& SuperTwist,
                                            Matrix<std::complex<RealType>>& Ctemp)
{
  Matrix<RealType> Creal(Ctemp.rows(), Ctemp.cols());
  Matrix<RealType> Ccmplx(Ctemp.rows(), Ctemp.cols());

  char name[72];
  std::string setname;
  sprintf(name, "%s%d", "/Super_Twist/eigenset_", setVal);
  setname = name;
  readRealMatrixFromH5(hin, setname, Creal);

  if (std::abs(SuperTwist[0] - 0.0) < 1e-6 && std::abs(SuperTwist[1] - 0.0) < 1e-6 &&
      std::abs(SuperTwist[2] - 0.0) < 1e-6)
  {
    for (int i = 0; i < Ccmplx.rows(); i++)
      for (int j = 0; j < Ccmplx.cols(); j++)
        Ccmplx[i][j] = 0.0;
  }
  else
  {
    setname = std::string(name) + "_imag";
    readRealMatrixFromH5(hin, setname, Ccmplx);
  }

  for (int i = 0; i < Ctemp.rows(); i++)
    for (int j = 0; j < Ctemp.cols(); j++)
      Ctemp[i][j] = std::complex<RealType>(Creal[i][j], Ccmplx[i][j]);
}

void LCAOrbitalBuilder::LoadFullCoefsFromH5(hdf_archive& hin, int setVal, PosType& SuperTwist, Matrix<RealType>& Creal)
{
  bool IsComplex = false;
  //FIXME: need to check the path to IsComplex in h5
  hin.read(IsComplex, "IsComplex");
  if (IsComplex &&
      (std::abs(SuperTwist[0]) >= 1e-6 || std::abs(SuperTwist[1]) >= 1e-6 || std::abs(SuperTwist[2]) >= 1e-6))
  {
    std::string setname("This Wavefunction is Complex and you are using the real version of QMCPACK. "
                        "Please re-run this job with the Complex build of QMCPACK.");
    APP_ABORT(setname.c_str());
  }

  char name[72];
  sprintf(name, "%s%d", "/Super_Twist/eigenset_", setVal);
  readRealMatrixFromH5(hin, name, Creal);
}

///Function Not yet called. Periodic Image Phase Factors computation to be determined
void LCAOrbitalBuilder::EvalPeriodicImagePhaseFactors(PosType SuperTwist)
{
//NEED Specialization
#if not defined(QMC_COMPLEX)
  PeriodicImagePhaseFactors.resize(3);
  PeriodicImagePhaseFactors[0] = 1.0;
  PeriodicImagePhaseFactors[1] = 1.0;
  PeriodicImagePhaseFactors[2] = 1.0;
#else
  ///Exp(ik.g) where i is imaginary, k is the supertwist and g is the translation vector PBCImage.
  int phase_idx = 0;
  int TransX, TransY, TransZ;
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
        RealType vec_scalar = (TransX * SuperTwist[0] + TransY * SuperTwist[1] + TransZ * SuperTwist[2]);
        sincos(-2 * RealType(M_PI) * vec_scalar, &s, &c);
        PeriodicImagePhaseFactors.emplace_back(c, s);
      }
    }
  }
#endif
}
} // namespace qmcplusplus

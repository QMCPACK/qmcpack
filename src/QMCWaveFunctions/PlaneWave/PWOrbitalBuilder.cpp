//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 * @brief Definition of a builder class for PWOrbitalSet
 */
#include "PWOrbitalBuilder.h"
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
PWOrbitalBuilder::PWOrbitalBuilder(Communicate* comm, ParticleSet& els, const PSetMap& psets)
    : WaveFunctionComponentBuilder(comm, els), ptclPool(psets), hfileID(-1), rootNode(NULL)
{
  myParam = new PWParameterSet(myComm);
}

PWOrbitalBuilder::~PWOrbitalBuilder() { delete myParam; }

//All data parsing is handled here, outside storage classes.
std::unique_ptr<WaveFunctionComponent> PWOrbitalBuilder::buildComponent(xmlNodePtr cur)
{
  std::unique_ptr<WaveFunctionComponent> slater_det;
  //save the parent
  rootNode = cur;
  //
  //Get wavefunction data and parameters from XML and HDF5
  //

  //close it if open
  if (hfileID > 0)
    H5Fclose(hfileID);
  //check the current href
  hfileID = getH5(cur, "href");
  //no file, check the root
  if (hfileID < 0)
    hfileID = getH5(rootNode, "href");
  //Move through the XML tree and read basis information
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "basisset")
    {
      const std::string a(getXMLAttributeValue(cur, "ecut"));
      if (!a.empty())
        myParam->Ecut = std::stod(a);
    }
    else if (cname == "coefficients")
    {
      //close
      if (hfileID > 0)
        H5Fclose(hfileID);
      hfileID = getH5(cur, "hdata");
    }
    else if (cname == sd_tag)
    {
      if (hfileID < 0)
        hfileID = getH5(cur, "href");
      if (hfileID < 0)
      {
        APP_ABORT("  Cannot create a SlaterDet due to missing h5 file\n");
        OHMMS::Controller->abort();
      }
      createPWBasis(cur);
      slater_det = putSlaterDet(cur);
    }
    cur = cur->next;
  }
  H5Fclose(hfileID);
  return slater_det;
}

std::unique_ptr<WaveFunctionComponent> PWOrbitalBuilder::putSlaterDet(xmlNodePtr cur)
{
  //catch parameters
  myParam->put(cur);

  std::vector<std::unique_ptr<DiracDeterminantBase>> dets;
  int spin_group = 0;
  cur            = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    //Which determinant?
    if (cname == "determinant")
    {
      std::string id("updet");
      std::string ref("0");
      OhmmsAttributeSet aAttrib;
      aAttrib.add(id, "id");
      aAttrib.add(ref, "ref");
      aAttrib.put(cur);
      if (ref == "0")
        ref = id;
      const int firstIndex = targetPtcl.first(spin_group);
      const int lastIndex  = targetPtcl.last(spin_group);
      std::map<std::string, SPOSetPtr>::iterator lit(spomap.find(ref));
      //int spin_group=0;
      if (lit == spomap.end())
      {
        app_log() << "  Create a PWOrbitalSet" << std::endl;
        std::unique_ptr<SPOSet> psi(createPW(cur, spin_group));
        spomap[ref] = psi.get();
        dets.push_back(std::make_unique<DiracDeterminant<>>(std::move(psi), firstIndex, lastIndex));
      }
      else
      {
        app_log() << "  Reuse a PWOrbitalSet" << std::endl;
        std::unique_ptr<SPOSet> psi((*lit).second->makeClone());
        dets.push_back(std::make_unique<DiracDeterminant<>>(std::move(psi), firstIndex, lastIndex));
      }
      app_log() << "    spin=" << spin_group << " id=" << id << " ref=" << ref << std::endl;
      spin_group++;
    }
    cur = cur->next;
  }

  if (spin_group)
    return std::make_unique<SlaterDet>(targetPtcl, std::move(dets));
  ;

  myComm->barrier_and_abort(" Failed to create a SlaterDet at PWOrbitalBuilder::putSlaterDet ");
  return nullptr;
}

/** The read routine - get data from XML and H5. Process it and build orbitals.
 *
 * - parameters
 *   -- num_tiwsts
 *   -- num_bands
 *   -- complex_coefficients
 *   -- maximum_ecut
 * - basis
 */
bool PWOrbitalBuilder::createPWBasis(xmlNodePtr cur)
{
  //recycle int and double reader
  int idata;
  double ddata;
  HDFAttribIO<int> hdfint(idata);
  HDFAttribIO<double> hdfdbl(ddata);
  //start of parameters
  hdfint.read(hfileID, "electrons/number_of_kpoints");
  int nkpts = idata;
  hdfint.read(hfileID, "electrons/number_of_spins");
  hdfint.read(hfileID, "electrons/kpoint_0/spin_0/number_of_states");
  int nbands        = idata;
  myParam->numBands = nbands;
  app_log() << "Number of bands = " << nbands << std::endl;
  // Cutoff no longer present in the HDF file
  RealType ecut = 0.0;
  //end of parameters
  //check if input parameters are valid
  int nup   = targetPtcl.last(0);
  int ndown = targetPtcl.getTotalNum() - nup;
  if (nbands < nup || nbands < ndown)
  {
    app_error() << "Not enough bands in h5 file" << std::endl;
    OHMMS::Controller->abort();
  }
  std::string tname = myParam->getTwistAngleName();
  TinyVector<double, OHMMS_DIM> TwistAngle_DP;
  HDFAttribIO<TinyVector<double, OHMMS_DIM>> hdfobj_twist(TwistAngle_DP);
  hdfobj_twist.read(hfileID, "/electrons/kpoint_0/reduced_k");
  TwistAngle = TwistAngle_DP;
  if (!myBasisSet)
  {
    myBasisSet = std::make_unique<PWBasis>(TwistAngle);
  }
  //Read the planewave basisset.
  //Note that the same data is opened here for each twist angle-avoids duplication in the
  //h5 file (which may become very large).
  //return the ecut to be used by the basis set
  RealType real_ecut = myParam->getEcut(ecut);
  //create at least one basis set but do resize the containers
  int nh5gvecs = myBasisSet->readbasis(hfileID, real_ecut, targetPtcl.getLattice(), myParam->pwTag, myParam->pwMultTag);
  app_log() << "  num_twist = " << nkpts << std::endl;
  app_log() << "  twist angle = " << TwistAngle << std::endl;
  app_log() << "  num_bands = " << nbands << std::endl;
  app_log() << "  input maximum_ecut = " << ecut << std::endl;
  app_log() << "  current maximum_ecut = " << real_ecut << std::endl;
  app_log() << "  num_planewaves = " << nh5gvecs << std::endl;
  return true;
}

SPOSet* PWOrbitalBuilder::createPW(xmlNodePtr cur, int spinIndex)
{
  int nb = targetPtcl.last(spinIndex) - targetPtcl.first(spinIndex);
  std::vector<int> occBand(nb);
  for (int i = 0; i < nb; i++)
    occBand[i] = i;
  using GIndex_t = PWBasis::GIndex_t;
  GIndex_t nG(1);
  bool transform2grid = false;
  cur                 = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "transform")
    {
      putContent(nG, cur);
      transform2grid = true;
    }
    else if (cname == "occupation")
    {
      std::string occMode("ground");
      int bandoffset(1);
      OhmmsAttributeSet aAttrib;
      aAttrib.add(spinIndex, "spindataset");
      aAttrib.add(occMode, "mode");
      aAttrib.add(bandoffset, "offset"); /* reserved for index offset */
      aAttrib.put(cur);
      if (occMode == "excited")
      {
        std::vector<int> occ;
        std::vector<int> deleted, added;
        putContent(occ, cur);
        for (int i = 0; i < occ.size(); i++)
        {
          if (occ[i] < 0)
            deleted.push_back(-occ[i]);
          else
            added.push_back(occ[i]);
        }
        if (deleted.size() != added.size())
        {
          app_error() << "  Numbers of deleted and added bands are not identical." << std::endl;
          OHMMS::Controller->abort();
        }
        for (int i = 0; i < deleted.size(); i++)
        {
          occBand[deleted[i] - bandoffset] = added[i] - bandoffset;
        }
        app_log() << "  mode=\"excited\" Occupied states: " << std::endl;
        copy(occBand.begin(), occBand.end(), std::ostream_iterator<int>(app_log(), " "));
        app_log() << std::endl;
      }
    }
    cur = cur->next;
  }
  //std::string tname=myParam->getTwistName();
  std::string tname = "kpoint_0";
  //hid_t es_grp_id = H5Gopen2(hfileID,myParam->eigTag.c_str(), H5P_DEFAULT);
  //hid_t es_grp_id = H5Gopen2(hfileID,"electrons/kpoint_0/spin_0/state_0", H5P_DEFAULT);
  hid_t es_grp_id = H5Gopen2(hfileID, "electrons", H5P_DEFAULT);

  //hid_t twist_grp_id = H5Gopen2(es_grp_id,tname.c_str(), H5P_DEFAULT);
  hid_t twist_grp_id = H5Gopen2(es_grp_id, "kpoint_0", H5P_DEFAULT);
  //create a single-particle orbital set
  SPOSetType* psi = new SPOSetType;
  if (transform2grid)
  {
    nb = myParam->numBands;
    occBand.resize(nb);
    for (int i = 0; i < nb; i++)
      occBand[i] = i;
  }
  //going to take care of occ
  psi->resize(new PWBasis(*myBasisSet), nb, true);
  if (myParam->hasComplexData(hfileID)) //input is complex
  {
    //app_log() << "  PW coefficients are complex." << std::endl;
    using TempVecType    = std::vector<std::complex<RealType>>;
    using TempVecType_DP = std::vector<std::complex<double>>;
    TempVecType_DP coefs_DP(myBasisSet->inputmap.size());
    HDFAttribIO<TempVecType_DP> hdfobj_coefs(coefs_DP);
    int ib = 0;
    while (ib < nb)
    {
      std::string bname(myParam->getBandName(occBand[ib], spinIndex));
      app_log() << "  Reading " << tname << "/" << bname << std::endl;
      hid_t band_grp_id = H5Gopen2(twist_grp_id, bname.c_str(), H5P_DEFAULT);
      hdfobj_coefs.read(band_grp_id, "psi_g");
      TempVecType coefs(coefs_DP.begin(), coefs_DP.end());
      psi->addVector(coefs, ib);
      H5Gclose(band_grp_id);
      ++ib;
    }
  }
  else
  {
    // It appears the coefficients are always stored as complex in the HDF file?
    //app_log() << "  PW coefficients are real." << std::endl;
    using ComplexTempVecType    = std::vector<std::complex<RealType>>;
    using ComplexTempVecType_DP = std::vector<std::complex<double>>;
    ComplexTempVecType_DP complex_coefs_DP(myBasisSet->inputmap.size());
    HDFAttribIO<ComplexTempVecType_DP> hdfobj_complex_coefs(complex_coefs_DP);
    int ib = 0;
    while (ib < nb)
    {
      std::string bname(myParam->getBandName(occBand[ib], spinIndex));
      app_log() << "  Reading " << tname << "/" << bname << std::endl;
      hid_t band_grp_id = H5Gopen2(twist_grp_id, bname.c_str(), H5P_DEFAULT);
      hdfobj_complex_coefs.read(band_grp_id, "psi_g");
      ComplexTempVecType complex_coefs(complex_coefs_DP.begin(), complex_coefs_DP.end());
      psi->addVector(complex_coefs, ib);
      H5Gclose(band_grp_id);
      ++ib;
    }
  }
  H5Gclose(twist_grp_id);
  H5Gclose(es_grp_id);
#if defined(QMC_COMPLEX)
  if (transform2grid)
  {
    app_warning() << "  Going to transform on grid " << std::endl;
    transform2GridData(nG, spinIndex, *psi);
  }
#endif
  return psi;
}

#if defined(QMC_COMPLEX)
void PWOrbitalBuilder::transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc)
{
  std::ostringstream splineTag;
  splineTag << "eigenstates_" << nG[0] << "_" << nG[1] << "_" << nG[2];
  herr_t status = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  app_log() << " splineTag " << splineTag.str() << std::endl;
  hid_t es_grp_id;
  if (H5Lexists(hfileID, splineTag.str().c_str(), H5P_DEFAULT) != true)
  {
    es_grp_id = H5Gcreate2(hfileID, splineTag.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDFAttribIO<PWBasis::GIndex_t> t(nG);
    t.write(es_grp_id, "grid");
  }
  else
  {
    es_grp_id = H5Gopen2(hfileID, splineTag.str().c_str(), H5P_DEFAULT);
  }
  std::string tname = myParam->getTwistName();
  hid_t twist_grp_id;
  if (H5Lexists(es_grp_id, tname.c_str(), H5P_DEFAULT) != true)
    twist_grp_id = H5Gcreate2(es_grp_id, tname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  else
    twist_grp_id = H5Gopen2(es_grp_id, tname.c_str(), H5P_DEFAULT);
  TinyVector<double, OHMMS_DIM> TwistAngle_DP;
  TwistAngle_DP = TwistAngle;
  HDFAttribIO<TinyVector<double, OHMMS_DIM>> hdfobj_twist(TwistAngle_DP);
  hdfobj_twist.write(twist_grp_id, "twist_angle");
  const ParticleSet::ParticleLayout& lattice(targetPtcl.getLattice());
  RealType dx = 1.0 / static_cast<RealType>(nG[0] - 1);
  RealType dy = 1.0 / static_cast<RealType>(nG[1] - 1);
  RealType dz = 1.0 / static_cast<RealType>(nG[2] - 1);
#if defined(VERYTINYMEMORY)
  using StorageType = Array<ParticleSet::SingleParticleValue, 3>;
  StorageType inData(nG[0], nG[1], nG[2]);
  int ib = 0;
  while (ib < myParam->numBands)
  {
    std::string bname(myParam->getBandName(ib));
    hid_t band_grp_id, spin_grp_id = -1;
    if (H5Lexists(twist_grp_id, bname.c_str(), H5P_DEFAULT) != true)
    {
      band_grp_id = H5Gcreate2(twist_grp_id, bname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
      band_grp_id = H5Gopen2(twist_grp_id, bname.c_str(), H5P_DEFAULT);
    }
    hid_t parent_id = band_grp_id;
    if (myParam->hasSpin)
    {
      bname = myParam->getSpinName(spinIndex);
      if (H5Lexists(band_grp_id, bname.c_str(), H5P_DEFAULT) != true)
      {
        spin_grp_id = H5Gcreate2(band_grp_id, bname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      }
      else
      {
        spin_grp_id = H5Gopen2(band_grp_id, bname.c_str(), H5P_DEFAULT);
      }
      parent_id = spin_grp_id;
    }
    for (int ig = 0; ig < nG[0]; ig++)
    {
      RealType x = ig * dx;
      for (int jg = 0; jg < nG[1]; jg++)
      {
        RealType y = jg * dy;
        for (int kg = 0; kg < nG[2]; kg++)
        {
          inData(ig, jg, kg) = pwFunc.evaluate(ib, lattice.toCart(PosType(x, y, kg * dz)));
        }
      }
    }
    app_log() << "  Add spline data " << ib << " h5path=" << tname << "/eigvector" << std::endl;
    HDFAttribIO<StorageType> t(inData);
    t.write(parent_id, myParam->eigvecTag.c_str());
    if (spin_grp_id >= 0)
      H5Gclose(spin_grp_id);
    H5Gclose(band_grp_id);
    ++ib;
  }
#else
  using StorageType = Array<ParticleSet::SingleParticleValue, 3>;
  std::vector<StorageType*> inData;
  int nb = myParam->numBands;
  for (int ib = 0; ib < nb; ib++)
    inData.push_back(new StorageType(nG[0], nG[1], nG[2]));
  PosType tAngle = targetPtcl.getLattice().k_cart(TwistAngle);
  PWOrbitalSet::ValueVector phi(nb);
  for (int ig = 0; ig < nG[0]; ig++)
  {
    RealType x = ig * dx;
    for (int jg = 0; jg < nG[1]; jg++)
    {
      RealType y = jg * dy;
      for (int kg = 0; kg < nG[2]; kg++)
      {
        targetPtcl.R[0] = lattice.toCart(PosType(x, y, kg * dz));
        pwFunc.evaluateValue(targetPtcl, 0, phi);
        RealType x(dot(targetPtcl.R[0], tAngle));
        ValueType phase(std::cos(x), -std::sin(x));
        for (int ib = 0; ib < nb; ib++)
          (*inData[ib])(ig, jg, kg) = phase * phi[ib];
      }
    }
  }
  for (int ib = 0; ib < nb; ib++)
  {
    std::string bname(myParam->getBandName(ib));
    hid_t band_grp_id, spin_grp_id = -1;
    if (H5Lexists(twist_grp_id, bname.c_str(), H5P_DEFAULT) != true)
    {
      band_grp_id = H5Gcreate2(twist_grp_id, bname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
      band_grp_id = H5Gopen2(twist_grp_id, bname.c_str(), H5P_DEFAULT);
    }
    hid_t parent_id = band_grp_id;
    if (myParam->hasSpin)
    {
      bname = myParam->getSpinName(spinIndex);
      if (H5Lexists(band_grp_id, bname.c_str(), H5P_DEFAULT) != true)
      {
        spin_grp_id = H5Gcreate2(band_grp_id, bname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      }
      else
      {
        spin_grp_id = H5Gopen2(band_grp_id, bname.c_str(), H5P_DEFAULT);
      }
      parent_id = spin_grp_id;
    }
    app_log() << "  Add spline data " << ib << " h5path=" << tname << "/eigvector" << std::endl;
    HDFAttribIO<StorageType> t(*(inData[ib]));
    t.write(parent_id, myParam->eigvecTag.c_str());
    if (spin_grp_id >= 0)
      H5Gclose(spin_grp_id);
    H5Gclose(band_grp_id);
  }
  for (int ib = 0; ib < nb; ib++)
    delete inData[ib];
#endif
  H5Gclose(twist_grp_id);
  H5Gclose(es_grp_id);
}
#endif

hid_t PWOrbitalBuilder::getH5(xmlNodePtr cur, const char* aname)
{
  const std::string a(getXMLAttributeValue(cur, aname));
  if (a.empty())
    return -1;
  hid_t h = H5Fopen(a.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h < 0)
  {
    app_error() << " Cannot open " << a << " file." << std::endl;
    OHMMS::Controller->abort();
  }
  myParam->checkVersion(h);
  //overwrite the parameters
  myParam->put(rootNode);
  return h;
}

} // namespace qmcplusplus

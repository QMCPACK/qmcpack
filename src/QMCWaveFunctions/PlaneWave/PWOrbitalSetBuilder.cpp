//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 * @brief Definition of a builder class for PWOrbitalSet
 */
#include "PWOrbitalSetBuilder.h"
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
PWOrbitalSetBuilder::PWOrbitalSetBuilder(const ParticleSet& p, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("Planewave", comm),
      targetPtcl(p),
      rootNode(cur),
      myParam{std::make_unique<PWParameterSet>(comm)},
      hfile{comm}
{
  //
  //Get wavefunction data and parameters from XML and HDF5
  //
  //catch parameters
  myParam->put(cur);
  //check the current href
  bool success = getH5(cur, "href");
  //Move through the XML tree and read basis information
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "basisset")
    {
      const std::string a(getXMLAttributeValue(element, "ecut"));
      if (!a.empty())
        myParam->Ecut = std::stod(a);
    }
    else if (cname == "coefficients")
    {
      //close
      if (success)
        hfile.close();
      success = getH5(element, "hdata");
    }
  });

  if (!success)
    throw std::runtime_error("h5 cannot be open for creating PW basis!");
  //create PW Basis
  createPWBasis();
}

PWOrbitalSetBuilder::~PWOrbitalSetBuilder() = default;

std::unique_ptr<SPOSet> PWOrbitalSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int spin_group = 0;
  std::string spo_object_name;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(spin_group, "spindataset");
  aAttrib.add(spo_object_name, "name");
  aAttrib.add(spo_object_name, "id");
  aAttrib.put(cur);
  return createPW(cur, spo_object_name, spin_group);
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
bool PWOrbitalSetBuilder::createPWBasis()
{
  //recycle int and double reader
  int idata = 0;
  //start of parameters
  hfile.read(idata, "electrons/number_of_kpoints");
  int nkpts = idata;
  hfile.read(idata, "electrons/number_of_spins");
  hfile.read(idata, "electrons/kpoint_0/spin_0/number_of_states");
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
  hfile.read(TwistAngle_DP, "/electrons/kpoint_0/reduced_k");
  TwistAngle = TwistAngle_DP;
  if (!myBasisSet)
    myBasisSet = std::make_unique<PWBasis>(TwistAngle);

  //Read the planewave basisset.
  //Note that the same data is opened here for each twist angle-avoids duplication in the
  //h5 file (which may become very large).
  //return the ecut to be used by the basis set
  RealType real_ecut = myParam->getEcut(ecut);
  //create at least one basis set but do resize the containers
  int nh5gvecs = myBasisSet->readbasis(hfile, real_ecut, targetPtcl.getLattice(), myParam->pwTag, myParam->pwMultTag);
  app_log() << "  num_twist = " << nkpts << std::endl;
  app_log() << "  twist angle = " << TwistAngle << std::endl;
  app_log() << "  num_bands = " << nbands << std::endl;
  app_log() << "  input maximum_ecut = " << ecut << std::endl;
  app_log() << "  current maximum_ecut = " << real_ecut << std::endl;
  app_log() << "  num_planewaves = " << nh5gvecs << std::endl;
  return true;
}

std::unique_ptr<SPOSet> PWOrbitalSetBuilder::createPW(xmlNodePtr cur, const std::string& objname, int spinIndex)
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
  std::string tname = "kpoint_0";
  hfile.push("electrons", false);
  hfile.push("kpoint_0", false);
  //create a single-particle orbital set
  auto psi = std::make_unique<SPOSetType>(objname);
  if (transform2grid)
  {
    nb = myParam->numBands;
    occBand.resize(nb);
    for (int i = 0; i < nb; i++)
      occBand[i] = i;
  }
  //going to take care of occ
  psi->resize(new PWBasis(*myBasisSet), nb, true);
  if (myParam->hasComplexData(hfile)) //input is complex
  {
    //app_log() << "  PW coefficients are complex." << std::endl;
    using TempVecType    = std::vector<std::complex<RealType>>;
    using TempVecType_DP = std::vector<std::complex<double>>;
    TempVecType_DP coefs_DP(myBasisSet->inputmap.size());
    int ib = 0;
    while (ib < nb)
    {
      std::string bname(myParam->getBandName(occBand[ib], spinIndex));
      app_log() << "  Reading " << tname << "/" << bname << std::endl;
      hfile.push(bname, false);
      hfile.read(coefs_DP, "psi_g");
      TempVecType coefs(coefs_DP.begin(), coefs_DP.end());
      psi->addVector(coefs, ib);
      hfile.pop();
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
    int ib = 0;
    while (ib < nb)
    {
      std::string bname(myParam->getBandName(occBand[ib], spinIndex));
      app_log() << "  Reading " << tname << "/" << bname << std::endl;
      hfile.push(bname, false);
      hfile.read(complex_coefs_DP, "psi_g");
      ComplexTempVecType complex_coefs(complex_coefs_DP.begin(), complex_coefs_DP.end());
      psi->addVector(complex_coefs, ib);
      hfile.pop();
      ++ib;
    }
  }
  hfile.pop();
  hfile.pop();
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
void PWOrbitalSetBuilder::transform2GridData(PWBasis::GIndex_t& nG, int spinIndex, PWOrbitalSet& pwFunc)
{
  std::ostringstream splineTag;
  splineTag << "eigenstates_" << nG[0] << "_" << nG[1] << "_" << nG[2];
  herr_t status            = H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  std::string splineTagStr = splineTag.str();
  app_log() << " splineTag " << splineTagStr << std::endl;
  if (!hfile.is_group(splineTagStr))
  {
    hfile.push(splineTagStr, true);
    hfile.write(nG, "grid");
  }
  else
  {
    hfile.push(splineTagStr, false);
  }
  std::string tname = myParam->getTwistName();
  hfile.push(tname, true);

  TinyVector<double, OHMMS_DIM> TwistAngle_DP;
  TwistAngle_DP = TwistAngle;
  hfile.write(TwistAngle_DP, "twist_angle");
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
    hfile.push(bname, true);
    if (myParam->hasSpin)
    {
      bname = myParam->getSpinName(spinIndex);
      hfile.push(bname, true);
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
    hfile.write(inData, myParam->eigvecTag);
    if (myParam->hasSpin)
      mfile.pop();
    mfile.pop();
    ++ib;
  }
#else
  using StorageType = Array<ParticleSet::SingleParticleValue, 3>;
  UPtrVector<StorageType> inData;
  int nb = myParam->numBands;
  for (int ib = 0; ib < nb; ib++)
    inData.push_back(std::make_unique<StorageType>(nG[0], nG[1], nG[2]));
  PosType tAngle = targetPtcl.getLattice().k_cart(TwistAngle);
  ParticleSet ptemp(targetPtcl.getSimulationCell());
  ptemp.create({1});
  PWOrbitalSet::ValueVector phi(nb);
  for (int ig = 0; ig < nG[0]; ig++)
  {
    RealType x = ig * dx;
    for (int jg = 0; jg < nG[1]; jg++)
    {
      RealType y = jg * dy;
      for (int kg = 0; kg < nG[2]; kg++)
      {
        ptemp.R[0] = lattice.toCart(PosType(x, y, kg * dz));
        pwFunc.evaluateValue(ptemp, 0, phi);
        RealType x(dot(ptemp.R[0], tAngle));
        ValueType phase(std::cos(x), -std::sin(x));
        for (int ib = 0; ib < nb; ib++)
          (*inData[ib])(ig, jg, kg) = phase * phi[ib];
      }
    }
  }
  for (int ib = 0; ib < nb; ib++)
  {
    std::string bname(myParam->getBandName(ib));
    hfile.push(bname, true);
    if (myParam->hasSpin)
    {
      bname = myParam->getSpinName(spinIndex);
      hfile.push(bname, true);
    }
    app_log() << "  Add spline data " << ib << " h5path=" << tname << "/eigvector" << std::endl;
    hfile.write(*inData[ib], myParam->eigvecTag);
    if (myParam->hasSpin)
      hfile.pop();
    hfile.pop();
  }
#endif
  hfile.pop();
  hfile.pop();
}
#endif

bool PWOrbitalSetBuilder::getH5(xmlNodePtr cur, const char* aname)
{
  const std::string a(getXMLAttributeValue(cur, aname));
  if (a.empty())
    return false;

  bool success = hfile.open(a, H5F_ACC_RDONLY | H5P_DEFAULT);
  if (!success)
  {
    app_error() << " Cannot open " << a << " file." << std::endl;
    OHMMS::Controller->abort();
  }
  myParam->checkVersion(hfile);
  //overwrite the parameters
  myParam->put(rootNode);
  return success;
}

} // namespace qmcplusplus

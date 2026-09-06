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
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include "QMCWaveFunctions/PlaneWave/PWBandInfo.h"
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/PlaneWave/PWTiling.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
PWOrbitalSetBuilder::PWOrbitalSetBuilder(ParticleSet& p, Communicate* comm, xmlNodePtr cur)
    : SPOSetBuilder("Planewave", comm),
      targetPtcl(p),
      rootNode(cur),
      myParam{std::make_unique<PWParameterSet>(comm)},
      hfile{comm}
{
  for (int i = 0; i < OHMMS_DIM; ++i)
    for (int j = 0; j < OHMMS_DIM; ++j)
      tileMatrix(i, j) = i == j ? 1 : 0;
  OhmmsAttributeSet attributes;
  attributes.add(tileMatrix, "tilematrix");
  attributes.add(sortBands, "sort");
  attributes.put(cur);

  //
  //Get wavefunction data and parameters from XML and HDF5
  //
  //catch parameters
  myParam->put(cur);
  const std::string twist_attribute(getXMLAttributeValue(cur, "twist"));
  if (!twist_attribute.empty())
  {
    std::istringstream twist_stream(twist_attribute);
    std::string trailing_value;
    if (twist_stream >> requestedTwist[0] >> requestedTwist[1] >> requestedTwist[2] &&
        !(twist_stream >> trailing_value))
      hasRequestedTwist = true;
    else
    {
      app_warning() << "Using the twist attribute as an HDF group tag is deprecated; use twist_tag instead.\n";
    }
  }
  if (sortBands < 0 || sortBands > 2)
    throw std::runtime_error("Plane-wave sort must be 0, 1, or 2");
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
  if (!twist_attribute.empty() && !hasRequestedTwist)
    myParam->twistTag = twist_attribute;
  //create PW Basis
  createPWBasis();
}

PWOrbitalSetBuilder::~PWOrbitalSetBuilder() = default;

std::unique_ptr<SPOSet> PWOrbitalSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int spin_group    = 0;
  int orbital_count = 0;
  std::string spo_object_name;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(spin_group, "spindataset");
  aAttrib.add(orbital_count, "size");
  aAttrib.add(spo_object_name, "name");
  aAttrib.add(spo_object_name, "id");
  aAttrib.put(cur);
  if (orbital_count < 1)
    throw std::runtime_error("Plane-wave SPO size must be positive");
  return createPW(cur, spo_object_name, spin_group, orbital_count);
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
  twistAngles.resize(nkpts);
  for (int kpoint_index = 0; kpoint_index < nkpts; ++kpoint_index)
  {
    TinyVector<double, OHMMS_DIM> twist_angle_dp;
    hfile.read(twist_angle_dp, "/electrons/kpoint_" + std::to_string(kpoint_index) + "/reduced_k");
    twistAngles[kpoint_index] = twist_angle_dp;
  }
  const pw::TwistGroups twist_groups = pw::groupPrimitiveTwists(tileMatrix, twistAngles);

  const pw::Lattice target_lattice = targetPtcl.getLattice().R;
  pw::Lattice primitive_lattice    = target_lattice;
  if (hfile.is_dataset("/supercell/primitive_vectors"))
    hfile.read(primitive_lattice, "/supercell/primitive_vectors");
  pw::validateLattice(tileMatrix, primitive_lattice, target_lattice);

  int selected_twist_index;
  if (hasRequestedTwist)
  {
    if (myParam->twistIndex >= 0)
      app_warning() << "twist attribute exists; twistnum/twistIndex is ignored.\n";
    selected_twist_index = pw::findTwistGroup(twist_groups, requestedTwist);
  }
  else if (myParam->twistIndex >= 0)
  {
    if (myParam->twistIndex >= twist_groups.primitive_indices.size())
      throw std::runtime_error("Requested plane-wave supercell twist index is outside the available range");
    selected_twist_index = myParam->twistIndex;
  }
  else
    selected_twist_index = pw::findTwistGroup(twist_groups, PosType(0));

  pw::validateTwistGroup(tileMatrix, twistAngles, twist_groups, selected_twist_index);
  superTwist          = twist_groups.super_twists[selected_twist_index];
  includedKPoints     = twist_groups.primitive_indices[selected_twist_index];
  myParam->twistIndex = selected_twist_index;
#if !defined(QMC_COMPLEX)
  if (!pw::isRealCompatible(superTwist))
    throw std::runtime_error("Requested plane-wave supercell twist requires a complex QMCPACK build");
#endif
  targetPtcl.setTwist(superTwist);

  TwistAngle = twistAngles[includedKPoints.front()];
  if (!myBasisSet)
    myBasisSet = std::make_unique<PWBasis>(TwistAngle);

  //Read the planewave basisset.
  //Note that the same data is opened here for each twist angle-avoids duplication in the
  //h5 file (which may become very large).
  //return the ecut to be used by the basis set
  RealType real_ecut = myParam->getEcut(ecut);
  //create at least one basis set but do resize the containers
  int nh5gvecs =
      myBasisSet->readbasis(hfile, real_ecut, targetPtcl.getLattice(), myParam->pwTag, myParam->pwMultTag, true);
  app_log() << "  num_twist = " << nkpts << std::endl;
  app_log() << "  included k-points = ";
  std::copy(includedKPoints.begin(), includedKPoints.end(), std::ostream_iterator<int>(app_log(), " "));
  app_log() << std::endl;
  app_log() << "  twist angle = " << TwistAngle << std::endl;
  app_log() << "  num_bands = " << nbands << std::endl;
  app_log() << "  input maximum_ecut = " << ecut << std::endl;
  app_log() << "  current maximum_ecut = " << real_ecut << std::endl;
  app_log() << "  num_planewaves = " << nh5gvecs << std::endl;
  return true;
}

std::unique_ptr<SPOSet> PWOrbitalSetBuilder::createPW(xmlNodePtr cur,
                                                      const std::string& objname,
                                                      int spinIndex,
                                                      int orbital_count)
{
  std::vector<int> occBand(orbital_count);
  for (int i = 0; i < orbital_count; i++)
    occBand[i] = i;
  using GIndex_t = PWBasis::GIndex_t;
  GIndex_t nG(1);
  bool transform2grid     = false;
  bool excited_occupation = false;
  cur                     = cur->children;
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
        excited_occupation = true;
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

  if (includedKPoints.size() > 1 && (transform2grid || excited_occupation))
    throw std::runtime_error("Multi-k-point plane-wave SPOs do not support transform or excited occupations");

  if (transform2grid)
  {
    orbital_count = myParam->numBands;
    occBand.resize(orbital_count);
    for (int i = 0; i < orbital_count; i++)
      occBand[i] = i;
  }

  if (includedKPoints.empty())
    throw std::runtime_error("No primitive k-points contribute to the requested plane-wave supercell twist");

  std::vector<pw::BandInfo> bands;
#if defined(QMC_COMPLEX)
  for (const int kpoint_index : includedKPoints)
#else
  for (const pw::DistinctTwist& distinct_twist : pw::findDistinctRealTwists(includedKPoints, twistAngles))
#endif
  {
#if defined(QMC_COMPLEX)
    constexpr bool make_two_copies = false;
#else
    const int kpoint_index     = distinct_twist.twist_index;
    const bool make_two_copies = distinct_twist.make_two_copies;
#endif
    std::vector<double> eigenvalues;
    hfile.read(eigenvalues,
               "/electrons/kpoint_" + std::to_string(kpoint_index) + "/spin_" + std::to_string(spinIndex) +
                   "/eigenvalues");
    for (int band_index = 0; band_index < eigenvalues.size(); ++band_index)
      bands.push_back({eigenvalues[band_index], kpoint_index, band_index, false, make_two_copies});
  }
  pw::sortBands(bands, sortBands);

  std::vector<pw::BandInfo> expanded_bands;
  for (const pw::BandInfo& band : bands)
  {
    expanded_bands.push_back(band);
    if (band.make_two_copies)
    {
      pw::BandInfo imaginary_band       = band;
      imaginary_band.use_imaginary_part = true;
      expanded_bands.push_back(imaginary_band);
    }
  }

  if (expanded_bands.size() < orbital_count)
    throw std::runtime_error("Not enough plane-wave bands for the requested tiled SPO");

  std::vector<pw::BandInfo> selected_bands;
  selected_bands.reserve(orbital_count);
  if (excited_occupation)
  {
    for (const int occupied_index : occBand)
    {
      if (occupied_index < 0 || occupied_index >= expanded_bands.size())
        throw std::runtime_error("Plane-wave excited occupation index is outside the available band range");
      selected_bands.push_back(expanded_bands[occupied_index]);
    }
  }
  else
    selected_bands.assign(expanded_bands.begin(), expanded_bands.begin() + orbital_count);

  struct BandRun
  {
    int kpoint_index;
    std::vector<pw::BandInfo> bands;
  };
  std::vector<BandRun> band_runs;
  for (const pw::BandInfo& band : selected_bands)
  {
    if (band_runs.empty() || band_runs.back().kpoint_index != band.twist_index)
      band_runs.push_back({band.twist_index, {}});
    band_runs.back().bands.push_back(band);
  }

  std::vector<std::unique_ptr<SPOSet>> components;
  components.reserve(band_runs.size());
  for (int run_index = 0; run_index < band_runs.size(); ++run_index)
  {
    const int kpoint_index                           = band_runs[run_index].kpoint_index;
    const std::vector<pw::BandInfo>& component_bands = band_runs[run_index].bands;

    std::unique_ptr<PWBasis> basis;
    if (kpoint_index == includedKPoints.front())
      basis = std::make_unique<PWBasis>(*myBasisSet);
    else
    {
      basis = std::make_unique<PWBasis>(twistAngles[kpoint_index]);
      basis->readbasis(hfile, myParam->getEcut(0.0), targetPtcl.getLattice(), myParam->pwTag, myParam->pwMultTag, true);
    }

    const std::string kpoint_name = "kpoint_" + std::to_string(kpoint_index);
    hfile.push("electrons", false);
    hfile.push(kpoint_name, false);
    int available_bands = 0;
    hfile.read(available_bands, "spin_" + std::to_string(spinIndex) + "/number_of_states");
    if (std::find_if(component_bands.begin(), component_bands.end(), [&](const pw::BandInfo& band) {
          return band.band_index >= available_bands;
        }) != component_bands.end())
      throw std::runtime_error("Not enough plane-wave bands for a folded primitive k-point");

    const std::string component_name = band_runs.size() == 1 ? objname : objname + "_run" + std::to_string(run_index);
    auto psi                         = std::make_unique<SPOSetType>(component_name, component_bands.size());
    const size_t coefficient_count   = basis->inputmap.size();
    psi->resize(basis.release(), true);

    using TempVecType    = std::vector<std::complex<RealType>>;
    using TempVecType_DP = std::vector<std::complex<double>>;
    TempVecType_DP coefs_DP(coefficient_count);
    for (int orbital_index = 0; orbital_index < component_bands.size(); ++orbital_index)
    {
      const int band_index = component_bands[orbital_index].band_index;
      std::string bname(myParam->getBandName(band_index, spinIndex));
      app_log() << "  Reading " << kpoint_name << "/" << bname << std::endl;
      hfile.push(bname, false);
      hfile.read(coefs_DP, "psi_g");
      TempVecType coefs(coefs_DP.begin(), coefs_DP.end());
#if defined(QMC_COMPLEX)
      psi->addVector(coefs, orbital_index);
#else
      psi->addVector(coefs, orbital_index, component_bands[orbital_index].use_imaginary_part);
#endif
      hfile.pop();
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
    components.emplace_back(std::move(psi));
  }

  if (components.size() == 1)
    return std::move(components.front());
  return std::make_unique<CompositeSPOSet<ValueType>>(objname, std::move(components));
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
  const Lattice& lattice(targetPtcl.getLattice());
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

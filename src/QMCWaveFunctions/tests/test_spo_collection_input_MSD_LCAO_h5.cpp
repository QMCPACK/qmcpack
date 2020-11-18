//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
void test_LiH_msd_xml_input(const std::string& spo_xml_string, const std::string& check_sponame, int check_spo_size, int check_basisset_size)
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSet ions_;
  ParticleSet elec_;

  ions_.setName("ion0");
  ions_.create({1,1});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {0.0, 0.0, 3.0139239693};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Li");
  int HIdx             = ispecies.addSpecies("H");

  elec_.setName("elec");
  elec_.create({2, 2});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 0.0, 1.0};
  elec_.R[2] = {0.0, 1.0, 0.0};
  elec_.R[3] = {1.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  // Need 1 electron and 1 proton, somehow
  //ParticleSet target = ParticleSet();
  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.addParticleSet(&elec_);
  ptcl.addParticleSet(&ions_);

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(ein_xml);

  SPOSet* spo_ptr(get_sposet(check_sponame));
  REQUIRE(spo_ptr != nullptr);
  REQUIRE(spo_ptr->getOrbitalSetSize() == check_spo_size);
  REQUIRE(spo_ptr->getBasisSetSize() == check_basisset_size);

  SPOSetBuilderFactory::clear();
}

TEST_CASE("SPO input spline from xml LiH_msd", "[wavefunction]")
{
  // capture 3 spo input styles, the 2nd and 3rd ones will be deprecated and removed eventually.
  // the first case should be simplified using SPOSetBuilderFactory instead of WaveFunctionFactory
  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd input style 1 using sposet_collection" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\"> \
        <detlist size=\"1487\" type=\"DETS\" nca=\"0\" ncb=\"0\" nea=\"2\" neb=\"2\" nstates=\"85\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input(spo_xml_string1, "spo-up", 85, 105);

  app_log() << "-----------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd input style 1 using sposet_collection with basisset added" << std::endl;
  app_log() << "-----------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_updated = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <basisset name=\"LCAOBSet\" key=\"GTO\" transform=\"yes\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
      </basisset> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\"> \
        <detlist size=\"1487\" type=\"DETS\" nca=\"0\" ncb=\"0\" nea=\"2\" neb=\"2\" nstates=\"85\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input(spo_xml_string1_updated, "spo-up", 85, 105);

  app_log() << "------------------------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd input style 1 using sposet_collection with basisset added no transform" << std::endl;
  app_log() << "------------------------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_updated_no_transform = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <basisset name=\"LCAOBSet\" key=\"GTO\" transform=\"no\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
      </basisset> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\"> \
        <detlist size=\"1487\" type=\"DETS\" nca=\"0\" ncb=\"0\" nea=\"2\" neb=\"2\" nstates=\"85\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input(spo_xml_string1_updated_no_transform, "spo-up", 85, 105);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd input style 2 using sposet_collection" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string2 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <sposet basisset=\"LCAOBSet\" name=\"spo\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_up=\"spo\" spo_dn=\"spo\"> \
        <detlist size=\"1487\" type=\"DETS\" nca=\"0\" ncb=\"0\" nea=\"2\" neb=\"2\" nstates=\"85\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input(spo_xml_string2, "spo", 85, 105);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd input style 3 sposet inside determinantset" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string3 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <determinantset type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"85\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"85\" spindataset=\"0\"/> \
      </sposet> \
      <multideterminant optimize=\"yes\" spo_up=\"spo-up\" spo_dn=\"spo-dn\"> \
        <detlist size=\"1487\" type=\"DETS\" nca=\"0\" ncb=\"0\" nea=\"2\" neb=\"2\" nstates=\"85\" cutoff=\"1e-20\" href=\"LiH.orbs.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input(spo_xml_string3, "spo-up", 85, 105);
}
} // namespace qmcplusplus

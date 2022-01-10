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
void test_LiH_msd_xml_input(const std::string& spo_xml_string,
                            const std::string& check_sponame,
                            int check_spo_size,
                            int check_basisset_size)
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {0.0, 0.0, 3.0139239693};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Li");
  int HIdx             = ispecies.addSpecies("H");

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2, 2});
  elec_.R[0] = {0.5, 0.5, 0.5};
  elec_.R[1] = {0.1, 0.1, 1.1};
  elec_.R[2] = {-0.5, -0.5, -0.5};
  elec_.R[3] = {-0.1, -0.1, 1.5};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(ein_xml);

  SPOSet* spo_ptr(wf_factory.getSPOSet(check_sponame));
  REQUIRE(spo_ptr != nullptr);
  REQUIRE(spo_ptr->getOrbitalSetSize() == check_spo_size);
  REQUIRE(spo_ptr->getBasisSetSize() == check_basisset_size);
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

void test_LiH_msd_xml_input_with_positron(const std::string& spo_xml_string,
                                          const std::string& check_sponame,
                                          int check_spo_size,
                                          int check_basisset_size)
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1, 1});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {0.0, 0.0, 3.0139239693};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int LiIdx            = ispecies.addSpecies("Li");
  int HIdx             = ispecies.addSpecies("H");

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2, 2, 1});
  elec_.R[0] = {0.5, 0.5, 0.5};
  elec_.R[1] = {0.1, 0.1, 1.1};
  elec_.R[2] = {-0.5, -0.5, -0.5};
  elec_.R[3] = {-0.1, -0.1, 1.5};
  elec_.R[4] = {0.0, -1.0, 2};

  SpeciesSet& tspecies             = elec_.getSpeciesSet();
  int upIdx                        = tspecies.addSpecies("u");
  int downIdx                      = tspecies.addSpecies("d");
  int positronIdx                  = tspecies.addSpecies("pos");
  int massIdx                      = tspecies.addAttribute("mass");
  int chargeIdx                    = tspecies.addAttribute("charge");
  tspecies(massIdx, upIdx)         = 1.0;
  tspecies(massIdx, downIdx)       = 1.0;
  tspecies(massIdx, positronIdx)   = 1.0;
  tspecies(chargeIdx, upIdx)       = -1.0;
  tspecies(chargeIdx, downIdx)     = -1.0;
  tspecies(chargeIdx, positronIdx) = 1.0;

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory("psi0", elec_, ptcl.getPool(), c);
  wf_factory.put(ein_xml);

  SPOSet* spo_ptr(wf_factory.getSPOSet(check_sponame));
  REQUIRE(spo_ptr != nullptr);
  REQUIRE(spo_ptr->getOrbitalSetSize() == check_spo_size);
  REQUIRE(spo_ptr->getBasisSetSize() == check_basisset_size);
}

TEST_CASE("SPO input spline from xml LiH_msd arbitrary species", "[wavefunction]")
{
  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd with positron xml input style" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-ps\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_0=\"spo-up\" spo_1=\"spo-dn\" spo_2=\"spo-ps\"> \
        <detlist size=\"2\" type=\"DETS\" nc0=\"0\" nc1=\"0\" nc2=\"0\" ne0=\"2\" ne1=\"2\" ne2=\"1\" nstates=\"5\" cutoff=\"1e-20\"> \
          <ci id=\"CIcoeff_0\" coeff=\"0.7071\" qchem_coeff=\"0.7071\" occ0=\"11000\" occ1=\"11000\" occ2=\"10000\"/> \
          <ci id=\"CIcoeff_1\" coeff=\"-0.7071\" qchem_coeff=\"-0.7071\" occ0=\"10100\" occ1=\"11000\" occ2=\"00100\" /> \
        </detlist> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input_with_positron(spo_xml_string1, "spo-ps", 5, 105);
}
TEST_CASE("SPO input spline from h5 LiH_msd arbitrary species", "[wavefunction]")
{
  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "LiH_msd with positron h5 input style" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = "<wavefunction name=\"psi0\" target=\"e\"> \
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\" href=\"LiH.orbs.h5\"> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-up\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-dn\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
      <sposet basisset=\"LCAOBSet\" name=\"spo-ps\" size=\"5\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"5\" spindataset=\"0\"/> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <multideterminant optimize=\"yes\" spo_0=\"spo-up\" spo_1=\"spo-dn\" spo_2=\"spo-ps\"> \
        <detlist size=\"2\" type=\"DETS\" cutoff=\"1e-20\" href=\"LiH.Multidet.h5\"/> \
      </multideterminant> \
    </determinantset> \
</wavefunction> \
";
  test_LiH_msd_xml_input_with_positron(spo_xml_string1, "spo-ps", 5, 105);
}
} // namespace qmcplusplus

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
void test_He_sto3g_xml_input(const std::string& spo_xml_string)
{
  Communicate* c = OHMMS::Controller;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& elec(*elec_uptr);
  ParticleSet& ions(*ions_uptr);

  elec.setName("e");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.create({1, 1});
  elec.R[0] = 0.0;

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  ions.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create({1});
  ions.R[0]            = 0.0;
  SpeciesSet& ispecies = ions.getSpeciesSet();
  int heIdx            = ispecies.addSpecies("He");
  ions.update();

  elec.addTable(ions);
  elec.update();

  Libxml2Document doc;
  bool okay = doc.parseFromString(spo_xml_string);
  REQUIRE(okay);

  xmlNodePtr ein_xml = doc.getRoot();

  WaveFunctionFactory wf_factory(elec, ptcl.getPool(), c);
  auto twf_ptr = wf_factory.buildTWF(ein_xml);

  std::unique_ptr<SPOSet> sposet(twf_ptr->getSPOSet("spo").makeClone());

  SPOSet::ValueVector values;
  SPOSet::GradVector dpsi;
  SPOSet::ValueVector d2psi;
  values.resize(1);
  dpsi.resize(1);
  d2psi.resize(1);

  // Call makeMove to compute the distances
  ParticleSet::SingleParticlePos newpos(0.0001, 0.0, 0.0);
  elec.makeMove(0, newpos);

  sposet->evaluateValue(elec, 0, values);

  // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
  REQUIRE(values[0] == ValueApprox(0.9996037001));

  sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);

  // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
  REQUIRE(values[0] == ValueApprox(0.9996037001));
  REQUIRE(dpsi[0][0] == ValueApprox(-0.0006678035459));
  REQUIRE(dpsi[0][1] == ValueApprox(0.0));
  REQUIRE(dpsi[0][2] == ValueApprox(0.0));
  REQUIRE(d2psi[0] == ValueApprox(-20.03410564));


  ParticleSet::SingleParticlePos disp(1.0, 0.0, 0.0);
  elec.makeMove(0, disp);

  sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
  // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
  REQUIRE(values[0] == ValueApprox(0.2315567641));
  REQUIRE(dpsi[0][0] == ValueApprox(-0.3805431885));
  REQUIRE(dpsi[0][1] == ValueApprox(0.0));
  REQUIRE(dpsi[0][2] == ValueApprox(0.0));
  REQUIRE(d2psi[0] == ValueApprox(-0.2618497452));
}

TEST_CASE("SPO input spline from xml He_sto3g", "[wavefunction]")
{
  // capture 3 spo input styles, the 2nd and 3rd ones will be deprecated and removed eventually.
  // the first case should be simplified using SPOSetBuilderFactory instead of WaveFunctionFactory
  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "He_sto3g input style 1 using sposet_collection" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1 = "<wavefunction name=\"psi0\" target=\"e\">\
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\"> \
      <basisset name=\"LCAOBSet\"> \
        <atomicBasisSet name=\"Gaussian\" angular=\"cartesian\" type=\"Gaussian\" elementType=\"He\" normalized=\"no\"> \
          <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
          <basisGroup rid=\"He00\" n=\"0\" l=\"0\" type=\"Gaussian\"> \
            <radfunc exponent=\"6.362421400000e+00\" contraction=\"1.543289672950e-01\"/> \
            <radfunc exponent=\"1.158923000000e+00\" contraction=\"5.353281422820e-01\"/> \
            <radfunc exponent=\"3.136498000000e-01\" contraction=\"4.446345421850e-01\"/> \
          </basisGroup> \
        </atomicBasisSet> \
      </basisset> \
      <sposet name=\"spo\" size=\"1\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"1\" id=\"updetC\"> \
          1.00000000000000e+00 \
        </coefficient> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <slaterdeterminant> \
        <determinant name=\"det_up\" sposet=\"spo\"/> \
        <determinant name=\"det_dn\" sposet=\"spo\"/> \
      </slaterdeterminant> \
    </determinantset> \
</wavefunction> \
";
  test_He_sto3g_xml_input(spo_xml_string1);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "He_sto3g input style 1 using sposet_collection updated" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string1_updated = "<wavefunction name=\"psi0\" target=\"e\">\
    <sposet_collection type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" cuspCorrection=\"no\"> \
      <basisset name=\"LCAOBSet\" transform=\"yes\"> \
        <atomicBasisSet name=\"Gaussian\" angular=\"cartesian\" type=\"Gaussian\" elementType=\"He\" normalized=\"no\"> \
          <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
          <basisGroup rid=\"He00\" n=\"0\" l=\"0\" type=\"Gaussian\"> \
            <radfunc exponent=\"6.362421400000e+00\" contraction=\"1.543289672950e-01\"/> \
            <radfunc exponent=\"1.158923000000e+00\" contraction=\"5.353281422820e-01\"/> \
            <radfunc exponent=\"3.136498000000e-01\" contraction=\"4.446345421850e-01\"/> \
          </basisGroup> \
        </atomicBasisSet> \
      </basisset> \
      <sposet name=\"spo\" basisset=\"LCAOBSet\" size=\"1\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"1\" id=\"updetC\"> \
          1.00000000000000e+00 \
        </coefficient> \
      </sposet> \
    </sposet_collection> \
    <determinantset> \
      <slaterdeterminant> \
        <determinant name=\"det_up\" sposet=\"spo\"/> \
        <determinant name=\"det_dn\" sposet=\"spo\"/> \
      </slaterdeterminant> \
    </determinantset> \
</wavefunction> \
";
  test_He_sto3g_xml_input(spo_xml_string1_updated);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "He_sto3g input style 2 sposet inside determinantset" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string2 = "<wavefunction name=\"psi0\" target=\"e\">\
    <determinantset type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\"> \
      <basisset name=\"LCAOBSet\"> \
        <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
        <atomicBasisSet name=\"Gaussian\" angular=\"cartesian\" type=\"Gaussian\" elementType=\"He\" normalized=\"no\"> \
          <basisGroup rid=\"He00\" n=\"0\" l=\"0\" type=\"Gaussian\"> \
            <radfunc exponent=\"6.362421400000e+00\" contraction=\"1.543289672950e-01\"/> \
            <radfunc exponent=\"1.158923000000e+00\" contraction=\"5.353281422820e-01\"/> \
            <radfunc exponent=\"3.136498000000e-01\" contraction=\"4.446345421850e-01\"/> \
          </basisGroup> \
        </atomicBasisSet> \
      </basisset> \
      <sposet name=\"spo\" size=\"1\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"1\" id=\"updetC\"> \
          1.00000000000000e+00 \
        </coefficient> \
      </sposet> \
      <sposet name=\"spo-down\" size=\"1\"> \
        <occupation mode=\"ground\"/> \
        <coefficient size=\"1\" id=\"downdetC\"> \
          1.00000000000000e+00 \
        </coefficient> \
      </sposet> \
      <slaterdeterminant> \
        <determinant name=\"det_up\" sposet=\"spo\"/> \
        <determinant name=\"det_dn\" sposet=\"spo-down\"/> \
      </slaterdeterminant> \
    </determinantset> \
</wavefunction> \
";
  test_He_sto3g_xml_input(spo_xml_string2);

  app_log() << "-------------------------------------------------------------" << std::endl;
  app_log() << "He_sto3g input style 3 sposet inside determinantset" << std::endl;
  app_log() << "-------------------------------------------------------------" << std::endl;
  const char* spo_xml_string3 = "<wavefunction name=\"psi0\" target=\"e\">\
    <determinantset type=\"MolecularOrbital\" name=\"LCAOBSet\" source=\"ion0\" transform=\"yes\" cuspCorrection=\"no\"> \
      <basisset name=\"LCAOBSet\"> \
        <atomicBasisSet name=\"Gaussian\" angular=\"cartesian\" type=\"Gaussian\" elementType=\"He\" normalized=\"no\"> \
          <grid type=\"log\" ri=\"1.e-6\" rf=\"1.e2\" npts=\"1001\"/> \
          <basisGroup rid=\"He00\" n=\"0\" l=\"0\" type=\"Gaussian\"> \
            <radfunc exponent=\"6.362421400000e+00\" contraction=\"1.543289672950e-01\"/> \
            <radfunc exponent=\"1.158923000000e+00\" contraction=\"5.353281422820e-01\"/> \
            <radfunc exponent=\"3.136498000000e-01\" contraction=\"4.446345421850e-01\"/> \
          </basisGroup> \
        </atomicBasisSet> \
      </basisset> \
      <slaterdeterminant> \
        <determinant name=\"spo\" size=\"1\"> \
          <occupation mode=\"ground\"/> \
          <coefficient size=\"1\" id=\"updetC\"> \
            1.00000000000000e+00 \
          </coefficient> \
        </determinant> \
        <determinant name=\"spo-down\" size=\"1\"> \
          <occupation mode=\"ground\"/> \
          <coefficient size=\"1\" id=\"downdetC\"> \
            1.00000000000000e+00 \
          </coefficient> \
        </determinant> \
      </slaterdeterminant> \
    </determinantset> \
</wavefunction> \
";
  test_He_sto3g_xml_input(spo_xml_string3);
}
} // namespace qmcplusplus

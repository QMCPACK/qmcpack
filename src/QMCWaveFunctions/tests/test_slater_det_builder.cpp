//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "OhmmsData/Libxml2Doc.h"
#include "WaveFunctionFactory.h"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{

TEST_CASE("SlaterDetBuilder throws", "[wavefunction]")
{
  Communicate& c(*OHMMS::Controller);
  RuntimeOptions runtime_options;

  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({1});
  ions_.R[0] = {0.0, 0.0, 0.0};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  ispecies.addSpecies("H");

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2, 2});
  elec_.R[0] = {0.1, 0.0, 0.0};
  elec_.R[1] = {0.2, 0.0, 0.0};
  elec_.R[2] = {0.3, 0.0, 0.0};
  elec_.R[3] = {0.4, 0.0, 0.0};
  SpeciesSet& tspecies = elec_.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  WaveFunctionFactory wf_factory(elec_, ptcl.getPool(), c);

  const char* base_spo_xml = R"(
    <sposet_collection type="MolecularOrbital" name="LCAOBSet" source="ion0">
      <basisset keyword="GTO" name="LCAOBSet">
        <atomicBasisSet name="Gaussian" angular="spherical" type="Gaussian" elementType="H" normalized="no">
          <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
          <basisGroup rid="H00" n="0" l="0" type="Gaussian">
            <radfunc exponent="1.0" contraction="1.0"/>
          </basisGroup>
          <basisGroup rid="H01" n="0" l="0" type="Gaussian">
            <radfunc exponent="0.5" contraction="1.0"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
      <sposet name="spo-up" size="2">
        <occupation mode="ground"/>
        <coefficient size="2" id="updetC">
          1.0 0.0
          0.0 1.0
        </coefficient>
      </sposet>
      <sposet name="spo-dn" size="2">
        <occupation mode="ground"/>
        <coefficient size="2" id="downdetC">
          1.0 0.0
          0.0 1.0
        </coefficient>
      </sposet>
      <sposet name="spo-small" size="1">
        <occupation mode="ground"/>
        <coefficient size="1" id="smalldetC">
          1.0 0.0
        </coefficient>
      </sposet>
    </sposet_collection>
  )";

  SECTION("Only single backflow block allowed")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <backflow/>
          <backflow/>
          <slaterdeterminant>
            <determinant sposet="spo-up"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("Only a single backflow block is allowed in the xml. Please collect all transformations into a single block."));
  }

  SECTION("Multiple slaterdeterminant entries not allowed")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <slaterdeterminant>
            <determinant sposet="spo-up"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
          <slaterdeterminant>
            <determinant sposet="spo-up"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("Only one slaterdeterminant or multideterminant entry allowed in XML"));
  }

  SECTION("Backflow not supported with MultiSlaterDeterminant table method")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <backflow/>
          <multideterminant optimize="yes" spo_up="spo-up" spo_dn="spo-dn">
            <detlist size="1" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="2" cutoff="1e-20">
              <ci id="CIcoeff_0" coeff="1.0" qchem_coeff="1.0" occ0="11" occ1="11"/>
            </detlist>
          </multideterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("Backflow is not supported by Multi-Slater determinants using the table method!"));
  }

  SECTION("MultiSlaterDeterminant missing SPOSet")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <multideterminant optimize="yes" spo_up="nonexistent_spo" spo_dn="spo-dn">
            <detlist size="1" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="2" cutoff="1e-20">
              <ci id="CIcoeff_0" coeff="1.0" qchem_coeff="1.0" occ0="11" occ1="11"/>
            </detlist>
          </multideterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("In SlaterDetBuilder: SPOSet \"nonexistent_spo\" is not found. Expected for MultiSlaterDeterminant.\n"));
  }

  SECTION("Single Slater determinant missing SPOSet")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <slaterdeterminant>
            <determinant sposet="nonexistent_spo"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("A sposet named \"nonexistent_spo\" cannot be found for constructing a Slater determinant! Please check the xml input file!\n"));
  }

  SECTION("Determinant SPOSet has fewer orbitals than particles")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <slaterdeterminant>
            <determinant sposet="spo-small"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("The SPOSet spo-small only has 1 orbitals but this determinant needs at least 2\n"));
  }

  SECTION("Delay rank invalid")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <slaterdeterminant delay_rank="-1">
            <determinant sposet="spo-up"/>
            <determinant sposet="spo-dn"/>
          </slaterdeterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("SlaterDetBuilder::putDeterminant delay_rank must be positive and no larger than the electron count within a determinant!\nAcceptable value [1,2], user input -1"));
  }

  SECTION("MSD occupancy particle count mismatch")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <multideterminant optimize="yes" spo_up="spo-up" spo_dn="spo-dn">
            <detlist size="1" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="2" cutoff="1e-20">
              <ci id="CIcoeff_0" coeff="1.0" qchem_coeff="1.0" occ0="10" occ1="11"/>
            </detlist>
          </multideterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("SlaterDetBuilder::createMSDFast Issues found in the particle group 0 :\nThe unique determinant configuration 0 contains 1 occupied orbitals not matching 2 particles.\n"));
  }

#ifndef QMC_COMPLEX
  SECTION("Complex CI expansion coefficients without QMC_COMPLEX")
  {
    std::string wf_xml = std::string("<wavefunction name=\"psi0\" target=\"elec\">") + base_spo_xml +
        R"(
        <determinantset>
          <multideterminant optimize="yes" spo_up="spo-up" spo_dn="spo-dn">
            <detlist size="1" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="2" cutoff="1e-20">
              <ci id="CIcoeff_0" coeff="1.0" coeff_imag="0.5" qchem_coeff="1.0" occ0="11" occ1="11"/>
            </detlist>
          </multideterminant>
        </determinantset>
      </wavefunction>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(wf_factory.buildTWF(doc.getRoot(), runtime_options), UniformCommunicateError,
                         Catch::Matchers::Message("SlaterDetBuilder::readDetList. Build with QMC_COMPLEX if using complex CI expansion coefficients."));
  }
#endif
}

} // namespace qmcplusplus

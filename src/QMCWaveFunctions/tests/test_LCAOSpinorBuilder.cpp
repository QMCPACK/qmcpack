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
#include "hdf/hdf_archive.h"
#include "LCAO/LCAOSpinorBuilder.h"
#include <filesystem>

namespace qmcplusplus
{

class TestLCAOSpinorBuilder : public LCAOSpinorBuilder
{
public:
  using LCAOSpinorBuilder::LCAOSpinorBuilder;
  void setH5Path(const std::string& p) { h5_path = p; }
};

TEST_CASE("LCAOSpinorBuilder throws", "[wavefunction]")
{
  Communicate& c(*OHMMS::Controller);
  ParticleSetPool ptcl = ParticleSetPool(c);
  auto ions_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr       = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create({1});
  SpeciesSet& ispecies = ions.getSpeciesSet();
  ispecies.addSpecies("H");

  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.create({1});

  SECTION("LCAOSpinorBuilder only works with href")
  {
    const char* wf_xml = R"(
      <tmp>
        <basisset keyword="GTO" name="LCAOBSet">
          <atomicBasisSet angular="spherical" elementType="H" type="Gaussian">
            <grid npts="1001" rf="1.e2" ri="1.e-6" type="log" />
            <basisGroup l="0" n="0" rid="H00" type="Gaussian">
              <radfunc contraction="1.0" exponent="1.0" />
            </basisGroup>
          </atomicBasisSet>
        </basisset>
      </tmp>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(wf_xml));
    CHECK_THROWS_MATCHES(LCAOSpinorBuilder(elec, ions, c, doc.getRoot()), UniformCommunicateError,
                         Catch::Matchers::Message("LCAOSpinorBuilder only works with href"));
  }

  SECTION("createSPOSetFromXML and loadMO exceptions")
  {
    const char* builder_xml = R"(
      <sposet_builder name="spinorbuilder" type="molecularorbital" href="lcao_spinor.h5" source="ion0" precision="float">
        <basisset transform="yes"/>
      </sposet_builder>
    )";
    Libxml2Document doc;
    REQUIRE(doc.parseFromString(builder_xml));
    TestLCAOSpinorBuilder builder(elec, ions, c, doc.getRoot());

    SECTION("missing basisset")
    {
      Libxml2Document spo_doc;
      REQUIRE(spo_doc.parseFromString(R"(<sposet name="myspo" basisset="nonexistent" size="1"/>)"));
      CHECK_THROWS_MATCHES(builder.createSPOSetFromXML(spo_doc.getRoot()), UniformCommunicateError,
                           Catch::Matchers::Message("basisset \"nonexistent\" cannot be found\n"));
    }

    SECTION("missing or incorrect path to H5 file")
    {
      builder.setH5Path("nonexistent_spinor.h5");
      Libxml2Document spo_doc;
      REQUIRE(spo_doc.parseFromString(R"(<sposet name="myspo" size="1"/>)"));
      CHECK_THROWS_MATCHES(builder.createSPOSetFromXML(spo_doc.getRoot()), UniformCommunicateError,
                           Catch::Matchers::Message("LCAOSpinorBuilder::loadMO missing or incorrect path to H5 file."));
    }

    SECTION("lcao spinors not implemented in PBC")
    {
      const std::string pbc_h5 = "test_spinor_pbc.h5";
      if (c.rank() == 0)
      {
        hdf_archive hout(c);
        hout.create(pbc_h5);
        hout.push("PBC", true);
        bool pbc = true;
        hout.write(pbc, "PBC");
        hout.pop();
        hout.close();
      }
      c.barrier();

      builder.setH5Path(pbc_h5);
      Libxml2Document spo_doc;
      REQUIRE(spo_doc.parseFromString(R"(<sposet name="myspo" size="1"/>)"));
      CHECK_THROWS_MATCHES(builder.createSPOSetFromXML(spo_doc.getRoot()), UniformCommunicateError,
                           Catch::Matchers::Message("LCAOSpinorBuilder::loadMO lcao spinors not implemented in PBC"));

      c.barrier();
      if (c.rank() == 0)
        std::filesystem::remove(pbc_h5);
    }
  }
}

} // namespace qmcplusplus

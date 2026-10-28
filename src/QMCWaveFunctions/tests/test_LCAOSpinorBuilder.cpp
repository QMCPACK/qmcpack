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

#ifdef QMC_COMPLEX
#include "LCAO/LCAOSpinorBuilder.h"

namespace qmcplusplus
{

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
  ispecies.addSpecies("C");

  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.create({1});

  SECTION("LCAOSpinorBuilder only works with href")
  {
    const char* wf_xml = R"(
      <tmp>
        <basisset keyword="GTO" name="LCAOBSet">
          <atomicBasisSet angular="spherical" elementType="C" type="Gaussian">
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
}

} // namespace qmcplusplus

#endif

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Particle/LongRange/StructFact.h"
#include "Particle/tests/MinimalParticlePool.h"

#include <cmath>
#include <limits>

namespace qmcplusplus
{
TEST_CASE("MinimalParticlePool factory reset computes ionic structure factors", "[qmcapp][particle]")
{
  auto pool = MinimalParticlePool::make_diamondC_1x1x1(*OHMMS::Controller);
  const ParticleSet* ions_ptr = pool.getParticleSet("ion");
  REQUIRE(ions_ptr != nullptr);
  const ParticleSet& ions = *ions_ptr;
  REQUIRE(ions.hasSK());
  const auto& sk = ions.getSK();
  const auto& kpts = ions.getSimulationCell().getKLists().getKptsCart();

  REQUIRE(ions.getTotalNum() == 2);
  REQUIRE(sk.rhok_r.size(0) == ions.groups());
  REQUIRE(sk.rhok_r.size(1) == kpts.size());
  REQUIRE(sk.rhok_i.size(0) == ions.groups());
  REQUIRE(sk.rhok_i.size(1) == kpts.size());
  const double tolerance = 100 * std::numeric_limits<StructFact::RealType>::epsilon();

  // The first nonzero k vector avoids passing on an uninitialized zero-k density.
  bool checked_nonzero_k = false;
  for (std::size_t ik = 0; ik < kpts.size(); ++ik)
  {
    const auto& k = kpts[ik];
    if (dot(k, k) == 0.0)
      continue;

    for (int species = 0; species < sk.rhok_r.size(0); ++species)
    {
      double expected_r = 0.0;
      double expected_i = 0.0;
      for (int particle = 0; particle < ions.getTotalNum(); ++particle)
        if (ions.GroupID[particle] == species)
        {
          const double phase = dot(k, ions.R[particle]);
          expected_r += std::cos(phase);
          expected_i += std::sin(phase);
        }
      CHECK(sk.rhok_r(species, ik) == Approx(expected_r).margin(tolerance));
      CHECK(sk.rhok_i(species, ik) == Approx(expected_i).margin(tolerance));
    }
    checked_nonzero_k = true;
    break;
  }
  REQUIRE(checked_nonzero_k);
}
} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "PerParticleHamiltonianLogger.h"

#include "Utilities/StdRandom.h"

namespace qmcplusplus
{

using QMCT = QMCTraits;
using Real = QMCT::RealType;

class MultiWalkerTalker
{
private:
  std::vector<ListenerVector<Real>> listener_vectors_;
  const std::string name_{"Talker"};
  int walkers_;

public:
  MultiWalkerTalker(const std::string& name, int walkers) : name_(name), walkers_(walkers){};
  void registerVector(ListenerVector<Real>& listener_vector) { listener_vectors_.push_back(listener_vector); }
  void reportVector()
  {
    Vector<Real> vec_part(4);

    for (auto& listener : listener_vectors_)
      for (int iw = 0; iw < walkers_; ++iw)
      {
        std::iota(vec_part.begin(), vec_part.end(), iw * walkers_);

        listener.report(iw, name_, vec_part);
      }
  }
};


TEST_CASE("PerParticleHamiltonianLogger_sum", "[estimators]")
{
  std::string_view xml{R"XML(
<PerParticleHamiltonianLogger to_stdout="true" validate_per_particle_sum="yes"/>
)XML"};

  Libxml2Document doc;
  bool okay = doc.parseFromString(xml);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  PerParticleHamiltonianLoggerInput pphli(node);

  CHECK(pphli.get_to_stdout());
  CHECK(pphli.get_validate_per_particle_sum());
  
  int rank = 0;
  PerParticleHamiltonianLogger rank_logger(std::move(pphli), rank);

  int ncrowds = 3;

  UPtrVector<OperatorEstBase> crowd_loggers;
  for (int ic = 0; ic < ncrowds; ++ic)
  {
    crowd_loggers.emplace_back(rank_logger.spawnCrowdClone());
  }

  int nwalkers = 3;

  const SimulationCell simulation_cell;
  std::vector<OperatorEstBase::MCPWalker> walkers;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(2);

  std::vector<ParticleSet> psets;
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    psets.emplace_back(simulation_cell);
    ParticleSet& pset = psets.back();
    pset.create({2});
    pset.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
    pset.R[1] = ParticleSet::PosType(0.68658058, 0.68658058, 0.68658058);
  }

  std::vector<TrialWaveFunction> wfns;

  auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);
  auto ref_wfns    = makeRefVector<TrialWaveFunction>(wfns);

  std::vector<MultiWalkerTalker> multi_walker_talkers{{"Talker1", nwalkers},
                                                      {"Talker2", nwalkers},
                                                      {"Talker3", nwalkers}};
  for (auto& crowd_oeb : crowd_loggers)
    for (auto& mwt : multi_walker_talkers)
    {
      auto& crowd_logger = dynamic_cast<PerParticleHamiltonianLogger&>(*crowd_oeb);
      ListenerVector<Real> listener("whatever", crowd_logger.getLogger());
      mwt.registerVector(listener);
    }

  for (auto& mwt : multi_walker_talkers)
    mwt.reportVector();

  RandomGenerator rng;

  int crowd_id = 0;
  for (auto& crowd_oeb : crowd_loggers)
  {
    crowd_oeb->accumulate(ref_walkers, ref_psets, ref_wfns, rng, crowd_id++);
  }
}

} // namespace qmcplusplus

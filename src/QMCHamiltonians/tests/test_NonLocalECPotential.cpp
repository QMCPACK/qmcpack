//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>

#include "Configuration.h"
#include "Numerics/Quadrature.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "TestListenerFunction.h"
#include "Utilities/StdRandom.h"
#include "Utilities/StlPrettyPrint.hpp"
#include "Utilities/RuntimeOptions.h"

namespace qmcplusplus
{
namespace testing
{
/** class to violate access control because evaluation of NonLocalECPotential uses RNG
 *  which we may not be able to control.
 */
class TestNonLocalECPotential
{
  using Real = QMCTraits::RealType;

public:
  static void copyGridUnrotatedForTest(NonLocalECPotential& nl_ecp)
  {
    nl_ecp.PPset[0]->rrotsgrid_m = nl_ecp.PPset[0]->sgridxyz_m;
  }

  static bool didGridChange(NonLocalECPotential& nl_ecp)
  {
    return nl_ecp.PPset[0]->rrotsgrid_m != nl_ecp.PPset[0]->sgridxyz_m;
  }

  static void evaluateImpl(NonLocalECPotential& nl_ecp,
                           TrialWaveFunction& psi,
                           ParticleSet& P,
                           bool compute_txy_all,
                           bool keep_grid)
  {
    nl_ecp.evaluateImpl(psi, P, compute_txy_all, keep_grid);
  }

  static void mw_evaluateImpl(NonLocalECPotential& nl_ecp,
                              const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& twf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              bool compute_txy_all,
                              const std::optional<ListenerOption<Real>> listener_opt,
                              bool keep_grid)
  {
    nl_ecp.mw_evaluateImpl(o_list, twf_list, p_list, compute_txy_all, listener_opt, keep_grid);
  }
};

} // namespace testing

TEST_CASE("NonLocalECPotential", "[hamiltonian]")
{
  using Real         = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::FullPrecRealType;
  using Position     = QMCTraits::PosType;
  using testing::getParticularListener;

  Lattice lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(20.0);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);

  ParticleSet ions(simulation_cell);

  ions.setName("ion");
  ions.create({2});
  ions.R[0] = {0.0, 1.0, 0.0};
  ions.R[1] = {0.0, -1.0, 0.0};

  SpeciesSet& ion_species                         = ions.getSpeciesSet();
  int index_species                               = ion_species.addSpecies("Na");
  int index_charge                                = ion_species.addAttribute("charge");
  int index_atomic_number                         = ion_species.addAttribute("atomic_number");
  ion_species(index_charge, index_species)        = 1;
  ion_species(index_atomic_number, index_species) = 1;
  ions.createSK();
  ions.resetGroups(); // test_ecp.cpp claims this is needed
  ions.update();      // elsewhere its implied this is needed

  ParticleSet ions2(ions);
  ions2.update();

  ParticleSet elec(simulation_cell);
  elec.setName("elec");
  elec.create({2, 1});
  elec.R[0] = {0.4, 0.0, 0.0};
  elec.R[1] = {1.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  int dnIdx                  = tspecies.addSpecies("d");
  chargeIdx                  = tspecies.addAttribute("charge");
  massIdx                    = tspecies.addAttribute("mass");
  tspecies(chargeIdx, dnIdx) = -1;
  tspecies(massIdx, dnIdx)   = 1.0;

  elec.createSK();
  elec.resetGroups();
  elec.addTable(ions);
  elec.update();

  ParticleSet elec2(elec);
  elec2.update();

  RefVectorWithLeader<ParticleSet> p_list(elec, {elec, elec2});

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  TrialWaveFunction psi2(runtime_options);
  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, psi2});

  NonLocalECPotential nl_ecp(ions, elec, false /*use_DLA*/, false /*use_VP*/);

  int num_walkers = 2;
  int max_values  = 10;
  Matrix<Real> local_pots(num_walkers, max_values);
  Matrix<Real> local_pots2(num_walkers, max_values);

  ResourceCollection pset_res("test_pset_res");
  elec.createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  std::vector<ListenerVector<Real>> listeners;
  listeners.emplace_back("nonlocalpotential", getParticularListener(local_pots));
  listeners.emplace_back("nonlocalpotential", getParticularListener(local_pots2));

  Matrix<Real> ion_pots(num_walkers, max_values);
  Matrix<Real> ion_pots2(num_walkers, max_values);

  std::vector<ListenerVector<Real>> ion_listeners;
  ion_listeners.emplace_back("nonlocalpotential", getParticularListener(ion_pots));
  ion_listeners.emplace_back("nonlocalpotential", getParticularListener(ion_pots2));


  // This took some time to sort out from the multistage mess of put and clones
  // but this accomplishes in a straight forward way what I interpret to be done by that code.
  Communicate* comm = OHMMS::Controller;
  ECPComponentBuilder ecp_comp_builder("test_read_ecp", comm, 4, 1);

  bool okay = ecp_comp_builder.read_pp_file("Na.BFD.xml");
  REQUIRE(okay);
  UPtr<NonLocalECPComponent> nl_ecp_comp = std::move(ecp_comp_builder.pp_nonloc);
  nl_ecp.addComponent(0, std::move(nl_ecp_comp));
  UPtr<OperatorBase> nl_ecp2_ptr = nl_ecp.makeClone(elec2, psi2);
  auto& nl_ecp2                  = dynamic_cast<NonLocalECPotential&>(*nl_ecp2_ptr);

  StdRandom<FullPrecReal> rng(10101);
  StdRandom<FullPrecReal> rng2(10201);
  nl_ecp.setRandomGenerator(&rng);
  nl_ecp2.setRandomGenerator(&rng2);

  RefVectorWithLeader<OperatorBase> o_list(nl_ecp, {nl_ecp, nl_ecp2});
  ResourceCollection nl_ecp_res("test_nl_ecp_res");
  nl_ecp.createResource(nl_ecp_res);
  ResourceCollectionTeamLock<OperatorBase> nl_ecp_lock(nl_ecp_res, o_list);

  // Despite what test_ecp.cpp says this does not need to be done.
  // I think because of the pp
  testing::TestNonLocalECPotential::copyGridUnrotatedForTest(nl_ecp);
  testing::TestNonLocalECPotential::copyGridUnrotatedForTest(nl_ecp2);

  CHECK(!testing::TestNonLocalECPotential::didGridChange(nl_ecp));

  ListenerOption<Real> listener_opt{listeners, ion_listeners};
  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, listener_opt, true);

  // for now we'll check against the single particle API
  testing::TestNonLocalECPotential::evaluateImpl(nl_ecp, psi, elec, false, true);
  const auto value = nl_ecp.getValue();

  double total_localpots = std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0);
  total_localpots += std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0);
  CHECK(total_localpots == Approx(value));
  double total_localpots2 = std::accumulate(local_pots[1], local_pots[1] + local_pots.cols(), 0.0);
  total_localpots2 += std::accumulate(ion_pots[1], ion_pots[1] + ion_pots.cols(), 0.0);
  CHECK(total_localpots2 == Approx(value));

  CHECK(!testing::TestNonLocalECPotential::didGridChange(nl_ecp));

  elec.R[0] = {0.5, 0.0, 2.0};
  elec.update();

  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, listener_opt, true);

  CHECK(!testing::TestNonLocalECPotential::didGridChange(nl_ecp));
  auto value_moved = o_list[0].evaluateDeterministic(psi, elec);

  total_localpots = std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0);
  total_localpots += std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0);
  CHECK(total_localpots == Approx(value_moved));
  // check the second walker which will be unchanged.
  total_localpots2 = std::accumulate(local_pots[1], local_pots[1] + local_pots.cols(), 0.0);
  total_localpots2 += std::accumulate(ion_pots[1], ion_pots[1] + ion_pots.cols(), 0.0);
  CHECK(total_localpots2 == Approx(value));

  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, listener_opt, false);
  auto value3     = o_list[0].evaluateDeterministic(twf_list[0], p_list[0]);
  total_localpots = std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0);
  total_localpots += std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0);

  CHECK(total_localpots == Approx(value3));
}

namespace
{
struct TmoveV1Result
{
  std::vector<int> accepts;
  std::vector<QMCTraits::PosType> final_R;
};

/** run a two-walker v1 T-move sweep with fixed quadrature grids and
 *  per-walker seeded RNGs; batched and serial paths must agree walker by
 *  walker in accepted counts and final electron positions.
 */
TmoveV1Result runTmoveV1(bool batched, bool use_VP)
{
  using FullPrecReal = QMCTraits::FullPrecRealType;

  Lattice lattice;
  lattice.BoxBConds = true;
  lattice.R.diagonal(20.0);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();
  const SimulationCell simulation_cell(lattice);

  ParticleSet ions(simulation_cell);
  ions.setName("ion");
  ions.create({2});
  ions.R[0]                                       = {0.0, 1.0, 0.0};
  ions.R[1]                                       = {0.0, -1.0, 0.0};
  SpeciesSet& ion_species                         = ions.getSpeciesSet();
  int index_species                               = ion_species.addSpecies("Na");
  int index_charge                                = ion_species.addAttribute("charge");
  int index_atomic_number                         = ion_species.addAttribute("atomic_number");
  ion_species(index_charge, index_species)        = 1;
  ion_species(index_atomic_number, index_species) = 1;
  ions.createSK();
  ions.resetGroups();
  ions.update();

  ParticleSet elec(simulation_cell);
  elec.setName("elec");
  elec.create({2, 1});
  elec.R[0]                  = {0.4, 0.0, 0.0};
  elec.R[1]                  = {1.0, 0.0, 0.0};
  elec.R[2]                  = {-0.4, 0.6, -0.3};
  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;
  int dnIdx                  = tspecies.addSpecies("d");
  chargeIdx                  = tspecies.addAttribute("charge");
  massIdx                    = tspecies.addAttribute("mass");
  tspecies(chargeIdx, dnIdx) = -1;
  tspecies(massIdx, dnIdx)   = 1.0;
  elec.createSK();
  elec.resetGroups();
  elec.addTable(ions);
  elec.update();

  ParticleSet elec2(elec);
  elec2.R[1] = {0.8, 0.3, 0.1};
  elec2.update();

  RefVectorWithLeader<ParticleSet> p_list(elec, {elec, elec2});
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  TrialWaveFunction psi2(runtime_options);
  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, psi2});

  NonLocalECPotential nl_ecp(ions, elec, false /*enable_DLA*/, use_VP);
  Communicate* comm = OHMMS::Controller;
  ECPComponentBuilder ecp_comp_builder("test_read_ecp", comm, 4, 1);
  bool okay = ecp_comp_builder.read_pp_file("Na.BFD.xml");
  REQUIRE(okay);
  UPtr<NonLocalECPComponent> nl_ecp_comp = std::move(ecp_comp_builder.pp_nonloc);
  nl_ecp.addComponent(0, std::move(nl_ecp_comp));
  UPtr<OperatorBase> nl_ecp2_ptr = nl_ecp.makeClone(elec2, psi2);
  auto& nl_ecp2                  = dynamic_cast<NonLocalECPotential&>(*nl_ecp2_ptr);

  StdRandom<FullPrecReal> rng(10101);
  StdRandom<FullPrecReal> rng2(10201);
  nl_ecp.setRandomGenerator(&rng);
  nl_ecp2.setRandomGenerator(&rng2);

  RefVectorWithLeader<OperatorBase> o_list(nl_ecp, {nl_ecp, nl_ecp2});
  ResourceCollection pset_res("test_pset_res");
  elec.createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);
  ResourceCollection nl_ecp_res("test_nl_ecp_res");
  nl_ecp.createResource(nl_ecp_res);
  ResourceCollectionTeamLock<OperatorBase> nl_ecp_lock(nl_ecp_res, o_list);

  testing::TestNonLocalECPotential::copyGridUnrotatedForTest(nl_ecp);
  testing::TestNonLocalECPotential::copyGridUnrotatedForTest(nl_ecp2);

  // the energy evaluation builds the neighbor lists the v1 sweep reads
  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, std::nullopt, true);

  NonLocalTOperator move_op(TmoveKind::V1, 0.5 /*tau*/, 0.0 /*alpha*/, 0.0 /*gamma*/);

  TmoveV1Result res;
  if (batched)
    res.accepts = NonLocalECPotential::mw_makeNonLocalMovesPbyP(o_list, twf_list, p_list, move_op);
  else
  {
    res.accepts.resize(2);
    res.accepts[0] = nl_ecp.makeNonLocalMovesPbyP(psi, elec, move_op);
    res.accepts[1] = nl_ecp2.makeNonLocalMovesPbyP(psi2, elec2, move_op);
  }
  for (int iat = 0; iat < elec.getTotalNum(); ++iat)
    res.final_R.push_back(elec.R[iat]);
  for (int iat = 0; iat < elec2.getTotalNum(); ++iat)
    res.final_R.push_back(elec2.R[iat]);
  return res;
}
} // namespace

TEST_CASE("NonLocalECPotential Tmove v1 batched matches serial", "[hamiltonian]")
{
  for (const bool use_VP : {false, true})
  {
    const auto serial  = runTmoveV1(false, use_VP);
    const auto batched = runTmoveV1(true, use_VP);

    REQUIRE(serial.accepts.size() == batched.accepts.size());
    CHECK(serial.accepts == batched.accepts);
    // a sweep with no accepted move would leave the accept path untested
    const int total_accepts = std::accumulate(serial.accepts.begin(), serial.accepts.end(), 0);
    CHECK(total_accepts > 0);

    REQUIRE(serial.final_R.size() == batched.final_R.size());
    for (size_t i = 0; i < serial.final_R.size(); ++i)
      for (int d = 0; d < 3; ++d)
        CHECK(serial.final_R[i][d] == Approx(batched.final_R[i][d]).epsilon(1e-12));
  }
}

} // namespace qmcplusplus

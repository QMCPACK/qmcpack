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

#include "catch.hpp"

#include "Configuration.h"
#include "Numerics/Quadrature.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "TestListenerFunction.h"
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
  static void mw_evaluateImpl(NonLocalECPotential& nl_ecp,
                              const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& twf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              bool Tmove,
                              const std::optional<ListenerOption<Real>> listener_opt,
                              bool keep_grid)
  {
    nl_ecp.mw_evaluateImpl(o_list, twf_list, p_list, Tmove, listener_opt, keep_grid);
  }
};

} // namespace testing

TEST_CASE("NonLocalECPotential", "[hamiltonian]")
{
  using Real         = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::FullPrecRealType;
  using Position     = QMCTraits::PosType;
  using testing::getParticularListener;

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
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
  const int ei_table_index = elec.addTable(ions);
  elec.update();

  ParticleSet elec2(elec);
  elec2.update();

  RefVector<ParticleSet> ptcls{elec, elec2};
  RefVectorWithLeader<ParticleSet> p_list(elec, ptcls);

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  TrialWaveFunction psi2(runtime_options);
  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, psi2});

  bool doForces = false;
  bool use_DLA  = false;

  NonLocalECPotential nl_ecp(ions, elec, psi, doForces, use_DLA);

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

  RefVector<OperatorBase> nl_ecps{nl_ecp, nl_ecp2};
  RefVectorWithLeader<OperatorBase> o_list(nl_ecp, nl_ecps);
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

  // I'd like to see this gone when legacy drivers are dropped but for now we'll check against
  // the single particle API
  auto value = o_list[0].evaluateDeterministic(p_list[0]);

  CHECK(std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0) == Approx(value));
  CHECK(std::accumulate(local_pots2.begin(), local_pots2.begin() + local_pots2.cols(), 0.0) == Approx(value));
  CHECK(std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0) == Approx(-2.1047365665));
  CHECK(std::accumulate(ion_pots2.begin(), ion_pots2.begin() + ion_pots2.cols(), 0.0) == Approx(-2.1047365665));

  CHECK(!testing::TestNonLocalECPotential::didGridChange(nl_ecp));

  elec.R[0] = {0.5, 0.0, 2.0};
  elec.update();

  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, listener_opt, true);

  CHECK(!testing::TestNonLocalECPotential::didGridChange(nl_ecp));
  auto value2 = o_list[0].evaluateDeterministic(elec);

  CHECK(std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0) == Approx(value2));
  // check the second walker which will be unchanged.
  CHECK(std::accumulate(local_pots2[1], local_pots2[1] + local_pots2.cols(), 0.0) == Approx(value));

  // Randomizing grid does nothing for Na pp
  testing::TestNonLocalECPotential::mw_evaluateImpl(nl_ecp, o_list, twf_list, p_list, false, listener_opt, false);
  auto value3 = o_list[0].evaluateDeterministic(p_list[0]);
  CHECK(std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0) == Approx(value3));
}

} // namespace qmcplusplus

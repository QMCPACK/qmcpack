//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia Nationaln Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "QMCHamiltonians/SOECPotential.h"
#include "QMCWaveFunctions/ElectronGas/FreeOrbital.h"
#include "QMCWaveFunctions/SpinorSet.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBatched.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "TestListenerFunction.h"
#include "Utilities/RuntimeOptions.h"
namespace qmcplusplus
{
namespace testing
{
class TestSOECPotential
{
  using Real = QMCTraits::RealType;

public:
  static void copyGridUnrotatedForTest(SOECPotential& so_ecp)
  {
    so_ecp.ppset_[0]->rrotsgrid_m_ = so_ecp.ppset_[0]->sgridxyz_m_;
  }
  static bool didGridChange(SOECPotential& so_ecp)
  {
    return so_ecp.ppset_[0]->rrotsgrid_m_ != so_ecp.ppset_[0]->sgridxyz_m_;
  }
  static void addVPs(const RefVectorWithLeader<OperatorBase>& o_list, const RefVectorWithLeader<ParticleSet>& p_list)
  {
    for (size_t iw = 0; iw < o_list.size(); iw++)
    {
      auto& sopp = o_list.getCastedElement<SOECPotential>(iw);
      auto& pset = p_list[iw];
      for (auto& uptr_comp : sopp.ppset_)
        uptr_comp.get()->initVirtualParticle(pset);
    }
  }
  static void mw_evaluateImpl(SOECPotential& so_ecp,
                              const RefVectorWithLeader<OperatorBase>& o_list,
                              const RefVectorWithLeader<TrialWaveFunction>& twf_list,
                              const RefVectorWithLeader<ParticleSet>& p_list,
                              const std::optional<ListenerOption<Real>> listener_opt,
                              bool keep_grid)
  {
    so_ecp.mw_evaluateImpl(o_list, twf_list, p_list, listener_opt, keep_grid);
  }
};
} // namespace testing

void doSOECPotentialTest(bool use_VPs)
{
  using Real         = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::FullPrecRealType;
  using Value        = QMCTraits::ValueType;
  using Pos          = QMCTraits::PosType;
  using testing::getParticularListener;

  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = false; // periodic
  lattice.R.diagonal(20);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_uptr = std::make_unique<ParticleSet>(simulation_cell);
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);
  ions.setName("ion0");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 0;
  ion_species(iatnumber, pIdx)  = 1;
  ions.createSK();

  elec.setName("e");
  elec.setSpinor(true);
  elec.create({2});
  elec.R[0]  = {0.138, -0.24, 0.216};
  elec.R[1]  = {-0.216, 0.24, -0.138};
  elec.spins = {0.2, 0.51};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();

  ions.resetGroups();
  elec.resetGroups();

  ParticleSet elec2(elec);
  elec2.update();

  RefVector<ParticleSet> ptcls{elec, elec2};
  RefVectorWithLeader<ParticleSet> p_list(elec, ptcls);
  ResourceCollection pset_res("test_pset_res");
  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, p_list);

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  std::vector<Pos> kup, kdn;
  std::vector<Real> k2up, k2dn;
  QMCTraits::IndexType nelec = elec.getTotalNum();
  REQUIRE(nelec == 2);

  kup.resize(nelec);
  kup[0] = Pos(1, 1, 1);
  kup[1] = Pos(2, 2, 2);

  kdn.resize(nelec);
  kdn[0] = Pos(2, 2, 2);
  kdn[1] = Pos(1, 1, 1);

  auto spo_up = std::make_unique<FreeOrbital>("free_orb_up", kup);
  auto spo_dn = std::make_unique<FreeOrbital>("free_orb_dn", kdn);

  auto spinor_set = std::make_unique<SpinorSet>("free_orb_spinor");
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));
  QMCTraits::IndexType norb = spinor_set->getOrbitalSetSize();
  REQUIRE(norb == 2);

  auto dd = std::make_unique<DiracDeterminantBatched<>>(std::move(spinor_set), 0, nelec);
  std::vector<std::unique_ptr<DiracDeterminantBase>> dirac_dets;
  dirac_dets.push_back(std::move(dd));
  auto sd = std::make_unique<SlaterDet>(elec, std::move(dirac_dets));
  psi.addComponent(std::move(sd));

  //Add the two body jastrow, parameters from test_J2_bspline
  //adding jastrow will allow for adding WF parameter derivatives since FreeOrbital doesn't
  //support that
  const char* particles = R"(<tmp>
  <jastrow name="J2" type="Two-Body" function="Bspline" print="yes" gpu="no">
      <correlation speciesA="u" speciesB="u" rcut="5" size="5">
      <coefficients id="uu" type="Array"> 0.02904699284 -0.1004179 -0.1752703883 -0.2232576505 -0.2728029201</coefficients>
        </correlation>
  </jastrow>
  </tmp>
  )";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr jas2 = xmlFirstElementChild(root);
  RadialJastrowBuilder jastrow(c, elec);
  psi.addComponent(jastrow.buildComponent(jas2));

  std::unique_ptr<TrialWaveFunction> psi_clone = psi.makeClone(elec2);
  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, *psi_clone});
  ResourceCollection twf_res("test_twf_res");
  psi.createResource(twf_res);
  ResourceCollectionTeamLock<TrialWaveFunction> mw_twf_lock(twf_res, twf_list);

  //Now we set up the SO ECP component.
  SOECPotential so_ecp(ions, elec, psi);
  ECPComponentBuilder ecp_comp_builder("test_read_soecp", c);
  okay = ecp_comp_builder.read_pp_file("so_ecp_test.xml");
  REQUIRE(okay);
  UPtr<SOECPComponent> so_ecp_comp = std::move(ecp_comp_builder.pp_so);
  so_ecp.addComponent(0, std::move(so_ecp_comp));
  UPtr<OperatorBase> so_ecp2_ptr = so_ecp.makeClone(elec2, *psi_clone);
  auto& so_ecp2                  = dynamic_cast<SOECPotential&>(*so_ecp2_ptr);

  StdRandom<FullPrecReal> rng(10101);
  StdRandom<FullPrecReal> rng2(10201);
  so_ecp.setRandomGenerator(&rng);
  so_ecp2.setRandomGenerator(&rng2);

  RefVector<OperatorBase> so_ecps{so_ecp, so_ecp2};
  RefVectorWithLeader<OperatorBase> o_list(so_ecp, so_ecps);
  if (use_VPs)
    testing::TestSOECPotential::addVPs(o_list, p_list);
  ResourceCollection so_ecp_res("test_so_ecp_res");
  so_ecp.createResource(so_ecp_res);
  ResourceCollectionTeamLock<OperatorBase> so_ecp_lock(so_ecp_res, o_list);

  testing::TestSOECPotential::copyGridUnrotatedForTest(so_ecp);
  testing::TestSOECPotential::copyGridUnrotatedForTest(so_ecp2);

  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp));
  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp2));

  int num_walkers = 2;
  int max_values  = 10;
  Matrix<Real> local_pots(num_walkers, max_values);
  Matrix<Real> local_pots2(num_walkers, max_values);
  std::vector<ListenerVector<Real>> listeners;
  listeners.emplace_back("sopotential", getParticularListener(local_pots));
  listeners.emplace_back("sopotential", getParticularListener(local_pots2));

  Matrix<Real> ion_pots(num_walkers, max_values);
  Matrix<Real> ion_pots2(num_walkers, max_values);
  std::vector<ListenerVector<Real>> ion_listeners;
  ion_listeners.emplace_back("sopotential", getParticularListener(ion_pots));
  ion_listeners.emplace_back("sopotential", getParticularListener(ion_pots2));

  ParticleSet::mw_update(p_list);
  TrialWaveFunction::mw_evaluateLog(twf_list, p_list);

  ListenerOption<Real> listener_opt{listeners, ion_listeners};
  testing::TestSOECPotential::mw_evaluateImpl(so_ecp, o_list, twf_list, p_list, listener_opt, true);

  //use single walker API to get reference value
  auto value = o_list[0].evaluateDeterministic(p_list[0]);
  CHECK(std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0) == Approx(value));
  CHECK(std::accumulate(local_pots2.begin(), local_pots2.begin() + local_pots2.cols(), 0.0) == Approx(value));
  CHECK(std::accumulate(ion_pots.begin(), ion_pots.begin() + ion_pots.cols(), 0.0) == Approx(value));
  CHECK(std::accumulate(ion_pots2.begin(), ion_pots2.begin() + ion_pots2.cols(), 0.0) == Approx(value));

  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp));
  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp2));

  elec.R[0] = {0.05, 0.0, -0.05};
  elec.update();
  testing::TestSOECPotential::mw_evaluateImpl(so_ecp, o_list, twf_list, p_list, listener_opt, true);

  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp));
  CHECK(!testing::TestSOECPotential::didGridChange(so_ecp2));
  auto value2 = o_list[0].evaluateDeterministic(elec);

  CHECK(std::accumulate(local_pots.begin(), local_pots.begin() + local_pots.cols(), 0.0) == Approx(value2));
  // check the second walker which will be unchanged.
  CHECK(std::accumulate(local_pots2[1], local_pots2[1] + local_pots2.cols(), 0.0) == Approx(value));

  //also check whether or not reference value from single_walker API is actually correct
  //this value comes directly from the reference code soecp_eval_reference.cpp
  CHECK(value == Approx(-3.530511241));
}

TEST_CASE("SOECPotential", "[hamiltonian]")
{
  //do test using VPs. This uses mw_ APIs for TWF in the mw_ SOECP APIs
  doSOECPotentialTest(true);
  //do test without VPs. This uses legacy APIs for TWF in the mw_ SOECP APIs
  doSOECPotentialTest(false);
}

} // namespace qmcplusplus

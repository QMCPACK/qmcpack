//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/tests/TestListenerFunction.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "Utilities/RuntimeOptions.h"

#include <functional>
#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

using QMCT = QMCTraits;
using Real = QMCT::RealType;

TEST_CASE("Bare Kinetic Energy", "[hamiltonian]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};
  elec.setName("elec");
  elec.create({2});
  elec.R[0]                = {0.0, 1.0, 0.0};
  elec.R[1]                = {1.0, 1.0, 0.0};
  SpeciesSet& tspecies     = elec.getSpeciesSet();
  int upIdx                = tspecies.addSpecies("u");
  int massIdx              = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;

  elec.addTable(ions);
  elec.update();

  const char* particles = R"(<tmp></tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr h1 = xmlFirstElementChild(root);

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  BareKineticEnergy bare_ke(elec, psi);
  bare_ke.put(h1);

  elec.L[0] = 1.0;
  elec.L[1] = 0.0;
  double v  = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);

  elec.L[0]    = 0.0;
  elec.L[1]    = 0.0;
  elec.G[0][0] = 1.0;
  elec.G[0][1] = 0.0;
  elec.G[0][2] = 0.0;
  elec.G[1][0] = 0.0;
  elec.G[1][1] = 0.0;
  elec.G[1][2] = 0.0;
  v            = bare_ke.evaluate(elec);
  REQUIRE(v == -0.5);
}

TEST_CASE("Bare KE Pulay PBC", "[hamiltonian]")
{
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(20);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion0");
  ions.create({2});
  ions.R[0]                     = {0.0, 0.0, 0.0};
  ions.R[1]                     = {6.0, 0.0, 0.0};
  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("Na");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(iatnumber, pIdx)  = 11;
  ions.createSK();

  elec.setName("e");
  std::vector<int> agroup(2, 1);
  elec.create(agroup);
  elec.R[0]                    = {2.0, 0.0, 0.0};
  elec.R[1]                    = {3.0, 0.0, 0.0};
  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

  elec.createSK();

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  //Cool.  Now to construct a wavefunction with 1 and 2 body jastrow (no determinant)
  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);

  //Add the two body jastrow
  const char* particles = R"(<tmp>
  <jastrow name="J2" type="Two-Body" function="Bspline" print="yes" gpu="no">
      <correlation speciesA="u" speciesB="d" rcut="10" size="8">
          <coefficients id="ud" type="Array"> 2.015599059 1.548994099 1.17959447 0.8769687661 0.6245736507 0.4133517767 0.2333851935 0.1035636904</coefficients>
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
  // Done with two body jastrow.

  //Add the one body jastrow.
  const char* particles2 = R"(<tmp>
  <jastrow name="J1" type="One-Body" function="Bspline" source="ion0" print="yes">
        <correlation elementType="Na" rcut="10" size="10" cusp="0">
          <coefficients id="eNa" type="Array"> 1.244201343 -1.188935609 -1.840397253 -1.803849126 -1.612058635 -1.35993202 -1.083353212 -0.8066295188 -0.5319252448 -0.3158819772</coefficients>
        </correlation>
      </jastrow>
  </tmp>
  )";
  bool okay3             = doc.parseFromString(particles2);
  REQUIRE(okay3);

  root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow1bdy(c, elec, ions);
  psi.addComponent(jastrow1bdy.buildComponent(jas1));

  const char* kexml = R"(<tmp></tmp>)";

  root = doc.getRoot();

  xmlNodePtr h1 = xmlFirstElementChild(root);

  BareKineticEnergy bare_ke(elec, psi);
  bare_ke.put(h1);

  // update all distance tables
  ions.update();
  elec.update();

  RealType logpsi = psi.evaluateLog(elec);

  RealType keval = bare_ke.evaluate(elec);

  //This is validated against an alternate code path (waveefunction tester for local energy).
  CHECK(keval == Approx(-0.147507745));

  ParticleSet::ParticlePos HFTerm, PulayTerm;
  HFTerm.resize(ions.getTotalNum());
  PulayTerm.resize(ions.getTotalNum());

  RealType keval2 = bare_ke.evaluateWithIonDerivs(elec, ions, psi, HFTerm, PulayTerm);

  CHECK(keval2 == Approx(-0.147507745));
  //These are validated against finite differences (delta=1e-6).
  CHECK(PulayTerm[0][0] == Approx(-0.13166));
  CHECK(PulayTerm[0][1] == Approx(0.0));
  CHECK(PulayTerm[0][2] == Approx(0.0));
  CHECK(PulayTerm[1][0] == Approx(-0.12145));
  CHECK(PulayTerm[1][1] == Approx(0.0));
  CHECK(PulayTerm[1][2] == Approx(0.0));
}

/** Provide a test scope parameterized on electron species mass that then can run a set of tests using
 *  its objects.
 *  I couldn't see a easier solution to deal with all the setup to get both coverage of this
 *  interesting feature considering various approaches to copy constructors the objects in the testing
 *  scope take.
 */
void testElecCase(double mass_up,
                  double mass_dn,
                  std::function<void(RefVectorWithLeader<OperatorBase>& o_list,
                                     RefVectorWithLeader<TrialWaveFunction>& twf_list,
                                     RefVectorWithLeader<ParticleSet>& p_list,
                                     Matrix<Real>& kinetic_energies,
                                     std::vector<ListenerVector<Real>>& listeners,
                                     std::vector<ListenerVector<Real>>& ion_listeners)> tests)
{
  using testing::getParticularListener;

  const SimulationCell simulation_cell;

  ParticleSet ions(simulation_cell);

  ions.setName("ion");
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};
  Matrix<Real> kinetic_energies(2);

  ParticleSet elec(simulation_cell);

  elec.setName("elec");
  elec.create({1, 1});
  elec.R[0]                    = {0.0, 1.0, 0.0};
  elec.R[1]                    = {1.0, 1.0, 0.0};
  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = mass_up;
  tspecies(massIdx, downIdx)   = mass_dn;

  elec.addTable(ions);
  elec.update();

  ParticleSet elec2(elec);
  elec2.R[1][2] = 1.0;
  elec2.update();

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  TrialWaveFunction psi_clone(runtime_options);

  RefVector<ParticleSet> ptcls{elec, elec2};
  RefVectorWithLeader<ParticleSet> p_list(elec, ptcls);

  BareKineticEnergy bare_ke(elec, psi);
  BareKineticEnergy bare_ke2(elec, psi_clone);

  RefVector<OperatorBase> bare_kes{bare_ke, bare_ke2};
  RefVectorWithLeader<OperatorBase> o_list(bare_ke, bare_kes);
  elec.L[0]  = 1.0;
  elec.L[1]  = 0.0;
  elec2.L[0] = 1.0;
  elec2.L[1] = 0.0;

  RefVectorWithLeader<TrialWaveFunction> twf_list(psi, {psi, psi_clone});

  ResourceCollection bare_ke_res("test_bare_ke_res");
  bare_ke.createResource(bare_ke_res);
  ResourceCollectionTeamLock<OperatorBase> bare_ke_lock(bare_ke_res, o_list);

  ResourceCollection pset_res("test_pset_res");
  elec.createResource(pset_res);
  ResourceCollectionTeamLock<ParticleSet> pset_lock(pset_res, p_list);

  std::vector<ListenerVector<Real>> listeners;
  listeners.emplace_back("kinetic", getParticularListener(kinetic_energies));

  std::vector<ListenerVector<Real>> ion_listeners;

  tests(o_list, twf_list, p_list, kinetic_energies, listeners, ion_listeners);
}

/** just set the values to test the electron weights */
auto getTestCaseForWeights(Real value1, Real value2, Real value3)
{
  return [value1, value2,
          value3](RefVectorWithLeader<OperatorBase>& o_list, RefVectorWithLeader<TrialWaveFunction>& twf_list,
                  RefVectorWithLeader<ParticleSet>& p_list, Matrix<Real>& kinetic_energies,
                  std::vector<ListenerVector<Real>>& listeners, std::vector<ListenerVector<Real>>& ion_listeners) {
    auto& bare_ke = o_list[0];
    bare_ke.mw_evaluatePerParticle(o_list, twf_list, p_list, listeners, ion_listeners);
    CHECK(bare_ke.getValue() == Approx(value1));
    auto& bare_ke2 = o_list[1];
    CHECK(bare_ke2.getValue() == Approx(value1));
    // Check that the sum of the particle energies is consistent with the total.
    CHECK(std::accumulate(kinetic_energies.begin(), kinetic_energies.begin() + kinetic_energies.cols(), 0.0) ==
          Approx(bare_ke.getValue()));
    // check a single particle
    CHECK(std::accumulate(kinetic_energies[1], kinetic_energies[1] + kinetic_energies.cols(), 0.0) == Approx(value1));

    auto& elec  = p_list[0];
    elec.L[0]   = 2.0;
    elec.L[1]   = 1.0;
    auto& elec2 = p_list[1];
    elec2.L[0]  = 2.0;
    elec2.L[1]  = 1.0;

    bare_ke.mw_evaluatePerParticle(o_list, twf_list, p_list, listeners, ion_listeners);
    CHECK(std::accumulate(kinetic_energies.begin(), kinetic_energies.begin() + kinetic_energies.cols(), 0.0) ==
          Approx(bare_ke.getValue()));
    CHECK(std::accumulate(kinetic_energies[1], kinetic_energies[1] + kinetic_energies.cols(), 0.0) == Approx(value2));

    // test that mw_evaluatePerParticleWithToperator decays to mw_evaluatePerParticle
    // If mw_evaluatePerParticle follows the standard behavior in OperatorBase the Listner will not receive a report
    // and tis values will remain as in the previous check.
    elec.L[0]  = 1.0;
    elec.L[1]  = 0.0;
    elec2.L[0] = 1.0;
    elec2.L[1] = 0.0;

    bare_ke.mw_evaluatePerParticleWithToperator(o_list, twf_list, p_list, listeners, ion_listeners);
    CHECK(std::accumulate(kinetic_energies[1], kinetic_energies[1] + kinetic_energies.cols(), 0.0) == Approx(value3));
  };
}

TEST_CASE("BareKineticEnergyListener", "[hamiltonian]")
{
  outputManager.pause();
  testElecCase(1.0, 1.0, getTestCaseForWeights(-0.5, -1.5, -0.5));
  testElecCase(0.5, 1.0, getTestCaseForWeights(-1.0, -2.5, -1.0));
  outputManager.resume();
}

} // namespace qmcplusplus

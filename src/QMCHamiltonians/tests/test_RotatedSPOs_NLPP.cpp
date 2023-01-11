//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//
// File created by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "type_traits/ConvertToReal.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/RotatedSPOs.h"

#include "QMCHamiltonians/HamiltonianFactory.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{



TEST_CASE("RotatedSPOs SplineR2R hcpBe values single det", "[wavefunction]")
{
  using RealType = QMCTraits::RealType;

  /*
    BEGIN Boilerplate stuff to make a simple SPOSet. Copied from test_einset.cpp
  */

  Communicate* c = OHMMS::Controller;

  ParticleSetPool pp(c);

  //setupParticleSetPoolHe(pp);

  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 4.32747284;
  lattice.R(0, 1) = 0.00000000;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = -2.16373642;
  lattice.R(1, 1) = 3.74770142;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = 0.00000000;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 6.78114995;

  lattice.BoxBConds = false;
  lattice.reset();




  const SimulationCell simulation_cell(lattice);
  pp.setSimulationCell(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  auto& elec(*elec_ptr);
  auto ions_uptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  ParticleSet& ions(*ions_uptr);
  //LRCoulombSingleton::CoulombHandler = std::make_unique<EwaldHandler3D>(ions);
  //LRCoulombSingleton::CoulombHandler->initBreakup(elec);

  elec.setName("e");
  elec.create({1,1});
  elec.R[0] = {0.1, 0.2, 0.3};
  elec.R[1] = {1.0, 1.0, 1.0};

  SpeciesSet& especies       = elec.getSpeciesSet();
  int upIdx                  = especies.addSpecies("u");
  int chargeIdx              = especies.addAttribute("charge");
  int massIdx                = especies.addAttribute("mass");
  especies(chargeIdx, upIdx) = -1;
  especies(massIdx, upIdx) = 1.0;
  int dnIdx                  = especies.addSpecies("d");
  especies(chargeIdx, dnIdx) = -1;
  especies(massIdx, dnIdx) = 1.0;
  elec.resetGroups(); // need to set Mass so


  elec.createSK();
  pp.addParticleSet(std::move(elec_ptr));

  ions.setName("ion0");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  SpeciesSet& tspecies       = ions.getSpeciesSet();
  int CIdx                   = tspecies.addSpecies("Be");
  int CchargeIdx             = tspecies.addAttribute("charge");
  int CatomicnumberIdx       = tspecies.addAttribute("atomicnumber");
  tspecies(CchargeIdx, CIdx) = 2;
  tspecies(CatomicnumberIdx, CIdx) = 4;

  ions.createSK();

  pp.addParticleSet(std::move(ions_uptr));

  WaveFunctionPool wp(pp, c);
  REQUIRE(wp.empty() == true);

  const char* wf_input = R"(
       <wavefunction name="psi0" target="e">
         <sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="double" truncate="no" save_coefs="yes">
            <sposet type="bspline" name="spo_up" size="2" spindataset="0" optimize="yes">
            </sposet>
            <sposet type="bspline" name="spo_down" size="2" spindataset="0" optimize="yes">
            </sposet>
         </sposet_builder>
         <determinantset>
            <slaterdeterminant>
               <determinant id="updet" group="u" sposet="spo_up" size="1"/>
               <determinant id="downdet" group="d" sposet="spo_down" size="1"/>
            </slaterdeterminant>
<!--
          <multideterminant optimize="no" spo_0="spo_up" spo_1="spo_down" algorithm="precomputed_table_method">
             <detlist size="1" type="DETS" nc0="0" ne0="1" nc1="0" ne1="1" nstates="2" cutoff="1e-20">
               <ci coeff="1.0" occ0="10" occ1="10"/>
             </detlist>
           </multideterminant>
-->


         </determinantset>
      </wavefunction>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  wp.put(root);
  TrialWaveFunction* psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != nullptr);
  REQUIRE(psi->getOrbitals().size() == 1);

  const char *ham_input = R"(
        <hamiltonian name="h0" type="generic" target="e">
<!--
         <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
         <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
-->
         <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml" algorithm="batched" pbc="no">
            <pseudo elementType="Be" href="Be.BFD.xml" nrule="2" disable_randomize_grid="yes"/>
         </pairpot>
      </hamiltonian>)";


  //HamiltonianFactory::PSetMap particle_set_map;
  //particle_set_map.emplace(ions.getName(), std::move(ions_uptr));
  //particle_set_map.emplace(elec.getName(), std::move(elec_ptr));

  //HamiltonianFactory::PsiPoolType psi_map;
  //psi_map.emplace("psi0", psi);

  //HamiltonianFactory hf("h0", elec, particle_set_map, psi_map, c);
  HamiltonianFactory hf("h0", elec, pp.getPool(), wp.getPool(), c);

  Libxml2Document doc2;
  bool okay2 = doc2.parseFromString(ham_input);
  REQUIRE(okay2);

  xmlNodePtr root2 = doc2.getRoot();
  hf.put(root2);

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  opt_vars.resetIndex();
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  elec.update();
  ions.update();

  double logval = psi->evaluateLog(elec);
  CHECK(logval == Approx(-1.2865633501081344));


  CHECK(elec.G[0][0] == ValueApprox(0.54752651));
  CHECK(elec.L[0] == ValueApprox(11.066512459947848));
  CHECK(elec.L[1] == ValueApprox(-0.4831061477045371));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(-2.97750823));
  CHECK(dlogpsi[1] == ValueApprox(-1.06146356));
  CHECK(dhpsioverpsi[0] == ValueApprox(-36.71707483));
  CHECK(dhpsioverpsi[1] == ValueApprox(-0.35274333));

  QMCHamiltonian* h = hf.getH();
  RandomGenerator myrng;
  h->setRandomGenerator(&myrng);

  app_log() << "H size = " << h->size() << std::endl;
  h->evaluate(elec);
  double loc_e = h->getLocalEnergy();
  double ke = h->getKineticEnergy();
  double loc_pot = h->getLocalPotential();
  std::cout << std::setprecision(16) << std::endl;
  std::cout << "local e = " << loc_e << std::endl;
  std::cout << "kinetic e = " << ke << std::endl;
  std::cout << "local potential = " << loc_pot << std::endl;
  // without e-n or pp
  CHECK(ke == Approx(-6.818620576308302));
  CHECK(loc_e == Approx(-3.562354739253797));

  Vector<ValueType> dlogpsi2(2);
  Vector<ValueType> dhpsioverpsi2(2);

  h->evaluateValueAndDerivatives(elec, opt_vars, dlogpsi2, dhpsioverpsi2);
  CHECK(dlogpsi2[0] == ValueApprox(-2.97750823));
  CHECK(dlogpsi2[1] == ValueApprox(-1.06146356));

  // No NLPP
  //CHECK(dhpsioverpsi2[0] == ValueApprox(-0.00043764));
  //CHECK(dhpsioverpsi2[1] == ValueApprox(-0.35274333));
  // With  NLPP
  CHECK(dhpsioverpsi2[0] == ValueApprox(-5.45054261));
  CHECK(dhpsioverpsi2[1] == ValueApprox(-0.34818307));
}

TEST_CASE("RotatedSPOs SplineR2R hcpBe values multi det", "[wavefunction]")
{
  using RealType = QMCTraits::RealType;

  /*
    BEGIN Boilerplate stuff to make a simple SPOSet. Copied from test_einset.cpp
  */

  Communicate* c = OHMMS::Controller;

  ParticleSetPool pp(c);

  //setupParticleSetPoolHe(pp);

  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 4.32747284;
  lattice.R(0, 1) = 0.00000000;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = -2.16373642;
  lattice.R(1, 1) = 3.74770142;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = 0.00000000;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 6.78114995;


  const SimulationCell simulation_cell(lattice);
  pp.setSimulationCell(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  auto& elec(*elec_ptr);
  auto ions_uptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  ParticleSet& ions(*ions_uptr);

  elec.setName("e");
  elec.create({1,1});
  elec.R[0] = {0.1, 0.2, 0.3};
  elec.R[1] = {1.0, 1.0, 1.0};

  SpeciesSet& especies       = elec.getSpeciesSet();
  int upIdx                  = especies.addSpecies("u");
  int chargeIdx              = especies.addAttribute("charge");
  int massIdx                = especies.addAttribute("mass");
  especies(chargeIdx, upIdx) = -1;
  especies(massIdx, upIdx) = 1.0;
  int dnIdx                  = especies.addSpecies("d");
  especies(chargeIdx, dnIdx) = -1;
  especies(massIdx, dnIdx) = 1.0;
  elec.resetGroups(); // need to set Mass so


  pp.addParticleSet(std::move(elec_ptr));

  ions.setName("ion0");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  SpeciesSet& tspecies       = ions.getSpeciesSet();
  int CIdx                   = tspecies.addSpecies("Be");
  int CchargeIdx             = tspecies.addAttribute("charge");
  int CatomicnumberIdx       = tspecies.addAttribute("atomicnumber");
  tspecies(CchargeIdx, CIdx) = 2;
  tspecies(CatomicnumberIdx, CIdx) = 4;


  pp.addParticleSet(std::move(ions_uptr));

  WaveFunctionPool wp(pp, c);
  REQUIRE(wp.empty() == true);

  const char* wf_input = R"(
       <wavefunction name="psi0" target="e">
         <sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="double" truncate="no" save_coefs="yes">
            <sposet type="bspline" name="spo_up" size="2" spindataset="0" optimize="yes">
            </sposet>
            <sposet type="bspline" name="spo_down" size="2" spindataset="0" optimize="yes">
            </sposet>
         </sposet_builder>
         <determinantset>
<!--
            <slaterdeterminant>
               <determinant id="updet" group="u" sposet="spo_up" size="1"/>
               <determinant id="downdet" group="d" sposet="spo_down" size="1"/>
            </slaterdeterminant>
-->
          <multideterminant optimize="no" spo_0="spo_up" spo_1="spo_down" algorithm="precomputed_table_method">
             <detlist size="1" type="DETS" nc0="0" ne0="1" nc1="0" ne1="1" nstates="2" cutoff="1e-20">
               <ci coeff="1.0" occ0="10" occ1="10"/>
             </detlist>
           </multideterminant>


         </determinantset>
      </wavefunction>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  wp.put(root);
  TrialWaveFunction* psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != nullptr);
  REQUIRE(psi->getOrbitals().size() == 1);


  const char *ham_input = R"(
        <hamiltonian name="h0" type="generic" target="e">
<!--
         <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
         <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
-->
         <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml" algorithm="non-batched" pbc="no">
            <pseudo elementType="Be" href="Be.BFD.xml" nrule="2" disable_randomize_grid="yes"/>
         </pairpot>
      </hamiltonian>)";


  HamiltonianFactory hf("h0", elec, pp.getPool(), wp.getPool(), c);

  Libxml2Document doc2;
  bool okay2 = doc2.parseFromString(ham_input);
  REQUIRE(okay2);

  xmlNodePtr root2 = doc2.getRoot();
  hf.put(root2);

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  opt_vars.resetIndex();
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  elec.update();

  double logval = psi->evaluateLog(elec);
  CHECK(logval == Approx(-1.2865633501081344));

  CHECK(elec.G[0][0] == ValueApprox(0.54752651));
  CHECK(elec.L[0] == ValueApprox(11.066512459947848));
  CHECK(elec.L[1] == ValueApprox(-0.4831061477045371));


  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(-2.97750823));
  CHECK(dlogpsi[1] == ValueApprox(-1.06146356));
  CHECK(dhpsioverpsi[0] == ValueApprox(-36.71707483));
  CHECK(dhpsioverpsi[1] == ValueApprox(-0.35274333));

  RefVectorWithLeader<TrialWaveFunction> wf_list(*psi, {*psi});
  RefVectorWithLeader<ParticleSet> p_list(elec, {elec});

  // Test list with one wavefunction

  int nparam = 2;
  int nentry = 1;
  RecordArray<ValueType> dlogpsi_list(nentry, nparam);
  RecordArray<ValueType> dhpsi_over_psi_list(nentry, nparam);

  TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, opt_vars, dlogpsi_list, dhpsi_over_psi_list);

  CHECK(dlogpsi_list[0][0] == dlogpsi[0]);
  CHECK(dlogpsi_list[0][1] == dlogpsi[1]);
  CHECK(dhpsi_over_psi_list[0][0] == dhpsioverpsi[0]);
  CHECK(dhpsi_over_psi_list[0][1] == dhpsioverpsi[1]);


  QMCHamiltonian* h = hf.getH();
  RandomGenerator myrng;
  h->setRandomGenerator(&myrng);

  app_log() << "H size = " << h->size() << std::endl;
  h->evaluate(elec);
  double loc_e = h->getLocalEnergy();
  double ke = h->getKineticEnergy();
  double loc_pot = h->getLocalPotential();
  std::cout << std::setprecision(16) << std::endl;
  std::cout << "local e = " << loc_e << std::endl;
  std::cout << "kinetic e = " << ke << std::endl;
  std::cout << "local potential = " << loc_pot << std::endl;
  // without e-n or pp
  CHECK(ke == Approx(-6.818620576308302));
  CHECK(loc_e == Approx(-3.562354739253797));

  auto *localECP_H = h->getHamiltonian("LocalECP");
  double local_pp = localECP_H->evaluate(elec);
  std::cout << "local pp = " << local_pp << std::endl;

  Vector<ValueType> dlogpsi2(2);
  Vector<ValueType> dhpsioverpsi2(2);

  h->evaluateValueAndDerivatives(elec, opt_vars, dlogpsi2, dhpsioverpsi2);
  //CHECK(dlogpsi2[0] == ValueApprox(-2.97750823));
  //CHECK(dlogpsi2[1] == ValueApprox(-1.06146356));
  CHECK(dlogpsi2[0] == dlogpsi[0]);
  CHECK(dlogpsi2[1] == dlogpsi[1]);

  // With  NLPP
  CHECK(dhpsioverpsi2[0] == ValueApprox(-5.45054261));
  CHECK(dhpsioverpsi2[1] == ValueApprox(-0.34818307));
  
  // batched interface
  RefVectorWithLeader<QMCHamiltonian> h_list(*h, {*h});

  RecordArray<ValueType> dlogpsi_list2(nentry, nparam);
  RecordArray<ValueType> dhpsi_over_psi_list2(nentry, nparam);

  h->mw_evaluateValueAndDerivatives(h_list, wf_list, p_list, opt_vars, dlogpsi_list2, dhpsi_over_psi_list2);

  CHECK(dlogpsi_list2[0][0] == ValueApprox(-2.97750823));
  CHECK(dlogpsi_list2[0][1] == ValueApprox(-1.06146356));

  CHECK(dhpsi_over_psi_list2[0][0] == ValueApprox(-5.45054261));
  CHECK(dhpsi_over_psi_list2[0][1] == ValueApprox(-0.34818307));


#if 0
  elec.R[0] = {2.0, 2.0, 4.0};

  elec.R[1] = {2.0, 4.0, 4.0};
  elec.update();
  h->evaluate(elec);
  ke = h->getKineticEnergy();
  std::cout << "second ke = " << ke << std::endl;
  local_pp = localECP_H->evaluate(elec);
  std::cout << "second local pp = " << local_pp << std::endl;
#endif
}

TEST_CASE("RotatedSPOs SplineR2R hcpBe rot1 values multi det", "[wavefunction]")
{
  using RealType = QMCTraits::RealType;

  /*
    BEGIN Boilerplate stuff to make a simple SPOSet. Copied from test_einset.cpp
  */

  Communicate* c = OHMMS::Controller;

  ParticleSetPool pp(c);

  //setupParticleSetPoolHe(pp);

  ParticleSet::ParticleLayout lattice;
  lattice.R(0, 0) = 4.32747284;
  lattice.R(0, 1) = 0.00000000;
  lattice.R(0, 2) = 0.00000000;
  lattice.R(1, 0) = -2.16373642;
  lattice.R(1, 1) = 3.74770142;
  lattice.R(1, 2) = 0.00000000;
  lattice.R(2, 0) = 0.00000000;
  lattice.R(2, 1) = 0.00000000;
  lattice.R(2, 2) = 6.78114995;


  const SimulationCell simulation_cell(lattice);
  pp.setSimulationCell(simulation_cell);
  auto elec_ptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  auto& elec(*elec_ptr);
  auto ions_uptr = std::make_unique<ParticleSet>(pp.getSimulationCell());
  ParticleSet& ions(*ions_uptr);

  elec.setName("e");
  elec.create({1,1});
  elec.R[0] = {0.1, 0.2, 0.3};
  elec.R[1] = {1.0, 1.0, 1.0};

  SpeciesSet& especies       = elec.getSpeciesSet();
  int upIdx                  = especies.addSpecies("u");
  int chargeIdx              = especies.addAttribute("charge");
  int massIdx                = especies.addAttribute("mass");
  especies(chargeIdx, upIdx) = -1;
  especies(massIdx, upIdx) = 1.0;
  int dnIdx                  = especies.addSpecies("d");
  especies(chargeIdx, dnIdx) = -1;
  especies(massIdx, dnIdx) = 1.0;
  elec.resetGroups(); // need to set Mass so


  pp.addParticleSet(std::move(elec_ptr));

  ions.setName("ion0");
  ions.create({1});
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  SpeciesSet& tspecies       = ions.getSpeciesSet();
  int CIdx                   = tspecies.addSpecies("Be");
  int CchargeIdx             = tspecies.addAttribute("charge");
  int CatomicnumberIdx       = tspecies.addAttribute("atomicnumber");
  tspecies(CchargeIdx, CIdx) = 2;
  tspecies(CatomicnumberIdx, CIdx) = 4;


  pp.addParticleSet(std::move(ions_uptr));

  WaveFunctionPool wp(pp, c);
  REQUIRE(wp.empty() == true);

  const char* wf_input = R"(
       <wavefunction name="psi0" target="e">
         <sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="double" truncate="no" save_coefs="yes">
            <sposet type="bspline" name="spo_up" size="2" spindataset="0" optimize="yes">
             <opt_vars>0.1</opt_vars>
            </sposet>
            <sposet type="bspline" name="spo_down" size="2" spindataset="0" optimize="yes">
             <opt_vars>0.1</opt_vars>
            </sposet>
         </sposet_builder>
         <determinantset>
<!--
            <slaterdeterminant>
               <determinant id="updet" group="u" sposet="spo_up" size="1"/>
               <determinant id="downdet" group="d" sposet="spo_down" size="1"/>
            </slaterdeterminant>
-->
          <multideterminant optimize="no" spo_0="spo_up" spo_1="spo_down" algorithm="precomputed_table_method">
             <detlist size="1" type="DETS" nc0="0" ne0="1" nc1="0" ne1="1" nstates="2" cutoff="1e-20">
               <ci coeff="1.0" occ0="10" occ1="10"/>
             </detlist>
           </multideterminant>


         </determinantset>
      </wavefunction>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  wp.put(root);
  TrialWaveFunction* psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != nullptr);
  REQUIRE(psi->getOrbitals().size() == 1);


  const char *ham_input = R"(
        <hamiltonian name="h0" type="generic" target="e">
<!--
         <pairpot type="coulomb" name="ElecElec" source="e" target="e"/>
         <pairpot type="coulomb" name="IonIon" source="ion0" target="ion0"/>
-->
         <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml" algorithm="non-batched" pbc="no">
            <pseudo elementType="Be" href="Be.BFD.xml" nrule="2" disable_randomize_grid="yes"/>
         </pairpot>
      </hamiltonian>)";


  HamiltonianFactory hf("h0", elec, pp.getPool(), wp.getPool(), c);

  Libxml2Document doc2;
  bool okay2 = doc2.parseFromString(ham_input);
  REQUIRE(okay2);

  xmlNodePtr root2 = doc2.getRoot();
  hf.put(root2);

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  opt_vars.resetIndex();
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  elec.update();

  double logval = psi->evaluateLog(elec);
  CHECK(logval == Approx(-1.7640777829734315));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(-4.38906396));
  CHECK(dlogpsi[1] == ValueApprox(-1.30028015));
  CHECK(dhpsioverpsi[0] == ValueApprox(-75.41699062082654));
  CHECK(dhpsioverpsi[1] == ValueApprox(-0.446294418102374));

  RefVectorWithLeader<TrialWaveFunction> wf_list(*psi, {*psi});
  RefVectorWithLeader<ParticleSet> p_list(elec, {elec});

  // Test list with one wavefunction

  int nparam = 2;
  int nentry = 1;
  RecordArray<ValueType> dlogpsi_list(nentry, nparam);
  RecordArray<ValueType> dhpsi_over_psi_list(nentry, nparam);

  TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, opt_vars, dlogpsi_list, dhpsi_over_psi_list);

  CHECK(dlogpsi_list[0][0] == Approx(dlogpsi[0]));
  CHECK(dlogpsi_list[0][1] == Approx(dlogpsi[1]));
  CHECK(dhpsi_over_psi_list[0][0] == Approx(dhpsioverpsi[0]));
  CHECK(dhpsi_over_psi_list[0][1] == Approx(dhpsioverpsi[1]));


  QMCHamiltonian* h = hf.getH();
  RandomGenerator myrng;
  h->setRandomGenerator(&myrng);

  app_log() << "H size = " << h->size() << std::endl;
  h->evaluate(elec);
  double loc_e = h->getLocalEnergy();
  double ke = h->getKineticEnergy();
  double loc_pot = h->getLocalPotential();
  std::cout << std::setprecision(16) << std::endl;
  std::cout << "local e = " << loc_e << std::endl;
  std::cout << "kinetic e = " << ke << std::endl;
  std::cout << "local potential = " << loc_pot << std::endl;
  // without e-n or pp
  CHECK(ke == Approx(-12.111681211640777));
  CHECK(loc_e == Approx(-4.381316707074211));

  auto *localECP_H = h->getHamiltonian("LocalECP");
  double local_pp = localECP_H->evaluate(elec);
  std::cout << "local pp = " << local_pp << std::endl;

  Vector<ValueType> dlogpsi2(2);
  Vector<ValueType> dhpsioverpsi2(2);

  h->evaluateValueAndDerivatives(elec, opt_vars, dlogpsi2, dhpsioverpsi2);
  CHECK(dlogpsi2[0] == Approx(dlogpsi[0]));
  CHECK(dlogpsi2[1] == Approx(dlogpsi[1]));
  CHECK(dhpsioverpsi2[0] == Approx(-11.195432172121736));
  CHECK(dhpsioverpsi2[1] == Approx(-0.440524726136756));
  
  // batched interface
  RefVectorWithLeader<QMCHamiltonian> h_list(*h, {*h});

  RecordArray<ValueType> dlogpsi_list2(nentry, nparam);
  RecordArray<ValueType> dhpsi_over_psi_list2(nentry, nparam);

  h->mw_evaluateValueAndDerivatives(h_list, wf_list, p_list, opt_vars, dlogpsi_list2, dhpsi_over_psi_list2);

  CHECK(dlogpsi_list2[0][0] == Approx(dlogpsi2[0]));
  CHECK(dlogpsi_list2[0][1] == Approx(dlogpsi2[1]));
  CHECK(dhpsi_over_psi_list2[0][0] == Approx(dhpsioverpsi2[0]));
  CHECK(dhpsi_over_psi_list2[0][1] == Approx(dhpsioverpsi2[1]));

}


} // namespace qmcplusplus

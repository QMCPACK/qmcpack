//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
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
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "Utilities/ProjectData.h"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
// Copied and extended from QMCWaveFunctions/tests/test_RotatedSPOs.cpp


void test_hcpBe_rotation(bool use_single_det, bool use_nlpp_batched)
{
  using RealType = QMCTraits::RealType;

  /*
    BEGIN Boilerplate stuff to make a simple SPOSet. Copied from test_einset.cpp
  */
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c = OHMMS::Controller;

  ParticleSetPool pp(c);

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
  elec.create({1, 1});
  elec.R[0] = {0.1, 0.2, 0.3};
  elec.R[1] = {1.0, 1.0, 1.0};

  SpeciesSet& especies       = elec.getSpeciesSet();
  int upIdx                  = especies.addSpecies("u");
  int chargeIdx              = especies.addAttribute("charge");
  int massIdx                = especies.addAttribute("mass");
  especies(chargeIdx, upIdx) = -1;
  especies(massIdx, upIdx)   = 1.0;
  int dnIdx                  = especies.addSpecies("d");
  especies(chargeIdx, dnIdx) = -1;
  especies(massIdx, dnIdx)   = 1.0;
  elec.resetGroups(); // need to set Mass so


  pp.addParticleSet(std::move(elec_ptr));

  ions.setName("ion0");
  ions.create({1});
  ions.R[0]                        = {0.0, 0.0, 0.0};
  SpeciesSet& tspecies             = ions.getSpeciesSet();
  int CIdx                         = tspecies.addSpecies("Be");
  int CchargeIdx                   = tspecies.addAttribute("charge");
  int CatomicnumberIdx             = tspecies.addAttribute("atomicnumber");
  tspecies(CchargeIdx, CIdx)       = 2;
  tspecies(CatomicnumberIdx, CIdx) = 4;


  pp.addParticleSet(std::move(ions_uptr));

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);
  REQUIRE(wp.empty() == true);

  const char* wf_input_multi_det = R"(
       <wavefunction name="psi0" target="e">
         <sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="double" truncate="no" save_coefs="no">
            <sposet type="bspline" name="spo_up" size="2" spindataset="0" optimize="yes">
            </sposet>
            <sposet type="bspline" name="spo_down" size="2" spindataset="0" optimize="yes">
            </sposet>
         </sposet_builder>
         <determinantset>
           <multideterminant optimize="no" spo_0="spo_up" spo_1="spo_down" algorithm="precomputed_table_method">
             <detlist size="1" type="DETS" nc0="0" ne0="1" nc1="0" ne1="1" nstates="2" cutoff="1e-20">
                <ci coeff="1.0" occ0="10" occ1="10"/>
             </detlist>
           </multideterminant>
         </determinantset>
      </wavefunction>)";

  const char* wf_input_single_det = R"(
       <wavefunction name="psi0" target="e">
         <sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion0" version="0.10" meshfactor="1.0" precision="double" truncate="no" save_coefs="no">
            <sposet type="bspline" name="spo_up" size="2" spindataset="0" optimize="yes">
            </sposet>
            <sposet type="bspline" name="spo_down" size="2" spindataset="0" optimize="yes">
            </sposet>
         </sposet_builder>
         <determinantset>
           <slaterdeterminant>
               <determinant id="spo_up" size="1">
                  <occupation mode="ground" spindataset="0"/>
               </determinant>
               <determinant id="spo_down" size="1">
                  <occupation mode="ground" spindataset="0"/>
               </determinant>
            </slaterdeterminant>
         </determinantset>
      </wavefunction>)";

  const char* wf_input = wf_input_multi_det;
  if (use_single_det)
    wf_input = wf_input_single_det;

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  wp.put(root);
  TrialWaveFunction* psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != nullptr);
  REQUIRE(psi->getOrbitals().size() == 1);


  // Note the pbc="no" setting to turn off long-range summation.
  const char* ham_input_nlpp_nonbatched = R"(
        <hamiltonian name="h0" type="generic" target="e">
         <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml" algorithm="non-batched" pbc="no">
            <pseudo elementType="Be" href="Be.BFD.xml" nrule="2" disable_randomize_grid="yes"/>
         </pairpot>
      </hamiltonian>)";

  const char* ham_input_nlpp_batched = R"(
        <hamiltonian name="h0" type="generic" target="e">
         <pairpot type="pseudo" name="PseudoPot" source="ion0" wavefunction="psi0" format="xml" algorithm="batched" pbc="no">
            <pseudo elementType="Be" href="Be.BFD.xml" nrule="2" disable_randomize_grid="yes"/>
         </pairpot>
      </hamiltonian>)";

  const char* ham_input = ham_input_nlpp_nonbatched;
  if (use_nlpp_batched)
    ham_input = ham_input_nlpp_batched;

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
#if defined(MIXED_PRECISION)
  CHECK(elec.L[1] == ValueApprox(-0.4830868071).epsilon(1e-3));
#else
  CHECK(elec.L[1] == ValueApprox(-0.4831061477045371));
#endif

  // Parameter derivatives of just the wavefunction
  // Values come from QMCWaveFunctions/tests/eval_bspline_spo.py
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

  CHECK(dlogpsi_list[0][0] == Approx(dlogpsi[0]));
  CHECK(dlogpsi_list[0][1] == Approx(dlogpsi[1]));
  CHECK(dhpsi_over_psi_list[0][0] == Approx(dhpsioverpsi[0]));
  CHECK(dhpsi_over_psi_list[0][1] == Approx(dhpsioverpsi[1]));


  QMCHamiltonian* h = hf.getH();
  RandomGenerator myrng;
  h->setRandomGenerator(&myrng);

  h->evaluate(elec);
  double loc_e = h->getLocalEnergy();
  double ke    = h->getKineticEnergy();
  CHECK(ke == Approx(-6.818620576308302));
  CHECK(loc_e == Approx(-3.562354739253797));

  auto* localECP_H = h->getHamiltonian("LocalECP");
  double local_pp  = localECP_H->evaluate(elec);

  Vector<ValueType> dlogpsi2(2);
  Vector<ValueType> dhpsioverpsi2(2);

  h->evaluateValueAndDerivatives(elec, opt_vars, dlogpsi2, dhpsioverpsi2);
  // Derivative the wavefunction is unchanged by NLPP
  CHECK(dlogpsi2[0] == Approx(dlogpsi[0]));
  CHECK(dlogpsi2[1] == Approx(dlogpsi[1]));

  // Derivative of H is different with NLPP included
  CHECK(dhpsioverpsi2[0] == ValueApprox(-5.45054261));
  CHECK(dhpsioverpsi2[1] == ValueApprox(-0.34818307));

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

TEST_CASE("RotatedSPOs SplineR2R hcpBe values", "[wavefunction]")
{
  SECTION("nlpp non-batched")
  {
    bool use_single_det   = GENERATE(true, false);
    bool use_nlpp_batched = false;
    test_hcpBe_rotation(use_single_det, use_nlpp_batched);
  }

  SECTION("nlpp batched")
  {
    bool use_single_det   = true;
    bool use_nlpp_batched = true;
    test_hcpBe_rotation(use_single_det, use_nlpp_batched);
  }
}

} // namespace qmcplusplus

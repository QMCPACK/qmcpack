//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"

namespace qmcplusplus
{
void setup_He_wavefunction(Communicate* c,
                           ParticleSet& elec,
                           ParticleSet& ions,
                           std::unique_ptr<WaveFunctionFactory>& wff,
                           WaveFunctionFactory::PtclPoolType& particle_set_map)
{
  std::vector<int> agroup(2);
  int nelec = 2;
  agroup[0] = 1;
  agroup[1] = 1;
  elec.setName("e");
  elec.create(agroup);
  elec.R[0][0] = 1.0;
  elec.R[0][1] = 2.0;
  elec.R[0][2] = 3.0;
  elec.R[1][0] = 1.0;
  elec.R[1][1] = 2.1;
  elec.R[1][2] = 2.2;

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  int downIdx          = tspecies.addSpecies("d");
  int massIdx          = tspecies.addAttribute("mass");
  // Define the charge so the Jastrow cusp is set automatically
  int e_chargeIdx = tspecies.addAttribute("charge");
  // Mass is set oddly so the parameter derivative code is tested properly
  tspecies(massIdx, upIdx)       = 3.0;
  tspecies(massIdx, downIdx)     = 3.0;
  tspecies(e_chargeIdx, upIdx)   = -1.0;
  tspecies(e_chargeIdx, downIdx) = -1.0;
  elec.resetGroups();

  particle_set_map["e"] = &elec;

  ions.setName("ion0");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  SpeciesSet& he_species      = ions.getSpeciesSet();
  int He_Idx                  = he_species.addSpecies("He");
  int chargeIdx               = he_species.addAttribute("charge");
  tspecies(chargeIdx, He_Idx) = 2.0;
  tspecies(massIdx, upIdx)    = 2.0;
  particle_set_map["ion0"]    = &ions;

  elec.addTable(ions);

  wff = std::make_unique<WaveFunctionFactory>(&elec, particle_set_map, c);

  const char* wavefunction_xml = "<wavefunction name=\"psi0\" target=\"e\">  \
     <jastrow name=\"Jee\" type=\"Two-Body\" function=\"pade\"> \
      <correlation speciesA=\"u\" speciesB=\"d\"> \
        <var id=\"jud_b\" name=\"B\">0.8</var> \
      </correlation> \
     </jastrow> \
     <determinantset type=\"MO\" key=\"STO\" transform=\"no\" source=\"ion0\"> \
      <basisset> \
        <atomicBasisSet type=\"STO\" elementType=\"He\"> \
          <basisGroup rid=\"R0\" n=\"1\" l=\"0\" m=\"0\" type=\"Slater\"> \
             <radfunc exponent=\"2.0\"/> \
          </basisGroup> \
        </atomicBasisSet> \
      </basisset> \
      <slaterdeterminant> \
        <determinant id=\"updet\" spin=\"1\" size=\"1\"> \
          <coefficient id=\"updetC\" type=\"Array\" size=\"1\"> \
            1.0 \
          </coefficient> \
        </determinant> \
        <determinant id=\"downdet\" spin=\"-1\" size=\"1\"> \
          <coefficient id=\"downdetC\" type=\"Array\" size=\"1\"> \
            1.0 \
          </coefficient> \
        </determinant> \
      </slaterdeterminant> \
      </determinantset> \
    </wavefunction>";


  Libxml2Document doc;
  bool okay = doc.parseFromString(wavefunction_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  wff->put(root);

  REQUIRE(wff->targetPsi != NULL);
  REQUIRE(wff->targetPsi->size() == 2);
}

TEST_CASE("TrialWaveFunction flex_evaluateParameterDerivatives", "[wavefunction]")
{
  using ValueType = QMCTraits::ValueType;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;
  ParticleSet elec;
  ParticleSet ions;
  std::unique_ptr<WaveFunctionFactory> wff;
  WaveFunctionFactory::PtclPoolType particle_set_map;
  setup_He_wavefunction(c, elec, ions, wff, particle_set_map);
  TrialWaveFunction& psi(*(wff->targetPsi));

  ions.update();
  elec.update();

  const int nparam = 1;
  optimize::VariableSet var_param;
  psi.checkInVariables(var_param);

  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);

  psi.evaluateDerivatives(elec, var_param, dlogpsi, dhpsioverpsi);


  RefVector<TrialWaveFunction> wf_list;
  RefVector<ParticleSet> p_list;

  // Test list with one wavefunction

  int nentry = 1;
  RecordArray<ValueType> dlogpsi_list(nparam, nentry);
  RecordArray<ValueType> dhpsi_over_psi_list(nparam, nentry);

  wf_list.push_back(psi);
  p_list.push_back(elec);
  psi.flex_evaluateParameterDerivatives(wf_list, p_list, var_param, dlogpsi_list, dhpsi_over_psi_list);

  CHECK(dlogpsi[0] == ValueApprox(dlogpsi_list.getValue(0, 0)));
  CHECK(dhpsioverpsi[0] == ValueApprox(dhpsi_over_psi_list.getValue(0, 0)));

  // Test list with two wavefunctions

  nentry = 2;
  dlogpsi_list.resize(nparam, nentry);
  dhpsi_over_psi_list.resize(nparam, nentry);

  ParticleSet elec2(elec);
  elec2.R[0][0] = 0.9;
  elec2.update();

  // Will re-using  the same TrialWaveFunction work, or should a new one be created.
  //  If a new one is needed, what is the easiest way to copy?
  wf_list.push_back(psi);
  p_list.push_back(elec2);
  psi.flex_evaluateParameterDerivatives(wf_list, p_list, var_param, dlogpsi_list, dhpsi_over_psi_list);

  std::vector<ValueType> dlogpsi2(nparam);
  std::vector<ValueType> dhpsioverpsi2(nparam);

  psi.evaluateDerivatives(elec2, var_param, dlogpsi2, dhpsioverpsi2);

  CHECK(dlogpsi[0] == ValueApprox(dlogpsi_list.getValue(0, 0)));
  CHECK(dhpsioverpsi[0] == ValueApprox(dhpsi_over_psi_list.getValue(0, 0)));

  CHECK(dlogpsi2[0] == ValueApprox(dlogpsi_list.getValue(0, 1)));
  CHECK(dhpsioverpsi2[0] == ValueApprox(dhpsi_over_psi_list.getValue(0, 1)));

  SPOSetBuilderFactory::clear();
}


UPtrVector<ParticleSet::ParticleGradient_t> create_particle_gradient(int nelec, int nentry)
{
  UPtrVector<ParticleSet::ParticleGradient_t> G_list;
  for (int i = 0; i < nentry; i++)
    G_list.emplace_back(std::make_unique<ParticleSet::ParticleGradient_t>(nelec));

  return G_list;
}

UPtrVector<ParticleSet::ParticleLaplacian_t> create_particle_laplacian(int nelec, int nentry)
{
  UPtrVector<ParticleSet::ParticleLaplacian_t> L_list;
  for (int i = 0; i < nentry; i++)
    L_list.emplace_back(std::make_unique<ParticleSet::ParticleLaplacian_t>(nelec));

  return L_list;
}

TEST_CASE("TrialWaveFunction flex_evaluateDeltaLogSetup", "[wavefunction]")
{
  using ValueType = QMCTraits::ValueType;
  using RealType  = QMCTraits::RealType;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;
  ParticleSet elec1;
  ParticleSet ions;
  std::unique_ptr<WaveFunctionFactory> wff;
  WaveFunctionFactory::PtclPoolType particle_set_map;
  // This He wavefunction has two components
  // The orbitals are fixed and have not optimizable parameters.
  // The Jastrow factor does have an optimizable parameter.
  setup_He_wavefunction(c, elec1, ions, wff, particle_set_map);
  TrialWaveFunction& psi(*(wff->targetPsi));
  ions.update();
  elec1.update();

  ParticleSet elec1b(elec1);
  elec1b.update();

  ParticleSet elec2(elec1);
  elec2.R[0][0] = 0.9;
  elec2.update();
  ParticleSet elec2b(elec2);
  elec2b.update();

  TrialWaveFunction psi2(c);
  WaveFunctionComponent* orb1 = psi.getOrbitals()[0]->makeClone(elec2);
  psi2.addComponent(orb1, "orb1");
  WaveFunctionComponent* orb2 = psi.getOrbitals()[1]->makeClone(elec2);
  psi2.addComponent(orb2, "orb2");


  // Prepare to compare using list with one wavefunction and particleset

  int nentry = 1;
  int nelec  = 2;

  RefVector<TrialWaveFunction> wf_list;
  RefVector<ParticleSet> p_list;
  wf_list.push_back(psi);
  p_list.push_back(elec1b);

  // Evaluate new flex_evaluateDeltaLogSetup

  std::vector<RealType> logpsi_fixed_list(nentry);
  std::vector<RealType> logpsi_opt_list(nentry);

  auto fixedG_list_ptr = create_particle_gradient(nelec, nentry);
  auto fixedL_list_ptr = create_particle_laplacian(nelec, nentry);
  auto fixedG_list     = convertUPtrToRefVector(fixedG_list_ptr);
  auto fixedL_list     = convertUPtrToRefVector(fixedL_list_ptr);

  psi.flex_evaluateDeltaLogSetup(wf_list, p_list, logpsi_fixed_list, logpsi_opt_list, fixedG_list, fixedL_list);


  // Evaluate old (single item) evaluateDeltaLog

  RealType logpsi_fixed_r1;
  RealType logpsi_opt_r1;
  ParticleSet::ParticleGradient_t fixedG1;
  ParticleSet::ParticleLaplacian_t fixedL1;
  fixedG1.resize(nelec);
  fixedL1.resize(nelec);

  psi.evaluateDeltaLog(elec1, logpsi_fixed_r1, logpsi_opt_r1, fixedG1, fixedL1);

  CHECK(logpsi_fixed_r1 == Approx(logpsi_fixed_list[0]));
  CHECK(logpsi_opt_r1 == Approx(logpsi_opt_list[0]));

  CHECK(fixedG1[0][0] == ValueApprox(fixedG_list[0].get()[0][0]));
  CHECK(fixedG1[0][1] == ValueApprox(fixedG_list[0].get()[0][1]));
  CHECK(fixedG1[0][2] == ValueApprox(fixedG_list[0].get()[0][2]));
  CHECK(fixedG1[1][0] == ValueApprox(fixedG_list[0].get()[1][0]));
  CHECK(fixedG1[1][1] == ValueApprox(fixedG_list[0].get()[1][1]));
  CHECK(fixedG1[1][2] == ValueApprox(fixedG_list[0].get()[1][2]));

  CHECK(fixedL1[0] == ValueApprox(fixedL_list[0].get()[0]));
  CHECK(fixedL1[1] == ValueApprox(fixedL_list[0].get()[1]));

  // Compare the ParticleSet gradient and laplacian storage
  // This should be temporary until these get removed from ParticleSet
  CHECK(elec1b.L[0] == ValueApprox(elec1.L[0]));
  CHECK(elec1b.L[1] == ValueApprox(elec1.L[1]));

  CHECK(elec1b.G[0][0] == ValueApprox(elec1.G[0][0]));
  CHECK(elec1b.G[1][1] == ValueApprox(elec1.G[1][1]));

  // Prepare to compare using list with two wavefunctions and particlesets

  nentry = 2;

  wf_list.push_back(psi2);
  p_list.push_back(elec2b);

  ParticleSet::ParticleGradient_t G2;
  ParticleSet::ParticleLaplacian_t L2;
  G2.resize(nelec);
  L2.resize(nelec);
  fixedG_list.push_back(G2);
  fixedL_list.push_back(L2);

  std::vector<RealType> logpsi_fixed_list2(nentry);
  std::vector<RealType> logpsi_opt_list2(nentry);

  RealType logpsi_fixed_r1b;
  RealType logpsi_opt_r1b;
  psi2.evaluateDeltaLog(elec1, logpsi_fixed_r1b, logpsi_opt_r1b, fixedG1, fixedL1);

  CHECK(logpsi_fixed_r1 == Approx(logpsi_fixed_r1b));
  CHECK(logpsi_opt_r1 == Approx(logpsi_opt_r1b));

  auto fixedG_list2_ptr = create_particle_gradient(nelec, nentry);
  auto fixedL_list2_ptr = create_particle_laplacian(nelec, nentry);
  auto fixedG_list2     = convertUPtrToRefVector(fixedG_list2_ptr);
  auto fixedL_list2     = convertUPtrToRefVector(fixedL_list2_ptr);

  psi2.flex_evaluateDeltaLogSetup(wf_list, p_list, logpsi_fixed_list2, logpsi_opt_list2, fixedG_list2, fixedL_list2);

  // Evaluate old (single item) evaluateDeltaLog corresponding to the second wavefunction/particleset

  RealType logpsi_fixed_r2;
  RealType logpsi_opt_r2;
  ParticleSet::ParticleGradient_t fixedG2;
  ParticleSet::ParticleLaplacian_t fixedL2;
  fixedG2.resize(nelec);
  fixedL2.resize(nelec);

  psi2.setLogPsi(0.0);
  psi2.evaluateDeltaLog(elec2, logpsi_fixed_r2, logpsi_opt_r2, fixedG2, fixedL2);


  CHECK(logpsi_fixed_r1 == Approx(logpsi_fixed_r1b));

  CHECK(logpsi_fixed_list[0] == Approx(logpsi_fixed_list2[0]));

  CHECK(logpsi_fixed_r1 == Approx(logpsi_fixed_list2[0]));
  CHECK(logpsi_opt_r1 == Approx(logpsi_opt_list2[0]));

  CHECK(logpsi_fixed_r2 == Approx(logpsi_fixed_list2[1]));
  CHECK(logpsi_opt_r2 == Approx(logpsi_opt_list2[1]));

  // Laplacian for first entry in the wavefunction/particleset list
  CHECK(fixedL1[0] == ValueApprox(fixedL_list2[0].get()[0]));
  CHECK(fixedL1[1] == ValueApprox(fixedL_list2[0].get()[1]));
  // Laplacian for second entry in the wavefunction/particleset list
  CHECK(fixedL2[0] == ValueApprox(fixedL_list2[1].get()[0]));
  CHECK(fixedL2[1] == ValueApprox(fixedL_list2[1].get()[1]));


  // First entry wavefunction/particleset list
  // Gradient for first electron
  CHECK(fixedG1[0][0] == ValueApprox(fixedG_list2[0].get()[0][0]));
  CHECK(fixedG1[0][1] == ValueApprox(fixedG_list2[0].get()[0][1]));
  CHECK(fixedG1[0][2] == ValueApprox(fixedG_list2[0].get()[0][2]));
  // Gradient for second electron
  CHECK(fixedG1[1][0] == ValueApprox(fixedG_list2[0].get()[1][0]));
  CHECK(fixedG1[1][1] == ValueApprox(fixedG_list2[0].get()[1][1]));
  CHECK(fixedG1[1][2] == ValueApprox(fixedG_list2[0].get()[1][2]));


  // Second entry wavefunction/particleset list
  // Gradient for first electron
  CHECK(fixedG2[0][0] == ValueApprox(fixedG_list2[1].get()[0][0]));
  CHECK(fixedG2[0][1] == ValueApprox(fixedG_list2[1].get()[0][1]));
  CHECK(fixedG2[0][2] == ValueApprox(fixedG_list2[1].get()[0][2]));
  // Gradient for second electron
  CHECK(fixedG2[1][0] == ValueApprox(fixedG_list2[1].get()[1][0]));
  CHECK(fixedG2[1][1] == ValueApprox(fixedG_list2[1].get()[1][1]));
  CHECK(fixedG2[1][2] == ValueApprox(fixedG_list2[1].get()[1][2]));

  // Compare the ParticleSet gradient and laplacian storage
  // This should be temporary until these get removed from ParticleSet
  CHECK(elec1b.L[0] == ValueApprox(elec1.L[0]));
  CHECK(elec1b.L[1] == ValueApprox(elec1.L[1]));
  CHECK(elec2b.L[0] == ValueApprox(elec2.L[0]));
  CHECK(elec2b.L[1] == ValueApprox(elec2.L[1]));

  CHECK(elec2b.G[0][0] == ValueApprox(elec2.G[0][0]));
  CHECK(elec2b.G[1][1] == ValueApprox(elec2.G[1][1]));


  // these lists not used if 'recompute' is false
  RefVector<ParticleSet::ParticleGradient_t> dummyG_list;
  RefVector<ParticleSet::ParticleLaplacian_t> dummyL_list;

  std::vector<RealType> logpsi_variable_list(nentry);
  TrialWaveFunction::flex_evaluateDeltaLog(wf_list, p_list, logpsi_variable_list, dummyG_list, dummyL_list, false);

  RealType logpsi1 = psi.evaluateDeltaLog(p_list[0], false);
  CHECK(logpsi1 == Approx(logpsi_variable_list[0]));

  RealType logpsi2 = psi2.evaluateDeltaLog(p_list[1], false);
  CHECK(logpsi2 == Approx(logpsi_variable_list[1]));


  // Now check with 'recompute = true'
  auto dummyG_list2_ptr = create_particle_gradient(nelec, nentry);
  auto dummyL_list2_ptr = create_particle_laplacian(nelec, nentry);
  auto dummyG_list2     = convertUPtrToRefVector(dummyG_list2_ptr);
  auto dummyL_list2     = convertUPtrToRefVector(dummyL_list2_ptr);

  std::vector<RealType> logpsi_variable_list2(nentry);

  TrialWaveFunction::flex_evaluateDeltaLog(wf_list, p_list, logpsi_variable_list2, dummyG_list2, dummyL_list2, true);

  RealType logpsi1b = psi.evaluateDeltaLog(p_list[0], true);
  CHECK(logpsi1b == Approx(logpsi_variable_list2[0]));

  RealType logpsi2b = psi2.evaluateDeltaLog(p_list[1], true);
  CHECK(logpsi2b == Approx(logpsi_variable_list2[1]));

  SPOSetBuilderFactory::clear();
}


} // namespace qmcplusplus

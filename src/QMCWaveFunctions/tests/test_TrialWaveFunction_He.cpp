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
                           std::unique_ptr<WaveFunctionFactory>& wff)
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

  WaveFunctionFactory::PtclPoolType particle_set_map;
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

  elec.addTable(ions, DT_SOA);

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
  setup_He_wavefunction(c, elec, ions, wff);
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
}

} // namespace qmcplusplus

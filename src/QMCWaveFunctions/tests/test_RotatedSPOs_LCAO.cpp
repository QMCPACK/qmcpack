//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"


#include <stdio.h>
#include <string>
#include <sstream>


// Reference values from gen_rotated_lcao_wf.py


namespace qmcplusplus
{
void setupParticleSetPool(ParticleSetPool& pp)
{
  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = R"(<tmp>
  <particleset name="ion0" size="1">
    <group name="He">
      <parameter name="charge">2</parameter>
    </group>
    <attrib name="position" datatype="posArray">
      0.0 0.0 0.0
    </attrib>
  </particleset>

  <particleset name="e" random="no">
    <group name="u" size="1">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
      1.0 2.0 3.0
    </attrib>
    </group>
    <group name="d" size="1">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
      0.0 1.1 2.2
    </attrib>
    </group>
  </particleset>
</tmp>)";
  Libxml2Document doc;

  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr part_ion = xmlFirstElementChild(root);
  pp.put(part_ion);
  xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
  pp.put(part_elec);
}


// No Jastrow, rotation angle theta1=0.1 and theta2=0.2
TEST_CASE("Rotated LCAO WF1", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(pp, c);

  REQUIRE(wp.empty() == true);


  const char* wf_input = R"(<wavefunction target='e'>
    <sposet_collection type="MolecularOrbital">
      <!-- Use a single Slater Type Orbital (STO) for the basis. Cusp condition is correct. -->
      <basisset keyword="STO" transform="no">
        <atomicBasisSet type="STO" elementType="He" normalized="no">
          <basisGroup rid="R0" l="0" m="0" type="Slater">
             <radfunc n="1" exponent="2.0"/>
          </basisGroup>
          <basisGroup rid="R1" l="0" m="0" type="Slater">
             <radfunc n="2" exponent="1.0"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="2" optimize="yes">
           <opt_vars>0.1</opt_vars>
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-down" size="2" optimize="yes">
           <opt_vars>0.2</opt_vars>
          <coefficient id="downdetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant id="spo-up" spin="1" size="2"/>
        <determinant id="spo-down" spin="-1" size="2"/>
      </slaterdeterminant>
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

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  ParticleSet* elec = pp.getParticleSet("e");
  elec->update();


  double logval = psi->evaluateLog(*elec);
  CHECK(logval == Approx(-9.26625670653773));

  CHECK(elec->G[0][0] == ValueApprox(-0.2758747113720909));
  CHECK(elec->L[0] == ValueApprox(-0.316459652026054));
  CHECK(elec->L[1] == ValueApprox(-0.6035591598540904));

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(2);
  std::vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);


  CHECK(dlogpsi[0] == ValueApprox(7.58753078998516));
  CHECK(dlogpsi[1] == ValueApprox(2.58896036829191));
  CHECK(dhpsioverpsi[0] == ValueApprox(2.59551625714144));
  CHECK(dhpsioverpsi[1] == ValueApprox(1.70071425070404));
}

// No Jastrow, rotation angle of 0
TEST_CASE("Rotated LCAO WF zero angle", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(pp, c);

  REQUIRE(wp.empty() == true);


  const char* wf_input = R"(<wavefunction target='e'>
    <sposet_collection type="MolecularOrbital">
      <!-- Use a single Slater Type Orbital (STO) for the basis. Cusp condition is correct. -->
      <basisset keyword="STO" transform="no">
        <atomicBasisSet type="STO" elementType="He" normalized="no">
          <basisGroup rid="R0" l="0" m="0" type="Slater">
             <radfunc n="1" exponent="2.0"/>
          </basisGroup>
          <basisGroup rid="R1" l="0" m="0" type="Slater">
             <radfunc n="2" exponent="1.0"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="2" optimize="yes">
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-down" size="2" optimize="yes">
          <coefficient id="downdetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant id="spo-up" spin="1" size="2"/>
        <determinant id="spo-down" spin="-1" size="2"/>
      </slaterdeterminant>
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

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  opt_vars.resetIndex();
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  ParticleSet* elec = pp.getParticleSet("e");
  elec->update();


  double logval = psi->evaluateLog(*elec);
  CHECK(logval == Approx(-11.467952668216984));

  CHECK(elec->G[0][0] == ValueApprox(-0.5345224838248487));
  CHECK(elec->L[0] == ValueApprox(-1.0690449676496974));
  CHECK(elec->L[1] == ValueApprox(-1.626231256363484));

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(2);
  std::vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(32.2062050179872));
  CHECK(dlogpsi[1] == ValueApprox(5.87482925187464));
  CHECK(dhpsioverpsi[0] == ValueApprox(46.0088643114103));
  CHECK(dhpsioverpsi[1] == ValueApprox(7.84119772047731));

  RefVectorWithLeader<TrialWaveFunction> wf_list(*psi, {*psi});
  RefVectorWithLeader<ParticleSet> p_list(*elec, {*elec});

  // Test list with one wavefunction

  int nparam = 2;
  int nentry = 1;
  RecordArray<ValueType> dlogpsi_list(nparam, nentry);
  RecordArray<ValueType> dhpsi_over_psi_list(nparam, nentry);

  TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, opt_vars, dlogpsi_list, dhpsi_over_psi_list);

  CHECK(dlogpsi_list.getValue(0, 0) == ValueApprox(32.2062050179872));
  CHECK(dlogpsi_list.getValue(1, 0) == ValueApprox(5.87482925187464));
  CHECK(dhpsi_over_psi_list.getValue(0, 0) == ValueApprox(46.0088643114103));
  CHECK(dhpsi_over_psi_list.getValue(1, 0) == ValueApprox(7.84119772047731));
}

// Rotation angle of 0 and add Jastrow factory
TEST_CASE("Rotated LCAO WF with jastrow", "[qmcapp]")
{
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(pp, c);

  REQUIRE(wp.empty() == true);


  const char* wf_input = R"(<wavefunction target='e'>

     <sposet_collection type="MolecularOrbital">
      <!-- Use a single Slater Type Orbital (STO) for the basis. Cusp condition is correct. -->
      <basisset keyword="STO" transform="no">
        <atomicBasisSet type="STO" elementType="He" normalized="no">
          <basisGroup rid="R0" l="0" m="0" type="Slater">
             <radfunc n="1" exponent="2.0"/>
          </basisGroup>
          <basisGroup rid="R1" l="0" m="0" type="Slater">
             <radfunc n="2" exponent="1.0"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="2" optimize="yes">
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-down" size="2" optimize="yes">
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
      </sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant id="spo-up" spin="1" size="2"/>
        <determinant id="spo-down" spin="-1" size="2"/>
      </slaterdeterminant>
    </determinantset>
    <jastrow name="Jee" type="Two-Body" function="pade">
      <correlation speciesA="u" speciesB="d">
        <var id="jud_b" name="B">0.1</var>
      </correlation>
    </jastrow>
   </wavefunction>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_input);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  wp.put(root);

  TrialWaveFunction* psi = wp.getWaveFunction("psi0");
  REQUIRE(psi != nullptr);
  REQUIRE(psi->getOrbitals().size() == 2);

  opt_variables_type opt_vars;
  psi->checkInVariables(opt_vars);
  opt_vars.resetIndex();
  psi->checkOutVariables(opt_vars);
  psi->resetParameters(opt_vars);

  ParticleSet* elec = pp.getParticleSet("e");
  elec->update();


  double logval = psi->evaluateLog(*elec);
  CHECK(logval == Approx(-15.791249652199634));

  CHECK(elec->G[0][0] == ValueApprox(-0.2956989647881321));
  CHECK(elec->L[0] == ValueApprox(-0.6560429678274734));
  CHECK(elec->L[1] == ValueApprox(-1.2132292565412577));

  using ValueType = QMCTraits::ValueType;
  std::vector<ValueType> dlogpsi(3);
  std::vector<ValueType> dhpsioverpsi(3);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(32.206205017987166));
  CHECK(dlogpsi[1] == ValueApprox(5.874829251874641));
  CHECK(dlogpsi[2] == ValueApprox(49.08414605622605));
  CHECK(dhpsioverpsi[0] == ValueApprox(32.462519534916666));
  CHECK(dhpsioverpsi[1] == ValueApprox(10.047601212881027));
  CHECK(dhpsioverpsi[2] == ValueApprox(2.0820644399551895));
}
} // namespace qmcplusplus

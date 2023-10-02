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
#include "QMCWaveFunctions/RotatedSPOs.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/LCAO/LCAOrbitalSet.h"
#include "checkMatrix.hpp"
#include "Utilities/ProjectData.h"

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

// Set particles for Be atom
void setupParticleSetPoolBe(ParticleSetPool& pp)
{
  // See ParticleIO/tests/test_xml_io.cpp for particle parsing
  const char* particles = R"(<qmcsystem>
  <particleset name="ion0" size="1">
    <group name="Be">
      <parameter name="charge">4</parameter>
      <parameter name="valence">4</parameter>
      <parameter name="atomicnumber">4</parameter>
    </group>
    <attrib name="position" datatype="posArray">
  0.0000000000e+00  0.0000000000e+00  0.0000000000e+00
</attrib>
    <attrib name="ionid" datatype="stringArray">
 Be
</attrib>
  </particleset>

  <particleset name="e" random="no">
    <group name="u" size="2">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
      0.7 2.0 3.0
      1.2 1.5 0.5
    </attrib>
    </group>
    <group name="d" size="2">
      <parameter name="charge">-1</parameter>
    <attrib name="position" datatype="posArray">
      1.5 1.6 1.5
      0.7 1.0 1.2
    </attrib>
    </group>
  </particleset>
</qmcsystem>)";

  Libxml2Document doc;

  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr part_ion = xmlFirstElementChild(root);
  pp.put(part_ion);
  xmlNodePtr part_elec = xmlNextElementSibling(part_ion);
  pp.put(part_elec);
}

std::string setupRotationXML(const std::string& rot_angle_up,
                             const std::string& rot_angle_down,
                             const std::string& coeff_up,
                             const std::string& coeff_down)
{
  // Replace with std::format when minimum standard is switched to C++20

  const std::string wf_input1 = R"(<wavefunction target='e'>
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
     <rotated_sposet name="rot-spo-up">
       <sposet basisset="LCAOBSet" name="spo-up">)";

  // Opt vars for up determinant
  //       <opt_vars>0.1</opt_vars>
  const std::string opt_vars_start_tag("<opt_vars>");
  const std::string opt_vars_end_tag("</opt_vars>");

  std::string rot_angle_up_element = opt_vars_start_tag + rot_angle_up + opt_vars_end_tag + "\n";

  // Construct the coefficient matrix XML element for the up determinant
  //       <coefficient id="updetC" type="Array" size="2">
  //         1.0 0.0
  //         0.0 1.0
  //       </coefficient>
  const std::string wf_input_coeff_up_start = R"(<coefficient id="updetC" type="Array" size="2">)";

  const std::string wf_input_coeff_up_end("</coefficient>");

  std::string coeff_up_element = wf_input_coeff_up_start + coeff_up + wf_input_coeff_up_end;

  const std::string sposet_end = R"(</sposet>)";

  // Middle part of XML input block
  const std::string wf_input2 = R"(
      </rotated_sposet>
      <rotated_sposet name="rot-spo-down">
        <sposet basisset="LCAOBSet" name="spo-down">)";

  // Opt vars for down determinant
  //   <opt_vars>0.2</opt_vars>
  std::string rot_angle_down_element = opt_vars_start_tag + rot_angle_down + opt_vars_end_tag + "\n";

  // Construct the coefficient matrix XML element for the down determinant
  //       <coefficient id="downdetC" type="Array" size="2">
  //         1.0 0.0
  //         0.0 1.0
  //       </coefficient>
  const std::string wf_input_coeff_down_start = R"(<coefficient id="downdetC" type="Array" size="2">)";

  const std::string wf_input_coeff_down_end("</coefficient>");

  std::string coeff_down_element = wf_input_coeff_down_start + coeff_down + wf_input_coeff_down_end;

  const std::string wf_input3 = R"(
      </rotated_sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant sposet="rot-spo-up"/>
        <determinant sposet="rot-spo-down"/>
      </slaterdeterminant>
    </determinantset>
   </wavefunction>)";


  // clang-format off
  std::string wf_input = std::string(wf_input1) + "\n" +
                         coeff_up_element + "\n" +
                         sposet_end + "\n" +
                         (rot_angle_up.empty() ? std::string() : rot_angle_up_element) +
                         wf_input2 + "\n" +
                         coeff_down_element + "\n" +
                         sposet_end + "\n" +
                         (rot_angle_down.empty() ? std::string() : rot_angle_down_element) +
                         std::string(wf_input3);
  // clang-format on

  return wf_input;
}

const std::string identity_coeff = R"(
            1.0 0.0
            0.0 1.0
          )";

// No Jastrow, rotation angle of 0. Identity coefficients.
TEST_CASE("Rotated LCAO WF0 zero angle", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);


  std::string wf_input = setupRotationXML("", "", identity_coeff, identity_coeff);

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
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
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
  RecordArray<ValueType> dlogpsi_list(nentry, nparam);
  RecordArray<ValueType> dhpsi_over_psi_list(nentry, nparam);

  TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, opt_vars, dlogpsi_list, dhpsi_over_psi_list);

  CHECK(dlogpsi_list[0][0] == ValueApprox(32.2062050179872));
  CHECK(dlogpsi_list[0][1] == ValueApprox(5.87482925187464));
  CHECK(dhpsi_over_psi_list[0][0] == ValueApprox(46.0088643114103));
  CHECK(dhpsi_over_psi_list[0][1] == ValueApprox(7.84119772047731));
}

// No Jastrow, rotation angle theta1=0.1 and theta2=0.2 from identity coefficients
TEST_CASE("Rotated LCAO WF1", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);


  std::string wf_input = setupRotationXML("0.1", "0.2", identity_coeff, identity_coeff);

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
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);


  CHECK(dlogpsi[0] == ValueApprox(7.58753078998516));
  CHECK(dlogpsi[1] == ValueApprox(2.58896036829191));
  CHECK(dhpsioverpsi[0] == ValueApprox(2.59551625714144));
  CHECK(dhpsioverpsi[1] == ValueApprox(1.70071425070404));
}

// Rotation angle of 0 and add Jastrow factory
TEST_CASE("Rotated LCAO WF2 with jastrow", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

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
      <rotated_sposet name="rot-spo-up">
        <sposet basisset="LCAOBSet" name="spo-up" method="history">
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
        </sposet>
      </rotated_sposet>
      <rotated_sposet name="rot-spo-down">
        <sposet basisset="LCAOBSet" name="spo-down" method="history">
          <coefficient id="updetC" type="Array" size="2">
            1.0 0.0
            0.0 1.0
          </coefficient>
        </sposet>
      </rotated_sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant sposet="rot-spo-up"/>
        <determinant sposet="rot-spo-down"/>
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
  Vector<ValueType> dlogpsi(3);
  Vector<ValueType> dhpsioverpsi(3);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(32.206205017987166));
  CHECK(dlogpsi[1] == ValueApprox(5.874829251874641));
  CHECK(dlogpsi[2] == ValueApprox(49.08414605622605));
  CHECK(dhpsioverpsi[0] == ValueApprox(32.462519534916666));
  CHECK(dhpsioverpsi[1] == ValueApprox(10.047601212881027));
  CHECK(dhpsioverpsi[2] == ValueApprox(2.0820644399551895));

  // Check for out-of-bounds array access when a variable is disabled.
  // When a variable is disabled, myVars.where() returns -1.
  opt_vars.insert("rot-spo-up_orb_rot_0000_0001", 0.0, false);
  psi->checkOutVariables(opt_vars);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);
}

// Test the case where the rotation has already been applied to
// the MO coefficients in the input file.
// Should give the same results as the "Rotated LCAO WF1 zero angle" test case

const std::string coeff_rot_by_point1 = R"(
            0.995004165278026 0.0998334166468282
           -0.0998334166468282 0.995004165278026
    )";

const std::string coeff_rot_by_point2 = R"(
             0.980066577841242  0.198669330795061
            -0.198669330795061  0.980066577841242
    )";

TEST_CASE("Rotated LCAO WF1, MO coeff rotated, zero angle", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);

  std::string wf_input = setupRotationXML("", "", coeff_rot_by_point1, coeff_rot_by_point2);

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
  CHECK(logval == Approx(-9.26625670653773));

  CHECK(elec->G[0][0] == ValueApprox(-0.2758747113720909));
  CHECK(elec->L[0] == ValueApprox(-0.316459652026054));
  CHECK(elec->L[1] == ValueApprox(-0.6035591598540904));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);


  CHECK(dlogpsi[0] == ValueApprox(7.58753078998516));
  CHECK(dlogpsi[1] == ValueApprox(2.58896036829191));
  CHECK(dhpsioverpsi[0] == ValueApprox(2.59551625714144));
  CHECK(dhpsioverpsi[1] == ValueApprox(1.70071425070404));
}
// Test the case where half the rotation has already been applied to
// the MO coefficients in the input file and half the rotation is
// applied through the input.
// Should give the same results as the "Rotated LCAO WF1 zero angle" test case

const std::string coeff_rot_by_point05 = R"(
             0.998750260394966 0.0499791692706783
            -0.0499791692706783 0.998750260394966
    )";

TEST_CASE("Rotated LCAO WF1 MO coeff rotated, half angle", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);

  std::string wf_input = setupRotationXML("0.05", "0.1", coeff_rot_by_point05, coeff_rot_by_point1);

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
  CHECK(logval == Approx(-9.26625670653773));

  CHECK(elec->G[0][0] == ValueApprox(-0.2758747113720909));
  CHECK(elec->L[0] == ValueApprox(-0.316459652026054));
  CHECK(elec->L[1] == ValueApprox(-0.6035591598540904));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(2);
  Vector<ValueType> dhpsioverpsi(2);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);


  CHECK(dlogpsi[0] == ValueApprox(7.58753078998516));
  CHECK(dlogpsi[1] == ValueApprox(2.58896036829191));
  CHECK(dhpsioverpsi[0] == ValueApprox(2.59551625714144));
  CHECK(dhpsioverpsi[1] == ValueApprox(1.70071425070404));
}

// Test rotation using stored coefficients
//  and test consistency between history list and global rotation
TEST_CASE("Rotated LCAO rotation consistency", "[qmcapp]")
{
  using RealType    = QMCTraits::RealType;
  using ValueType   = QMCTraits::ValueType;
  using ValueMatrix = SPOSet::ValueMatrix;

  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPool(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);
  REQUIRE(wp.empty() == true);

  // Only care that this wavefunction has 3 SPOs and a 3x3 coefficient matrix
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
          <basisGroup rid="R2" l="0" m="0" type="Slater">
             <radfunc n="3" exponent="1.0"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
      <rotated_sposet name="rot-spo-up">
        <sposet basisset="LCAOBSet" name="spo-up">
          <coefficient id="updetC" type="Array" size="3">
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
          </coefficient>
        </sposet>
      </rotated_sposet>
      <rotated_sposet name="rot-spo-down">
        <sposet basisset="LCAOBSet" name="spo-down">
          <coefficient id="updetC" type="Array" size="3">
            1.0 0.0 0.0
            0.0 1.0 0.0
            0.0 0.0 1.0
          </coefficient>
        </sposet>
      </rotated_sposet>
    </sposet_collection>
    <determinantset type="MO" key="STO" transform="no" source="ion0">
      <slaterdeterminant>
        <determinant id="rot-spo-up"/>
        <determinant id="rot-spo-down"/>
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

  // Type should be pointer to SlaterDet
  auto orb1       = psi->getOrbitals()[0].get();
  SlaterDet* sdet = dynamic_cast<SlaterDet*>(orb1);
  REQUIRE(sdet != nullptr);

  // Use the SPOs from different spins to separately track rotation applied using the stored coefficients
  // versus the regular coefficients.
  // The coefficient matrices should match after the same rotations are applied to each.
  SPOSetPtr spoptr     = sdet->getPhi(0);
  RotatedSPOs* rot_spo = dynamic_cast<RotatedSPOs*>(spoptr);
  REQUIRE(rot_spo != nullptr);

  SPOSetPtr spoptr1           = sdet->getPhi(1);
  RotatedSPOs* global_rot_spo = dynamic_cast<RotatedSPOs*>(spoptr1);
  REQUIRE(global_rot_spo != nullptr);

  std::vector<RealType> params1 = {0.1, 0.2};

  // Apply against the existing coefficients
  rot_spo->apply_rotation(params1, false);

  // Apply against the stored coefficients
  global_rot_spo->apply_rotation(params1, true);

  // Value after first rotation, computed from gen_matrix_ops.py
  // clang-format off
  std::vector<ValueType> rot_data0 =
        {  0.975103993210479,   0.0991687475215628,  0.198337495043126,
          -0.0991687475215628,  0.995020798642096,  -0.00995840271580824,
          -0.198337495043126,  -0.00995840271580824, 0.980083194568384 };
  // clang-format on

  ValueMatrix new_rot_m0(rot_data0.data(), 3, 3);

  LCAOrbitalSet* lcao1       = dynamic_cast<LCAOrbitalSet*>(rot_spo->Phi.get());
  LCAOrbitalSet* lcao_global = dynamic_cast<LCAOrbitalSet*>(global_rot_spo->Phi.get());

  CheckMatrixResult check_matrix_result = checkMatrix(*lcao1->C, *lcao_global->C, true);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  CheckMatrixResult check_matrix_result0 = checkMatrix(*lcao_global->C, new_rot_m0, true);
  CHECKED_ELSE(check_matrix_result0.result) { FAIL(check_matrix_result0.result_message); }

  std::vector<RealType> old_params = {0.0, 0.0, 0.0};
  std::vector<RealType> new_params(3);
  global_rot_spo->applyDeltaRotation(params1, old_params, new_params);

  std::vector<RealType> params2 = {0.3, 0.15};
  rot_spo->apply_rotation(params2, false);

  std::vector<RealType> new_params2(3);
  global_rot_spo->applyDeltaRotation(params2, new_params, new_params2);
  CheckMatrixResult check_matrix_result2 = checkMatrix(*lcao1->C, *lcao_global->C, true);
  CHECKED_ELSE(check_matrix_result2.result) { FAIL(check_matrix_result2.result_message); }

  // Value after two rotations, computed from gen_matrix_ops.py
  // clang-format off
  std::vector<ValueType> rot_data3 =
    {  0.862374825309137,  0.38511734273482,   0.328624851461217,
      -0.377403929117215,  0.921689108007811, -0.0897522281988318,
      -0.337455085840952, -0.046624248032951,  0.940186281826872 };
  // clang-format on

  ValueMatrix new_rot_m3(rot_data3.data(), 3, 3);

  CheckMatrixResult check_matrix_result3 = checkMatrix(*lcao1->C, new_rot_m3, true);
  CHECKED_ELSE(check_matrix_result3.result) { FAIL(check_matrix_result3.result_message); }

  // Need to flip the sign on the first two entries to match the output from gen_matrix_ops.py
  std::vector<ValueType> expected_param = {0.3998099017676912, 0.34924318065960236, -0.02261313113492491};
  for (int i = 0; i < expected_param.size(); i++)
    CHECK(new_params2[i] == Approx(expected_param[i]));
}

// Reference values from rot_be_sto_wf.py
// Uses single determinant code path
TEST_CASE("Rotated LCAO Be single determinant", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPoolBe(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);

  Libxml2Document doc;
  bool okay = doc.parse("rot_Be_STO.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  wp.put(xmlFirstElementChild(root));


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
  CHECK(logval == Approx(-17.768474132175342));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(10);
  Vector<ValueType> dhpsioverpsi(10);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(0.24797938203759148));
  CHECK(dlogpsi[1] == ValueApprox(0.41454059122930453));
  CHECK(dlogpsi[2] == ValueApprox(0.7539626161586822));
  CHECK(dlogpsi[3] == ValueApprox(3.13489394217799));
  CHECK(dlogpsi[4] == ValueApprox(8.47176722646749));
  CHECK(dlogpsi[5] == ValueApprox(-0.48182453464906033));
  CHECK(dlogpsi[6] == ValueApprox(2.269206401396164));
  CHECK(dlogpsi[7] == ValueApprox(-1.883221269688377));
  CHECK(dlogpsi[8] == ValueApprox(-19.450964163527598));
  CHECK(dlogpsi[9] == ValueApprox(-47.28198556252034));

  CHECK(dhpsioverpsi[0] == ValueApprox(0.3662586398420111));
  CHECK(dhpsioverpsi[1] == ValueApprox(-5.544323554018982));
  CHECK(dhpsioverpsi[2] == ValueApprox(-0.7790656028274846));
  CHECK(dhpsioverpsi[3] == ValueApprox(24.930187483208087));
  CHECK(dhpsioverpsi[4] == ValueApprox(71.30301022344871));
  CHECK(dhpsioverpsi[5] == ValueApprox(-1.1614358798793771));
  CHECK(dhpsioverpsi[6] == ValueApprox(17.678711245652913));
  CHECK(dhpsioverpsi[7] == ValueApprox(2.491238469662668));
  CHECK(dhpsioverpsi[8] == ValueApprox(-79.37464297365679));
  CHECK(dhpsioverpsi[9] == ValueApprox(-227.0976672502695));
}

// Reference values from rot_be_sto_wf.py
// Uses multi-determinant code path with one determinant
TEST_CASE("Rotated LCAO Be multi determinant with one determinant", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPoolBe(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);

  Libxml2Document doc;
  bool okay = doc.parse("rot_multi_1det_Be_STO.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  wp.put(xmlFirstElementChild(root));

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
  CHECK(logval == Approx(-17.768474132175342));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(10);
  Vector<ValueType> dhpsioverpsi(10);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(0.24797938203759148));
  CHECK(dlogpsi[1] == ValueApprox(0.41454059122930453));
  CHECK(dlogpsi[2] == ValueApprox(0.7539626161586822));
  CHECK(dlogpsi[3] == ValueApprox(3.13489394217799));
  CHECK(dlogpsi[4] == ValueApprox(8.47176722646749));
  CHECK(dlogpsi[5] == ValueApprox(-0.48182453464906033));
  CHECK(dlogpsi[6] == ValueApprox(2.269206401396164));
  CHECK(dlogpsi[7] == ValueApprox(-1.883221269688377));
  CHECK(dlogpsi[8] == ValueApprox(-19.450964163527598));
  CHECK(dlogpsi[9] == ValueApprox(-47.28198556252034));

  CHECK(dhpsioverpsi[0] == ValueApprox(0.3662586398420111));
  CHECK(dhpsioverpsi[1] == ValueApprox(-5.544323554018982));
  CHECK(dhpsioverpsi[2] == ValueApprox(-0.7790656028274846));
  CHECK(dhpsioverpsi[3] == ValueApprox(24.930187483208087));
  CHECK(dhpsioverpsi[4] == ValueApprox(71.30301022344871));
  CHECK(dhpsioverpsi[5] == ValueApprox(-1.1614358798793771));
  CHECK(dhpsioverpsi[6] == ValueApprox(17.678711245652913));
  CHECK(dhpsioverpsi[7] == ValueApprox(2.491238469662668));
  CHECK(dhpsioverpsi[8] == ValueApprox(-79.37464297365679));
  CHECK(dhpsioverpsi[9] == ValueApprox(-227.0976672502695));
}

// Reference values from rot_multi_be_sto_wf.py
// Uses multi-determinant code path with two determinants
TEST_CASE("Rotated LCAO Be two determinant", "[qmcapp]")
{
  ProjectData test_project("test", ProjectData::DriverVersion::BATCH);
  Communicate* c;
  c = OHMMS::Controller;

  ParticleSetPool pp(c);
  setupParticleSetPoolBe(pp);

  WaveFunctionPool wp(test_project.getRuntimeOptions(), pp, c);

  REQUIRE(wp.empty() == true);

  Libxml2Document doc;
  bool okay = doc.parse("rot_multi_2det_Be_STO.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  wp.put(xmlFirstElementChild(root));

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
  CHECK(logval == Approx(-17.762687110866413));

  using ValueType = QMCTraits::ValueType;
  Vector<ValueType> dlogpsi(16);
  Vector<ValueType> dhpsioverpsi(16);
  psi->evaluateDerivatives(*elec, opt_vars, dlogpsi, dhpsioverpsi);

  CHECK(dlogpsi[0] == ValueApprox(0.05770308755290168));
  CHECK(dlogpsi[1] == ValueApprox(0.00593995768443123));
  CHECK(dlogpsi[2] == ValueApprox(0.24654846443828843));
  CHECK(dlogpsi[3] == ValueApprox(0.4214539468865001));
  CHECK(dlogpsi[4] == ValueApprox(0.7484015451192123));
  CHECK(dlogpsi[5] == ValueApprox(3.076586144487743));
  CHECK(dlogpsi[6] == ValueApprox(8.329621106110908));
  CHECK(dlogpsi[7] == ValueApprox(-0.4311398324864351));
  CHECK(dlogpsi[8] == ValueApprox(2.2561123798306273));
  CHECK(dlogpsi[9] == ValueApprox(-1.8723545015077454));
  CHECK(dlogpsi[10] == ValueApprox(-19.33872609471596));
  CHECK(dlogpsi[11] == ValueApprox(-47.00915390726143));
  CHECK(dlogpsi[12] == ValueApprox(-0.05463186141658209));
  CHECK(dlogpsi[13] == ValueApprox(0.045055811131004785));
  CHECK(dlogpsi[14] == ValueApprox(0.46675941272234));
  CHECK(dlogpsi[15] == ValueApprox(1.1352711502777513));


  CHECK(dhpsioverpsi[0] == ValueApprox(0.2761674423047662));
  CHECK(dhpsioverpsi[1] == ValueApprox(0.022999975062422046));
  CHECK(dhpsioverpsi[2] == ValueApprox(0.3572968312376671));
  CHECK(dhpsioverpsi[3] == ValueApprox(-5.459873357259045));
  CHECK(dhpsioverpsi[4] == ValueApprox(-0.792225084691375));
  CHECK(dhpsioverpsi[5] == ValueApprox(24.453138754349123));
  CHECK(dhpsioverpsi[6] == ValueApprox(70.0280297306038));
  CHECK(dhpsioverpsi[7] == ValueApprox(-1.0272848501840672));
  CHECK(dhpsioverpsi[8] == ValueApprox(17.514031530576368));
  CHECK(dhpsioverpsi[9] == ValueApprox(2.52887169464403));
  CHECK(dhpsioverpsi[10] == ValueApprox(-78.37945447401765));
  CHECK(dhpsioverpsi[11] == ValueApprox(-224.4814690906403));
  CHECK(dhpsioverpsi[12] == ValueApprox(-0.6346957697642424));
  CHECK(dhpsioverpsi[13] == ValueApprox(0.03270289146243591));
  CHECK(dhpsioverpsi[14] == ValueApprox(3.263830358386392));
  CHECK(dhpsioverpsi[15] == ValueApprox(8.944714289946793));
}

} // namespace qmcplusplus

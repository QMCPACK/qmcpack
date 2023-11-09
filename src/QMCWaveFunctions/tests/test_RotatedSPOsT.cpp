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

#include "FakeSPOT.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSetPoolT.h"
#include "Particle/ParticleSetT.h"
#include "QMCWaveFunctions/EinsplineSetBuilderT.h"
#include "QMCWaveFunctions/RotatedSPOsT.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "catch.hpp"
#include "checkMatrix.hpp"
#include "type_traits/ConvertToReal.h"
#include "type_traits/template_types.hpp"
#include <ResourceCollection.h>
#include <stdio.h>

#include <limits>
#include <string>
#include <tuple>

using std::string;

namespace qmcplusplus
{
template<typename T>
struct ValueApproxHelper
{
  using Type = Catch::Detail::Approx;
};
template<typename T>
struct ValueApproxHelper<std::complex<T>>
{
  using Type = Catch::Detail::ComplexApprox;
};

template<typename T>
using ValueApprox = typename ValueApproxHelper<T>::Type;

namespace testing
{
OptVariablesTypeT<float>& getMyVars(SPOSetT<float>& rot) { return rot.myVars; }
OptVariablesTypeT<double>& getMyVars(SPOSetT<double>& rot) { return rot.myVars; }
OptVariablesTypeT<float>& getMyVarsFull(RotatedSPOsT<float>& rot) { return rot.myVarsFull; }
OptVariablesTypeT<double>& getMyVarsFull(RotatedSPOsT<double>& rot) { return rot.myVarsFull; }
std::vector<std::vector<float>>& getHistoryParams(RotatedSPOsT<float>& rot) { return rot.history_params_; }

std::vector<std::vector<double>>& getHistoryParams(RotatedSPOsT<double>& rot) { return rot.history_params_; }
} // namespace testing

#ifndef QMC_COMPLEX
#ifndef MIXED_PRECISION
using TestTypeList = std::tuple<double>;
#else
using TestTypeList = std::tuple<float>;
#endif
#else
using TestTypeList = std::tuple<>;
#endif

/*
  JPT 04.01.2022: Adapted from test_einset.cpp
  Test the spline rotated machinery for SplineR2R (extend to others later).
*/
TEMPLATE_LIST_TEST_CASE("RotatedSPOs via SplineR2R", "[wavefunction][template]", TestTypeList)
{
  using RealType = typename SPOSetT<TestType>::RealType;

  /*
      BEGIN Boilerplate stuff to make a simple SPOSet. Copied from
      test_einset.cpp
    */

  Communicate* c = OHMMS::Controller;

  // We get a "Mismatched supercell lattices" error due to default ctor?
  typename ParticleSetT<TestType>::ParticleLayout lattice;

  // diamondC_1x1x1
  lattice.R = {3.37316115, 3.37316115, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPoolT<TestType> ptcl = ParticleSetPoolT<TestType>(c);
  ptcl.setSimulationCell(lattice);
  // LAttice seems fine after this point...

  auto ions_uptr = std::make_unique<ParticleSetT<TestType>>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSetT<TestType>>(ptcl.getSimulationCell());
  ParticleSetT<TestType>& ions_(*ions_uptr);
  ParticleSetT<TestType>& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};
  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0]                 = {0.0, 0.0, 0.0};
  elec_.R[1]                 = {0.0, 1.0, 0.0};
  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // diamondC_1x1x1 - 8 bands available
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="float" size="8"/>
</tmp>
)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilderT<TestType> einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  /*
      END Boilerplate stuff. Now we have a SplineR2R wavefunction
      ready for rotation. What follows is the actual test.
    */

  // SplineR2R only for the moment, so skip if QMC_COMPLEX is set
#if !defined(QMC_COMPLEX)

  spo->storeParamsBeforeRotation();
  // 1.) Make a RotatedSPOs object so that we can use the rotation routines
  auto rot_spo = std::make_unique<RotatedSPOsT<TestType>>("one_rotated_set", std::move(spo));

  // Sanity check for orbs. Expect 2 electrons, 8 orbitals, & 79507 coefs/orb.
  const auto orbitalsetsize = rot_spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 8);

  // 2.) Get data for unrotated orbitals. Check that there's no rotation
  rot_spo->buildOptVariables(elec_.R.size());
  typename SPOSetT<TestType>::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  rot_spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // This stuff checks that no rotation was applied. Copied from
  // test_einset.cpp. value
  CHECK(std::real(psiM_bare[1][0]) == ValueApprox<TestType>(-0.8886948824));
  CHECK(std::real(psiM_bare[1][1]) == ValueApprox<TestType>(1.4194120169));
  // grad
  CHECK(std::real(dpsiM_bare[1][0][0]) == ValueApprox<TestType>(-0.0000183403));
  CHECK(std::real(dpsiM_bare[1][0][1]) == ValueApprox<TestType>(0.1655139178));
  CHECK(std::real(dpsiM_bare[1][0][2]) == ValueApprox<TestType>(-0.0000193077));
  CHECK(std::real(dpsiM_bare[1][1][0]) == ValueApprox<TestType>(-1.3131694794));
  CHECK(std::real(dpsiM_bare[1][1][1]) == ValueApprox<TestType>(-1.1174004078));
  CHECK(std::real(dpsiM_bare[1][1][2]) == ValueApprox<TestType>(-0.8462534547));
  // lapl
  CHECK(std::real(d2psiM_bare[1][0]) == ValueApprox<TestType>(1.3313053846));
  CHECK(std::real(d2psiM_bare[1][1]) == ValueApprox<TestType>(-4.712583065));

  /*
       3.) Apply a rotation to the orbitals
           To do this, construct a params vector and call the
       RotatedSPOs::apply_rotation(params) method. That should do the
       right thing for this particular spline class.

       For 2 electrons in 8 orbs, we expect 2*(8-2) = 12 params.
    */
  const auto rot_size = rot_spo->m_act_rot_inds.size();
  REQUIRE(rot_size == 12); // = Nelec*(Norbs - Nelec) = 2*(8-2) = 12
  std::vector<RealType> param(rot_size);
  for (auto i = 0; i < rot_size; i++)
  {
    param[i] = 0.01 * static_cast<RealType>(i);
  }
  rot_spo->apply_rotation(param, false); // Expect this to call SplineR2R::applyRotation()

  // 4.) Get data for rotated orbitals.
  typename SPOSetT<TestType>::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  rot_spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  /*
       Manually encode the unitary transformation. Ugly, but it works.
       @TODO: Use the total rotation machinery when it's implemented

       NB: This is truncated to 5 sig-figs, so there is some slop here as
           compared to what is done in the splines via apply_rotation().
       So below we reduce the threshold for comparison. This can
       probably be ditched once we have a way to grab the actual
       rotation matrix...
    */
  typename SPOSetT<TestType>::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat[0][0] = 0.99726;
  rot_mat[0][1] = -0.00722;
  rot_mat[0][2] = 0.00014;
  rot_mat[0][3] = -0.00982;
  rot_mat[0][4] = -0.01979;
  rot_mat[0][5] = -0.02976;
  rot_mat[0][6] = -0.03972;
  rot_mat[0][7] = -0.04969;
  rot_mat[1][0] = -0.00722;
  rot_mat[1][1] = 0.97754;
  rot_mat[1][2] = -0.05955;
  rot_mat[1][3] = -0.06945;
  rot_mat[1][4] = -0.07935;
  rot_mat[1][5] = -0.08925;
  rot_mat[1][6] = -0.09915;
  rot_mat[1][7] = -0.10905;
  rot_mat[2][0] = -0.00014;
  rot_mat[2][1] = 0.05955;
  rot_mat[2][2] = 0.99821;
  rot_mat[2][3] = -0.00209;
  rot_mat[2][4] = -0.00239;
  rot_mat[2][5] = -0.00269;
  rot_mat[2][6] = -0.00299;
  rot_mat[2][7] = -0.00329;
  rot_mat[3][0] = 0.00982;
  rot_mat[3][1] = 0.06945;
  rot_mat[3][2] = -0.00209;
  rot_mat[3][3] = 0.99751;
  rot_mat[3][4] = -0.00289;
  rot_mat[3][5] = -0.00329;
  rot_mat[3][6] = -0.00368;
  rot_mat[3][7] = -0.00408;
  rot_mat[4][0] = 0.01979;
  rot_mat[4][1] = 0.07935;
  rot_mat[4][2] = -0.00239;
  rot_mat[4][3] = -0.00289;
  rot_mat[4][4] = 0.99661;
  rot_mat[4][5] = -0.00388;
  rot_mat[4][6] = -0.00438;
  rot_mat[4][7] = -0.00488;
  rot_mat[5][0] = 0.02976;
  rot_mat[5][1] = 0.08925;
  rot_mat[5][2] = -0.00269;
  rot_mat[5][3] = -0.00329;
  rot_mat[5][4] = -0.00388;
  rot_mat[5][5] = 0.99552;
  rot_mat[5][6] = -0.00508;
  rot_mat[5][7] = -0.00568;
  rot_mat[6][0] = 0.03972;
  rot_mat[6][1] = 0.09915;
  rot_mat[6][2] = -0.00299;
  rot_mat[6][3] = -0.00368;
  rot_mat[6][4] = -0.00438;
  rot_mat[6][5] = -0.00508;
  rot_mat[6][6] = 0.99422;
  rot_mat[6][7] = -0.00647;
  rot_mat[7][0] = 0.04969;
  rot_mat[7][1] = 0.10905;
  rot_mat[7][2] = -0.00329;
  rot_mat[7][3] = -0.00408;
  rot_mat[7][4] = -0.00488;
  rot_mat[7][5] = -0.00568;
  rot_mat[7][6] = -0.00647;
  rot_mat[7][7] = 0.99273;

  // Now compute the expected values by hand using the transformation above
  double val1 = 0.;
  double val2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    val1 += psiM_bare[0][i] * rot_mat[i][0];
    val2 += psiM_bare[1][i] * rot_mat[i][0];
  }

  // value
  CHECK(std::real(psiM_rot[0][0]) == ValueApprox<TestType>(val1));
  CHECK(std::real(psiM_rot[1][0]) == ValueApprox<TestType>(val2));

  std::vector<double> grad1(3);
  std::vector<double> grad2(3);
  for (auto j = 0; j < grad1.size(); j++)
  {
    for (auto i = 0; i < rot_mat.size1(); i++)
    {
      grad1[j] += dpsiM_bare[0][i][j] * rot_mat[i][0];
      grad2[j] += dpsiM_bare[1][i][j] * rot_mat[i][0];
    }
  }

  // grad
  CHECK(dpsiM_rot[0][0][0] == ValueApprox<TestType>(grad1[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][1] == ValueApprox<TestType>(grad1[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][2] == ValueApprox<TestType>(grad1[2]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][0] == ValueApprox<TestType>(grad2[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][1] == ValueApprox<TestType>(grad2[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][2] == ValueApprox<TestType>(grad2[2]).epsilon(0.0001));

  double lap1 = 0.;
  double lap2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    lap1 += d2psiM_bare[0][i] * rot_mat[i][0];
    lap2 += d2psiM_bare[1][i] * rot_mat[i][0];
  }

  // Lapl
  CHECK(std::real(d2psiM_rot[0][0]) == ValueApprox<TestType>(lap1).epsilon(0.0001));
  CHECK(std::real(d2psiM_rot[1][0]) == ValueApprox<TestType>(lap2).epsilon(0.0001));

#endif
}

TEMPLATE_LIST_TEST_CASE("RotatedSPOs createRotationIndices", "[wavefunction][template]", TestTypeList)
{
  // No active-active or virtual-virtual rotations
  // Only active-virtual
  typename RotatedSPOsT<TestType>::RotationIndices rot_ind;
  int nel = 1;
  int nmo = 3;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind);
  CHECK(rot_ind.size() == 2);

  // Full rotation contains all rotations
  // Size should be number of pairs of orbitals: nmo*(nmo-1)/2
  typename RotatedSPOsT<TestType>::RotationIndices full_rot_ind;
  RotatedSPOsT<TestType>::createRotationIndicesFull(nel, nmo, full_rot_ind);
  CHECK(full_rot_ind.size() == 3);

  nel = 2;
  typename RotatedSPOsT<TestType>::RotationIndices rot_ind2;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind2);
  CHECK(rot_ind2.size() == 2);

  typename RotatedSPOsT<TestType>::RotationIndices full_rot_ind2;
  RotatedSPOsT<TestType>::createRotationIndicesFull(nel, nmo, full_rot_ind2);
  CHECK(full_rot_ind2.size() == 3);

  nmo = 4;
  typename RotatedSPOsT<TestType>::RotationIndices rot_ind3;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind3);
  CHECK(rot_ind3.size() == 4);

  typename RotatedSPOsT<TestType>::RotationIndices full_rot_ind3;
  RotatedSPOsT<TestType>::createRotationIndicesFull(nel, nmo, full_rot_ind3);
  CHECK(full_rot_ind3.size() == 6);
}

TEMPLATE_LIST_TEST_CASE("RotatedSPOs constructAntiSymmetricMatrix", "[wavefunction][template]", TestTypeList)
{
  using ValueType   = typename SPOSetT<TestType>::ValueType;
  using ValueMatrix = typename SPOSetT<TestType>::ValueMatrix;

  typename RotatedSPOsT<TestType>::RotationIndices rot_ind;
  int nel = 1;
  int nmo = 3;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind);

  ValueMatrix m3(nmo, nmo);
  m3                            = ValueType(0);
  std::vector<ValueType> params = {0.1, 0.2};

  RotatedSPOsT<TestType>::constructAntiSymmetricMatrix(rot_ind, params, m3);

  // clang-format off
  std::vector<ValueType> expected_data = { 0.0,  -0.1, -0.2,
                                           0.1,   0.0,  0.0,
                                           0.2,   0.0,  0.0 };
  // clang-format on

  ValueMatrix expected_m3(expected_data.data(), 3, 3);

  CheckMatrixResult check_matrix_result = checkMatrix(m3, expected_m3, true);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }

  std::vector<ValueType> params_out(2);
  RotatedSPOsT<TestType>::extractParamsFromAntiSymmetricMatrix(rot_ind, m3, params_out);
  CHECK(params_out[0] == ValueApprox<TestType>(0.1));
  CHECK(params_out[1] == ValueApprox<TestType>(0.2));
}

// Expected values of the matrix exponential come from gen_matrix_ops.py
TEMPLATE_LIST_TEST_CASE("RotatedSPOs exponentiate matrix", "[wavefunction][template]", TestTypeList)
{
  using ValueType   = typename SPOSetT<TestType>::ValueType;
  using ValueMatrix = typename SPOSetT<TestType>::ValueMatrix;

  std::vector<typename SPOSetT<TestType>::ValueType> mat1_data = {0.0};
  typename SPOSetT<TestType>::ValueMatrix m1(mat1_data.data(), 1, 1);
  RotatedSPOsT<TestType>::exponentiate_antisym_matrix(m1);
  // Always return 1.0 (the only possible anti-symmetric 1x1 matrix is 0)
  CHECK(m1(0, 0) == ValueApprox<TestType>(1.0));

  // clang-format off
  std::vector<typename SPOSetT<TestType>::ValueType> mat2_data = { 0.0, -0.1,
                                               0.1,  0.0 };
  // clang-format on

  typename SPOSetT<TestType>::ValueMatrix m2(mat2_data.data(), 2, 2);
  RotatedSPOsT<TestType>::exponentiate_antisym_matrix(m2);

  // clang-format off
  std::vector<ValueType> expected_rot2 = {  0.995004165278026,  -0.0998334166468282,
                                            0.0998334166468282,  0.995004165278026 };
  // clang-format on

  ValueMatrix expected_m2(expected_rot2.data(), 2, 2);
  CheckMatrixResult check_matrix_result2 = checkMatrix(m2, expected_m2, true);
  CHECKED_ELSE(check_matrix_result2.result) { FAIL(check_matrix_result2.result_message); }

  // clang-format off
  std::vector<ValueType> m3_input_data = { 0.0,  -0.3, -0.1,
                                           0.3,   0.0, -0.2,
                                           0.1,   0.2,  0.0 };


  std::vector<ValueType> expected_rot3 = {  0.950580617906092, -0.302932713402637, -0.0680313164049401,
                                            0.283164960565074,  0.935754803277919, -0.210191705950743,
                                            0.127334574917630,  0.180540076694398,  0.975290308953046 };

  // clang-format on

  ValueMatrix m3(m3_input_data.data(), 3, 3);
  ValueMatrix expected_m3(expected_rot3.data(), 3, 3);

  RotatedSPOsT<TestType>::exponentiate_antisym_matrix(m3);

  CheckMatrixResult check_matrix_result3 = checkMatrix(m3, expected_m3, true);
  CHECKED_ELSE(check_matrix_result3.result) { FAIL(check_matrix_result3.result_message); }
}

TEMPLATE_LIST_TEST_CASE("RotatedSPOs log matrix", "[wavefunction][template]", TestTypeList)
{
  using ValueType   = typename SPOSetT<TestType>::ValueType;
  using ValueMatrix = typename SPOSetT<TestType>::ValueMatrix;

  std::vector<typename SPOSetT<TestType>::ValueType> mat1_data = {1.0};
  typename SPOSetT<TestType>::ValueMatrix m1(mat1_data.data(), 1, 1);
  typename SPOSetT<TestType>::ValueMatrix out_m1(1, 1);
  RotatedSPOsT<TestType>::log_antisym_matrix(m1, out_m1);
  // Should always be 1.0 (the only possible anti-symmetric 1x1 matrix is 0)
  CHECK(out_m1(0, 0) == ValueApprox<TestType>(0.0));

  // clang-format off
  std::vector<ValueType> start_rot2 = {  0.995004165278026,  -0.0998334166468282,
                                         0.0998334166468282,  0.995004165278026 };

  std::vector<typename SPOSetT<TestType>::ValueType> mat2_data = { 0.0, -0.1,
                                               0.1,  0.0 };
  // clang-format on

  ValueMatrix rot_m2(start_rot2.data(), 2, 2);
  ValueMatrix out_m2(2, 2);
  RotatedSPOsT<TestType>::log_antisym_matrix(rot_m2, out_m2);

  typename SPOSetT<TestType>::ValueMatrix m2(mat2_data.data(), 2, 2);
  CheckMatrixResult check_matrix_result2 = checkMatrix(m2, out_m2, true);
  CHECKED_ELSE(check_matrix_result2.result) { FAIL(check_matrix_result2.result_message); }

  // clang-format off
  std::vector<ValueType> start_rot3 = {  0.950580617906092, -0.302932713402637, -0.0680313164049401,
                                         0.283164960565074,  0.935754803277919, -0.210191705950743,
                                         0.127334574917630,  0.180540076694398,  0.975290308953046 };

  std::vector<ValueType> m3_input_data = { 0.0,  -0.3, -0.1,
                                           0.3,   0.0, -0.2,
                                           0.1,   0.2,  0.0 };
  // clang-format on
  ValueMatrix rot_m3(start_rot3.data(), 3, 3);
  ValueMatrix out_m3(3, 3);
  RotatedSPOsT<TestType>::log_antisym_matrix(rot_m3, out_m3);

  typename SPOSetT<TestType>::ValueMatrix m3(m3_input_data.data(), 3, 3);
  CheckMatrixResult check_matrix_result3 = checkMatrix(m3, out_m3, true);
  CHECKED_ELSE(check_matrix_result3.result) { FAIL(check_matrix_result3.result_message); }
}

// Test round trip A -> exp(A) -> log(exp(A))
// The log is multi-valued so this test may fail if the rotation parameters are
// too large. The exponentials will be the same, though
//   exp(log(exp(A))) == exp(A)
TEMPLATE_LIST_TEST_CASE("RotatedSPOs exp-log matrix", "[wavefunction][template]", TestTypeList)
{
  using ValueType   = typename SPOSetT<TestType>::ValueType;
  using ValueMatrix = typename SPOSetT<TestType>::ValueMatrix;

  typename RotatedSPOsT<TestType>::RotationIndices rot_ind;
  int nel = 2;
  int nmo = 4;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind);

  ValueMatrix rot_m4(nmo, nmo);
  rot_m4 = ValueType(0);

  std::vector<ValueType> params4 = {-1.1, 1.5, 0.2, -0.15};

  RotatedSPOsT<TestType>::constructAntiSymmetricMatrix(rot_ind, params4, rot_m4);
  ValueMatrix orig_rot_m4 = rot_m4;
  ValueMatrix out_m4(nmo, nmo);

  RotatedSPOsT<TestType>::exponentiate_antisym_matrix(rot_m4);

  RotatedSPOsT<TestType>::log_antisym_matrix(rot_m4, out_m4);

  CheckMatrixResult check_matrix_result4 = checkMatrix(out_m4, orig_rot_m4, true);
  CHECKED_ELSE(check_matrix_result4.result) { FAIL(check_matrix_result4.result_message); }

  std::vector<ValueType> params4out(4);
  RotatedSPOsT<TestType>::extractParamsFromAntiSymmetricMatrix(rot_ind, out_m4, params4out);
  for (int i = 0; i < params4.size(); i++)
  {
    CHECK(params4[i] == ValueApprox<TestType>(params4out[i]));
  }
}

TEMPLATE_LIST_TEST_CASE("RotatedSPOs hcpBe", "[wavefunction][template]", TestTypeList)
{
  using RealType = typename OrbitalSetTraits<TestType>::RealType;
  Communicate* c = OHMMS::Controller;

  typename ParticleSetT<TestType>::ParticleLayout lattice;
  lattice.R = {4.32747284, 0.00000000, 0.00000000, -2.16373642, 3.74770142,
               0.00000000, 0.00000000, 0.00000000, 6.78114995};

  ParticleSetPoolT<TestType> ptcl = ParticleSetPoolT<TestType>(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSetT<TestType>>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSetT<TestType>>(ptcl.getSimulationCell());
  ParticleSetT<TestType>& ions(*ions_uptr);
  ParticleSetT<TestType>& elec(*elec_uptr);

  ions.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions.create({1});
  ions.R[0] = {0.0, 0.0, 0.0};

  elec.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // Add the attribute save_coefs="yes" to the sposet_builder tag to generate
  // the spline file for use in eval_bspline_spo.py

  const char* particles = R"(<tmp>
<sposet_builder type="bspline" href="hcpBe.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="double">
      <sposet type="bspline" name="spo_ud" spindataset="0" size="2"/>
</sposet_builder>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr sposet_builder = xmlFirstElementChild(root);
  xmlNodePtr sposet_ptr     = xmlFirstElementChild(sposet_builder);

  EinsplineSetBuilderT<TestType> einSet(elec, ptcl.getPool(), c, sposet_builder);
  auto spo = einSet.createSPOSetFromXML(sposet_ptr);
  REQUIRE(spo);

  spo->storeParamsBeforeRotation();
  auto rot_spo = std::make_unique<RotatedSPOsT<TestType>>("one_rotated_set", std::move(spo));

  // Sanity check for orbs. Expect 1 electron, 2 orbitals
  const auto orbitalsetsize = rot_spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 2);

  rot_spo->buildOptVariables(elec.R.size());

  typename SPOSetT<TestType>::ValueMatrix psiM_bare(elec.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::GradMatrix dpsiM_bare(elec.R.size(), orbitalsetsize);
  typename SPOSetT<TestType>::ValueMatrix d2psiM_bare(elec.R.size(), orbitalsetsize);
  rot_spo->evaluate_notranspose(elec, 0, elec.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Values generated from eval_bspline_spo.py, the
  // generate_point_values_hcpBe function
  CHECK(std::real(psiM_bare[0][0]) == ValueApprox<TestType>(0.210221765375514));
  CHECK(std::real(psiM_bare[0][1]) == ValueApprox<TestType>(-2.984345024542937e-06));

  CHECK(std::real(d2psiM_bare[0][0]) == ValueApprox<TestType>(5.303848362116568));

  OptVariablesTypeT<TestType> opt_vars;
  rot_spo->checkInVariablesExclusive(opt_vars);
  opt_vars.resetIndex();
  rot_spo->checkOutVariables(opt_vars);
  rot_spo->resetParametersExclusive(opt_vars);

  using ValueType = TestType;
  Vector<ValueType> dlogpsi(1);
  Vector<ValueType> dhpsioverpsi(1);
  rot_spo->evaluateDerivatives(elec, opt_vars, dlogpsi, dhpsioverpsi, 0, 1);

  CHECK(dlogpsi[0] == ValueApprox<TestType>(-1.41961753e-05));
  CHECK(dhpsioverpsi[0] == ValueApprox<TestType>(-0.00060853));

  std::vector<RealType> params = {0.1};
  rot_spo->apply_rotation(params, false);

  rot_spo->evaluate_notranspose(elec, 0, elec.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);
  CHECK(std::real(psiM_bare[0][0]) == ValueApprox<TestType>(0.20917123424337608));
  CHECK(std::real(psiM_bare[0][1]) == ValueApprox<TestType>(-0.02099012652669549));

  CHECK(std::real(d2psiM_bare[0][0]) == ValueApprox<TestType>(5.277362065087747));

  dlogpsi[0]      = 0.0;
  dhpsioverpsi[0] = 0.0;

  rot_spo->evaluateDerivatives(elec, opt_vars, dlogpsi, dhpsioverpsi, 0, 1);
  CHECK(dlogpsi[0] == ValueApprox<TestType>(-0.10034901119468914));
  CHECK(dhpsioverpsi[0] == ValueApprox<TestType>(32.96939041498753));
}

// Test construction of delta rotation
TEMPLATE_LIST_TEST_CASE("RotatedSPOs construct delta matrix", "[wavefunction][template]", TestTypeList)
{
  using ValueType   = typename SPOSetT<TestType>::ValueType;
  using ValueMatrix = typename SPOSetT<TestType>::ValueMatrix;

  int nel = 2;
  int nmo = 4;
  typename RotatedSPOsT<TestType>::RotationIndices rot_ind;
  RotatedSPOsT<TestType>::createRotationIndices(nel, nmo, rot_ind);
  typename RotatedSPOsT<TestType>::RotationIndices full_rot_ind;
  RotatedSPOsT<TestType>::createRotationIndicesFull(nel, nmo, full_rot_ind);
  // rot_ind size is 4 and full rot_ind size is 6

  ValueMatrix rot_m4(nmo, nmo);
  rot_m4 = ValueType(0);

  // When comparing with gen_matrix_ops.py, be aware of the order of indices
  // in full_rot
  // rot_ind is (0,2) (0,3) (1,2) (1,3)
  // full_rot_ind is (0,2) (0,3) (1,2) (1,3) (0,1) (2,3)
  // The extra indices go at the back
  std::vector<ValueType> old_params   = {1.5, 0.2, -0.15, 0.03, -1.1, 0.05};
  std::vector<ValueType> delta_params = {0.1, 0.3, 0.2, -0.1};
  std::vector<ValueType> new_params(6);

  RotatedSPOsT<TestType>::constructDeltaRotation(delta_params, old_params, rot_ind, full_rot_ind, new_params, rot_m4);

  // clang-format off
  std::vector<ValueType> rot_data4 =
    { -0.371126931484737,  0.491586564957393,   -0.784780958819798,   0.0687480658200083,
      -0.373372784561548,  0.66111547793048,     0.610450337985578,   0.225542620014052,
       0.751270334458895,  0.566737323353515,   -0.0297901110611425, -0.336918744155143,
       0.398058348785074,  0.00881931472604944, -0.102867783149713,   0.911531672428406 };
  // clang-format on

  ValueMatrix new_rot_m4(rot_data4.data(), 4, 4);

  CheckMatrixResult check_matrix_result4 = checkMatrix(rot_m4, new_rot_m4, true);
  CHECKED_ELSE(check_matrix_result4.result) { FAIL(check_matrix_result4.result_message); }

  // Reminder: Ordering!
  std::vector<ValueType> expected_new_param = {1.6813965019790489,   0.3623564254653294,  -0.05486544454559908,
                                               -0.20574472941408453, -0.9542513302873077, 0.27497788909911774};
  for (int i = 0; i < new_params.size(); i++)
    CHECK(new_params[i] == ValueApprox<TestType>(expected_new_param[i]));

  // Rotated back to original position

  std::vector<ValueType> new_params2(6);
  std::vector<ValueType> reverse_delta_params = {-0.1, -0.3, -0.2, 0.1};
  RotatedSPOsT<TestType>::constructDeltaRotation(reverse_delta_params, new_params, rot_ind, full_rot_ind, new_params2,
                                                 rot_m4);
  for (int i = 0; i < new_params2.size(); i++)
    CHECK(new_params2[i] == ValueApprox<TestType>(old_params[i]));
}

// Test using global rotation
TEMPLATE_LIST_TEST_CASE("RotatedSPOs read and write parameters", "[wavefunction][template]", TestTypeList)
{
  auto fake_spo = std::make_unique<FakeSPOT<TestType>>();
  fake_spo->setOrbitalSetSize(4);
  RotatedSPOsT<TestType> rot("fake_rot", std::move(fake_spo));
  int nel = 2;
  rot.buildOptVariables(nel);

  optimize::VariableSetT<TestType> vs;
  rot.checkInVariablesExclusive(vs);
  vs[0] = 0.1;
  vs[1] = 0.15;
  vs[2] = 0.2;
  vs[3] = 0.25;
  rot.resetParametersExclusive(vs);

  {
    hdf_archive hout;
    vs.writeToHDF("rot_vp.h5", hout);

    rot.writeVariationalParameters(hout);
  }

  auto fake_spo2 = std::make_unique<FakeSPOT<TestType>>();
  fake_spo2->setOrbitalSetSize(4);

  RotatedSPOsT<TestType> rot2("fake_rot", std::move(fake_spo2));
  rot2.buildOptVariables(nel);

  optimize::VariableSetT<TestType> vs2;
  rot2.checkInVariablesExclusive(vs2);

  hdf_archive hin;
  vs2.readFromHDF("rot_vp.h5", hin);
  rot2.readVariationalParameters(hin);

  auto& var = testing::getMyVars(rot2);
  CHECK(var[0] == ValueApprox<TestType>(vs[0]));
  CHECK(var[1] == ValueApprox<TestType>(vs[1]));
  CHECK(var[2] == ValueApprox<TestType>(vs[2]));
  CHECK(var[3] == ValueApprox<TestType>(vs[3]));

  auto& full_var = testing::getMyVarsFull(rot2);
  CHECK(full_var[0] == ValueApprox<TestType>(vs[0]));
  CHECK(full_var[1] == ValueApprox<TestType>(vs[1]));
  CHECK(full_var[2] == ValueApprox<TestType>(vs[2]));
  CHECK(full_var[3] == ValueApprox<TestType>(vs[3]));
  CHECK(full_var[4] == ValueApprox<TestType>(0.0));
  CHECK(full_var[5] == ValueApprox<TestType>(0.0));
}

// Test using history list.
TEMPLATE_LIST_TEST_CASE("RotatedSPOs read and write parameters history", "[wavefunction][template]", TestTypeList)
{
  auto fake_spo = std::make_unique<FakeSPOT<TestType>>();
  fake_spo->setOrbitalSetSize(4);
  RotatedSPOsT<TestType> rot("fake_rot", std::move(fake_spo));
  rot.set_use_global_rotation(false);
  int nel = 2;
  rot.buildOptVariables(nel);

  optimize::VariableSetT<TestType> vs;
  rot.checkInVariablesExclusive(vs);
  vs[0] = 0.1;
  vs[1] = 0.15;
  vs[2] = 0.2;
  vs[3] = 0.25;
  rot.resetParametersExclusive(vs);

  {
    hdf_archive hout;
    vs.writeToHDF("rot_vp_hist.h5", hout);

    rot.writeVariationalParameters(hout);
  }

  auto fake_spo2 = std::make_unique<FakeSPOT<TestType>>();
  fake_spo2->setOrbitalSetSize(4);

  RotatedSPOsT<TestType> rot2("fake_rot", std::move(fake_spo2));
  rot2.buildOptVariables(nel);

  optimize::VariableSetT<TestType> vs2;
  rot2.checkInVariablesExclusive(vs2);

  hdf_archive hin;
  vs2.readFromHDF("rot_vp_hist.h5", hin);
  rot2.readVariationalParameters(hin);

  auto& var = testing::getMyVars(rot2);
  CHECK(var[0] == ValueApprox<TestType>(vs[0]));
  CHECK(var[1] == ValueApprox<TestType>(vs[1]));
  CHECK(var[2] == ValueApprox<TestType>(vs[2]));
  CHECK(var[3] == ValueApprox<TestType>(vs[3]));

  auto hist = testing::getHistoryParams(rot2);
  REQUIRE(hist.size() == 1);
  REQUIRE(hist[0].size() == 4);
}

template<typename T>
class DummySPOSetWithoutMWT : public SPOSetT<T>
{
public:
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;

  DummySPOSetWithoutMWT(const std::string& my_name) : SPOSetT<T>(my_name) {}
  void setOrbitalSetSize(int norbs) override {}
  void evaluateValue(const ParticleSetT<T>& P, int iat, typename SPOSetT<T>::ValueVector& psi) override
  {
    assert(psi.size() == 3);
    psi[0] = 123;
    psi[1] = 456;
    psi[2] = 789;
  }
  void evaluateVGL(const ParticleSetT<T>& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {}
  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override
  {}
  std::string getClassName() const override { return this->my_name_; }
};

template<typename T>
class DummySPOSetWithMWT : public DummySPOSetWithoutMWT<T>
{
public:
  using ValueVector = typename DummySPOSetWithoutMWT<T>::ValueVector;

  DummySPOSetWithMWT(const std::string& my_name) : DummySPOSetWithoutMWT<T>(my_name) {}
  void mw_evaluateValue(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
                        const RefVectorWithLeader<ParticleSetT<T>>& P_list,
                        int iat,
                        const RefVector<ValueVector>& psi_v_list) const override
  {
    for (auto& psi : psi_v_list)
    {
      assert(psi.get().size() == 3);
      psi.get()[0] = 321;
      psi.get()[1] = 654;
      psi.get()[2] = 987;
    }
  }
};

TEMPLATE_LIST_TEST_CASE("RotatedSPOs mw_ APIs", "[wavefunction][template]", TestTypeList)
{
  // checking that mw_ API works in RotatedSPOs and is not defaulting to
  // SPOSet default implementation
  {
    // First check calling the mw_ APIs for RotatedSPOs, for which the
    // underlying implementation just calls the underlying SPOSet mw_ API
    // In the case that the underlying SPOSet doesn't specialize the mw_
    // API, the underlying SPOSet will fall back to the default SPOSet mw_,
    // which is just a loop over the single walker API.
    RotatedSPOsT<TestType> rot_spo0("rotated0", std::make_unique<DummySPOSetWithoutMWT<TestType>>("no mw 0"));
    RotatedSPOsT<TestType> rot_spo1("rotated1", std::make_unique<DummySPOSetWithoutMWT<TestType>>("no mw 1"));
    RefVectorWithLeader<SPOSetT<TestType>> spo_list(rot_spo0, {rot_spo0, rot_spo1});

    ResourceCollection spo_res("test_rot_res");
    rot_spo0.createResource(spo_res);
    ResourceCollectionTeamLock<SPOSetT<TestType>> mw_sposet_lock(spo_res, spo_list);

    const SimulationCellT<TestType> simulation_cell;
    ParticleSetT<TestType> elec0(simulation_cell);
    ParticleSetT<TestType> elec1(simulation_cell);
    RefVectorWithLeader<ParticleSetT<TestType>> p_list(elec0, {elec0, elec1});

    typename SPOSetT<TestType>::ValueVector psi0(3);
    typename SPOSetT<TestType>::ValueVector psi1(3);
    RefVector<typename SPOSetT<TestType>::ValueVector> psi_v_list{psi0, psi1};

    rot_spo0.mw_evaluateValue(spo_list, p_list, 0, psi_v_list);
    for (int iw = 0; iw < spo_list.size(); iw++)
    {
      CHECK(psi_v_list[iw].get()[0] == ValueApprox<TestType>(123));
      CHECK(psi_v_list[iw].get()[1] == ValueApprox<TestType>(456));
      CHECK(psi_v_list[iw].get()[2] == ValueApprox<TestType>(789));
    }
  }
  {
    // In the case that the underlying SPOSet DOES have mw_ specializations,
    // we want to make sure that RotatedSPOs are triggering that
    // appropriately This will mean that the underlying SPOSets will do the
    // appropriate offloading To check this, DummySPOSetWithMW has an
    // explicit mw_evaluateValue which sets different values than what gets
    // set in evaluateValue. By doing this, we are ensuring that
    // RotatedSPOs->mw_evaluaeValue is calling the specialization in the
    // underlying SPO and not using the default SPOSet implementation which
    // loops over single walker APIs (which have different values enforced
    // in
    //  DummySPOSetWithoutMW

    RotatedSPOsT<TestType> rot_spo0("rotated0", std::make_unique<DummySPOSetWithMWT<TestType>>("mw 0"));
    RotatedSPOsT<TestType> rot_spo1("rotated1", std::make_unique<DummySPOSetWithMWT<TestType>>("mw 1"));
    RefVectorWithLeader<SPOSetT<TestType>> spo_list(rot_spo0, {rot_spo0, rot_spo1});

    ResourceCollection spo_res("test_rot_res");
    rot_spo0.createResource(spo_res);
    ResourceCollectionTeamLock<SPOSetT<TestType>> mw_sposet_lock(spo_res, spo_list);

    const SimulationCellT<TestType> simulation_cell;
    ParticleSetT<TestType> elec0(simulation_cell);
    ParticleSetT<TestType> elec1(simulation_cell);
    RefVectorWithLeader<ParticleSetT<TestType>> p_list(elec0, {elec0, elec1});

    typename SPOSetT<TestType>::ValueVector psi0(3);
    typename SPOSetT<TestType>::ValueVector psi1(3);
    RefVector<typename SPOSetT<TestType>::ValueVector> psi_v_list{psi0, psi1};

    rot_spo0.mw_evaluateValue(spo_list, p_list, 0, psi_v_list);
    for (int iw = 0; iw < spo_list.size(); iw++)
    {
      CHECK(psi_v_list[iw].get()[0] == ValueApprox<TestType>(321));
      CHECK(psi_v_list[iw].get()[1] == ValueApprox<TestType>(654));
      CHECK(psi_v_list[iw].get()[2] == ValueApprox<TestType>(987));
    }
  }
}

} // namespace qmcplusplus
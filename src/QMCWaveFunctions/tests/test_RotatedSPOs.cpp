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
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/RotatedSPOs.h"
#include "checkMatrix.hpp"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
/*
  JPT 04.01.2022: Adapted from test_einset.cpp
  Test the spline rotated machinery for SplineR2R (extend to others later).
*/
TEST_CASE("RotatedSPOs via SplineR2R", "[wavefunction]")
{
  using RealType = QMCTraits::RealType;

  /*
    BEGIN Boilerplate stuff to make a simple SPOSet. Copied from test_einset.cpp
  */

  Communicate* c = OHMMS::Controller;

  // We get a "Mismatched supercell lattices" error due to default ctor?
  ParticleSet::ParticleLayout lattice;

  // diamondC_1x1x1
  lattice.R(0, 0) = 3.37316115;
  lattice.R(0, 1) = 3.37316115;
  lattice.R(0, 2) = 0.0;
  lattice.R(1, 0) = 0.0;
  lattice.R(1, 1) = 3.37316115;
  lattice.R(1, 2) = 3.37316115;
  lattice.R(2, 0) = 3.37316115;
  lattice.R(2, 1) = 0.0;
  lattice.R(2, 2) = 3.37316115;

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  // LAttice seems fine after this point...

  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0][0] = 0.0;
  ions_.R[0][1] = 0.0;
  ions_.R[0][2] = 0.0;
  ions_.R[1][0] = 1.68658058;
  ions_.R[1][1] = 1.68658058;
  ions_.R[1][2] = 1.68658058;


  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0][0] = 0.0;
  elec_.R[0][1] = 0.0;
  elec_.R[0][2] = 0.0;
  elec_.R[1][0] = 0.0;
  elec_.R[1][1] = 1.0;
  elec_.R[1][2] = 0.0;

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  //diamondC_1x1x1 - 8 bands available
  const char* particles = "<tmp> \
<determinantset type=\"einspline\" href=\"diamondC_1x1x1.pwscf.h5\" tilematrix=\"1 0 0 0 1 0 0 0 1\" twistnum=\"0\" source=\"ion\" meshfactor=\"1.0\" precision=\"float\" size=\"8\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr ein1 = xmlFirstElementChild(root);

  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  /*
    END Boilerplate stuff. Now we have a SplineR2R wavefunction 
    ready for rotation. What follows is the actual test.
  */

  // SplineR2R only for the moment, so skip if QMC_COMPLEX is set
#if !defined(QMC_CUDA) && !defined(QMC_COMPLEX)

  // 1.) Make a RotatedSPOs object so that we can use the rotation routines
  auto rot_spo = std::make_unique<RotatedSPOs>(std::move(spo));

  // Sanity check for orbs. Expect 2 electrons, 8 orbitals, & 79507 coefs/orb.
  const auto orbitalsetsize = rot_spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 8);

  // 2.) Get data for unrotated orbitals. Check that there's no rotation
  rot_spo->buildOptVariables(elec_.R.size());
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  rot_spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // This stuff checks that no rotation was applied. Copied from test_einset.cpp.
  // value
  CHECK(std::real(psiM_bare[1][0]) == Approx(-0.8886948824));
  CHECK(std::real(psiM_bare[1][1]) == Approx(1.4194120169));
  // grad
  CHECK(std::real(dpsiM_bare[1][0][0]) == Approx(-0.0000183403));
  CHECK(std::real(dpsiM_bare[1][0][1]) == Approx(0.1655139178));
  CHECK(std::real(dpsiM_bare[1][0][2]) == Approx(-0.0000193077));
  CHECK(std::real(dpsiM_bare[1][1][0]) == Approx(-1.3131694794));
  CHECK(std::real(dpsiM_bare[1][1][1]) == Approx(-1.1174004078));
  CHECK(std::real(dpsiM_bare[1][1][2]) == Approx(-0.8462534547));
  // lapl
  CHECK(std::real(d2psiM_bare[1][0]) == Approx(1.3313053846));
  CHECK(std::real(d2psiM_bare[1][1]) == Approx(-4.712583065));

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
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
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
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
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
  CHECK(std::real(psiM_rot[0][0]) == Approx(val1));
  CHECK(std::real(psiM_rot[1][0]) == Approx(val2));

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
  CHECK(dpsiM_rot[0][0][0] == Approx(grad1[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][1] == Approx(grad1[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[0][0][2] == Approx(grad1[2]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][0] == Approx(grad2[0]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][1] == Approx(grad2[1]).epsilon(0.0001));
  CHECK(dpsiM_rot[1][0][2] == Approx(grad2[2]).epsilon(0.0001));

  double lap1 = 0.;
  double lap2 = 0.;
  for (auto i = 0; i < rot_mat.size1(); i++)
  {
    lap1 += d2psiM_bare[0][i] * rot_mat[i][0];
    lap2 += d2psiM_bare[1][i] * rot_mat[i][0];
  }

  // Lapl
  CHECK(std::real(d2psiM_rot[0][0]) == Approx(lap1).epsilon(0.0001));
  CHECK(std::real(d2psiM_rot[1][0]) == Approx(lap2).epsilon(0.0001));

#endif
}

TEST_CASE("RotatedSPOs createRotationIndices", "[wavefunction]")
{
  // No active-active or virtual-virtual rotations
  // Only active-virtual
  RotatedSPOs::RotationIndices rot_ind;
  int nel = 1;
  int nmo = 3;
  RotatedSPOs::createRotationIndices(nel, nmo, rot_ind);

  CHECK(rot_ind.size() == 2);

  nel = 2;
  RotatedSPOs::RotationIndices rot_ind2;
  RotatedSPOs::createRotationIndices(nel, nmo, rot_ind2);
  CHECK(rot_ind2.size() == 2);

  nmo = 4;
  RotatedSPOs::RotationIndices rot_ind3;
  RotatedSPOs::createRotationIndices(nel, nmo, rot_ind3);
  CHECK(rot_ind3.size() == 4);
}

TEST_CASE("RotatedSPOs constructAntiSymmetricMatrix", "[wavefunction]")
{
  using ValueType   = SPOSet::ValueType;
  using ValueMatrix = SPOSet::ValueMatrix;

  RotatedSPOs::RotationIndices rot_ind;
  int nel = 1;
  int nmo = 3;
  RotatedSPOs::createRotationIndices(nel, nmo, rot_ind);

  ValueMatrix m3(nmo, nmo);
  m3                            = ValueType(0);
  std::vector<ValueType> params = {0.1, 0.2};

  RotatedSPOs::constructAntiSymmetricMatrix(rot_ind, params, m3);

  // clang-format off
  std::vector<ValueType> expected_data = { 0.0,  -0.1, -0.2,
                                           0.1,   0.0,  0.0,
                                           0.2,   0.0,  0.0 };
  // clang-format on

  ValueMatrix expected_m3(expected_data.data(), 3, 3);

  CheckMatrixResult check_matrix_result = checkMatrix(m3, expected_m3, true);
  CHECKED_ELSE(check_matrix_result.result) { FAIL(check_matrix_result.result_message); }
}

} // namespace qmcplusplus

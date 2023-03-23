//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//
// File created by: Joshua Townsend, jptowns@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/SplineC2C.h"
#include "Utilities/for_testing/checkMatrix.hpp"

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{
  TEST_CASE("Spline applyRotation one rotation", "[wavefunction]")
{
  // How to get rid of all this annoying boilerplate?
  Communicate* c = OHMMS::Controller;

  // diamondC_1x1x1
  ParticleSet::ParticleLayout lattice;
  lattice.R = {3.37316115, 3.37316115, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};

  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({2});
  ions_.R[0] = {0.0, 0.0, 0.0};
  ions_.R[1] = {1.68658058, 1.68658058, 1.68658058};

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({2});
  elec_.R[0] = {0.0, 0.0, 0.0};
  elec_.R[1] = {0.0, 1.0, 0.0};

  SpeciesSet& tspecies       = elec_.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx) = -1;

  // Load diamondC_1x1x1 wfn and explicitly construct a SplineC2C object with 7 orbitals
  // This results in padding of the spline coefs table and thus is a more stringent test.
  const char* particles = R"(<tmp>
<determinantset type="einspline" href="diamondC_1x1x1.pwscf.h5" tilematrix="1 0 0 0 1 0 0 0 1" twistnum="0" source="ion" meshfactor="1.0" precision="double" size="7"/>
</tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr ein1 = xmlFirstElementChild(root);
  EinsplineSetBuilder einSet(elec_, ptcl.getPool(), c, ein1);
  auto spo = einSet.createSPOSetFromXML(ein1);
  REQUIRE(spo);

  /*
      Here we test the low-level spline implementation of applyRotation().
      Manually encode the unitary transformation. Ugly, but it works.
    */
  const auto orbitalsetsize = spo->getOrbitalSetSize();
  REQUIRE(orbitalsetsize == 7);
  SPOSet::ValueMatrix psiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_bare(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_bare(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_bare, dpsiM_bare, d2psiM_bare);

  // Check before rotation is applied. From 'test_einset.cpp'
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

  /* Apply a rotation, U = exp(-K) as shown below, to the splines.
     Because K = -K.T, U is unitary.
     K=
      0.00000000e+00  -2.33449314e-01   5.21735139e-01   1.67744276e-01  -1.68565994e-01   1.26632312e-02   2.29272040e-01
      2.33449314e-01   0.00000000e+00  -3.81982621e-01  -6.76995896e-01   7.15727948e-02   9.10926961e-01  -6.94864205e-01
     -5.21735139e-01   3.81982621e-01   0.00000000e+00   2.37433139e-01  -7.09878105e-01  -6.33841289e-01  -7.49009582e-01
     -1.67744276e-01   6.76995896e-01  -2.37433139e-01   0.00000000e+00  -6.72002172e-01   4.27152921e-01  -4.82983743e-03
      1.68565994e-01  -7.15727948e-02   7.09878105e-01   6.72002172e-01   0.00000000e+00  -7.93860012e-01  -7.76624731e-01
     -1.26632312e-02  -9.10926961e-01   6.33841289e-01  -4.27152921e-01   7.93860012e-01   0.00000000e+00  -2.74931052e-01
     -2.29272040e-01   6.94864205e-01   7.49009582e-01   4.82983743e-03   7.76624731e-01   2.74931052e-01   0.00000000e+00
     U=
      8.06061880e-01   3.54921598e-01  -3.40706426e-01  -5.90163619e-02   1.11650454e-01  -1.99768450e-01  -2.28818375e-01
      1.58069821e-02   2.10363421e-01   2.74922448e-01   2.12581764e-01   2.64602356e-01  -6.34971914e-01   6.01265560e-01
      4.43696646e-01  -6.06912539e-02   2.61413193e-01  -1.98368802e-01  -6.43234645e-02   5.75880430e-01   5.96646495e-01
      1.45865363e-01  -5.16577220e-01  -2.09866367e-01   5.55699395e-01   5.73062990e-01   1.74778224e-01   8.77170506e-03
     -3.24609748e-01   2.89664179e-01  -7.50613752e-01  -1.63060005e-01   1.70803377e-01   1.63784167e-01   4.05850414e-01
      1.32570771e-03   3.80299846e-01  -5.08439810e-02   7.59141791e-01  -4.77844928e-01   2.12149087e-01   5.60882349e-02
      1.62781208e-01  -5.75073150e-01  -3.60485665e-01  -1.70070331e-02  -5.72251258e-01  -3.50549638e-01   2.49394158e-01
     UU.T=
      1.00000000e+00  -2.77555756e-17   5.55111512e-17  -8.67361738e-18  -1.11022302e-16  -6.07153217e-17   6.93889390e-17
     -2.77555756e-17   1.00000000e+00   5.55111512e-17  -1.84748050e-16   5.55111512e-17  -6.93889390e-18   8.32667268e-17
      5.55111512e-17   5.55111512e-17   1.00000000e+00   1.45716772e-16  -5.55111512e-17   1.38777878e-16   5.55111512e-17
     -8.67361738e-18  -1.84748050e-16   1.45716772e-16   1.00000000e+00  -1.95156391e-17   1.80411242e-16   8.02309608e-17
     -1.11022302e-16   5.55111512e-17  -5.55111512e-17  -1.95156391e-17   1.00000000e+00   7.28583860e-17   1.11022302e-16
     -6.07153217e-17  -6.93889390e-18   1.38777878e-16   1.80411242e-16   7.28583860e-17   1.00000000e+00  -1.12757026e-16
      6.93889390e-17   8.32667268e-17   5.55111512e-17   8.02309608e-17   1.11022302e-16  -1.12757026e-16   1.00000000e+00

      NB: There's nothing special about this rotation. I purposefully made something 'ugly' to make a worst case scenario...
    */
  SPOSet::ValueMatrix rot_mat(orbitalsetsize, orbitalsetsize);
  rot_mat[0][0] =  8.06061880e-01;
  rot_mat[0][1] =  3.54921598e-01;
  rot_mat[0][2] = -3.40706426e-01;
  rot_mat[0][3] = -5.90163619e-02;
  rot_mat[0][4] =  1.11650454e-01;
  rot_mat[0][5] = -1.99768450e-01;
  rot_mat[0][6] = -2.28818375e-01;
  rot_mat[1][0] =  1.58069821e-02;
  rot_mat[1][1] =  2.10363421e-01;
  rot_mat[1][2] =  2.74922448e-01;
  rot_mat[1][3] =  2.12581764e-01;
  rot_mat[1][4] =  2.64602356e-01;
  rot_mat[1][5] = -6.34971914e-01;
  rot_mat[1][6] =  6.01265560e-01;
  rot_mat[2][0] =  4.43696646e-01;
  rot_mat[2][1] = -6.06912539e-02;
  rot_mat[2][2] =  2.61413193e-01;
  rot_mat[2][3] = -1.98368802e-01;
  rot_mat[2][4] = -6.43234645e-02;
  rot_mat[2][5] =  5.75880430e-01;
  rot_mat[2][6] =  5.96646495e-01;
  rot_mat[3][0] =  1.45865363e-01;
  rot_mat[3][1] = -5.16577220e-01;
  rot_mat[3][2] = -2.09866367e-01;
  rot_mat[3][3] =  5.55699395e-01;
  rot_mat[3][4] =  5.73062990e-01;
  rot_mat[3][5] =  1.74778224e-01;
  rot_mat[3][6] =  8.77170506e-03;
  rot_mat[4][0] = -3.24609748e-01;
  rot_mat[4][1] =  2.89664179e-01;
  rot_mat[4][2] = -7.50613752e-01;
  rot_mat[4][3] = -1.63060005e-01;
  rot_mat[4][4] =  1.70803377e-01;
  rot_mat[4][5] =  1.63784167e-01;
  rot_mat[4][6] =  4.05850414e-01;
  rot_mat[5][0] =  1.32570771e-03;
  rot_mat[5][1] =  3.80299846e-01;
  rot_mat[5][2] = -5.08439810e-02;
  rot_mat[5][3] =  7.59141791e-01;
  rot_mat[5][4] = -4.77844928e-01;
  rot_mat[5][5] =  2.12149087e-01;
  rot_mat[5][6] =  5.60882349e-02;
  rot_mat[6][0] =  1.62781208e-01;
  rot_mat[6][1] = -5.75073150e-01;
  rot_mat[6][2] = -3.60485665e-01;
  rot_mat[6][3] = -1.70070331e-02;
  rot_mat[6][4] = -5.72251258e-01;
  rot_mat[6][5] = -3.50549638e-01;
  rot_mat[6][6] =  2.49394158e-01;
  spo->storeParamsBeforeRotation(); // avoids coefs_copy_ nullptr segfault
  spo->applyRotation(rot_mat, false);

  // Get the data for rotated orbitals
  SPOSet::ValueMatrix psiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::GradMatrix dpsiM_rot(elec_.R.size(), orbitalsetsize);
  SPOSet::ValueMatrix d2psiM_rot(elec_.R.size(), orbitalsetsize);
  spo->evaluate_notranspose(elec_, 0, elec_.R.size(), psiM_rot, dpsiM_rot, d2psiM_rot);

  // Compute the expected value, gradient, and laplacian by hand using the transformation above
  // Check value
  SPOSet::ValueMatrix psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (auto i = 0; i < elec_.R.size(); i++)
    {
      for (auto j = 0; j < orbitalsetsize; j++)
        {
          psiM_rot_manual[i][j] = 0.;
          for (auto k = 0; k < orbitalsetsize; k++)
            {
              psiM_rot_manual[i][j] += psiM_bare[i][k] * rot_mat[k][j];
            }
        }
    }
  checkMatrix(psiM_rot_manual, psiM_rot);

  // Check grad
  SPOSet::GradMatrix dpsiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (auto i = 0; i < elec_.R.size(); i++)
  {
    for (auto j = 0; j < orbitalsetsize; j++)
      {
        dpsiM_rot_manual[i][j][0] = 0.;
        dpsiM_rot_manual[i][j][1] = 0.;
        dpsiM_rot_manual[i][j][2] = 0.;
        for (auto k = 0; k < orbitalsetsize; k++)
          {
            for (auto l = 0; l < 3; l++)
              {
                dpsiM_rot_manual[i][j][l] += dpsiM_bare[i][k][l] * rot_mat[k][j];
              }
          }
      }
  }

  // No checkMatrix for tensors? Gotta do it by hand...
  double res = 0.;
  for (auto i=0; i<elec_.R.size(); i++ )
    {
      for (auto j=0; j<orbitalsetsize; j++ )
        {
          for (auto k=0; k<3; k++ )
            {
              res += std::sqrt(std::norm(dpsiM_rot_manual[i][j][k] - dpsiM_rot[i][j][k]));
            }
        }
    }
  CHECK(res == Approx(0.));

  // Check laplacian
  SPOSet::ValueMatrix d2psiM_rot_manual(elec_.R.size(), orbitalsetsize);
  for (auto i = 0; i < elec_.R.size(); i++)
    {
      for (auto j = 0; j < orbitalsetsize; j++)
        {
          d2psiM_rot_manual[i][j] = 0.;
          for (auto k = 0; k < orbitalsetsize; k++)
            {
              d2psiM_rot_manual[i][j] += d2psiM_bare[i][k] * rot_mat[k][j];
            }
        }
    }
  checkMatrix(d2psiM_rot_manual, d2psiM_rot);

} // TEST_CASE


} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Shiv Upadhyay, shivnupadhyay@gmail.com, University of Washington
//
// File created by: Shiv Upadhyay, shivnupadhyay@gmail.com, University of Washington
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>

#include "Numerics/SoaSphericalTensor.h"
#include "Utilities/for_testing/Catch2Approx.h"

#include <array>
#include <vector>

namespace qmcplusplus
{

// The hard-coded reference values below are the analytic real solid harmonics
// r^l S_l^m and their Cartesian derivatives, in the addsign=true (Gaussian)
// convention. They can be regenerated for every (l,m) rather than the spot
// checks used here by running gen_gto.py, which fills in
// test_full_soa_spherical_tensor.cpp.in from SymPy.

TEST_CASE("SoA Spherical Tensor evaluateVGL subset", "[numerics]")
{
  SoaSphericalTensor<double> st(5, true);

  const double x = 1.3;
  const double y = 1.2;
  const double z = -0.5;
  st.evaluateVGL(x, y, z);

  const double* Ylm = st[0];
  const double* gr0 = st[1];
  const double* gr1 = st[2];
  const double* gr2 = st[3];
  const double* lap = st[4];

  CHECK(Ylm[0] == Approx(0.282094791773878));
  CHECK(gr0[0] == Approx(0.0));
  CHECK(gr1[0] == Approx(0.0));
  CHECK(gr2[0] == Approx(0.0));
  CHECK(lap[0] == Approx(0.0));

  CHECK(Ylm[1] == Approx(0.586323014283504));
  CHECK(gr0[1] == Approx(0.0));
  CHECK(gr1[1] == Approx(0.488602511902920));
  CHECK(gr2[1] == Approx(0.0));
  CHECK(lap[1] == Approx(0.0));

  CHECK(Ylm[6] == Approx(-0.829479816614128));
  CHECK(gr0[6] == Approx(-0.820018069656552));
  CHECK(gr1[6] == Approx(-0.756939756606048));
  CHECK(gr2[6] == Approx(-0.630783130505040));
  CHECK(lap[6] == Approx(0.0));

  CHECK(Ylm[10] == Approx(-2.25467692525963));
  CHECK(gr0[10] == Approx(-1.73436686558433));
  CHECK(gr1[10] == Approx(-1.87889743771636));
  CHECK(gr2[10] == Approx(4.50935385051926));
  CHECK(lap[10] == Approx(0.0));

  CHECK(Ylm[23] == Approx(3.02603855093879));
  CHECK(gr0[23] == Approx(-0.663799038667474));
  CHECK(gr1[23] == Approx(8.28421200257007));
  CHECK(gr2[23] == Approx(-6.05207710187758));
  CHECK(lap[23] == Approx(0.0));

  CHECK(Ylm[26] == Approx(-1.61901660560721));
  CHECK(gr0[26] == Approx(-18.0831700872436));
  CHECK(gr1[26] == Approx(14.1933789091566));
  CHECK(gr2[26] == Approx(3.23803321121443));
  CHECK(lap[26] == Approx(0.0));
}

TEST_CASE("SoA Spherical Tensor evaluateVGH subset", "[numerics]")
{
  SoaSphericalTensor<double> st(5, true);

  const double x = 1.3;
  const double y = 1.2;
  const double z = -0.5;
  st.evaluateVGH(x, y, z);

  const double* Ylm = st[0];
  const double* gr0 = st[1];
  const double* gr1 = st[2];
  const double* gr2 = st[3];
  const double* h00 = st[4];
  const double* h01 = st[5];
  const double* h02 = st[6];
  const double* h11 = st[7];
  const double* h12 = st[8];
  const double* h22 = st[9];

  CHECK(Ylm[0] == Approx(0.282094791773878));
  CHECK(gr0[0] == Approx(0.0));
  CHECK(gr1[0] == Approx(0.0));
  CHECK(gr2[0] == Approx(0.0));
  CHECK(h00[0] == Approx(0.0));
  CHECK(h01[0] == Approx(0.0));
  CHECK(h02[0] == Approx(0.0));
  CHECK(h11[0] == Approx(0.0));
  CHECK(h12[0] == Approx(0.0));
  CHECK(h22[0] == Approx(0.0));

  CHECK(Ylm[3] == Approx(0.635183265473796));
  CHECK(gr0[3] == Approx(0.488602511902920));
  CHECK(gr1[3] == Approx(0.0));
  CHECK(gr2[3] == Approx(0.0));
  CHECK(h00[3] == Approx(0.0));
  CHECK(h01[3] == Approx(0.0));
  CHECK(h02[3] == Approx(0.0));
  CHECK(h11[3] == Approx(0.0));
  CHECK(h12[3] == Approx(0.0));
  CHECK(h22[3] == Approx(0.0));

  CHECK(Ylm[6] == Approx(-0.829479816614128));
  CHECK(gr0[6] == Approx(-0.820018069656552));
  CHECK(gr1[6] == Approx(-0.756939756606048));
  CHECK(gr2[6] == Approx(-0.630783130505040));
  CHECK(h00[6] == Approx(-0.630783130505040));
  CHECK(h01[6] == Approx(0.0));
  CHECK(h02[6] == Approx(0.0));
  CHECK(h11[6] == Approx(-0.630783130505040));
  CHECK(h12[6] == Approx(0.0));
  CHECK(h22[6] == Approx(1.26156626101008));

  CHECK(Ylm[13] == Approx(-1.26555981871711));
  CHECK(gr0[13] == Approx(-2.51832235504921));
  CHECK(gr1[13] == Approx(-1.42598289432913));
  CHECK(gr2[13] == Approx(-2.37663815721522));
  CHECK(h00[13] == Approx(-3.56495723582283));
  CHECK(h01[13] == Approx(-1.09690991871472));
  CHECK(h02[13] == Approx(-1.82818319785786));
  CHECK(h11[13] == Approx(-1.18831907860761));
  CHECK(h12[13] == Approx(0.0));
  CHECK(h22[13] == Approx(4.75327631443044));

  CHECK(Ylm[16] == Approx(0.976303747300715));
  CHECK(gr0[16] == Approx(10.9045618544664));
  CHECK(gr1[16] == Approx(-8.55892951800293));
  CHECK(gr2[16] == Approx(0.0));
  CHECK(h00[16] == Approx(23.4312899352172));
  CHECK(h01[16] == Approx(1.87750720634753));
  CHECK(h02[16] == Approx(0.0));
  CHECK(h11[16] == Approx(-23.4312899352172));
  CHECK(h12[16] == Approx(0.0));
  CHECK(h22[16] == Approx(0.0));

  CHECK(Ylm[35] == Approx(-9.48174731062297));
  CHECK(gr0[35] == Approx(-31.7423080777622));
  CHECK(gr1[35] == Approx(-5.11978004335332));
  CHECK(gr2[35] == Approx(0.0));
  CHECK(h00[35] == Approx(-44.8834050467308));
  CHECK(h01[35] == Approx(-57.1840047919156));
  CHECK(h02[35] == Approx(0.0));
  CHECK(h11[35] == Approx(44.8834050467308));
  CHECK(h12[35] == Approx(0.0));
  CHECK(h22[35] == Approx(0.0));
}

namespace
{
using ValueTable    = std::array<std::vector<double>, 4>;
using GradientTable = std::array<std::vector<double>, 3>;

ValueTable evaluateVGL(const int lmax, const bool addsign, const std::array<double, 3>& position)
{
  SoaSphericalTensor<double> st(lmax, addsign);
  st.evaluateVGL(position[0], position[1], position[2]);

  ValueTable values;
  for (size_t component = 0; component < values.size(); ++component)
    values[component].assign(st[component], st[component] + st.size());
  return values;
}

GradientTable evaluateGradients(const int lmax, const bool addsign, const std::array<double, 3>& position)
{
  const ValueTable vgl = evaluateVGL(lmax, addsign, position);
  return {vgl[1], vgl[2], vgl[3]};
}
} // namespace

TEST_CASE("SoA Spherical Tensor evaluateVGH finite difference", "[numerics]")
{
  constexpr int lmax                = 6;
  constexpr double h                = 1e-2;
  const int hessian_component[3][3] = {{4, 5, 6}, {5, 7, 8}, {6, 8, 9}};

  // The last two positions sit on the z axis and in the xy plane, which drive
  // the r2xy < eps2 branch of evaluate_bare and the sin(theta) == 0 case that
  // the generic points never reach.
  const std::array<std::array<double, 3>, 3> positions{{{0.37, -0.29, 0.43}, {0.0, 0.0, 1.7}, {2.1, 0.0, 0.0}}};

  for (const bool addsign : {false, true})
    for (const auto& position : positions)
    {
      CAPTURE(addsign, position[0], position[1], position[2]);
      const ValueTable reference_vgl = evaluateVGL(lmax, addsign, position);

      SoaSphericalTensor<double> st(lmax, addsign);
      st.evaluateVGH(position[0], position[1], position[2]);

      for (size_t component = 0; component < reference_vgl.size(); ++component)
        for (size_t lm = 0; lm < st.size(); ++lm)
        {
          CAPTURE(component, lm);
          CHECK(st[component][lm] == Approx(reference_vgl[component][lm]).margin(1e-13));
        }

      // For l <= 6 the gradients are polynomials of degree at most five, so the
      // six-point derivative stencil is exact apart from floating-point roundoff.
      for (int derivative_direction = 0; derivative_direction < 3; ++derivative_direction)
      {
        std::array<GradientTable, 6> displaced_gradients;
        constexpr std::array<int, 6> displacements{-3, -2, -1, 1, 2, 3};
        for (size_t sample = 0; sample < displacements.size(); ++sample)
        {
          auto displaced_position = position;
          displaced_position[derivative_direction] += displacements[sample] * h;
          displaced_gradients[sample] = evaluateGradients(lmax, addsign, displaced_position);
        }

        for (int gradient_component = 0; gradient_component < 3; ++gradient_component)
          for (size_t lm = 0; lm < st.size(); ++lm)
          {
            CAPTURE(derivative_direction, gradient_component, lm);
            const double finite_difference = (-displaced_gradients[0][gradient_component][lm] +
                                              9.0 * displaced_gradients[1][gradient_component][lm] -
                                              45.0 * displaced_gradients[2][gradient_component][lm] +
                                              45.0 * displaced_gradients[3][gradient_component][lm] -
                                              9.0 * displaced_gradients[4][gradient_component][lm] +
                                              displaced_gradients[5][gradient_component][lm]) /
                (60.0 * h);
            CHECK(st[hessian_component[gradient_component][derivative_direction]][lm] ==
                  Approx(finite_difference).margin(2e-10).epsilon(2e-10));
          }
      }

      // r^l S_l^m is harmonic, so the Hessian is traceless.
      for (size_t lm = 0; lm < st.size(); ++lm)
      {
        CAPTURE(lm);
        CHECK(st[4][lm] + st[7][lm] + st[9][lm] == Approx(0.0).margin(2e-12));
      }
    }
}

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Numerics/Ylm.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("Legendre", "[numerics]")
{
  // l=0 should be 1.0 for all values
  double v = LegendrePll(0, 0.5);
  REQUIRE(v == Approx(1.0));

  // check endpoints
  for (int i = 0; i < 10; i++)
  {
    double vp1 = LegendrePll(0, 2.0);
    REQUIRE(vp1 == Approx(1.0));
    double vm1 = LegendrePll(0, 2.0);
    REQUIRE(std::abs(vm1) == Approx(1.0));
  }
}

TEST_CASE("Spherical Harmonics", "[numerics]")
{
  // l=0 should be 1.0 for all values
  using vec_t = TinyVector<double, 3>;
  vec_t v;
  v[0]                     = 1.0;
  v[1]                     = 0.0;
  v[2]                     = 0.0;
  std::complex<double> out = Ylm(0, 0, v);
  REQUIRE(out.real() == Approx(0.28209479177387814)); // Y00 = 1/2 1/sqrt(pi)
}

TEST_CASE("Spherical Harmonics Many", "[numerics]")
{
  struct Point
  {
    double x;
    double y;
    double z;
  };

  struct YlmValue
  {
    Point p;
    int l;
    int m;
    double y_re;
    double y_im;
  };

#if 0
  // Use gen_ylm.py to create this file to test more values
#include "ylm.inc"
#else
  // Test a small subset of values
  const int N      = 11;
  YlmValue Vals[N] = {
      {{0, 0, 1}, 1, -1, 0, 0},
      {{0.587785, 0, 0.809017}, 1, -1, 0.203076, 0},
      {{0.587785, 0, 0.809017}, 1, 0, 0.395288, 0},
      {{0.951057, 0, 0.309017}, 1, 1, -0.328584, 0},
      {{0.587785, 0, 0.809017}, 2, -2, 0.133454, 0},
      {{0.293893, 0.904508, 0.309017}, 2, -1, 0.0701612, -0.215934},
      {{0.587785, 0, -0.809017}, 2, 0, 0.303888, 0},
      {{0.293893, 0.904508, -0.309017}, 2, 1, 0.0701612, 0.215934},
      {{0.293893, 0.904508, 0.309017}, 2, 2, -0.282661, 0.205365},
      {{-0.475528, -0.345492, -0.809017}, 3, -3, 0.0261823, 0.0805808},
      {{-0.475528, -0.345492, 0.809017}, 4, 4, -0.0427344, 0.0310484},
  };
#endif

  using vec_t = TinyVector<double, 3>;
  for (int i = 0; i < N; i++)
  {
    YlmValue& v = Vals[i];
    vec_t w;
    w[0] = v.p.z; // first component appears to be aligned along the z-axis
    w[1] = v.p.x;
    w[2] = v.p.y;

    std::complex<double> out = Ylm(v.l, v.m, w);
    //printf("%d %d  expected %g %g   actual %g %g\n",v.l,v.m,v.y_re,v.y_im, out.real(), out.imag());
    REQUIRE(v.y_re == Approx(out.real()));
    REQUIRE(v.y_im == Approx(out.imag()));
  }
}

TEST_CASE("Derivatives of Spherical Harmonics", "[numerics]")
{
  struct Point
  {
    double x;
    double y;
    double z;
  };

  struct YlmDerivValue
  {
    Point p;
    int l;
    int m;
    double th_re;
    double th_im;
    double ph_re;
    double ph_im;
  };

  //reference values from sympy
  //d/dtheta Ylm and d/dphi Ylm
  const int N           = 10;
  YlmDerivValue Vals[N] = {
      {{0.587785, 0, 0.809017}, 1, -1, 0.27951007, 0.0, 0.0, -0.2030763},
      {{0.587785, 0, 0.809017}, 1, 0, -0.2871933, 0.0, 0.0, 0.0},
      {{0.951057, 0, 0.309017}, 1, 1, -0.1067635, 0.0, 0.0, -0.3285845},
      {{0.587785, 0, 0.809017}, 2, -2, 0.3673685, 0.0, 0.0, -0.2669088},
      {{0.293893, 0.904508, 0.309017}, 2, -1, -0.1931373, 0.5944147, -0.2159338, -0.0701613},
      {{0.587785, 0, -0.809017}, 2, 0, 0.8998655, 0.0, 0.0, 0.0},
      {{0.293893, 0.904508, -0.309017}, 2, 1, 0.1931373, 0.5944146, -0.2159339, 0.0701613},
      {{0.293893, 0.904508, 0.309017}, 2, 2, -0.1836842, 0.1334547, -0.4107311, -0.5653217},
      {{-0.475528, -0.345492, -0.809017}, 3, -3, -0.1081114, -0.3327295, 0.2417422, -0.0785476},
      {{-0.475528, -0.345492, 0.809017}, 4, 4, -0.2352762, 0.1709368, -0.1241928, -0.1709382},
  };

  using vec_t = TinyVector<double, 3>;
  for (int i = 0; i < N; i++)
  {
    YlmDerivValue& v = Vals[i];
    vec_t w;
    w[0] = v.p.z; // first component appears to be aligned along the z-axis
    w[1] = v.p.x;
    w[2] = v.p.y;

    std::complex<double> theta, phi;
    derivYlmSpherical(v.l, v.m, w, theta, phi, false);
    //printf("%d %d  expected %g %g %g %g  actual %g %g %g %g\n",v.l,v.m,v.th_re,v.th_im, v.ph_re, v.ph_im, theta.real(), theta.imag(), phi.real(), phi.imag());
    REQUIRE(std::real(theta) == Approx(v.th_re));
    REQUIRE(std::imag(theta) == Approx(v.th_im));
    REQUIRE(std::real(phi) == Approx(v.ph_re));
    REQUIRE(std::imag(phi) == Approx(v.ph_im));
  }
}

TEST_CASE("Spherical Harmonics Wrapper", "[numerics]")
{
  struct Point
  {
    double x;
    double y;
    double z;
  };

  struct YlmValue
  {
    Point p;
    int l;
    int m;
    double y_re;
    double y_im;
  };

  // Test a small subset of values, same test as Spherical Harmonics except some are scaled to not be unit vectors
  const int N      = 11;
  YlmValue Vals[N] = {
      {{0, 0, 5}, 1, -1, 0, 0},
      {{2.938925, 0, 4.045085}, 1, -1, 0.203076, 0},
      {{2.938925, 0, 4.045085}, 1, 0, 0.395288, 0},
      {{0.951057, 0, 0.309017}, 1, 1, -0.328584, 0},
      {{0.587785, 0, 0.809017}, 2, -2, 0.133454, 0},
      {{1.469465, 4.52254, 1.545085}, 2, -1, 0.0701612, -0.215934},
      {{0.587785, 0, -0.809017}, 2, 0, 0.303888, 0},
      {{0.293893, 0.904508, -0.309017}, 2, 1, 0.0701612, 0.215934},
      {{1.469465, 4.52254, 1.545085}, 2, 2, -0.282661, 0.205365},
      {{-0.475528, -0.345492, -0.809017}, 3, -3, 0.0261823, 0.0805808},
      {{-2.37764, -1.72746, 4.045085}, 4, 4, -0.0427344, 0.0310484},
  };

  using vec_t = TinyVector<double, 3>;
  for (int i = 0; i < N; i++)
  {
    YlmValue& v = Vals[i];
    vec_t w;
    w[0] = v.p.x;
    w[1] = v.p.y;
    w[2] = v.p.z;

    std::complex<double> out = sphericalHarmonic(v.l, v.m, w);
    //printf("%d %d  expected %g %g   actual %g %g\n",v.l,v.m,v.y_re,v.y_im, out.real(), out.imag());
    REQUIRE(v.y_re == Approx(out.real()));
    REQUIRE(v.y_im == Approx(out.imag()));
  }
}

} // namespace qmcplusplus

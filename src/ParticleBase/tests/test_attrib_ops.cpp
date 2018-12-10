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

#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "ParticleBase/ParticleAttrib.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/ParticleUtility.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

template<unsigned int D> void double_test_case()
{
  TinyVector<double, D> v1;
  v1 = 2.0;
  double val1 = dot(v1, v1);
  REQUIRE(val1 == Approx(4.0*D));


  ParticleAttrib<TinyVector<double, D> > PA1;
  PA1.resize(3);
  REQUIRE(PA1.size() == 3);

  // Whole array operation
  PA1 = 1.0;

  double val = Dot(PA1, PA1);
  REQUIRE(val == Approx(3*1.0*D));
}


TEST_CASE("particle_attrib_ops_double", "[particle_base]")
{
  OHMMS::Controller->initialize(0, NULL);

  SECTION("dim = 1") {
    double_test_case<1>();
  }
  SECTION("dim = 2") {
    double_test_case<2>();
  }
  SECTION("dim = 3") {
    double_test_case<3>();
  }
  SECTION("dim = 4") {
    double_test_case<4>();
  }
}

template <unsigned int D> void complex_test_case()
{
  TinyVector<std::complex<double>, D> v1;

  v1 = std::complex<double> (2.0, 1.0);
  double val1 = OTCDot<double, double, D>::apply(v1, v1);
  REQUIRE(val1 == Approx(3.0*D));

  double val1_cc = OTCDot_CC<double, double, D>::apply(v1, v1);
  REQUIRE(val1_cc == Approx(5.0*D));


  ParticleAttrib<TinyVector<std::complex<double>, D> > PA1;
  PA1.resize(3);
  REQUIRE(PA1.size() == 3);

  // Whole array operation
  PA1 = std::complex<double> (1.0, 2.0);

  double val = Dot(PA1, PA1);
  REQUIRE(val == Approx(-3.0*3*D));

  double val_cc = Dot_CC(PA1, PA1);
  REQUIRE(val_cc == Approx(5.0*3*D));
}

TEST_CASE("particle_attrib_ops_complex", "[particle_base]")
{
  OHMMS::Controller->initialize(0, NULL);

  SECTION("dim = 1") {
    complex_test_case<1>();
  }
  SECTION("dim = 2") {
    complex_test_case<2>();
  }
  SECTION("dim = 3") {
    complex_test_case<3>();
  }
  SECTION("dim = 4") {
    complex_test_case<4>();
  }
}

}

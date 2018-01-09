//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "Numerics/CartesianTensor.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

typedef OHMMS_PRECISION real_type;


// Use gen_gto.py to generate the checks

TEST_CASE("Cartesian Tensor", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  //TinyVector<double, 3> pt(1.3, 0.0, 0.0);

  ct.evaluate(pt);

  //for (int i = 0; i < 35; i++) {
  //std::cout << "XYZ = " << i << " " << ct.getYlm(i) << std::endl;
  //}
  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getYlm(2) == Approx(0.586323014284));
  REQUIRE(ct.getYlm(3) == Approx(-0.244301255951));
  REQUIRE(ct.getYlm(4) == Approx(1.06602349055));
  REQUIRE(ct.getYlm(5) == Approx(0.908327707927));
  REQUIRE(ct.getYlm(6) == Approx(0.157695782626));
  REQUIRE(ct.getYlm(7) == Approx(1.70437555172));
  REQUIRE(ct.getYlm(8) == Approx(-0.710156479885));
  REQUIRE(ct.getYlm(9) == Approx(-0.655529058355));
  REQUIRE(ct.getYlm(10) == Approx(1.6397368054));
  REQUIRE(ct.getYlm(11) == Approx(1.28969740543));
  REQUIRE(ct.getYlm(12) == Approx(-0.0932940831475));
  REQUIRE(ct.getYlm(13) == Approx(3.38451965731));
  REQUIRE(ct.getYlm(14) == Approx(-1.41021652388));
  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getYlm(16) == Approx(-1.20160461206));
  REQUIRE(ct.getYlm(17) == Approx(0.542390970723));
  REQUIRE(ct.getYlm(18) == Approx(0.500668588359));
  REQUIRE(ct.getYlm(19) == Approx(-2.25467692526));
  REQUIRE(ct.getYlm(20) == Approx(2.41707280436));
  REQUIRE(ct.getYlm(21) == Approx(1.75485528067));
  REQUIRE(ct.getYlm(22) == Approx(0.0528927734576));
  REQUIRE(ct.getYlm(23) == Approx(5.90305249944));
  REQUIRE(ct.getYlm(24) == Approx(-2.4596052081));
  REQUIRE(ct.getYlm(25) == Approx(5.02981988118));
  REQUIRE(ct.getYlm(26) == Approx(-1.93454610815));
  REQUIRE(ct.getYlm(27) == Approx(-0.363846924275));
  REQUIRE(ct.getYlm(28) == Approx(-0.335858699331));
  REQUIRE(ct.getYlm(29) == Approx(7.03459200681));
  REQUIRE(ct.getYlm(30) == Approx(1.22128333452));
  REQUIRE(ct.getYlm(31) == Approx(1.04062011935));
  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getYlm(33) == Approx(-4.68625798704));
  REQUIRE(ct.getYlm(34) == Approx(1.9526074946));
}

TEST_CASE("Cartesian Tensor evaluateAll subset", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateAll(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGradYlm(i) << "  " << ct.getLaplYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(0) == Approx(0));

  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(1) == Approx(0));

  REQUIRE(ct.getYlm(8) == Approx(-0.710156479885));
  REQUIRE(ct.getGradYlm(8)[0] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(8)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(8)[2] == Approx(1.42031295977));
  REQUIRE(ct.getLaplYlm(8) == Approx(0));

  REQUIRE(ct.getYlm(23) == Approx(5.90305249944));
  REQUIRE(ct.getGradYlm(23)[0] == Approx(13.6224288449));
  REQUIRE(ct.getGradYlm(23)[1] == Approx(4.9192104162));
  REQUIRE(ct.getGradYlm(23)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(23) == Approx(20.9575828383));

  REQUIRE(ct.getYlm(34) == Approx(1.9526074946));
  REQUIRE(ct.getGradYlm(34)[0] == Approx(1.50200576508));
  REQUIRE(ct.getGradYlm(34)[1] == Approx(1.62717291217));
  REQUIRE(ct.getGradYlm(34)[2] == Approx(-7.81042997841));
  REQUIRE(ct.getLaplYlm(34) == Approx(15.6208599568));
}

TEST_CASE("Cartesian Tensor evaluateWithHessian subset", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithHessian(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getHessYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(0)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(1)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  REQUIRE(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(15)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(15)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(0,1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(1,0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(1,1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(15)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  REQUIRE(ct.getGradYlm(32)[2] == Approx(10.1535589719));

  REQUIRE(ct.getHessYlm(32)(0,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(32)(0,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(0,2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(1,0) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(32)(1,2) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,0) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(2,1) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,2) == Approx(0));

}

TEST_CASE("Cartesian Tensor evaluateWithThirdDeriv subset", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithThirdDeriv(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGGGYlm(i) <<  std::endl;
  //}


  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));


  REQUIRE(ct.getYlm(27) == Approx(-0.363846924275));
  REQUIRE(ct.getGradYlm(27)[0] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(27)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(27)[2] == Approx(2.18308154565));

  REQUIRE(ct.getHessYlm(27)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,0) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,2) == Approx(-8.73232618261));


  REQUIRE(ct.getGGGYlm(27)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,2) == Approx(17.4646523652));
}


TEST_CASE("Cartesian Tensor evaluateThirdDerivOnly subset", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateThirdDerivOnly(pt);

  REQUIRE(ct.getGGGYlm(14)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(33)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[0](1,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](2,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](1,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](1,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,2) == Approx(0));

}

// Activate for a more complete test of CartesianTensor
#if 0

TEST_CASE("Cartesian Tensor evaluateAll", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateAll(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGradYlm(i) << "  " << ct.getLaplYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(0) == Approx(0));
  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(1) == Approx(0));
  REQUIRE(ct.getYlm(2) == Approx(0.586323014284));
  REQUIRE(ct.getGradYlm(2)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(2)[1] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(2)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(2) == Approx(0));
  REQUIRE(ct.getYlm(3) == Approx(-0.244301255951));
  REQUIRE(ct.getGradYlm(3)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[2] == Approx(0.488602511903));
  REQUIRE(ct.getLaplYlm(3) == Approx(0));
  REQUIRE(ct.getYlm(4) == Approx(1.06602349055));
  REQUIRE(ct.getGradYlm(4)[0] == Approx(1.64003613931));
  REQUIRE(ct.getGradYlm(4)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(4)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(4) == Approx(1.26156626101));
  REQUIRE(ct.getYlm(5) == Approx(0.908327707927));
  REQUIRE(ct.getGradYlm(5)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(5)[1] == Approx(1.51387951321));
  REQUIRE(ct.getGradYlm(5)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(5) == Approx(1.26156626101));
  REQUIRE(ct.getYlm(6) == Approx(0.157695782626));
  REQUIRE(ct.getGradYlm(6)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[2] == Approx(-0.630783130505));
  REQUIRE(ct.getLaplYlm(6) == Approx(1.26156626101));
  REQUIRE(ct.getYlm(7) == Approx(1.70437555172));
  REQUIRE(ct.getGradYlm(7)[0] == Approx(1.31105811671));
  REQUIRE(ct.getGradYlm(7)[1] == Approx(1.42031295977));
  REQUIRE(ct.getGradYlm(7)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(7) == Approx(0));
  REQUIRE(ct.getYlm(8) == Approx(-0.710156479885));
  REQUIRE(ct.getGradYlm(8)[0] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(8)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(8)[2] == Approx(1.42031295977));
  REQUIRE(ct.getLaplYlm(8) == Approx(0));
  REQUIRE(ct.getYlm(9) == Approx(-0.655529058355));
  REQUIRE(ct.getGradYlm(9)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(9)[1] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(9)[2] == Approx(1.31105811671));
  REQUIRE(ct.getLaplYlm(9) == Approx(0));
  REQUIRE(ct.getYlm(10) == Approx(1.6397368054));
  REQUIRE(ct.getGradYlm(10)[0] == Approx(3.78400801246));
  REQUIRE(ct.getGradYlm(10)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(10)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(10) == Approx(5.82155078841));
  REQUIRE(ct.getYlm(11) == Approx(1.28969740543));
  REQUIRE(ct.getGradYlm(11)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(11)[1] == Approx(3.22424351358));
  REQUIRE(ct.getGradYlm(11)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(11) == Approx(5.3737391893));
  REQUIRE(ct.getYlm(12) == Approx(-0.0932940831475));
  REQUIRE(ct.getGradYlm(12)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[2] == Approx(0.559764498885));
  REQUIRE(ct.getLaplYlm(12) == Approx(-2.23905799554));
  REQUIRE(ct.getYlm(13) == Approx(3.38451965731));
  REQUIRE(ct.getGradYlm(13)[0] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(13)[1] == Approx(2.82043304776));
  REQUIRE(ct.getGradYlm(13)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(13) == Approx(4.00534870687));
  REQUIRE(ct.getYlm(14) == Approx(-1.41021652388));
  REQUIRE(ct.getGradYlm(14)[0] == Approx(-2.16956388289));
  REQUIRE(ct.getGradYlm(14)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(14)[2] == Approx(2.82043304776));
  REQUIRE(ct.getLaplYlm(14) == Approx(-1.66889529453));
  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  REQUIRE(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(15)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(15) == Approx(4.33912776578));
  REQUIRE(ct.getYlm(16) == Approx(-1.20160461206));
  REQUIRE(ct.getGradYlm(16)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(16)[1] == Approx(-2.00267435344));
  REQUIRE(ct.getGradYlm(16)[2] == Approx(2.40320922412));
  REQUIRE(ct.getLaplYlm(16) == Approx(-1.66889529453));
  REQUIRE(ct.getYlm(17) == Approx(0.542390970723));
  REQUIRE(ct.getGradYlm(17)[0] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(17)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(17)[2] == Approx(-2.16956388289));
  REQUIRE(ct.getLaplYlm(17) == Approx(4.33912776578));
  REQUIRE(ct.getYlm(18) == Approx(0.500668588359));
  REQUIRE(ct.getGradYlm(18)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(18)[1] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(18)[2] == Approx(-2.00267435344));
  REQUIRE(ct.getLaplYlm(18) == Approx(4.00534870687));
  REQUIRE(ct.getYlm(19) == Approx(-2.25467692526));
  REQUIRE(ct.getGradYlm(19)[0] == Approx(-1.73436686558));
  REQUIRE(ct.getGradYlm(19)[1] == Approx(-1.87889743772));
  REQUIRE(ct.getGradYlm(19)[2] == Approx(4.50935385052));
  REQUIRE(ct.getLaplYlm(19) == Approx(0));
  REQUIRE(ct.getYlm(20) == Approx(2.41707280436));
  REQUIRE(ct.getGradYlm(20)[0] == Approx(7.43714709033));
  REQUIRE(ct.getGradYlm(20)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(20)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(20) == Approx(17.1626471315));
  REQUIRE(ct.getYlm(21) == Approx(1.75485528067));
  REQUIRE(ct.getGradYlm(21)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(21)[1] == Approx(5.84951760222));
  REQUIRE(ct.getGradYlm(21)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(21) == Approx(14.6237940056));
  REQUIRE(ct.getYlm(22) == Approx(0.0528927734576));
  REQUIRE(ct.getGradYlm(22)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[2] == Approx(-0.423142187661));
  REQUIRE(ct.getLaplYlm(22) == Approx(2.53885312596));
  REQUIRE(ct.getYlm(23) == Approx(5.90305249944));
  REQUIRE(ct.getGradYlm(23)[0] == Approx(13.6224288449));
  REQUIRE(ct.getGradYlm(23)[1] == Approx(4.9192104162));
  REQUIRE(ct.getGradYlm(23)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(23) == Approx(20.9575828383));
  REQUIRE(ct.getYlm(24) == Approx(-2.4596052081));
  REQUIRE(ct.getGradYlm(24)[0] == Approx(-5.6760120187));
  REQUIRE(ct.getGradYlm(24)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(24)[2] == Approx(4.9192104162));
  REQUIRE(ct.getLaplYlm(24) == Approx(-8.73232618261));
  REQUIRE(ct.getYlm(25) == Approx(5.02981988118));
  REQUIRE(ct.getGradYlm(25)[0] == Approx(3.86909221629));
  REQUIRE(ct.getGradYlm(25)[1] == Approx(12.574549703));
  REQUIRE(ct.getGradYlm(25)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(25) == Approx(20.9575828383));
  REQUIRE(ct.getYlm(26) == Approx(-1.93454610815));
  REQUIRE(ct.getGradYlm(26)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(26)[1] == Approx(-4.83636527037));
  REQUIRE(ct.getGradYlm(26)[2] == Approx(3.86909221629));
  REQUIRE(ct.getLaplYlm(26) == Approx(-8.06060878395));
  REQUIRE(ct.getYlm(27) == Approx(-0.363846924275));
  REQUIRE(ct.getGradYlm(27)[0] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(27)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(27)[2] == Approx(2.18308154565));
  REQUIRE(ct.getLaplYlm(27) == Approx(-8.73232618261));
  REQUIRE(ct.getYlm(28) == Approx(-0.335858699331));
  REQUIRE(ct.getGradYlm(28)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(28)[1] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(28)[2] == Approx(2.01515219599));
  REQUIRE(ct.getLaplYlm(28) == Approx(-8.06060878395));
  REQUIRE(ct.getYlm(29) == Approx(7.03459200681));
  REQUIRE(ct.getGradYlm(29)[0] == Approx(10.8224492412));
  REQUIRE(ct.getGradYlm(29)[1] == Approx(11.7243200114));
  REQUIRE(ct.getGradYlm(29)[2] == Approx(0));
  REQUIRE(ct.getLaplYlm(29) == Approx(18.0952276309));
  REQUIRE(ct.getYlm(30) == Approx(1.22128333452));
  REQUIRE(ct.getGradYlm(30)[0] == Approx(1.87889743772));
  REQUIRE(ct.getGradYlm(30)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(30)[2] == Approx(-4.88513333806));
  REQUIRE(ct.getLaplYlm(30) == Approx(11.2155723974));
  REQUIRE(ct.getYlm(31) == Approx(1.04062011935));
  REQUIRE(ct.getGradYlm(31)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(31)[1] == Approx(1.73436686558));
  REQUIRE(ct.getGradYlm(31)[2] == Approx(-4.1624804774));
  REQUIRE(ct.getLaplYlm(31) == Approx(9.77026667613));
  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  REQUIRE(ct.getGradYlm(32)[2] == Approx(10.1535589719));
  REQUIRE(ct.getLaplYlm(32) == Approx(-6.00802306031));
  REQUIRE(ct.getYlm(33) == Approx(-4.68625798704));
  REQUIRE(ct.getGradYlm(33)[0] == Approx(-3.60481383619));
  REQUIRE(ct.getGradYlm(33)[1] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(33)[2] == Approx(9.37251597409));
  REQUIRE(ct.getLaplYlm(33) == Approx(-6.50869164867));
  REQUIRE(ct.getYlm(34) == Approx(1.9526074946));
  REQUIRE(ct.getGradYlm(34)[0] == Approx(1.50200576508));
  REQUIRE(ct.getGradYlm(34)[1] == Approx(1.62717291217));
  REQUIRE(ct.getGradYlm(34)[2] == Approx(-7.81042997841));
  REQUIRE(ct.getLaplYlm(34) == Approx(15.6208599568));
}

TEST_CASE("Cartesian Tensor evaluateWithHessian", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithHessian(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getHessYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(0)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(1)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(2) == Approx(0.586323014284));
  REQUIRE(ct.getGradYlm(2)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(2)[1] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(2)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(2)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(3) == Approx(-0.244301255951));
  REQUIRE(ct.getGradYlm(3)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[2] == Approx(0.488602511903));

  REQUIRE(ct.getHessYlm(3)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(4) == Approx(1.06602349055));
  REQUIRE(ct.getGradYlm(4)[0] == Approx(1.64003613931));
  REQUIRE(ct.getGradYlm(4)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(4)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(4)(0,0) == Approx(1.26156626101));
  REQUIRE(ct.getHessYlm(4)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(5) == Approx(0.908327707927));
  REQUIRE(ct.getGradYlm(5)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(5)[1] == Approx(1.51387951321));
  REQUIRE(ct.getGradYlm(5)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(5)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(1,1) == Approx(1.26156626101));
  REQUIRE(ct.getHessYlm(5)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(6) == Approx(0.157695782626));
  REQUIRE(ct.getGradYlm(6)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[2] == Approx(-0.630783130505));

  REQUIRE(ct.getHessYlm(6)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,2) == Approx(1.26156626101));


  REQUIRE(ct.getYlm(7) == Approx(1.70437555172));
  REQUIRE(ct.getGradYlm(7)[0] == Approx(1.31105811671));
  REQUIRE(ct.getGradYlm(7)[1] == Approx(1.42031295977));
  REQUIRE(ct.getGradYlm(7)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(7)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(0,1) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(7)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(1,0) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(7)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(8) == Approx(-0.710156479885));
  REQUIRE(ct.getGradYlm(8)[0] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(8)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(8)[2] == Approx(1.42031295977));

  REQUIRE(ct.getHessYlm(8)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(0,2) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(8)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(2,0) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(8)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(9) == Approx(-0.655529058355));
  REQUIRE(ct.getGradYlm(9)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(9)[1] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(9)[2] == Approx(1.31105811671));

  REQUIRE(ct.getHessYlm(9)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,2) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(9)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(2,1) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(9)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(10) == Approx(1.6397368054));
  REQUIRE(ct.getGradYlm(10)[0] == Approx(3.78400801246));
  REQUIRE(ct.getGradYlm(10)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(10)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(10)(0,0) == Approx(5.82155078841));
  REQUIRE(ct.getHessYlm(10)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(11) == Approx(1.28969740543));
  REQUIRE(ct.getGradYlm(11)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(11)[1] == Approx(3.22424351358));
  REQUIRE(ct.getGradYlm(11)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(11)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(1,1) == Approx(5.3737391893));
  REQUIRE(ct.getHessYlm(11)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(12) == Approx(-0.0932940831475));
  REQUIRE(ct.getGradYlm(12)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[2] == Approx(0.559764498885));

  REQUIRE(ct.getHessYlm(12)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,2) == Approx(-2.23905799554));


  REQUIRE(ct.getYlm(13) == Approx(3.38451965731));
  REQUIRE(ct.getGradYlm(13)[0] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(13)[1] == Approx(2.82043304776));
  REQUIRE(ct.getGradYlm(13)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(13)(0,0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(13)(0,1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(13)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(1,0) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(13)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(14) == Approx(-1.41021652388));
  REQUIRE(ct.getGradYlm(14)[0] == Approx(-2.16956388289));
  REQUIRE(ct.getGradYlm(14)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(14)[2] == Approx(2.82043304776));

  REQUIRE(ct.getHessYlm(14)(0,0) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(14)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(0,2) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(14)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(2,0) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(14)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  REQUIRE(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(15)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(15)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(0,1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(1,0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(1,1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(15)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(16) == Approx(-1.20160461206));
  REQUIRE(ct.getGradYlm(16)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(16)[1] == Approx(-2.00267435344));
  REQUIRE(ct.getGradYlm(16)[2] == Approx(2.40320922412));

  REQUIRE(ct.getHessYlm(16)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(1,1) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(16)(1,2) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(16)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(2,1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(16)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(17) == Approx(0.542390970723));
  REQUIRE(ct.getGradYlm(17)[0] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(17)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(17)[2] == Approx(-2.16956388289));

  REQUIRE(ct.getHessYlm(17)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(0,2) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(17)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(2,0) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(17)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(2,2) == Approx(4.33912776578));


  REQUIRE(ct.getYlm(18) == Approx(0.500668588359));
  REQUIRE(ct.getGradYlm(18)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(18)[1] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(18)[2] == Approx(-2.00267435344));

  REQUIRE(ct.getHessYlm(18)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,2) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(18)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(2,1) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(18)(2,2) == Approx(4.00534870687));


  REQUIRE(ct.getYlm(19) == Approx(-2.25467692526));
  REQUIRE(ct.getGradYlm(19)[0] == Approx(-1.73436686558));
  REQUIRE(ct.getGradYlm(19)[1] == Approx(-1.87889743772));
  REQUIRE(ct.getGradYlm(19)[2] == Approx(4.50935385052));

  REQUIRE(ct.getHessYlm(19)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(19)(0,1) == Approx(-1.44530572132));
  REQUIRE(ct.getHessYlm(19)(0,2) == Approx(3.46873373117));
  REQUIRE(ct.getHessYlm(19)(1,0) == Approx(-1.44530572132));
  REQUIRE(ct.getHessYlm(19)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(19)(1,2) == Approx(3.75779487543));
  REQUIRE(ct.getHessYlm(19)(2,0) == Approx(3.46873373117));
  REQUIRE(ct.getHessYlm(19)(2,1) == Approx(3.75779487543));
  REQUIRE(ct.getHessYlm(19)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(20) == Approx(2.41707280436));
  REQUIRE(ct.getGradYlm(20)[0] == Approx(7.43714709033));
  REQUIRE(ct.getGradYlm(20)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(20)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(20)(0,0) == Approx(17.1626471315));
  REQUIRE(ct.getHessYlm(20)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(21) == Approx(1.75485528067));
  REQUIRE(ct.getGradYlm(21)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(21)[1] == Approx(5.84951760222));
  REQUIRE(ct.getGradYlm(21)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(21)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(1,1) == Approx(14.6237940056));
  REQUIRE(ct.getHessYlm(21)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(22) == Approx(0.0528927734576));
  REQUIRE(ct.getGradYlm(22)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[2] == Approx(-0.423142187661));

  REQUIRE(ct.getHessYlm(22)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,2) == Approx(2.53885312596));


  REQUIRE(ct.getYlm(23) == Approx(5.90305249944));
  REQUIRE(ct.getGradYlm(23)[0] == Approx(13.6224288449));
  REQUIRE(ct.getGradYlm(23)[1] == Approx(4.9192104162));
  REQUIRE(ct.getGradYlm(23)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(23)(0,0) == Approx(20.9575828383));
  REQUIRE(ct.getHessYlm(23)(0,1) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(23)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(1,0) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(23)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(24) == Approx(-2.4596052081));
  REQUIRE(ct.getGradYlm(24)[0] == Approx(-5.6760120187));
  REQUIRE(ct.getGradYlm(24)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(24)[2] == Approx(4.9192104162));

  REQUIRE(ct.getHessYlm(24)(0,0) == Approx(-8.73232618261));
  REQUIRE(ct.getHessYlm(24)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(0,2) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(24)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(2,0) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(24)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(25) == Approx(5.02981988118));
  REQUIRE(ct.getGradYlm(25)[0] == Approx(3.86909221629));
  REQUIRE(ct.getGradYlm(25)[1] == Approx(12.574549703));
  REQUIRE(ct.getGradYlm(25)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(25)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(0,1) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(25)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(1,0) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(25)(1,1) == Approx(20.9575828383));
  REQUIRE(ct.getHessYlm(25)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(26) == Approx(-1.93454610815));
  REQUIRE(ct.getGradYlm(26)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(26)[1] == Approx(-4.83636527037));
  REQUIRE(ct.getGradYlm(26)[2] == Approx(3.86909221629));

  REQUIRE(ct.getHessYlm(26)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(1,1) == Approx(-8.06060878395));
  REQUIRE(ct.getHessYlm(26)(1,2) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(26)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(2,1) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(26)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(27) == Approx(-0.363846924275));
  REQUIRE(ct.getGradYlm(27)[0] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(27)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(27)[2] == Approx(2.18308154565));

  REQUIRE(ct.getHessYlm(27)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,0) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,2) == Approx(-8.73232618261));


  REQUIRE(ct.getYlm(28) == Approx(-0.335858699331));
  REQUIRE(ct.getGradYlm(28)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(28)[1] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(28)[2] == Approx(2.01515219599));

  REQUIRE(ct.getHessYlm(28)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(28)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(2,1) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(28)(2,2) == Approx(-8.06060878395));


  REQUIRE(ct.getYlm(29) == Approx(7.03459200681));
  REQUIRE(ct.getGradYlm(29)[0] == Approx(10.8224492412));
  REQUIRE(ct.getGradYlm(29)[1] == Approx(11.7243200114));
  REQUIRE(ct.getGradYlm(29)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(29)(0,0) == Approx(8.3249609548));
  REQUIRE(ct.getHessYlm(29)(0,1) == Approx(18.0374154021));
  REQUIRE(ct.getHessYlm(29)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(1,0) == Approx(18.0374154021));
  REQUIRE(ct.getHessYlm(29)(1,1) == Approx(9.77026667613));
  REQUIRE(ct.getHessYlm(29)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(30) == Approx(1.22128333452));
  REQUIRE(ct.getGradYlm(30)[0] == Approx(1.87889743772));
  REQUIRE(ct.getGradYlm(30)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(30)[2] == Approx(-4.88513333806));

  REQUIRE(ct.getHessYlm(30)(0,0) == Approx(1.44530572132));
  REQUIRE(ct.getHessYlm(30)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(0,2) == Approx(-7.51558975087));
  REQUIRE(ct.getHessYlm(30)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(2,0) == Approx(-7.51558975087));
  REQUIRE(ct.getHessYlm(30)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(2,2) == Approx(9.77026667613));


  REQUIRE(ct.getYlm(31) == Approx(1.04062011935));
  REQUIRE(ct.getGradYlm(31)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(31)[1] == Approx(1.73436686558));
  REQUIRE(ct.getGradYlm(31)[2] == Approx(-4.1624804774));

  REQUIRE(ct.getHessYlm(31)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(1,1) == Approx(1.44530572132));
  REQUIRE(ct.getHessYlm(31)(1,2) == Approx(-6.93746746234));
  REQUIRE(ct.getHessYlm(31)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(2,1) == Approx(-6.93746746234));
  REQUIRE(ct.getHessYlm(31)(2,2) == Approx(8.3249609548));


  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  REQUIRE(ct.getGradYlm(32)[2] == Approx(10.1535589719));

  REQUIRE(ct.getHessYlm(32)(0,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(32)(0,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(0,2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(1,0) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(32)(1,2) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,0) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(2,1) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(33) == Approx(-4.68625798704));
  REQUIRE(ct.getGradYlm(33)[0] == Approx(-3.60481383619));
  REQUIRE(ct.getGradYlm(33)[1] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(33)[2] == Approx(9.37251597409));

  REQUIRE(ct.getHessYlm(33)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(33)(0,1) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(33)(0,2) == Approx(7.20962767237));
  REQUIRE(ct.getHessYlm(33)(1,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(33)(1,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(33)(1,2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(33)(2,0) == Approx(7.20962767237));
  REQUIRE(ct.getHessYlm(33)(2,1) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(33)(2,2) == Approx(0));


  REQUIRE(ct.getYlm(34) == Approx(1.9526074946));
  REQUIRE(ct.getGradYlm(34)[0] == Approx(1.50200576508));
  REQUIRE(ct.getGradYlm(34)[1] == Approx(1.62717291217));
  REQUIRE(ct.getGradYlm(34)[2] == Approx(-7.81042997841));

  REQUIRE(ct.getHessYlm(34)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(34)(0,1) == Approx(1.2516714709));
  REQUIRE(ct.getHessYlm(34)(0,2) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(34)(1,0) == Approx(1.2516714709));
  REQUIRE(ct.getHessYlm(34)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(34)(1,2) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(34)(2,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(34)(2,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(34)(2,2) == Approx(15.6208599568));

}

TEST_CASE("Cartesian Tensor evaluateWithThirdDeriv", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithThirdDeriv(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGGGYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(0)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(0)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(1)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(1)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(2) == Approx(0.586323014284));
  REQUIRE(ct.getGradYlm(2)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(2)[1] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(2)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(2)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(2)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(2)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(3) == Approx(-0.244301255951));
  REQUIRE(ct.getGradYlm(3)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(3)[2] == Approx(0.488602511903));

  REQUIRE(ct.getHessYlm(3)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(3)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(3)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(4) == Approx(1.06602349055));
  REQUIRE(ct.getGradYlm(4)[0] == Approx(1.64003613931));
  REQUIRE(ct.getGradYlm(4)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(4)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(4)(0,0) == Approx(1.26156626101));
  REQUIRE(ct.getHessYlm(4)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(4)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(4)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(5) == Approx(0.908327707927));
  REQUIRE(ct.getGradYlm(5)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(5)[1] == Approx(1.51387951321));
  REQUIRE(ct.getGradYlm(5)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(5)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(1,1) == Approx(1.26156626101));
  REQUIRE(ct.getHessYlm(5)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(5)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(5)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(6) == Approx(0.157695782626));
  REQUIRE(ct.getGradYlm(6)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(6)[2] == Approx(-0.630783130505));

  REQUIRE(ct.getHessYlm(6)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(6)(2,2) == Approx(1.26156626101));


  REQUIRE(ct.getGGGYlm(6)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(7) == Approx(1.70437555172));
  REQUIRE(ct.getGradYlm(7)[0] == Approx(1.31105811671));
  REQUIRE(ct.getGradYlm(7)[1] == Approx(1.42031295977));
  REQUIRE(ct.getGradYlm(7)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(7)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(0,1) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(7)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(1,0) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(7)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(7)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(7)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(8) == Approx(-0.710156479885));
  REQUIRE(ct.getGradYlm(8)[0] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(8)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(8)[2] == Approx(1.42031295977));

  REQUIRE(ct.getHessYlm(8)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(0,2) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(8)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(2,0) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(8)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(8)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(8)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(9) == Approx(-0.655529058355));
  REQUIRE(ct.getGradYlm(9)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(9)[1] == Approx(-0.546274215296));
  REQUIRE(ct.getGradYlm(9)[2] == Approx(1.31105811671));

  REQUIRE(ct.getHessYlm(9)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(1,2) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(9)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(9)(2,1) == Approx(1.09254843059));
  REQUIRE(ct.getHessYlm(9)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(9)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(10) == Approx(1.6397368054));
  REQUIRE(ct.getGradYlm(10)[0] == Approx(3.78400801246));
  REQUIRE(ct.getGradYlm(10)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(10)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(10)(0,0) == Approx(5.82155078841));
  REQUIRE(ct.getHessYlm(10)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(10)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(10)[0](0,0) == Approx(4.47811599108));
  REQUIRE(ct.getGGGYlm(10)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(11) == Approx(1.28969740543));
  REQUIRE(ct.getGradYlm(11)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(11)[1] == Approx(3.22424351358));
  REQUIRE(ct.getGradYlm(11)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(11)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(1,1) == Approx(5.3737391893));
  REQUIRE(ct.getHessYlm(11)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(11)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(11)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](1,1) == Approx(4.47811599108));
  REQUIRE(ct.getGGGYlm(11)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(12) == Approx(-0.0932940831475));
  REQUIRE(ct.getGradYlm(12)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(12)[2] == Approx(0.559764498885));

  REQUIRE(ct.getHessYlm(12)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(12)(2,2) == Approx(-2.23905799554));


  REQUIRE(ct.getGGGYlm(12)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,2) == Approx(4.47811599108));


  REQUIRE(ct.getYlm(13) == Approx(3.38451965731));
  REQUIRE(ct.getGradYlm(13)[0] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(13)[1] == Approx(2.82043304776));
  REQUIRE(ct.getGradYlm(13)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(13)(0,0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(13)(0,1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(13)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(1,0) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(13)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(13)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(13)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](0,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](1,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](0,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(14) == Approx(-1.41021652388));
  REQUIRE(ct.getGradYlm(14)[0] == Approx(-2.16956388289));
  REQUIRE(ct.getGradYlm(14)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(14)[2] == Approx(2.82043304776));

  REQUIRE(ct.getHessYlm(14)(0,0) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(14)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(0,2) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(14)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(2,0) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(14)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(14)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(14)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  REQUIRE(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(15)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(15)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(0,1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(1,0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(1,1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(15)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(15)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](1,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](0,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](1,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(16) == Approx(-1.20160461206));
  REQUIRE(ct.getGradYlm(16)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(16)[1] == Approx(-2.00267435344));
  REQUIRE(ct.getGradYlm(16)[2] == Approx(2.40320922412));

  REQUIRE(ct.getHessYlm(16)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(1,1) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(16)(1,2) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(16)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(16)(2,1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(16)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(16)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](2,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](1,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(17) == Approx(0.542390970723));
  REQUIRE(ct.getGradYlm(17)[0] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(17)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(17)[2] == Approx(-2.16956388289));

  REQUIRE(ct.getHessYlm(17)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(0,2) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(17)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(2,0) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(17)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(17)(2,2) == Approx(4.33912776578));


  REQUIRE(ct.getGGGYlm(17)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](2,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(18) == Approx(0.500668588359));
  REQUIRE(ct.getGradYlm(18)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(18)[1] == Approx(0.417223823633));
  REQUIRE(ct.getGradYlm(18)[2] == Approx(-2.00267435344));

  REQUIRE(ct.getHessYlm(18)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(1,2) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(18)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(18)(2,1) == Approx(-1.66889529453));
  REQUIRE(ct.getHessYlm(18)(2,2) == Approx(4.00534870687));


  REQUIRE(ct.getGGGYlm(18)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](2,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(19) == Approx(-2.25467692526));
  REQUIRE(ct.getGradYlm(19)[0] == Approx(-1.73436686558));
  REQUIRE(ct.getGradYlm(19)[1] == Approx(-1.87889743772));
  REQUIRE(ct.getGradYlm(19)[2] == Approx(4.50935385052));

  REQUIRE(ct.getHessYlm(19)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(19)(0,1) == Approx(-1.44530572132));
  REQUIRE(ct.getHessYlm(19)(0,2) == Approx(3.46873373117));
  REQUIRE(ct.getHessYlm(19)(1,0) == Approx(-1.44530572132));
  REQUIRE(ct.getHessYlm(19)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(19)(1,2) == Approx(3.75779487543));
  REQUIRE(ct.getHessYlm(19)(2,0) == Approx(3.46873373117));
  REQUIRE(ct.getHessYlm(19)(2,1) == Approx(3.75779487543));
  REQUIRE(ct.getHessYlm(19)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(19)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,2) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](2,1) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,2) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](2,0) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](0,1) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](1,0) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(20) == Approx(2.41707280436));
  REQUIRE(ct.getGradYlm(20)[0] == Approx(7.43714709033));
  REQUIRE(ct.getGradYlm(20)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(20)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(20)(0,0) == Approx(17.1626471315));
  REQUIRE(ct.getHessYlm(20)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(20)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(20)[0](0,0) == Approx(26.40407251));
  REQUIRE(ct.getGGGYlm(20)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(21) == Approx(1.75485528067));
  REQUIRE(ct.getGradYlm(21)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(21)[1] == Approx(5.84951760222));
  REQUIRE(ct.getGradYlm(21)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(21)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(1,1) == Approx(14.6237940056));
  REQUIRE(ct.getHessYlm(21)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(21)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(21)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](1,1) == Approx(24.3729900093));
  REQUIRE(ct.getGGGYlm(21)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(22) == Approx(0.0528927734576));
  REQUIRE(ct.getGradYlm(22)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(22)[2] == Approx(-0.423142187661));

  REQUIRE(ct.getHessYlm(22)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(22)(2,2) == Approx(2.53885312596));


  REQUIRE(ct.getGGGYlm(22)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,2) == Approx(-10.1554125039));


  REQUIRE(ct.getYlm(23) == Approx(5.90305249944));
  REQUIRE(ct.getGradYlm(23)[0] == Approx(13.6224288449));
  REQUIRE(ct.getGradYlm(23)[1] == Approx(4.9192104162));
  REQUIRE(ct.getGradYlm(23)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(23)(0,0) == Approx(20.9575828383));
  REQUIRE(ct.getHessYlm(23)(0,1) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(23)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(1,0) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(23)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(23)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(23)[0](0,0) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(23)[0](0,1) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](1,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](0,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(24) == Approx(-2.4596052081));
  REQUIRE(ct.getGradYlm(24)[0] == Approx(-5.6760120187));
  REQUIRE(ct.getGradYlm(24)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(24)[2] == Approx(4.9192104162));

  REQUIRE(ct.getHessYlm(24)(0,0) == Approx(-8.73232618261));
  REQUIRE(ct.getHessYlm(24)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(0,2) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(24)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(2,0) == Approx(11.3520240374));
  REQUIRE(ct.getHessYlm(24)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(24)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(24)[0](0,0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(24)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](0,2) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](2,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](0,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(25) == Approx(5.02981988118));
  REQUIRE(ct.getGradYlm(25)[0] == Approx(3.86909221629));
  REQUIRE(ct.getGradYlm(25)[1] == Approx(12.574549703));
  REQUIRE(ct.getGradYlm(25)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(25)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(0,1) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(25)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(1,0) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(25)(1,1) == Approx(20.9575828383));
  REQUIRE(ct.getHessYlm(25)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(25)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(25)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](1,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](0,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](1,0) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[1](1,1) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(25)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(26) == Approx(-1.93454610815));
  REQUIRE(ct.getGradYlm(26)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(26)[1] == Approx(-4.83636527037));
  REQUIRE(ct.getGradYlm(26)[2] == Approx(3.86909221629));

  REQUIRE(ct.getHessYlm(26)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(1,1) == Approx(-8.06060878395));
  REQUIRE(ct.getHessYlm(26)(1,2) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(26)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(26)(2,1) == Approx(9.67273054074));
  REQUIRE(ct.getHessYlm(26)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(26)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](1,1) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(26)[1](1,2) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](2,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](1,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(27) == Approx(-0.363846924275));
  REQUIRE(ct.getGradYlm(27)[0] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(27)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(27)[2] == Approx(2.18308154565));

  REQUIRE(ct.getHessYlm(27)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0,2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,0) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2,2) == Approx(-8.73232618261));


  REQUIRE(ct.getGGGYlm(27)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,2) == Approx(17.4646523652));


  REQUIRE(ct.getYlm(28) == Approx(-0.335858699331));
  REQUIRE(ct.getGradYlm(28)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(28)[1] == Approx(-0.279882249443));
  REQUIRE(ct.getGradYlm(28)[2] == Approx(2.01515219599));

  REQUIRE(ct.getHessYlm(28)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(1,2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(28)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(28)(2,1) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(28)(2,2) == Approx(-8.06060878395));


  REQUIRE(ct.getGGGYlm(28)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](2,1) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](2,2) == Approx(16.1212175679));


  REQUIRE(ct.getYlm(29) == Approx(7.03459200681));
  REQUIRE(ct.getGradYlm(29)[0] == Approx(10.8224492412));
  REQUIRE(ct.getGradYlm(29)[1] == Approx(11.7243200114));
  REQUIRE(ct.getGradYlm(29)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(29)(0,0) == Approx(8.3249609548));
  REQUIRE(ct.getHessYlm(29)(0,1) == Approx(18.0374154021));
  REQUIRE(ct.getHessYlm(29)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(1,0) == Approx(18.0374154021));
  REQUIRE(ct.getHessYlm(29)(1,1) == Approx(9.77026667613));
  REQUIRE(ct.getHessYlm(29)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(29)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(29)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](0,1) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](1,0) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[0](1,1) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](0,0) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[1](0,1) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](1,0) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(30) == Approx(1.22128333452));
  REQUIRE(ct.getGradYlm(30)[0] == Approx(1.87889743772));
  REQUIRE(ct.getGradYlm(30)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(30)[2] == Approx(-4.88513333806));

  REQUIRE(ct.getHessYlm(30)(0,0) == Approx(1.44530572132));
  REQUIRE(ct.getHessYlm(30)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(0,2) == Approx(-7.51558975087));
  REQUIRE(ct.getHessYlm(30)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(1,2) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(2,0) == Approx(-7.51558975087));
  REQUIRE(ct.getHessYlm(30)(2,1) == Approx(0));
  REQUIRE(ct.getHessYlm(30)(2,2) == Approx(9.77026667613));


  REQUIRE(ct.getGGGYlm(30)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](0,2) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](2,0) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](2,2) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](0,0) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](0,2) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](2,0) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(31) == Approx(1.04062011935));
  REQUIRE(ct.getGradYlm(31)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(31)[1] == Approx(1.73436686558));
  REQUIRE(ct.getGradYlm(31)[2] == Approx(-4.1624804774));

  REQUIRE(ct.getHessYlm(31)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(0,1) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(0,2) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(1,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(1,1) == Approx(1.44530572132));
  REQUIRE(ct.getHessYlm(31)(1,2) == Approx(-6.93746746234));
  REQUIRE(ct.getHessYlm(31)(2,0) == Approx(0));
  REQUIRE(ct.getHessYlm(31)(2,1) == Approx(-6.93746746234));
  REQUIRE(ct.getHessYlm(31)(2,2) == Approx(8.3249609548));


  REQUIRE(ct.getGGGYlm(31)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,2) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](2,1) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[1](2,2) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](1,1) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[2](1,2) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](2,1) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  REQUIRE(ct.getGradYlm(32)[2] == Approx(10.1535589719));

  REQUIRE(ct.getHessYlm(32)(0,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(32)(0,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(0,2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(1,0) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(32)(1,2) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,0) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(2,1) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(32)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[0](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[0](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[0](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[0](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[0](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[0](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](0,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](0,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](2,0) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](0,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[2](0,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](1,0) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(33) == Approx(-4.68625798704));
  REQUIRE(ct.getGradYlm(33)[0] == Approx(-3.60481383619));
  REQUIRE(ct.getGradYlm(33)[1] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(33)[2] == Approx(9.37251597409));

  REQUIRE(ct.getHessYlm(33)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(33)(0,1) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(33)(0,2) == Approx(7.20962767237));
  REQUIRE(ct.getHessYlm(33)(1,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(33)(1,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(33)(1,2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(33)(2,0) == Approx(7.20962767237));
  REQUIRE(ct.getHessYlm(33)(2,1) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(33)(2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(33)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[0](1,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](2,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](1,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](1,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,2) == Approx(0));


  REQUIRE(ct.getYlm(34) == Approx(1.9526074946));
  REQUIRE(ct.getGradYlm(34)[0] == Approx(1.50200576508));
  REQUIRE(ct.getGradYlm(34)[1] == Approx(1.62717291217));
  REQUIRE(ct.getGradYlm(34)[2] == Approx(-7.81042997841));

  REQUIRE(ct.getHessYlm(34)(0,0) == Approx(0));
  REQUIRE(ct.getHessYlm(34)(0,1) == Approx(1.2516714709));
  REQUIRE(ct.getHessYlm(34)(0,2) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(34)(1,0) == Approx(1.2516714709));
  REQUIRE(ct.getHessYlm(34)(1,1) == Approx(0));
  REQUIRE(ct.getHessYlm(34)(1,2) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(34)(2,0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(34)(2,1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(34)(2,2) == Approx(15.6208599568));


  REQUIRE(ct.getGGGYlm(34)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,2) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](2,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[0](2,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](0,2) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](2,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](2,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[2](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[2](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[2](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[2](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[2](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](2,2) == Approx(0));

}


TEST_CASE("Cartesian Tensor evaluateThirdDerivOnly", "[numerics]") {
  CartesianTensor<double, TinyVector<double, 3>> ct(4);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateThirdDerivOnly(pt);

  REQUIRE(ct.getGGGYlm(0)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(0)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(1)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(1)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(2)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(2)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(3)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(3)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(4)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(4)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(5)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(5)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(6)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(6)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(7)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(7)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(8)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(8)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(9)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(9)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(10)[0](0,0) == Approx(4.47811599108));
  REQUIRE(ct.getGGGYlm(10)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(10)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(11)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](1,1) == Approx(4.47811599108));
  REQUIRE(ct.getGGGYlm(11)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(11)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(12)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(12)[2](2,2) == Approx(4.47811599108));


  REQUIRE(ct.getGGGYlm(13)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](0,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](1,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](0,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(13)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(13)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(14)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(15)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](1,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](0,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](1,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(15)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(15)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(16)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](1,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[1](2,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](1,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(16)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(16)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(17)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[0](2,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](0,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](2,0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(17)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(17)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(18)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[1](2,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](1,2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(18)[2](2,1) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(18)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(19)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](1,2) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[0](2,1) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](0,2) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](2,0) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](0,1) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](1,0) == Approx(2.89061144264));
  REQUIRE(ct.getGGGYlm(19)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(19)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(20)[0](0,0) == Approx(26.40407251));
  REQUIRE(ct.getGGGYlm(20)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(20)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(21)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](1,1) == Approx(24.3729900093));
  REQUIRE(ct.getGGGYlm(21)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(21)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(22)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(22)[2](2,2) == Approx(-10.1554125039));


  REQUIRE(ct.getGGGYlm(23)[0](0,0) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(23)[0](0,1) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](1,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](0,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(23)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(23)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(24)[0](0,0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(24)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](0,2) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](2,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](0,0) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(24)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(24)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(25)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](1,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](0,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](1,0) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(25)[1](1,1) == Approx(17.4646523652));
  REQUIRE(ct.getGGGYlm(25)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(25)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(26)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](1,1) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(26)[1](1,2) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[1](2,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](1,1) == Approx(16.1212175679));
  REQUIRE(ct.getGGGYlm(26)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(26)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(27)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2,2) == Approx(17.4646523652));


  REQUIRE(ct.getGGGYlm(28)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[1](2,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](1,2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(28)[2](2,1) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(28)[2](2,2) == Approx(16.1212175679));


  REQUIRE(ct.getGGGYlm(29)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](0,1) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](1,0) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[0](1,1) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](0,0) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(29)[1](0,1) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](1,0) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(29)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(29)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(30)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](0,2) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](2,0) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[0](2,2) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](0,0) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(30)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](0,2) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](2,0) == Approx(15.0311795017));
  REQUIRE(ct.getGGGYlm(30)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(30)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(31)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](1,2) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[1](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[1](2,1) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[1](2,2) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](1,1) == Approx(-5.78122288528));
  REQUIRE(ct.getGGGYlm(31)[2](1,2) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(31)[2](2,1) == Approx(13.8749349247));
  REQUIRE(ct.getGGGYlm(31)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(32)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[0](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[0](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[0](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[0](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[0](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[0](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](0,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(32)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](0,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](2,0) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](0,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(32)[2](0,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](1,0) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(32)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(32)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(33)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[0](1,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](2,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0,1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](1,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](1,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[2](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2,2) == Approx(0));


  REQUIRE(ct.getGGGYlm(34)[0](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](0,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](1,2) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[0](2,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[0](2,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[0](2,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[1](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](0,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](0,2) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[1](1,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](1,2) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](2,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[1](2,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[1](2,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](0,0) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[2](0,1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[2](0,2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[2](1,0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(34)[2](1,1) == Approx(0));
  REQUIRE(ct.getGGGYlm(34)[2](1,2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](2,0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(34)[2](2,1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(34)[2](2,2) == Approx(0));


}
#endif

}


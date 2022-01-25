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
using real_type = OHMMS_PRECISION;


// Use gen_gto.py to generate the checks

TEST_CASE("Cartesian Tensor", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

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
  REQUIRE(ct.getYlm(35) == Approx(3.47382688598));
  REQUIRE(ct.getYlm(36) == Approx(2.32807861094));
  REQUIRE(ct.getYlm(37) == Approx(-0.0292375806134));
  REQUIRE(ct.getYlm(38) == Approx(9.61982829963));
  REQUIRE(ct.getYlm(39) == Approx(-4.00826179151));
  REQUIRE(ct.getYlm(40) == Approx(7.56625548555));
  REQUIRE(ct.getYlm(41) == Approx(-2.91009826367));
  REQUIRE(ct.getYlm(42) == Approx(0.228053128784));
  REQUIRE(ct.getYlm(43) == Approx(0.210510580416));
  REQUIRE(ct.getYlm(44) == Approx(13.5641819555));
  REQUIRE(ct.getYlm(45) == Approx(2.35489270061));
  REQUIRE(ct.getYlm(46) == Approx(12.5207833436));
  REQUIRE(ct.getYlm(47) == Approx(1.85218688514));
  REQUIRE(ct.getYlm(48) == Approx(-0.905727961775));
  REQUIRE(ct.getYlm(49) == Approx(-0.771744535477));
  REQUIRE(ct.getYlm(50) == Approx(-9.78910512921));
  REQUIRE(ct.getYlm(51) == Approx(-8.34101265448));
  REQUIRE(ct.getYlm(52) == Approx(-1.44809247474));
  REQUIRE(ct.getYlm(53) == Approx(-11.6655511199));
  REQUIRE(ct.getYlm(54) == Approx(4.86064629996));
  REQUIRE(ct.getYlm(55) == Approx(4.48675043074));
  REQUIRE(ct.getYlm(56) == Approx(4.90938236205));
  REQUIRE(ct.getYlm(57) == Approx(3.03706593382));
  REQUIRE(ct.getYlm(58) == Approx(0.0158923005669));
  REQUIRE(ct.getYlm(59) == Approx(15.0300731514));
  REQUIRE(ct.getYlm(60) == Approx(-6.26253047974));
  REQUIRE(ct.getYlm(61) == Approx(10.9122088466));
  REQUIRE(ct.getYlm(62) == Approx(-4.19700340252));
  REQUIRE(ct.getYlm(63) == Approx(-0.137042874894));
  REQUIRE(ct.getYlm(64) == Approx(-0.126501115286));
  REQUIRE(ct.getYlm(65) == Approx(24.0303233904));
  REQUIRE(ct.getYlm(66) == Approx(4.17193114417));
  REQUIRE(ct.getYlm(67) == Approx(20.4755418238));
  REQUIRE(ct.getYlm(68) == Approx(3.0289263053));
  REQUIRE(ct.getYlm(69) == Approx(0.61714957754));
  REQUIRE(ct.getYlm(70) == Approx(0.525855261336));
  REQUIRE(ct.getYlm(71) == Approx(-17.3423920977));
  REQUIRE(ct.getYlm(72) == Approx(-13.6402610582));
  REQUIRE(ct.getYlm(73) == Approx(0.986708699234));
  REQUIRE(ct.getYlm(74) == Approx(26.2459034569));
  REQUIRE(ct.getYlm(75) == Approx(-1.89857519219));
  REQUIRE(ct.getYlm(76) == Approx(-1.49328080661));
  REQUIRE(ct.getYlm(77) == Approx(-24.4531767752));
  REQUIRE(ct.getYlm(78) == Approx(10.1888236563));
  REQUIRE(ct.getYlm(79) == Approx(-22.5721631771));
  REQUIRE(ct.getYlm(80) == Approx(8.68160122197));
  REQUIRE(ct.getYlm(81) == Approx(-3.91877832936));
  REQUIRE(ct.getYlm(82) == Approx(-3.61733384249));
  REQUIRE(ct.getYlm(83) == Approx(12.1418905657));
}

TEST_CASE("Cartesian Tensor evaluateAll subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

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

  REQUIRE(ct.getYlm(50) == Approx(-9.78910512921));
  REQUIRE(ct.getGradYlm(50)[0] == Approx(-22.5902426059));
  REQUIRE(ct.getGradYlm(50)[1] == Approx(-8.15758760768));
  REQUIRE(ct.getGradYlm(50)[2] == Approx(19.5782102584));
  REQUIRE(ct.getLaplYlm(50) == Approx(-34.7542193937));

  REQUIRE(ct.getYlm(83) == Approx(12.1418905657));
  REQUIRE(ct.getGradYlm(83)[0] == Approx(18.6798316395));
  REQUIRE(ct.getGradYlm(83)[1] == Approx(20.2364842761));
  REQUIRE(ct.getGradYlm(83)[2] == Approx(-48.5675622627));
  REQUIRE(ct.getLaplYlm(83) == Approx(128.367962683));
}

TEST_CASE("Cartesian Tensor evaluateWithHessian subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithHessian(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getHessYlm(i) <<  std::endl;
  //}

  REQUIRE(ct.getYlm(0) == Approx(0.282094791774));
  REQUIRE(ct.getGradYlm(0)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(0)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(0)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(0, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(1, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(0)(2, 2) == Approx(0));


  REQUIRE(ct.getYlm(1) == Approx(0.635183265474));
  REQUIRE(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  REQUIRE(ct.getGradYlm(1)[1] == Approx(0));
  REQUIRE(ct.getGradYlm(1)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(1)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(0, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(1, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(1)(2, 2) == Approx(0));


  REQUIRE(ct.getYlm(15) == Approx(3.12417199136));
  REQUIRE(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  REQUIRE(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  REQUIRE(ct.getGradYlm(15)[2] == Approx(0));

  REQUIRE(ct.getHessYlm(15)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(0, 1) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(0, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(1, 0) == Approx(4.00534870687));
  REQUIRE(ct.getHessYlm(15)(1, 1) == Approx(4.33912776578));
  REQUIRE(ct.getHessYlm(15)(1, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(15)(2, 2) == Approx(0));


  REQUIRE(ct.getYlm(32) == Approx(-5.07677948596));
  REQUIRE(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  REQUIRE(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  REQUIRE(ct.getGradYlm(32)[2] == Approx(10.1535589719));

  REQUIRE(ct.getHessYlm(32)(0, 0) == Approx(-6.00802306031));
  REQUIRE(ct.getHessYlm(32)(0, 1) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(0, 2) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(1, 0) == Approx(-6.50869164867));
  REQUIRE(ct.getHessYlm(32)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(32)(1, 2) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2, 0) == Approx(15.6208599568));
  REQUIRE(ct.getHessYlm(32)(2, 1) == Approx(8.46129914327));
  REQUIRE(ct.getHessYlm(32)(2, 2) == Approx(0));


  REQUIRE(ct.getYlm(52) == Approx(-1.44809247474));
  REQUIRE(ct.getGradYlm(52)[0] == Approx(-1.11391728826));
  REQUIRE(ct.getGradYlm(52)[1] == Approx(-1.20674372895));
  REQUIRE(ct.getGradYlm(52)[2] == Approx(8.68855484841));

  REQUIRE(ct.getHessYlm(52)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(52)(0, 1) == Approx(-0.928264406882));
  REQUIRE(ct.getHessYlm(52)(0, 2) == Approx(6.68350372955));
  REQUIRE(ct.getHessYlm(52)(1, 0) == Approx(-0.928264406882));
  REQUIRE(ct.getHessYlm(52)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(52)(1, 2) == Approx(7.24046237368));
  REQUIRE(ct.getHessYlm(52)(2, 0) == Approx(6.68350372955));
  REQUIRE(ct.getHessYlm(52)(2, 1) == Approx(7.24046237368));
  REQUIRE(ct.getHessYlm(52)(2, 2) == Approx(-34.7542193937));


  REQUIRE(ct.getYlm(71) == Approx(-17.3423920977));
  REQUIRE(ct.getGradYlm(71)[0] == Approx(-53.3612064546));
  REQUIRE(ct.getGradYlm(71)[1] == Approx(-14.4519934148));
  REQUIRE(ct.getGradYlm(71)[2] == Approx(34.6847841955));

  REQUIRE(ct.getHessYlm(71)(0, 0) == Approx(-123.141245664));
  REQUIRE(ct.getHessYlm(71)(0, 1) == Approx(-44.4676720455));
  REQUIRE(ct.getHessYlm(71)(0, 2) == Approx(106.722412909));
  REQUIRE(ct.getHessYlm(71)(1, 0) == Approx(-44.4676720455));
  REQUIRE(ct.getHessYlm(71)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(71)(1, 2) == Approx(28.9039868295));
  REQUIRE(ct.getHessYlm(71)(2, 0) == Approx(106.722412909));
  REQUIRE(ct.getHessYlm(71)(2, 1) == Approx(28.9039868295));
  REQUIRE(ct.getHessYlm(71)(2, 2) == Approx(0));
}

TEST_CASE("Cartesian Tensor evaluateWithThirdDeriv subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

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

  REQUIRE(ct.getHessYlm(27)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(0, 2) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(1, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(1, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2, 0) == Approx(1.67929349666));
  REQUIRE(ct.getHessYlm(27)(2, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(27)(2, 2) == Approx(-8.73232618261));


  REQUIRE(ct.getGGGYlm(27)[0](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[0](2, 2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[1](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[1](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](0, 2) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2, 0) == Approx(-6.71717398662));
  REQUIRE(ct.getGGGYlm(27)[2](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(27)[2](2, 2) == Approx(17.4646523652));


  REQUIRE(ct.getYlm(62) == Approx(-4.19700340252));
  REQUIRE(ct.getGradYlm(62)[0] == Approx(0));
  REQUIRE(ct.getGradYlm(62)[1] == Approx(-17.4875141772));
  REQUIRE(ct.getGradYlm(62)[2] == Approx(8.39400680505));

  REQUIRE(ct.getHessYlm(62)(0, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(62)(0, 1) == Approx(0));
  REQUIRE(ct.getHessYlm(62)(0, 2) == Approx(0));
  REQUIRE(ct.getHessYlm(62)(1, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(62)(1, 1) == Approx(-58.291713924));
  REQUIRE(ct.getHessYlm(62)(1, 2) == Approx(34.9750283544));
  REQUIRE(ct.getHessYlm(62)(2, 0) == Approx(0));
  REQUIRE(ct.getHessYlm(62)(2, 1) == Approx(34.9750283544));
  REQUIRE(ct.getHessYlm(62)(2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[0](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](1, 1) == Approx(-145.72928481));
  REQUIRE(ct.getGGGYlm(62)[1](1, 2) == Approx(116.583427848));
  REQUIRE(ct.getGGGYlm(62)[1](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[1](2, 1) == Approx(116.583427848));
  REQUIRE(ct.getGGGYlm(62)[1](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](1, 1) == Approx(116.583427848));
  REQUIRE(ct.getGGGYlm(62)[2](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(62)[2](2, 2) == Approx(0));
}


TEST_CASE("Cartesian Tensor evaluateThirdDerivOnly subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateThirdDerivOnly(pt);

  REQUIRE(ct.getGGGYlm(14)[0](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](0, 2) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2, 0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[0](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[0](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[1](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0, 0) == Approx(3.33779058906));
  REQUIRE(ct.getGGGYlm(14)[2](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(14)[2](2, 2) == Approx(0));


  REQUIRE(ct.getGGGYlm(33)[0](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](1, 1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[0](1, 2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[0](2, 1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[0](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](0, 1) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](0, 2) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](1, 0) == Approx(-5.00668588359));
  REQUIRE(ct.getGGGYlm(33)[1](1, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[1](1, 2) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2, 0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[1](2, 1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[1](2, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](0, 1) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](1, 0) == Approx(12.0160461206));
  REQUIRE(ct.getGGGYlm(33)[2](1, 1) == Approx(13.0173832973));
  REQUIRE(ct.getGGGYlm(33)[2](1, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(33)[2](2, 2) == Approx(0));


  REQUIRE(ct.getGGGYlm(80)[0](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[0](0, 1) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[0](0, 2) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[0](1, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[0](1, 1) == Approx(27.8256449422));
  REQUIRE(ct.getGGGYlm(80)[0](1, 2) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[0](2, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[0](2, 1) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[0](2, 2) == Approx(53.425238289));
  REQUIRE(ct.getGGGYlm(80)[1](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[1](0, 1) == Approx(27.8256449422));
  REQUIRE(ct.getGGGYlm(80)[1](0, 2) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[1](1, 0) == Approx(27.8256449422));
  REQUIRE(ct.getGGGYlm(80)[1](1, 1) == Approx(30.1444486874));
  REQUIRE(ct.getGGGYlm(80)[1](1, 2) == Approx(-144.6933537));
  REQUIRE(ct.getGGGYlm(80)[1](2, 0) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[1](2, 1) == Approx(-144.6933537));
  REQUIRE(ct.getGGGYlm(80)[1](2, 2) == Approx(173.632024439));
  REQUIRE(ct.getGGGYlm(80)[2](0, 0) == Approx(0));
  REQUIRE(ct.getGGGYlm(80)[2](0, 1) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[2](0, 2) == Approx(53.425238289));
  REQUIRE(ct.getGGGYlm(80)[2](1, 0) == Approx(-66.7815478613));
  REQUIRE(ct.getGGGYlm(80)[2](1, 1) == Approx(-144.6933537));
  REQUIRE(ct.getGGGYlm(80)[2](1, 2) == Approx(173.632024439));
  REQUIRE(ct.getGGGYlm(80)[2](2, 0) == Approx(53.425238289));
  REQUIRE(ct.getGGGYlm(80)[2](2, 1) == Approx(173.632024439));
  REQUIRE(ct.getGGGYlm(80)[2](2, 2) == Approx(0));
}
} // namespace qmcplusplus

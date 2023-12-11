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
  CHECK(ct.getYlm(0) == Approx(0.282094791774));
  CHECK(ct.getYlm(1) == Approx(0.635183265474));
  CHECK(ct.getYlm(2) == Approx(0.586323014284));
  CHECK(ct.getYlm(3) == Approx(-0.244301255951));
  CHECK(ct.getYlm(4) == Approx(1.06602349055));
  CHECK(ct.getYlm(5) == Approx(0.908327707927));
  CHECK(ct.getYlm(6) == Approx(0.157695782626));
  CHECK(ct.getYlm(7) == Approx(1.70437555172));
  CHECK(ct.getYlm(8) == Approx(-0.710156479885));
  CHECK(ct.getYlm(9) == Approx(-0.655529058355));
  CHECK(ct.getYlm(10) == Approx(1.6397368054));
  CHECK(ct.getYlm(11) == Approx(1.28969740543));
  CHECK(ct.getYlm(12) == Approx(-0.0932940831475));
  CHECK(ct.getYlm(13) == Approx(3.38451965731));
  CHECK(ct.getYlm(14) == Approx(-1.41021652388));
  CHECK(ct.getYlm(15) == Approx(3.12417199136));
  CHECK(ct.getYlm(16) == Approx(-1.20160461206));
  CHECK(ct.getYlm(17) == Approx(0.542390970723));
  CHECK(ct.getYlm(18) == Approx(0.500668588359));
  CHECK(ct.getYlm(19) == Approx(-2.25467692526));
  CHECK(ct.getYlm(20) == Approx(2.41707280436));
  CHECK(ct.getYlm(21) == Approx(1.75485528067));
  CHECK(ct.getYlm(22) == Approx(0.0528927734576));
  CHECK(ct.getYlm(23) == Approx(5.90305249944));
  CHECK(ct.getYlm(24) == Approx(-2.4596052081));
  CHECK(ct.getYlm(25) == Approx(5.02981988118));
  CHECK(ct.getYlm(26) == Approx(-1.93454610815));
  CHECK(ct.getYlm(27) == Approx(-0.363846924275));
  CHECK(ct.getYlm(28) == Approx(-0.335858699331));
  CHECK(ct.getYlm(29) == Approx(7.03459200681));
  CHECK(ct.getYlm(30) == Approx(1.22128333452));
  CHECK(ct.getYlm(31) == Approx(1.04062011935));
  CHECK(ct.getYlm(32) == Approx(-5.07677948596));
  CHECK(ct.getYlm(33) == Approx(-4.68625798704));
  CHECK(ct.getYlm(34) == Approx(1.9526074946));
  CHECK(ct.getYlm(35) == Approx(3.47382688598));
  CHECK(ct.getYlm(36) == Approx(2.32807861094));
  CHECK(ct.getYlm(37) == Approx(-0.0292375806134));
  CHECK(ct.getYlm(38) == Approx(9.61982829963));
  CHECK(ct.getYlm(39) == Approx(-4.00826179151));
  CHECK(ct.getYlm(40) == Approx(7.56625548555));
  CHECK(ct.getYlm(41) == Approx(-2.91009826367));
  CHECK(ct.getYlm(42) == Approx(0.228053128784));
  CHECK(ct.getYlm(43) == Approx(0.210510580416));
  CHECK(ct.getYlm(44) == Approx(13.5641819555));
  CHECK(ct.getYlm(45) == Approx(2.35489270061));
  CHECK(ct.getYlm(46) == Approx(12.5207833436));
  CHECK(ct.getYlm(47) == Approx(1.85218688514));
  CHECK(ct.getYlm(48) == Approx(-0.905727961775));
  CHECK(ct.getYlm(49) == Approx(-0.771744535477));
  CHECK(ct.getYlm(50) == Approx(-9.78910512921));
  CHECK(ct.getYlm(51) == Approx(-8.34101265448));
  CHECK(ct.getYlm(52) == Approx(-1.44809247474));
  CHECK(ct.getYlm(53) == Approx(-11.6655511199));
  CHECK(ct.getYlm(54) == Approx(4.86064629996));
  CHECK(ct.getYlm(55) == Approx(4.48675043074));
  CHECK(ct.getYlm(56) == Approx(4.90938236205));
  CHECK(ct.getYlm(57) == Approx(3.03706593382));
  CHECK(ct.getYlm(58) == Approx(0.0158923005669));
  CHECK(ct.getYlm(59) == Approx(15.0300731514));
  CHECK(ct.getYlm(60) == Approx(-6.26253047974));
  CHECK(ct.getYlm(61) == Approx(10.9122088466));
  CHECK(ct.getYlm(62) == Approx(-4.19700340252));
  CHECK(ct.getYlm(63) == Approx(-0.137042874894));
  CHECK(ct.getYlm(64) == Approx(-0.126501115286));
  CHECK(ct.getYlm(65) == Approx(24.0303233904));
  CHECK(ct.getYlm(66) == Approx(4.17193114417));
  CHECK(ct.getYlm(67) == Approx(20.4755418238));
  CHECK(ct.getYlm(68) == Approx(3.0289263053));
  CHECK(ct.getYlm(69) == Approx(0.61714957754));
  CHECK(ct.getYlm(70) == Approx(0.525855261336));
  CHECK(ct.getYlm(71) == Approx(-17.3423920977));
  CHECK(ct.getYlm(72) == Approx(-13.6402610582));
  CHECK(ct.getYlm(73) == Approx(0.986708699234));
  CHECK(ct.getYlm(74) == Approx(26.2459034569));
  CHECK(ct.getYlm(75) == Approx(-1.89857519219));
  CHECK(ct.getYlm(76) == Approx(-1.49328080661));
  CHECK(ct.getYlm(77) == Approx(-24.4531767752));
  CHECK(ct.getYlm(78) == Approx(10.1888236563));
  CHECK(ct.getYlm(79) == Approx(-22.5721631771));
  CHECK(ct.getYlm(80) == Approx(8.68160122197));
  CHECK(ct.getYlm(81) == Approx(-3.91877832936));
  CHECK(ct.getYlm(82) == Approx(-3.61733384249));
  CHECK(ct.getYlm(83) == Approx(12.1418905657));
}

TEST_CASE("Cartesian Tensor evaluateAll subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateAll(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGradYlm(i) << "  " << ct.getLaplYlm(i) <<  std::endl;
  //}

  CHECK(ct.getYlm(0) == Approx(0.282094791774));
  CHECK(ct.getGradYlm(0)[0] == Approx(0));
  CHECK(ct.getGradYlm(0)[1] == Approx(0));
  CHECK(ct.getGradYlm(0)[2] == Approx(0));
  CHECK(ct.getLaplYlm(0) == Approx(0));

  CHECK(ct.getYlm(1) == Approx(0.635183265474));
  CHECK(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  CHECK(ct.getGradYlm(1)[1] == Approx(0));
  CHECK(ct.getGradYlm(1)[2] == Approx(0));
  CHECK(ct.getLaplYlm(1) == Approx(0));

  CHECK(ct.getYlm(8) == Approx(-0.710156479885));
  CHECK(ct.getGradYlm(8)[0] == Approx(-0.546274215296));
  CHECK(ct.getGradYlm(8)[1] == Approx(0));
  CHECK(ct.getGradYlm(8)[2] == Approx(1.42031295977));
  CHECK(ct.getLaplYlm(8) == Approx(0));

  CHECK(ct.getYlm(23) == Approx(5.90305249944));
  CHECK(ct.getGradYlm(23)[0] == Approx(13.6224288449));
  CHECK(ct.getGradYlm(23)[1] == Approx(4.9192104162));
  CHECK(ct.getGradYlm(23)[2] == Approx(0));
  CHECK(ct.getLaplYlm(23) == Approx(20.9575828383));

  CHECK(ct.getYlm(34) == Approx(1.9526074946));
  CHECK(ct.getGradYlm(34)[0] == Approx(1.50200576508));
  CHECK(ct.getGradYlm(34)[1] == Approx(1.62717291217));
  CHECK(ct.getGradYlm(34)[2] == Approx(-7.81042997841));
  CHECK(ct.getLaplYlm(34) == Approx(15.6208599568));

  CHECK(ct.getYlm(50) == Approx(-9.78910512921));
  CHECK(ct.getGradYlm(50)[0] == Approx(-22.5902426059));
  CHECK(ct.getGradYlm(50)[1] == Approx(-8.15758760768));
  CHECK(ct.getGradYlm(50)[2] == Approx(19.5782102584));
  CHECK(ct.getLaplYlm(50) == Approx(-34.7542193937));

  CHECK(ct.getYlm(83) == Approx(12.1418905657));
  CHECK(ct.getGradYlm(83)[0] == Approx(18.6798316395));
  CHECK(ct.getGradYlm(83)[1] == Approx(20.2364842761));
  CHECK(ct.getGradYlm(83)[2] == Approx(-48.5675622627));
  CHECK(ct.getLaplYlm(83) == Approx(128.367962683));
}

TEST_CASE("Cartesian Tensor evaluateWithHessian subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithHessian(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getHessYlm(i) <<  std::endl;
  //}

  CHECK(ct.getYlm(0) == Approx(0.282094791774));
  CHECK(ct.getGradYlm(0)[0] == Approx(0));
  CHECK(ct.getGradYlm(0)[1] == Approx(0));
  CHECK(ct.getGradYlm(0)[2] == Approx(0));

  CHECK(ct.getHessYlm(0)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(0)(0, 1) == Approx(0));
  CHECK(ct.getHessYlm(0)(0, 2) == Approx(0));
  CHECK(ct.getHessYlm(0)(1, 0) == Approx(0));
  CHECK(ct.getHessYlm(0)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(0)(1, 2) == Approx(0));
  CHECK(ct.getHessYlm(0)(2, 0) == Approx(0));
  CHECK(ct.getHessYlm(0)(2, 1) == Approx(0));
  CHECK(ct.getHessYlm(0)(2, 2) == Approx(0));


  CHECK(ct.getYlm(1) == Approx(0.635183265474));
  CHECK(ct.getGradYlm(1)[0] == Approx(0.488602511903));
  CHECK(ct.getGradYlm(1)[1] == Approx(0));
  CHECK(ct.getGradYlm(1)[2] == Approx(0));

  CHECK(ct.getHessYlm(1)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(1)(0, 1) == Approx(0));
  CHECK(ct.getHessYlm(1)(0, 2) == Approx(0));
  CHECK(ct.getHessYlm(1)(1, 0) == Approx(0));
  CHECK(ct.getHessYlm(1)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(1)(1, 2) == Approx(0));
  CHECK(ct.getHessYlm(1)(2, 0) == Approx(0));
  CHECK(ct.getHessYlm(1)(2, 1) == Approx(0));
  CHECK(ct.getHessYlm(1)(2, 2) == Approx(0));


  CHECK(ct.getYlm(15) == Approx(3.12417199136));
  CHECK(ct.getGradYlm(15)[0] == Approx(2.40320922412));
  CHECK(ct.getGradYlm(15)[1] == Approx(5.20695331894));
  CHECK(ct.getGradYlm(15)[2] == Approx(0));

  CHECK(ct.getHessYlm(15)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(15)(0, 1) == Approx(4.00534870687));
  CHECK(ct.getHessYlm(15)(0, 2) == Approx(0));
  CHECK(ct.getHessYlm(15)(1, 0) == Approx(4.00534870687));
  CHECK(ct.getHessYlm(15)(1, 1) == Approx(4.33912776578));
  CHECK(ct.getHessYlm(15)(1, 2) == Approx(0));
  CHECK(ct.getHessYlm(15)(2, 0) == Approx(0));
  CHECK(ct.getHessYlm(15)(2, 1) == Approx(0));
  CHECK(ct.getHessYlm(15)(2, 2) == Approx(0));


  CHECK(ct.getYlm(32) == Approx(-5.07677948596));
  CHECK(ct.getGradYlm(32)[0] == Approx(-7.81042997841));
  CHECK(ct.getGradYlm(32)[1] == Approx(-4.23064957164));
  CHECK(ct.getGradYlm(32)[2] == Approx(10.1535589719));

  CHECK(ct.getHessYlm(32)(0, 0) == Approx(-6.00802306031));
  CHECK(ct.getHessYlm(32)(0, 1) == Approx(-6.50869164867));
  CHECK(ct.getHessYlm(32)(0, 2) == Approx(15.6208599568));
  CHECK(ct.getHessYlm(32)(1, 0) == Approx(-6.50869164867));
  CHECK(ct.getHessYlm(32)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(32)(1, 2) == Approx(8.46129914327));
  CHECK(ct.getHessYlm(32)(2, 0) == Approx(15.6208599568));
  CHECK(ct.getHessYlm(32)(2, 1) == Approx(8.46129914327));
  CHECK(ct.getHessYlm(32)(2, 2) == Approx(0));


  CHECK(ct.getYlm(52) == Approx(-1.44809247474));
  CHECK(ct.getGradYlm(52)[0] == Approx(-1.11391728826));
  CHECK(ct.getGradYlm(52)[1] == Approx(-1.20674372895));
  CHECK(ct.getGradYlm(52)[2] == Approx(8.68855484841));

  CHECK(ct.getHessYlm(52)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(52)(0, 1) == Approx(-0.928264406882));
  CHECK(ct.getHessYlm(52)(0, 2) == Approx(6.68350372955));
  CHECK(ct.getHessYlm(52)(1, 0) == Approx(-0.928264406882));
  CHECK(ct.getHessYlm(52)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(52)(1, 2) == Approx(7.24046237368));
  CHECK(ct.getHessYlm(52)(2, 0) == Approx(6.68350372955));
  CHECK(ct.getHessYlm(52)(2, 1) == Approx(7.24046237368));
  CHECK(ct.getHessYlm(52)(2, 2) == Approx(-34.7542193937));


  CHECK(ct.getYlm(71) == Approx(-17.3423920977));
  CHECK(ct.getGradYlm(71)[0] == Approx(-53.3612064546));
  CHECK(ct.getGradYlm(71)[1] == Approx(-14.4519934148));
  CHECK(ct.getGradYlm(71)[2] == Approx(34.6847841955));

  CHECK(ct.getHessYlm(71)(0, 0) == Approx(-123.141245664));
  CHECK(ct.getHessYlm(71)(0, 1) == Approx(-44.4676720455));
  CHECK(ct.getHessYlm(71)(0, 2) == Approx(106.722412909));
  CHECK(ct.getHessYlm(71)(1, 0) == Approx(-44.4676720455));
  CHECK(ct.getHessYlm(71)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(71)(1, 2) == Approx(28.9039868295));
  CHECK(ct.getHessYlm(71)(2, 0) == Approx(106.722412909));
  CHECK(ct.getHessYlm(71)(2, 1) == Approx(28.9039868295));
  CHECK(ct.getHessYlm(71)(2, 2) == Approx(0));
}

TEST_CASE("Cartesian Tensor evaluateWithThirdDeriv subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateWithThirdDeriv(pt);

  //for (int i = 0; i < 35; i++) {
  //  std::cout << "XYZ = " << i << " " << ct.getYlm(i) << " " << ct.getGGGYlm(i) <<  std::endl;
  //}


  CHECK(ct.getYlm(0) == Approx(0.282094791774));


  CHECK(ct.getYlm(27) == Approx(-0.363846924275));
  CHECK(ct.getGradYlm(27)[0] == Approx(-0.279882249443));
  CHECK(ct.getGradYlm(27)[1] == Approx(0));
  CHECK(ct.getGradYlm(27)[2] == Approx(2.18308154565));

  CHECK(ct.getHessYlm(27)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(27)(0, 1) == Approx(0));
  CHECK(ct.getHessYlm(27)(0, 2) == Approx(1.67929349666));
  CHECK(ct.getHessYlm(27)(1, 0) == Approx(0));
  CHECK(ct.getHessYlm(27)(1, 1) == Approx(0));
  CHECK(ct.getHessYlm(27)(1, 2) == Approx(0));
  CHECK(ct.getHessYlm(27)(2, 0) == Approx(1.67929349666));
  CHECK(ct.getHessYlm(27)(2, 1) == Approx(0));
  CHECK(ct.getHessYlm(27)(2, 2) == Approx(-8.73232618261));


  CHECK(ct.getGGGYlm(27)[0](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[0](2, 2) == Approx(-6.71717398662));
  CHECK(ct.getGGGYlm(27)[1](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[1](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](0, 2) == Approx(-6.71717398662));
  CHECK(ct.getGGGYlm(27)[2](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](2, 0) == Approx(-6.71717398662));
  CHECK(ct.getGGGYlm(27)[2](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(27)[2](2, 2) == Approx(17.4646523652));


  CHECK(ct.getYlm(62) == Approx(-4.19700340252));
  CHECK(ct.getGradYlm(62)[0] == Approx(0));
  CHECK(ct.getGradYlm(62)[1] == Approx(-17.4875141772));
  CHECK(ct.getGradYlm(62)[2] == Approx(8.39400680505));

  CHECK(ct.getHessYlm(62)(0, 0) == Approx(0));
  CHECK(ct.getHessYlm(62)(0, 1) == Approx(0));
  CHECK(ct.getHessYlm(62)(0, 2) == Approx(0));
  CHECK(ct.getHessYlm(62)(1, 0) == Approx(0));
  CHECK(ct.getHessYlm(62)(1, 1) == Approx(-58.291713924));
  CHECK(ct.getHessYlm(62)(1, 2) == Approx(34.9750283544));
  CHECK(ct.getHessYlm(62)(2, 0) == Approx(0));
  CHECK(ct.getHessYlm(62)(2, 1) == Approx(34.9750283544));
  CHECK(ct.getHessYlm(62)(2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[0](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](1, 1) == Approx(-145.72928481));
  CHECK(ct.getGGGYlm(62)[1](1, 2) == Approx(116.583427848));
  CHECK(ct.getGGGYlm(62)[1](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[1](2, 1) == Approx(116.583427848));
  CHECK(ct.getGGGYlm(62)[1](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](1, 1) == Approx(116.583427848));
  CHECK(ct.getGGGYlm(62)[2](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(62)[2](2, 2) == Approx(0));
}


TEST_CASE("Cartesian Tensor evaluateThirdDerivOnly subset", "[numerics]")
{
  CartesianTensor<double, TinyVector<double, 3>> ct(6);

  TinyVector<double, 3> pt(1.3, 1.2, -0.5);
  ct.evaluateThirdDerivOnly(pt);

  CHECK(ct.getGGGYlm(14)[0](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](0, 2) == Approx(3.33779058906));
  CHECK(ct.getGGGYlm(14)[0](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](2, 0) == Approx(3.33779058906));
  CHECK(ct.getGGGYlm(14)[0](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[0](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[1](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](0, 0) == Approx(3.33779058906));
  CHECK(ct.getGGGYlm(14)[2](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(14)[2](2, 2) == Approx(0));


  CHECK(ct.getGGGYlm(33)[0](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[0](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(33)[0](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(33)[0](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[0](1, 1) == Approx(-5.00668588359));
  CHECK(ct.getGGGYlm(33)[0](1, 2) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[0](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[0](2, 1) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[0](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(33)[1](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[1](0, 1) == Approx(-5.00668588359));
  CHECK(ct.getGGGYlm(33)[1](0, 2) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[1](1, 0) == Approx(-5.00668588359));
  CHECK(ct.getGGGYlm(33)[1](1, 1) == Approx(0));
  CHECK(ct.getGGGYlm(33)[1](1, 2) == Approx(13.0173832973));
  CHECK(ct.getGGGYlm(33)[1](2, 0) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[1](2, 1) == Approx(13.0173832973));
  CHECK(ct.getGGGYlm(33)[1](2, 2) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](0, 1) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[2](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](1, 0) == Approx(12.0160461206));
  CHECK(ct.getGGGYlm(33)[2](1, 1) == Approx(13.0173832973));
  CHECK(ct.getGGGYlm(33)[2](1, 2) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](2, 1) == Approx(0));
  CHECK(ct.getGGGYlm(33)[2](2, 2) == Approx(0));


  CHECK(ct.getGGGYlm(80)[0](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(80)[0](0, 1) == Approx(0));
  CHECK(ct.getGGGYlm(80)[0](0, 2) == Approx(0));
  CHECK(ct.getGGGYlm(80)[0](1, 0) == Approx(0));
  CHECK(ct.getGGGYlm(80)[0](1, 1) == Approx(27.8256449422));
  CHECK(ct.getGGGYlm(80)[0](1, 2) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[0](2, 0) == Approx(0));
  CHECK(ct.getGGGYlm(80)[0](2, 1) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[0](2, 2) == Approx(53.425238289));
  CHECK(ct.getGGGYlm(80)[1](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(80)[1](0, 1) == Approx(27.8256449422));
  CHECK(ct.getGGGYlm(80)[1](0, 2) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[1](1, 0) == Approx(27.8256449422));
  CHECK(ct.getGGGYlm(80)[1](1, 1) == Approx(30.1444486874));
  CHECK(ct.getGGGYlm(80)[1](1, 2) == Approx(-144.6933537));
  CHECK(ct.getGGGYlm(80)[1](2, 0) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[1](2, 1) == Approx(-144.6933537));
  CHECK(ct.getGGGYlm(80)[1](2, 2) == Approx(173.632024439));
  CHECK(ct.getGGGYlm(80)[2](0, 0) == Approx(0));
  CHECK(ct.getGGGYlm(80)[2](0, 1) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[2](0, 2) == Approx(53.425238289));
  CHECK(ct.getGGGYlm(80)[2](1, 0) == Approx(-66.7815478613));
  CHECK(ct.getGGGYlm(80)[2](1, 1) == Approx(-144.6933537));
  CHECK(ct.getGGGYlm(80)[2](1, 2) == Approx(173.632024439));
  CHECK(ct.getGGGYlm(80)[2](2, 0) == Approx(53.425238289));
  CHECK(ct.getGGGYlm(80)[2](2, 1) == Approx(173.632024439));
  CHECK(ct.getGGGYlm(80)[2](2, 2) == Approx(0));
}
} // namespace qmcplusplus

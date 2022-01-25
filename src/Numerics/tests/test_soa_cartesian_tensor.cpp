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


// This test is technically in the wrong location.  To be testing code in QMCWaveFunctions, it should
// be in QMCWaveFunctions/test.   However, it seems better to keep the two CartesianTensor tests together.

#include "catch.hpp"
#include "Message/Communicate.h"
#include "QMCWaveFunctions/LCAO/SoaCartesianTensor.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
using real_type = OHMMS_PRECISION;


// Use gen_gto.py to generate the checks

TEST_CASE("SoA Cartesian Tensor", "[numerics]")
{
  SoaCartesianTensor<double> ct(6);

  double x = 1.3;
  double y = 1.2;
  double z = -0.5;
  ct.evaluateV(x, y, z);

  double* XYZ = ct.cXYZ.data(0);

  REQUIRE(XYZ[0] == Approx(0.282094791774));
  REQUIRE(XYZ[1] == Approx(0.635183265474));
  REQUIRE(XYZ[2] == Approx(0.586323014284));
  REQUIRE(XYZ[3] == Approx(-0.244301255951));
  REQUIRE(XYZ[4] == Approx(1.06602349055));
  REQUIRE(XYZ[5] == Approx(0.908327707927));
  REQUIRE(XYZ[6] == Approx(0.157695782626));
  REQUIRE(XYZ[7] == Approx(1.70437555172));
  REQUIRE(XYZ[8] == Approx(-0.710156479885));
  REQUIRE(XYZ[9] == Approx(-0.655529058355));
  REQUIRE(XYZ[10] == Approx(1.6397368054));
  REQUIRE(XYZ[11] == Approx(1.28969740543));
  REQUIRE(XYZ[12] == Approx(-0.0932940831475));
  REQUIRE(XYZ[13] == Approx(3.38451965731));
  REQUIRE(XYZ[14] == Approx(-1.41021652388));
  REQUIRE(XYZ[15] == Approx(3.12417199136));
  REQUIRE(XYZ[16] == Approx(-1.20160461206));
  REQUIRE(XYZ[17] == Approx(0.542390970723));
  REQUIRE(XYZ[18] == Approx(0.500668588359));
  REQUIRE(XYZ[19] == Approx(-2.25467692526));
  REQUIRE(XYZ[20] == Approx(2.41707280436));
  REQUIRE(XYZ[21] == Approx(1.75485528067));
  REQUIRE(XYZ[22] == Approx(0.0528927734576));
  REQUIRE(XYZ[23] == Approx(5.90305249944));
  REQUIRE(XYZ[24] == Approx(-2.4596052081));
  REQUIRE(XYZ[25] == Approx(5.02981988118));
  REQUIRE(XYZ[26] == Approx(-1.93454610815));
  REQUIRE(XYZ[27] == Approx(-0.363846924275));
  REQUIRE(XYZ[28] == Approx(-0.335858699331));
  REQUIRE(XYZ[29] == Approx(7.03459200681));
  REQUIRE(XYZ[30] == Approx(1.22128333452));
  REQUIRE(XYZ[31] == Approx(1.04062011935));
  REQUIRE(XYZ[32] == Approx(-5.07677948596));
  REQUIRE(XYZ[33] == Approx(-4.68625798704));
  REQUIRE(XYZ[34] == Approx(1.9526074946));
  REQUIRE(XYZ[35] == Approx(3.47382688598));
  REQUIRE(XYZ[36] == Approx(2.32807861094));
  REQUIRE(XYZ[37] == Approx(-0.0292375806134));
  REQUIRE(XYZ[38] == Approx(9.61982829963));
  REQUIRE(XYZ[39] == Approx(-4.00826179151));
  REQUIRE(XYZ[40] == Approx(7.56625548555));
  REQUIRE(XYZ[41] == Approx(-2.91009826367));
  REQUIRE(XYZ[42] == Approx(0.228053128784));
  REQUIRE(XYZ[43] == Approx(0.210510580416));
  REQUIRE(XYZ[44] == Approx(13.5641819555));
  REQUIRE(XYZ[45] == Approx(2.35489270061));
  REQUIRE(XYZ[46] == Approx(12.5207833436));
  REQUIRE(XYZ[47] == Approx(1.85218688514));
  REQUIRE(XYZ[48] == Approx(-0.905727961775));
  REQUIRE(XYZ[49] == Approx(-0.771744535477));
  REQUIRE(XYZ[50] == Approx(-9.78910512921));
  REQUIRE(XYZ[51] == Approx(-8.34101265448));
  REQUIRE(XYZ[52] == Approx(-1.44809247474));
  REQUIRE(XYZ[53] == Approx(-11.6655511199));
  REQUIRE(XYZ[54] == Approx(4.86064629996));
  REQUIRE(XYZ[55] == Approx(4.48675043074));
  REQUIRE(XYZ[56] == Approx(4.90938236205));
  REQUIRE(XYZ[57] == Approx(3.03706593382));
  REQUIRE(XYZ[58] == Approx(0.0158923005669));
  REQUIRE(XYZ[59] == Approx(15.0300731514));
  REQUIRE(XYZ[60] == Approx(-6.26253047974));
  REQUIRE(XYZ[61] == Approx(10.9122088466));
  REQUIRE(XYZ[62] == Approx(-4.19700340252));
  REQUIRE(XYZ[63] == Approx(-0.137042874894));
  REQUIRE(XYZ[64] == Approx(-0.126501115286));
  REQUIRE(XYZ[65] == Approx(24.0303233904));
  REQUIRE(XYZ[66] == Approx(4.17193114417));
  REQUIRE(XYZ[67] == Approx(20.4755418238));
  REQUIRE(XYZ[68] == Approx(3.0289263053));
  REQUIRE(XYZ[69] == Approx(0.61714957754));
  REQUIRE(XYZ[70] == Approx(0.525855261336));
  REQUIRE(XYZ[71] == Approx(-17.3423920977));
  REQUIRE(XYZ[72] == Approx(-13.6402610582));
  REQUIRE(XYZ[73] == Approx(0.986708699234));
  REQUIRE(XYZ[74] == Approx(26.2459034569));
  REQUIRE(XYZ[75] == Approx(-1.89857519219));
  REQUIRE(XYZ[76] == Approx(-1.49328080661));
  REQUIRE(XYZ[77] == Approx(-24.4531767752));
  REQUIRE(XYZ[78] == Approx(10.1888236563));
  REQUIRE(XYZ[79] == Approx(-22.5721631771));
  REQUIRE(XYZ[80] == Approx(8.68160122197));
  REQUIRE(XYZ[81] == Approx(-3.91877832936));
  REQUIRE(XYZ[82] == Approx(-3.61733384249));
  REQUIRE(XYZ[83] == Approx(12.1418905657));
}

TEST_CASE("SoA Cartesian Tensor evaluateVGL subset", "[numerics]")
{
  SoaCartesianTensor<double> ct(6);

  double x = 1.3;
  double y = 1.2;
  double z = -0.5;
  ct.evaluateVGL(x, y, z);

  double* XYZ = ct.cXYZ.data(0);
  double* gr0 = ct.cXYZ.data(1);
  double* gr1 = ct.cXYZ.data(2);
  double* gr2 = ct.cXYZ.data(3);
  double* lap = ct.cXYZ.data(4);

  REQUIRE(XYZ[0] == Approx(0.282094791774));
  REQUIRE(gr0[0] == Approx(0));
  REQUIRE(gr1[0] == Approx(0));
  REQUIRE(gr2[0] == Approx(0));
  REQUIRE(lap[0] == Approx(0));

  REQUIRE(XYZ[1] == Approx(0.635183265474));
  REQUIRE(gr0[1] == Approx(0.488602511903));
  REQUIRE(gr1[1] == Approx(0));
  REQUIRE(gr2[1] == Approx(0));
  REQUIRE(lap[1] == Approx(0));

  REQUIRE(XYZ[8] == Approx(-0.710156479885));
  REQUIRE(gr0[8] == Approx(-0.546274215296));
  REQUIRE(gr1[8] == Approx(0));
  REQUIRE(gr2[8] == Approx(1.42031295977));
  REQUIRE(lap[8] == Approx(0));

  REQUIRE(XYZ[23] == Approx(5.90305249944));
  REQUIRE(gr0[23] == Approx(13.6224288449));
  REQUIRE(gr1[23] == Approx(4.9192104162));
  REQUIRE(gr2[23] == Approx(0));
  REQUIRE(lap[23] == Approx(20.9575828383));

  REQUIRE(XYZ[34] == Approx(1.9526074946));
  REQUIRE(gr0[34] == Approx(1.50200576508));
  REQUIRE(gr1[34] == Approx(1.62717291217));
  REQUIRE(gr2[34] == Approx(-7.81042997841));
  REQUIRE(lap[34] == Approx(15.6208599568));

  REQUIRE(XYZ[50] == Approx(-9.78910512921));
  REQUIRE(gr0[50] == Approx(-22.5902426059));
  REQUIRE(gr1[50] == Approx(-8.15758760768));
  REQUIRE(gr2[50] == Approx(19.5782102584));
  REQUIRE(lap[50] == Approx(-34.7542193937));

  REQUIRE(XYZ[83] == Approx(12.1418905657));
  REQUIRE(gr0[83] == Approx(18.6798316395));
  REQUIRE(gr1[83] == Approx(20.2364842761));
  REQUIRE(gr2[83] == Approx(-48.5675622627));
  REQUIRE(lap[83] == Approx(128.367962683));
}

TEST_CASE("SoA Cartesian Tensor evaluateVGH subset", "[numerics]")
{
  SoaCartesianTensor<double> ct(6);

  double x = 1.3;
  double y = 1.2;
  double z = -0.5;
  ct.evaluateVGH(x, y, z);

  double* XYZ = ct.cXYZ.data(0);
  double* gr0 = ct.cXYZ.data(1);
  double* gr1 = ct.cXYZ.data(2);
  double* gr2 = ct.cXYZ.data(3);
  double* h00 = ct.cXYZ.data(4);
  double* h01 = ct.cXYZ.data(5);
  double* h02 = ct.cXYZ.data(6);
  double* h11 = ct.cXYZ.data(7);
  double* h12 = ct.cXYZ.data(8);
  double* h22 = ct.cXYZ.data(9);

  REQUIRE(XYZ[0] == Approx(0.282094791774));
  REQUIRE(gr0[0] == Approx(0));
  REQUIRE(gr1[0] == Approx(0));
  REQUIRE(gr2[0] == Approx(0));

  REQUIRE(h00[0] == Approx(0));
  REQUIRE(h01[0] == Approx(0));
  REQUIRE(h02[0] == Approx(0));
  REQUIRE(h11[0] == Approx(0));
  REQUIRE(h12[0] == Approx(0));
  REQUIRE(h22[0] == Approx(0));


  REQUIRE(XYZ[1] == Approx(0.635183265474));
  REQUIRE(gr0[1] == Approx(0.488602511903));
  REQUIRE(gr1[1] == Approx(0));
  REQUIRE(gr2[1] == Approx(0));

  REQUIRE(h00[1] == Approx(0));
  REQUIRE(h01[1] == Approx(0));
  REQUIRE(h02[1] == Approx(0));
  REQUIRE(h11[1] == Approx(0));
  REQUIRE(h12[1] == Approx(0));
  REQUIRE(h22[1] == Approx(0));


  REQUIRE(XYZ[15] == Approx(3.12417199136));
  REQUIRE(gr0[15] == Approx(2.40320922412));
  REQUIRE(gr1[15] == Approx(5.20695331894));
  REQUIRE(gr2[15] == Approx(0));

  REQUIRE(h00[15] == Approx(0));
  REQUIRE(h01[15] == Approx(4.00534870687));
  REQUIRE(h02[15] == Approx(0));
  REQUIRE(h11[15] == Approx(4.33912776578));
  REQUIRE(h12[15] == Approx(0));
  REQUIRE(h22[15] == Approx(0));


  REQUIRE(XYZ[32] == Approx(-5.07677948596));
  REQUIRE(gr0[32] == Approx(-7.81042997841));
  REQUIRE(gr1[32] == Approx(-4.23064957164));
  REQUIRE(gr2[32] == Approx(10.1535589719));

  REQUIRE(h00[32] == Approx(-6.00802306031));
  REQUIRE(h01[32] == Approx(-6.50869164867));
  REQUIRE(h02[32] == Approx(15.6208599568));
  REQUIRE(h11[32] == Approx(0));
  REQUIRE(h12[32] == Approx(8.46129914327));
  REQUIRE(h22[32] == Approx(0));


  REQUIRE(XYZ[52] == Approx(-1.44809247474));
  REQUIRE(gr0[52] == Approx(-1.11391728826));
  REQUIRE(gr1[52] == Approx(-1.20674372895));
  REQUIRE(gr2[52] == Approx(8.68855484841));

  REQUIRE(h00[52] == Approx(0));
  REQUIRE(h01[52] == Approx(-0.928264406882));
  REQUIRE(h02[52] == Approx(6.68350372955));
  REQUIRE(h11[52] == Approx(0));
  REQUIRE(h12[52] == Approx(7.24046237368));
  REQUIRE(h22[52] == Approx(-34.7542193937));


  REQUIRE(XYZ[71] == Approx(-17.3423920977));
  REQUIRE(gr0[71] == Approx(-53.3612064546));
  REQUIRE(gr1[71] == Approx(-14.4519934148));
  REQUIRE(gr2[71] == Approx(34.6847841955));

  REQUIRE(h00[71] == Approx(-123.141245664));
  REQUIRE(h01[71] == Approx(-44.4676720455));
  REQUIRE(h02[71] == Approx(106.722412909));
  REQUIRE(h11[71] == Approx(0));
  REQUIRE(h12[71] == Approx(28.9039868295));
  REQUIRE(h22[71] == Approx(0));
}

} // namespace qmcplusplus

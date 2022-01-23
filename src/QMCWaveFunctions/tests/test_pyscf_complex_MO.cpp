//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Chandler Bennett, bennettcc@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/*  This is a deterministic unit test that verifies the molecular orbital values from QMCPACK
 *  match those from PySCF in the specific case of a complex periodic gaussian basis. 
 *  For off gamma carbon diamond, the real and imaginary parts of all molecular orbitals along a 
 *  real-space path that spans regions both inside and outside of the cell are required to agree with
 *  the reference values obtained from PySCF. The reference values were generated from Carbon1x1x1-tw1_gen_mos.py.   
*/

#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "Numerics/OneDimGridBase.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Numerics/GaussianBasisSet.h"
#include "QMCWaveFunctions/LCAO/LCAOrbitalBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
void test_C_diamond()
{
  std::ostringstream section_name;
  section_name << "Carbon diamond off gamma unit test: ";

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    Libxml2Document doc;
    bool okay = doc.parse("C_diamond-twist-third.structure.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    ParticleSet::ParticleLayout lattice;
    // BCC H
    lattice.R(0, 0) = 3.37316115;
    lattice.R(0, 1) = 3.37316115;
    lattice.R(0, 2) = 0.0;
    lattice.R(1, 0) = 0.0;
    lattice.R(1, 1) = 3.37316115;
    lattice.R(1, 2) = 3.37316115;
    lattice.R(2, 0) = 3.37316115;
    lattice.R(2, 1) = 0.0;
    lattice.R(2, 2) = 3.37316115;
    lattice.reset();

    const SimulationCell simulation_cell(lattice);
    ParticleSet ions(simulation_cell);
    XMLParticleParser parse_ions(ions);
    OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
    REQUIRE(particleset_ion.size() == 1);
    parse_ions.put(particleset_ion[0]);

    REQUIRE(ions.groups() == 1);
    REQUIRE(ions.R.size() == 2);
    ions.update();

    ParticleSet elec(simulation_cell);
    XMLParticleParser parse_elec(elec);
    OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
    REQUIRE(particleset_elec.size() == 1);
    parse_elec.put(particleset_elec[0]);

    REQUIRE(elec.groups() == 2);
    REQUIRE(elec.R.size() == 8);

    elec.R = 0.0;

    elec.addTable(ions);
    elec.update();

    Libxml2Document doc2;

    okay = doc2.parse("C_diamond-twist-third.wfj.xml");
    REQUIRE(okay);
    xmlNodePtr root2 = doc2.getRoot();

    WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
    particle_set_map["e"]    = &elec;
    particle_set_map["ion0"] = &ions;

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    auto& bb = bf.createSPOSetBuilder(MO_base[0]);

    OhmmsXPathObject slater_base("//sposet", doc2.getXPathContext());
    SPOSet* sposet = bb.createSPOSet(slater_base[0]);

    SPOSet::ValueVector values;
    values.resize(26);

    // BEGIN generated C++ input from Carbon1x1x1-tw1_gen_mos.py (pyscf version 1.6.2) on 2019-11-19 15:08:42.652893

    //Move electron 0 to position [[-10. -10. -10.]] a.u.:
    elec.R[0] = {-10.0, -10.0, -10.0};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[-10. -10. -10.]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(-0.060195105580765275));
    REQUIRE(std::imag(values[0]) == Approx(-0.011833831357235357));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.0486805973426202));
    REQUIRE(std::imag(values[1]) == Approx(0.02535990721099494));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(2.0679269928380872e-13));
    REQUIRE(std::imag(values[2]) == Approx(-7.649005669303763e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-8.497271896529268e-13));
    REQUIRE(std::imag(values[3]) == Approx(1.7299795493710973e-13));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.014722495006984714));
    REQUIRE(std::imag(values[4]) == Approx(-0.07164030105855403));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(2.5059967122265113e-13));
    REQUIRE(std::imag(values[5]) == Approx(2.6819557859270038e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(-1.4149372320505777e-13));
    REQUIRE(std::imag(values[6]) == Approx(-6.166027974666144e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.1625931972310664));
    REQUIRE(std::imag(values[7]) == Approx(-0.06588438085509145));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(0.14211538487895298));
    REQUIRE(std::imag(values[8]) == Approx(-0.10262731977374147));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(6.853515151228338e-13));
    REQUIRE(std::imag(values[9]) == Approx(-6.336955699765e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(-2.242625356252415e-12));
    REQUIRE(std::imag(values[10]) == Approx(1.2973476459787747e-12));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(0.16369199334396695));
    REQUIRE(std::imag(values[11]) == Approx(-0.10523147038857521));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(-6.491920716278332e-14));
    REQUIRE(std::imag(values[12]) == Approx(-2.5135275805165893e-14));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(2.477731561589818e-13));
    REQUIRE(std::imag(values[13]) == Approx(-6.782768791069354e-14));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-7.389267149549019e-13));
    REQUIRE(std::imag(values[14]) == Approx(1.8454929424926995e-12));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(1.551052689063858e-12));
    REQUIRE(std::imag(values[15]) == Approx(-5.080043156968638e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(0.02643126908254878));
    REQUIRE(std::imag(values[16]) == Approx(-0.10382334171056015));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(-0.008572223551589144));
    REQUIRE(std::imag(values[17]) == Approx(0.1360121319699977));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-8.83004606933017e-14));
    REQUIRE(std::imag(values[18]) == Approx(6.17399794483608e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(-1.2893852652240422e-13));
    REQUIRE(std::imag(values[19]) == Approx(3.7305644684515455e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.02873591847713064));
    REQUIRE(std::imag(values[20]) == Approx(-0.153030888627405));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-5.410008378781134e-13));
    REQUIRE(std::imag(values[21]) == Approx(-1.9873859502526956e-14));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(1.3256132302963402e-12));
    REQUIRE(std::imag(values[22]) == Approx(9.89605099255277e-13));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(1.1945652800271454e-13));
    REQUIRE(std::imag(values[23]) == Approx(-2.047390035286864e-13));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(6.814895869844398e-13));
    REQUIRE(std::imag(values[24]) == Approx(-9.42856903662901e-14));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(-0.015289434969691054));
    REQUIRE(std::imag(values[25]) == Approx(-0.03615194381469279));

    //Move electron 0 to position [[-6.666667 -6.666667 -6.666667]] a.u.:
    elec.R[0] = {-6.666667, -6.666667, -6.666667};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[-6.666667 -6.666667 -6.666667]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(0.0956158275805242));
    REQUIRE(std::imag(values[0]) == Approx(0.03287266037323013));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.11506570134621005));
    REQUIRE(std::imag(values[1]) == Approx(-0.0900472787152517));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(3.6262524017215473e-13));
    REQUIRE(std::imag(values[2]) == Approx(-8.165796056348355e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-1.499448351923057e-12));
    REQUIRE(std::imag(values[3]) == Approx(1.207191898465219e-13));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.09133965128853623));
    REQUIRE(std::imag(values[4]) == Approx(-0.07359368062940903));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(-2.1379426007743085e-12));
    REQUIRE(std::imag(values[5]) == Approx(2.7089146897775157e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(6.268878374154817e-13));
    REQUIRE(std::imag(values[6]) == Approx(-7.328572427840399e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.1450364317288277));
    REQUIRE(std::imag(values[7]) == Approx(0.1499844328343919));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(0.1533421330529442));
    REQUIRE(std::imag(values[8]) == Approx(0.09647361404443654));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(-8.187964194921162e-13));
    REQUIRE(std::imag(values[9]) == Approx(4.4221050437064805e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(4.14457530803798e-12));
    REQUIRE(std::imag(values[10]) == Approx(-3.1273247862003916e-13));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(-0.23130832819832292));
    REQUIRE(std::imag(values[11]) == Approx(-0.000920632077534666));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(2.3089689879702785e-13));
    REQUIRE(std::imag(values[12]) == Approx(8.097805194661028e-13));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(-1.2650332172739602e-12));
    REQUIRE(std::imag(values[13]) == Approx(-1.7846488183430846e-13));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-1.4239386578211e-12));
    REQUIRE(std::imag(values[14]) == Approx(1.6467786336155512e-12));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(3.445370725198817e-12));
    REQUIRE(std::imag(values[15]) == Approx(-2.88321103094349e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(0.06850124380096241));
    REQUIRE(std::imag(values[16]) == Approx(-0.004634947034247484));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(0.35713027737573216));
    REQUIRE(std::imag(values[17]) == Approx(0.09047102606218091));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(2.2383223715827679e-13));
    REQUIRE(std::imag(values[18]) == Approx(9.967200661795276e-14));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(9.707838699186186e-12));
    REQUIRE(std::imag(values[19]) == Approx(5.315903969755244e-13));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.24626463443252364));
    REQUIRE(std::imag(values[20]) == Approx(-0.07802697773073676));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-6.978536668841118e-13));
    REQUIRE(std::imag(values[21]) == Approx(2.91436579306753e-13));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(1.3880337858309117e-12));
    REQUIRE(std::imag(values[22]) == Approx(-2.2989839200496397e-12));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(-1.4148543446437055e-12));
    REQUIRE(std::imag(values[23]) == Approx(3.6835638711982286e-13));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(-5.500253025100187e-13));
    REQUIRE(std::imag(values[24]) == Approx(-1.235497815173411e-12));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(0.7850948237148717));
    REQUIRE(std::imag(values[25]) == Approx(-0.6716776823047774));

    //Move electron 0 to position [[-3.333334 -3.333334 -3.333334]] a.u.:
    elec.R[0] = {-3.333334, -3.333334, -3.333334};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[-3.333334 -3.333334 -3.333334]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(-0.06005230809177282));
    REQUIRE(std::imag(values[0]) == Approx(-0.005290783996832645));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.04593960387952771));
    REQUIRE(std::imag(values[1]) == Approx(0.03337244273082512));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(2.1968109839728915e-13));
    REQUIRE(std::imag(values[2]) == Approx(-8.417851376747011e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-8.999305207334697e-13));
    REQUIRE(std::imag(values[3]) == Approx(2.0359793834210627e-13));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.018164912111658126));
    REQUIRE(std::imag(values[4]) == Approx(-0.07922559625400716));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(1.4349914486982305e-13));
    REQUIRE(std::imag(values[5]) == Approx(2.9226239484116607e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(-1.1288736058805674e-13));
    REQUIRE(std::imag(values[6]) == Approx(-6.832621491132385e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.1651619894114805));
    REQUIRE(std::imag(values[7]) == Approx(-0.06134476438049227));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(0.14028043757888803));
    REQUIRE(std::imag(values[8]) == Approx(-0.10727459702035165));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(7.288215170092265e-13));
    REQUIRE(std::imag(values[9]) == Approx(-6.500147642486826e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(-2.3690801924617754e-12));
    REQUIRE(std::imag(values[10]) == Approx(1.2809042021000265e-12));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(0.17134778976483186));
    REQUIRE(std::imag(values[11]) == Approx(-0.09873043450786173));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(-4.721136675333712e-14));
    REQUIRE(std::imag(values[12]) == Approx(1.479762495375543e-14));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(1.827221642785641e-13));
    REQUIRE(std::imag(values[13]) == Approx(-3.08052194664315e-14));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-4.352638042034029e-13));
    REQUIRE(std::imag(values[14]) == Approx(1.8446439037660834e-12));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(6.46854097962043e-13));
    REQUIRE(std::imag(values[15]) == Approx(-4.927302489657481e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(0.007153805199742103));
    REQUIRE(std::imag(values[16]) == Approx(-0.09671197376753325));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(0.0127598393467023));
    REQUIRE(std::imag(values[17]) == Approx(0.16334397744406476));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-6.968514298788527e-14));
    REQUIRE(std::imag(values[18]) == Approx(6.594685735382753e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(4.700545509464493e-13));
    REQUIRE(std::imag(values[19]) == Approx(4.414564200230628e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.04159487544411769));
    REQUIRE(std::imag(values[20]) == Approx(-0.1385550610888857));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-5.91423394723593e-13));
    REQUIRE(std::imag(values[21]) == Approx(-4.260524213468368e-14));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(1.554485706631173e-12));
    REQUIRE(std::imag(values[22]) == Approx(1.0536316828100446e-12));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(1.2028225633590926e-13));
    REQUIRE(std::imag(values[23]) == Approx(-2.8023416921979087e-13));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(7.943472267278781e-13));
    REQUIRE(std::imag(values[24]) == Approx(-1.860109288800894e-13));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(-0.02037392105518354));
    REQUIRE(std::imag(values[25]) == Approx(-0.057434079996872056));

    //Move electron 0 to position [[-9.99999999e-07 -9.99999999e-07 -9.99999999e-07]] a.u.:
    elec.R[0] = {-9.999999992515995e-07, -9.999999992515995e-07, -9.999999992515995e-07};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[-9.99999999e-07 -9.99999999e-07 -9.99999999e-07]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(0.10167888698741122));
    REQUIRE(std::imag(values[0]) == Approx(-0.005941230707015811));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.06770487558800753));
    REQUIRE(std::imag(values[1]) == Approx(0.008637814352396151));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(3.60922270259699e-13));
    REQUIRE(std::imag(values[2]) == Approx(-5.756254305688054e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-1.5147241100318388e-12));
    REQUIRE(std::imag(values[3]) == Approx(4.370202116763758e-14));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.1285247939510426));
    REQUIRE(std::imag(values[4]) == Approx(-0.010293562806552263));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(-3.1382097441498023e-12));
    REQUIRE(std::imag(values[5]) == Approx(5.680356358934007e-13));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(8.714109078422965e-13));
    REQUIRE(std::imag(values[6]) == Approx(-1.8866876919070565e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.0787451415811546));
    REQUIRE(std::imag(values[7]) == Approx(0.00141109567072648));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(-0.0387351183536941));
    REQUIRE(std::imag(values[8]) == Approx(0.012263153973258889));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(-4.3274845513970603e-14));
    REQUIRE(std::imag(values[9]) == Approx(-9.283134158055235e-14));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(2.4316486325254445e-13));
    REQUIRE(std::imag(values[10]) == Approx(1.6584216636762436e-13));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(-0.010875903895641862));
    REQUIRE(std::imag(values[11]) == Approx(0.04714572895763401));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(2.0172275308035739e-13));
    REQUIRE(std::imag(values[12]) == Approx(5.497217874676079e-13));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(-8.496293946065942e-13));
    REQUIRE(std::imag(values[13]) == Approx(-1.036705443495618e-13));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-2.263782260614978e-12));
    REQUIRE(std::imag(values[14]) == Approx(5.795082295856767e-13));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(6.377089828402728e-12));
    REQUIRE(std::imag(values[15]) == Approx(-6.716103368246984e-13));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(0.12847308968746957));
    REQUIRE(std::imag(values[16]) == Approx(0.02035282496639674));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(0.1737553721002816));
    REQUIRE(std::imag(values[17]) == Approx(-0.0376579821006186));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-1.4876446427583875e-13));
    REQUIRE(std::imag(values[18]) == Approx(-3.7704344851099533e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(4.382987028924654e-12));
    REQUIRE(std::imag(values[19]) == Approx(-2.5937429287598488e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.15399434198150766));
    REQUIRE(std::imag(values[20]) == Approx(-0.02568058890861615));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-1.869012757649718e-13));
    REQUIRE(std::imag(values[21]) == Approx(2.032548120769792e-13));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(-9.61637020079959e-13));
    REQUIRE(std::imag(values[22]) == Approx(-1.8784349076636863e-12));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(-1.0113611337362861e-12));
    REQUIRE(std::imag(values[23]) == Approx(1.0432885024476444e-12));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(-8.348460811998449e-13));
    REQUIRE(std::imag(values[24]) == Approx(-1.0209888494617845e-13));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(0.8037234583780705));
    REQUIRE(std::imag(values[25]) == Approx(-0.8579924081605256));

    //Move electron 0 to position [[3.333332 3.333332 3.333332]] a.u.:
    elec.R[0] = {3.3333320000000004, 3.3333320000000004, 3.3333320000000004};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[3.333332 3.333332 3.333332]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(-0.059806620261567495));
    REQUIRE(std::imag(values[0]) == Approx(0.0014605026493625755));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.042130656010843565));
    REQUIRE(std::imag(values[1]) == Approx(0.041154352141271694));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(2.299461077958894e-13));
    REQUIRE(std::imag(values[2]) == Approx(-9.228330396886522e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-9.411403950013916e-13));
    REQUIRE(std::imag(values[3]) == Approx(2.373896446652552e-13));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.022449104164198187));
    REQUIRE(std::imag(values[4]) == Approx(-0.08741397581582212));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(5.8607649950285695e-15));
    REQUIRE(std::imag(values[5]) == Approx(3.172440606911535e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(-7.517818458469839e-14));
    REQUIRE(std::imag(values[6]) == Approx(-7.526353676741086e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.16660016621038506));
    REQUIRE(std::imag(values[7]) == Approx(-0.05698359576596997));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(0.13740625808310925));
    REQUIRE(std::imag(values[8]) == Approx(-0.10946922834461466));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(7.572904391658913e-13));
    REQUIRE(std::imag(values[9]) == Approx(-6.606785411995795e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(-2.4508546356377402e-12));
    REQUIRE(std::imag(values[10]) == Approx(1.2663325239134368e-12));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(0.17681058034809968));
    REQUIRE(std::imag(values[11]) == Approx(-0.09275809690899624));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(-2.9056613057169995e-14));
    REQUIRE(std::imag(values[12]) == Approx(4.9557581445141706e-14));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(1.2658321585137323e-13));
    REQUIRE(std::imag(values[13]) == Approx(1.070843942851782e-14));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-8.521829595450246e-14));
    REQUIRE(std::imag(values[14]) == Approx(1.77439224452637e-12));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(-3.6150944959863987e-13));
    REQUIRE(std::imag(values[15]) == Approx(-4.5876243176725545e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(-0.014036906825284585));
    REQUIRE(std::imag(values[16]) == Approx(-0.08614468966547723));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(0.033216162124121054));
    REQUIRE(std::imag(values[17]) == Approx(0.18639424345387834));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-4.552953546266951e-14));
    REQUIRE(std::imag(values[18]) == Approx(6.769691296399057e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(1.0413475613943049e-12));
    REQUIRE(std::imag(values[19]) == Approx(4.9922947429389975e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.049807632374415926));
    REQUIRE(std::imag(values[20]) == Approx(-0.11666665209986808));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-6.182341843173099e-13));
    REQUIRE(std::imag(values[21]) == Approx(-6.314824100072268e-14));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(1.7345326427769266e-12));
    REQUIRE(std::imag(values[22]) == Approx(1.094659994231126e-12));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(1.290911743878159e-13));
    REQUIRE(std::imag(values[23]) == Approx(-3.4848512356417114e-13));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(8.886571862605712e-13));
    REQUIRE(std::imag(values[24]) == Approx(-2.685212860435276e-13));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(-0.02637454702199341));
    REQUIRE(std::imag(values[25]) == Approx(-0.08028180462518811));

    //Move electron 0 to position [[6.666665 6.666665 6.666665]] a.u.:
    elec.R[0] = {6.666665000000002, 6.666665000000002, 6.666665000000002};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[6.666665 6.666665 6.666665]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(0.11380752148567085));
    REQUIRE(std::imag(values[0]) == Approx(-0.04390587226150721));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.02588030482486917));
    REQUIRE(std::imag(values[1]) == Approx(0.10687921643358252));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(3.197341334549546e-13));
    REQUIRE(std::imag(values[2]) == Approx(-2.3414238599544394e-14));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-1.3703920844638083e-12));
    REQUIRE(std::imag(values[3]) == Approx(-4.923708702683163e-14));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.18010843596428783));
    REQUIRE(std::imag(values[4]) == Approx(0.04962376866516952));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(-4.672976839154819e-12));
    REQUIRE(std::imag(values[5]) == Approx(-1.544795977090896e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(1.2509722576569074e-12));
    REQUIRE(std::imag(values[6]) == Approx(3.563501771870871e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.010394738060449391));
    REQUIRE(std::imag(values[7]) == Approx(-0.15080653874439975));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(-0.22742833831292314));
    REQUIRE(std::imag(values[8]) == Approx(-0.08681115895058802));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(7.017990606811553e-13));
    REQUIRE(std::imag(values[9]) == Approx(-6.006243721193688e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(-3.5441754061629804e-12));
    REQUIRE(std::imag(values[10]) == Approx(5.510674056680197e-13));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(0.19602654856168172));
    REQUIRE(std::imag(values[11]) == Approx(0.08793147752366054));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(1.1100410142953716e-13));
    REQUIRE(std::imag(values[12]) == Approx(2.1529784776147745e-13));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(-3.1113994995362083e-13));
    REQUIRE(std::imag(values[13]) == Approx(-9.16124845162143e-14));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(-2.8425823215546873e-12));
    REQUIRE(std::imag(values[14]) == Approx(-5.967361937077563e-13));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(8.500061674827252e-12));
    REQUIRE(std::imag(values[15]) == Approx(1.6117350650343773e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(0.1685048554426682));
    REQUIRE(std::imag(values[16]) == Approx(0.040751301285559566));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(-0.05992911952275292));
    REQUIRE(std::imag(values[17]) == Approx(-0.147726013001218));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-6.090732565369773e-13));
    REQUIRE(std::imag(values[18]) == Approx(-7.443355801535639e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(-2.2154083501205196e-12));
    REQUIRE(std::imag(values[19]) == Approx(-5.02729118133428e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.061780347627032536));
    REQUIRE(std::imag(values[20]) == Approx(0.01686579925949514));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(3.0787695380514886e-13));
    REQUIRE(std::imag(values[21]) == Approx(1.6918178698834226e-14));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(-3.221671335678293e-12));
    REQUIRE(std::imag(values[22]) == Approx(-9.562454010396891e-13));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(-4.3583890429925625e-13));
    REQUIRE(std::imag(values[23]) == Approx(1.578955716802118e-12));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(-9.202570529310635e-13));
    REQUIRE(std::imag(values[24]) == Approx(1.13342674089293e-12));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(0.7098778186417222));
    REQUIRE(std::imag(values[25]) == Approx(-0.928839302822032));

    //Move electron 0 to position [[9.999998 9.999998 9.999998]] a.u.:
    elec.R[0] = {9.999998000000001, 9.999998000000001, 9.999998000000001};
    elec.update();
    sposet->evaluateValue(elec, 0, values);

    // Position [[9.999998 9.999998 9.999998]] a.u.
    // Verifying values of SPO 0
    REQUIRE(std::real(values[0]) == Approx(-0.05950698756399922));
    REQUIRE(std::imag(values[0]) == Approx(0.008451159090256224));

    // Verifying values of SPO 1
    REQUIRE(std::real(values[1]) == Approx(0.0372927266605709));
    REQUIRE(std::imag(values[1]) == Approx(0.04865086992861883));

    // Verifying values of SPO 2
    REQUIRE(std::real(values[2]) == Approx(2.3700317966038143e-13));
    REQUIRE(std::imag(values[2]) == Approx(-1.0025983406577301e-13));

    // Verifying values of SPO 3
    REQUIRE(std::real(values[3]) == Approx(-9.70793757903483e-13));
    REQUIRE(std::imag(values[3]) == Approx(2.724608733149554e-13));

    // Verifying values of SPO 4
    REQUIRE(std::real(values[4]) == Approx(0.02755063738588984));
    REQUIRE(std::imag(values[4]) == Approx(-0.09640122180554679));

    // Verifying values of SPO 5
    REQUIRE(std::real(values[5]) == Approx(-1.6269971476154967e-13));
    REQUIRE(std::imag(values[5]) == Approx(3.438238409269704e-12));

    // Verifying values of SPO 6
    REQUIRE(std::real(values[6]) == Approx(-2.855643508474324e-14));
    REQUIRE(std::imag(values[6]) == Approx(-8.265646196873323e-13));

    // Verifying values of SPO 7
    REQUIRE(std::real(values[7]) == Approx(0.16726281762463338));
    REQUIRE(std::imag(values[7]) == Approx(-0.052898517004949845));

    // Verifying values of SPO 8
    REQUIRE(std::real(values[8]) == Approx(0.13334645793969357));
    REQUIRE(std::imag(values[8]) == Approx(-0.10880327269811665));

    // Verifying values of SPO 9
    REQUIRE(std::real(values[9]) == Approx(7.692422002437098e-13));
    REQUIRE(std::imag(values[9]) == Approx(-6.640972494690012e-13));

    // Verifying values of SPO 10
    REQUIRE(std::real(values[10]) == Approx(-2.4821733211548425e-12));
    REQUIRE(std::imag(values[10]) == Approx(1.2520175866069816e-12));

    // Verifying values of SPO 11
    REQUIRE(std::real(values[11]) == Approx(0.1798235956119672));
    REQUIRE(std::imag(values[11]) == Approx(-0.0874643862941136));

    // Verifying values of SPO 12
    REQUIRE(std::real(values[12]) == Approx(-1.1709383431389981e-14));
    REQUIRE(std::imag(values[12]) == Approx(7.639071676987064e-14));

    // Verifying values of SPO 13
    REQUIRE(std::real(values[13]) == Approx(8.413148674931816e-14));
    REQUIRE(std::imag(values[13]) == Approx(5.638892139842643e-14));

    // Verifying values of SPO 14
    REQUIRE(std::real(values[14]) == Approx(3.033437215039551e-13));
    REQUIRE(std::imag(values[14]) == Approx(1.6158520918629484e-12));

    // Verifying values of SPO 15
    REQUIRE(std::real(values[15]) == Approx(-1.4472936742422837e-12));
    REQUIRE(std::imag(values[15]) == Approx(-4.016332839335127e-12));

    // Verifying values of SPO 16
    REQUIRE(std::real(values[16]) == Approx(-0.03653573483382949));
    REQUIRE(std::imag(values[16]) == Approx(-0.07139149065641964));

    // Verifying values of SPO 17
    REQUIRE(std::real(values[17]) == Approx(0.05195911753492378));
    REQUIRE(std::imag(values[17]) == Approx(0.20377702124340888));

    // Verifying values of SPO 18
    REQUIRE(std::real(values[18]) == Approx(-1.3913349264311425e-14));
    REQUIRE(std::imag(values[18]) == Approx(6.650238087623492e-13));

    // Verifying values of SPO 19
    REQUIRE(std::real(values[19]) == Approx(1.5612684760915313e-12));
    REQUIRE(std::imag(values[19]) == Approx(5.4326976941529405e-12));

    // Verifying values of SPO 20
    REQUIRE(std::real(values[20]) == Approx(0.05183727962592069));
    REQUIRE(std::imag(values[20]) == Approx(-0.08693036223348635));

    // Verifying values of SPO 21
    REQUIRE(std::real(values[21]) == Approx(-6.168889188205947e-13));
    REQUIRE(std::imag(values[21]) == Approx(-7.934625127721342e-14));

    // Verifying values of SPO 22
    REQUIRE(std::real(values[22]) == Approx(1.854426333865271e-12));
    REQUIRE(std::imag(values[22]) == Approx(1.1091023934392614e-12));

    // Verifying values of SPO 23
    REQUIRE(std::real(values[23]) == Approx(1.4450940425974687e-13));
    REQUIRE(std::imag(values[23]) == Approx(-4.0532507912929545e-13));

    // Verifying values of SPO 24
    REQUIRE(std::real(values[24]) == Approx(9.63345721944599e-13));
    REQUIRE(std::imag(values[24]) == Approx(-3.381600555043806e-13));

    // Verifying values of SPO 25
    REQUIRE(std::real(values[25]) == Approx(-0.03316664551198499));
    REQUIRE(std::imag(values[25]) == Approx(-0.10326287192308459));

    // END generated C++ input from Carbon1x1x1-tw1_gen_mos.py (pyscf version 1.6.2) on 2019-11-19 15:08:42.847403
  }
}

TEST_CASE("ReadMolecularOrbital GTO Carbon Diamond", "[wavefunction]") { test_C_diamond(); }

} // namespace qmcplusplus

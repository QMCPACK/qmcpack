//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


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
void test_He(bool transform)
{
  std::ostringstream section_name;
  section_name << "He, transform orbitals to grid: " << (transform ? "T" : "F");

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    const SimulationCell simulation_cell;
    ParticleSet elec(simulation_cell);
    std::vector<int> agroup(2);
    agroup[0] = 1;
    agroup[1] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0] = 0.0;

    SpeciesSet& tspecies       = elec.getSpeciesSet();
    int upIdx                  = tspecies.addSpecies("u");
    int downIdx                = tspecies.addSpecies("d");
    int massIdx                = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx)   = 1.0;
    tspecies(massIdx, downIdx) = 1.0;

    ParticleSet ions(simulation_cell);
    ions.setName("ion0");
    ions.create(1);
    ions.R[0]            = 0.0;
    SpeciesSet& ispecies = ions.getSpeciesSet();
    int heIdx            = ispecies.addSpecies("He");
    ions.update();


    elec.addTable(ions);
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("he_sto3g.wfj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
    particle_set_map["e"]    = &elec;
    particle_set_map["ion0"] = &ions;


    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);
    if (transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
    }
    else
    {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar*)"transform", (const xmlChar*)"no");
      xmlSetProp(MO_base[0], (const xmlChar*)"key", (const xmlChar*)"GTO");
    }

    auto& bb = bf.createSPOSetBuilder(MO_base[0]);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    SPOSet* sposet = bb.createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSet::ValueVector_t values;
    SPOSet::GradVector_t dpsi;
    SPOSet::ValueVector_t d2psi;
    values.resize(1);
    dpsi.resize(1);
    d2psi.resize(1);

    // Call makeMove to compute the distances
    ParticleSet::SingleParticlePos_t newpos(0.0001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    REQUIRE(values[0] == Approx(0.9996037001));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    REQUIRE(values[0] == Approx(0.9996037001));
    REQUIRE(dpsi[0][0] == Approx(-0.0006678035459));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));
    REQUIRE(d2psi[0] == Approx(-20.03410564));


    ParticleSet::SingleParticlePos_t disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(0.2315567641));
    REQUIRE(dpsi[0][0] == Approx(-0.3805431885));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));
    REQUIRE(d2psi[0] == Approx(-0.2618497452));
  }
}

TEST_CASE("ReadMolecularOrbital GTO He", "[wavefunction]") { test_He(false); }

TEST_CASE("ReadMolecularOrbital Numerical He", "[wavefunction]") { test_He(true); }


void test_Ne(bool transform)
{
  std::ostringstream section_name;
  section_name << "Neon, transform orbitals to grid: " << (transform ? "T" : "F");

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    const SimulationCell simulation_cell;
    ParticleSet elec(simulation_cell);

    std::vector<int> agroup(2);
    agroup[0] = 1;
    agroup[1] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0] = 0.0;

    SpeciesSet& tspecies       = elec.getSpeciesSet();
    int upIdx                  = tspecies.addSpecies("u");
    int downIdx                = tspecies.addSpecies("d");
    int massIdx                = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx)   = 1.0;
    tspecies(massIdx, downIdx) = 1.0;

    ParticleSet ions(simulation_cell);
    ions.setName("ion0");
    ions.create(1);
    ions.R[0][0]         = 0.0;
    ions.R[0][1]         = 0.0;
    ions.R[0][2]         = 0.0;
    SpeciesSet& ispecies = ions.getSpeciesSet();
    int heIdx            = ispecies.addSpecies("Ne");
    ions.update();


    elec.addTable(ions);
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("ne_def2_svp.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
    particle_set_map["e"]    = &elec;
    particle_set_map["ion0"] = &ions;


    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
    }
    else
    {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar*)"transform", (const xmlChar*)"no");
      xmlSetProp(MO_base[0], (const xmlChar*)"key", (const xmlChar*)"GTO");
    }

    auto& bb = bf.createSPOSetBuilder(MO_base[0]);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    SPOSet* sposet = bb.createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSet::ValueVector_t values;
    SPOSet::GradVector_t dpsi;
    SPOSet::ValueVector_t d2psi;
    values.resize(5);
    dpsi.resize(5);
    d2psi.resize(5);

    ParticleSet::SingleParticlePos_t newpos(0.00001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    std::cout << "values = " << values << std::endl;

    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-16.11819042));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    std::cout << "values = " << values << std::endl;
    std::cout << "dpsi = " << dpsi << std::endl;
    std::cout << "d2psi = " << d2psi << std::endl;
    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-16.11819042));
    REQUIRE(dpsi[0][0] == Approx(0.1747261458));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));


    ParticleSet::SingleParticlePos_t disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);
    sposet->evaluateValue(elec, 0, values);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-0.005041631374));


    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-0.005041631374));
    REQUIRE(dpsi[0][0] == Approx(0.01862216578));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));
    REQUIRE(d2psi[0] == Approx(-0.01551755818));
  }
}

TEST_CASE("ReadMolecularOrbital GTO Ne", "[wavefunction]") { test_Ne(false); }

TEST_CASE("ReadMolecularOrbital Numerical Ne", "[wavefunction]") { test_Ne(true); }

TEST_CASE("ReadMolecularOrbital HCN", "[wavefunction]") {}

void test_HCN(bool transform)
{
  std::ostringstream section_name;
  section_name << "HCN, transform orbitals to grid: " << (transform ? "T" : "F");

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    Libxml2Document doc;
    bool okay = doc.parse("hcn.structure.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    const SimulationCell simulation_cell;
    ParticleSet ions(simulation_cell);
    XMLParticleParser parse_ions(ions);
    OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
    REQUIRE(particleset_ion.size() == 1);
    parse_ions.put(particleset_ion[0]);

    REQUIRE(ions.groups() == 3);
    REQUIRE(ions.R.size() == 3);
    ions.update();

    ParticleSet elec(simulation_cell);
    XMLParticleParser parse_elec(elec);
    OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
    REQUIRE(particleset_elec.size() == 1);
    parse_elec.put(particleset_elec[0]);

    REQUIRE(elec.groups() == 2);
    REQUIRE(elec.R.size() == 14);

    elec.R = 0.0;

    elec.addTable(ions);
    elec.update();

    Libxml2Document doc2;
    okay = doc2.parse("hcn.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root2 = doc2.getRoot();

    WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
    particle_set_map["e"]    = &elec;
    particle_set_map["ion0"] = &ions;

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
    }
    else
    {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar*)"transform", (const xmlChar*)"no");
      xmlSetProp(MO_base[0], (const xmlChar*)"key", (const xmlChar*)"GTO");
    }

    auto& bb = bf.createSPOSetBuilder(MO_base[0]);

    OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
    SPOSet* sposet = bb.createSPOSet(slater_base[0]);


    SPOSet::ValueVector_t values;
    SPOSet::GradVector_t dpsi;
    SPOSet::ValueVector_t d2psi;
    values.resize(7);
    dpsi.resize(7);
    d2psi.resize(7);

    ParticleSet::SingleParticlePos_t newpos;
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    REQUIRE(values[0] == Approx(0.009452265234));
    REQUIRE(values[1] == Approx(0.02008357407));
    REQUIRE(values[2] == Approx(0.4163749594));
    REQUIRE(values[3] == Approx(-0.08854428343));
    REQUIRE(values[4] == Approx(0.273158705));
    REQUIRE(values[5] == Approx(0));
    REQUIRE(values[6] == Approx(0));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);


    // Generated form gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(0.009452265234));
    REQUIRE(dpsi[0][0] == Approx(-0.05400764372));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));
    REQUIRE(d2psi[0] == Approx(0.2532157143));

    REQUIRE(values[1] == Approx(0.02008357407));
    REQUIRE(dpsi[1][0] == Approx(0.1009262252));
    REQUIRE(dpsi[1][1] == Approx(0));
    REQUIRE(dpsi[1][2] == Approx(0));
    REQUIRE(d2psi[1] == Approx(0.3423520138));

    REQUIRE(values[2] == Approx(0.4163749594));
    REQUIRE(dpsi[2][0] == Approx(-0.1202256419));
    REQUIRE(dpsi[2][1] == Approx(0));
    REQUIRE(dpsi[2][2] == Approx(0));
    REQUIRE(d2psi[2] == Approx(-1.178149899));

    REQUIRE(values[3] == Approx(-0.08854428343));
    REQUIRE(dpsi[3][0] == Approx(-0.004505552544));
    REQUIRE(dpsi[3][1] == Approx(0));
    REQUIRE(dpsi[3][2] == Approx(0));
    REQUIRE(d2psi[3] == Approx(0.2838238091));

    REQUIRE(values[4] == Approx(0.273158705));
    REQUIRE(dpsi[4][0] == Approx(-0.01125044248));
    REQUIRE(dpsi[4][1] == Approx(0));
    REQUIRE(dpsi[4][2] == Approx(0));
    REQUIRE(d2psi[4] == Approx(-0.9173261582));

    REQUIRE(values[5] == Approx(0));
    REQUIRE(dpsi[5][0] == Approx(0));
    REQUIRE(dpsi[5][1] == Approx(0.4221165864));
    REQUIRE(dpsi[5][2] == Approx(-0.08191634629));
    REQUIRE(d2psi[5] == Approx(0));

    REQUIRE(values[6] == Approx(0));
    REQUIRE(dpsi[6][0] == Approx(0));
    REQUIRE(dpsi[6][1] == Approx(0.08191634629));
    REQUIRE(dpsi[6][2] == Approx(0.4221165864));
    REQUIRE(d2psi[6] == Approx(0));

    //==========Hessian and Grad Hessian Test==========
    SPOSet::HessVector_t dhpsi;
    dhpsi.resize(7);

    sposet->evaluateVGH(elec, 0, values, dpsi, dhpsi);

    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(0.009452265234));
    REQUIRE(dpsi[0][0] == Approx(-0.05400764372));
    REQUIRE(dpsi[0][1] == Approx(0));
    REQUIRE(dpsi[0][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    REQUIRE(dhpsi[0](0, 0) == Approx(0.3523924743));
    REQUIRE(dhpsi[0](0, 1) == Approx(0));
    REQUIRE(dhpsi[0](0, 2) == Approx(0));
    REQUIRE(dhpsi[0](1, 1) == Approx(-0.04958838002));
    REQUIRE(dhpsi[0](1, 2) == Approx(0));
    REQUIRE(dhpsi[0](2, 2) == Approx(-0.04958838002));

    /////////////////////////////////////////////////////////////////////////////
    //Now we're going to test the API's called by SPOSet for higher derivatives//
    /////////////////////////////////////////////////////////////////////////////
    SPOSet::ValueMatrix_t psiM(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GradMatrix_t dpsiM(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::HessMatrix_t hesspsiV(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GGGMatrix_t d3psiV(elec.R.size(), sposet->getOrbitalSetSize());

    sposet->evaluate_notranspose(elec, 0, elec.R.size(), psiM, dpsiM, hesspsiV, d3psiV);


    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(psiM[0][0] == Approx(0.009452265234));
    REQUIRE(dpsiM[0][0][0] == Approx(-0.05400764372));
    REQUIRE(dpsiM[0][0][1] == Approx(0));
    REQUIRE(dpsiM[0][0][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    REQUIRE(hesspsiV[0][0](0, 0) == Approx(0.3523924743));
    REQUIRE(hesspsiV[0][0](0, 1) == Approx(0));
    REQUIRE(hesspsiV[0][0](0, 2) == Approx(0));
    REQUIRE(hesspsiV[0][0](1, 1) == Approx(-0.04958838002));
    REQUIRE(hesspsiV[0][0](1, 2) == Approx(0));
    REQUIRE(hesspsiV[0][0](2, 2) == Approx(-0.04958838002));

    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    REQUIRE(d3psiV[0][0][0](0, 0) == Approx(-2.241965465));
    REQUIRE(d3psiV[0][0][0](0, 1) == Approx(0));
    REQUIRE(d3psiV[0][0][0](0, 2) == Approx(0));
    REQUIRE(d3psiV[0][0][0](1, 1) == Approx(0.3714481861));
    REQUIRE(d3psiV[0][0][0](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][0][0](2, 2) == Approx(0.3714481861));
    REQUIRE(d3psiV[0][0][1](1, 1) == Approx(0));
    REQUIRE(d3psiV[0][0][1](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][0][1](2, 2) == Approx(0));
    REQUIRE(d3psiV[0][0][2](2, 2) == Approx(0));

    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(psiM[0][1] == Approx(0.02008357407));
    REQUIRE(dpsiM[0][1][0] == Approx(0.1009262252));
    REQUIRE(dpsiM[0][1][1] == Approx(0));
    REQUIRE(dpsiM[0][1][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    REQUIRE(hesspsiV[0][1](0, 0) == Approx(0.5298289497));
    REQUIRE(hesspsiV[0][1](0, 1) == Approx(0));
    REQUIRE(hesspsiV[0][1](0, 2) == Approx(0));
    REQUIRE(hesspsiV[0][1](1, 1) == Approx(-0.09373846794));
    REQUIRE(hesspsiV[0][1](1, 2) == Approx(0));
    REQUIRE(hesspsiV[0][1](2, 2) == Approx(-0.09373846794));

    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    REQUIRE(d3psiV[0][1][0](0, 0) == Approx(2.594787656));
    REQUIRE(d3psiV[0][1][0](0, 1) == Approx(0));
    REQUIRE(d3psiV[0][1][0](0, 2) == Approx(0));
    REQUIRE(d3psiV[0][1][0](1, 1) == Approx(-0.5720485625));
    REQUIRE(d3psiV[0][1][0](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][1][0](2, 2) == Approx(-0.5720485625));
    REQUIRE(d3psiV[0][1][1](1, 1) == Approx(0));
    REQUIRE(d3psiV[0][1][1](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][1][1](2, 2) == Approx(0));
    REQUIRE(d3psiV[0][1][2](2, 2) == Approx(0));
    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(psiM[0][2] == Approx(0.4163749594));
    REQUIRE(dpsiM[0][2][0] == Approx(-0.1202256419));
    REQUIRE(dpsiM[0][2][1] == Approx(0));
    REQUIRE(dpsiM[0][2][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    REQUIRE(hesspsiV[0][2](0, 0) == Approx(-0.02607695984));
    REQUIRE(hesspsiV[0][2](0, 1) == Approx(0));
    REQUIRE(hesspsiV[0][2](0, 2) == Approx(0));
    REQUIRE(hesspsiV[0][2](1, 1) == Approx(-0.5760364698));
    REQUIRE(hesspsiV[0][2](1, 2) == Approx(0));
    REQUIRE(hesspsiV[0][2](2, 2) == Approx(-0.5760364698));
    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    REQUIRE(d3psiV[0][2][0](0, 0) == Approx(-0.227147312));
    REQUIRE(d3psiV[0][2][0](0, 1) == Approx(0));
    REQUIRE(d3psiV[0][2][0](0, 2) == Approx(0));
    REQUIRE(d3psiV[0][2][0](1, 1) == Approx(0.2992015499));
    REQUIRE(d3psiV[0][2][0](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][2][0](2, 2) == Approx(0.2992015499));
    REQUIRE(d3psiV[0][2][1](1, 1) == Approx(0));
    REQUIRE(d3psiV[0][2][1](1, 2) == Approx(0));
    REQUIRE(d3psiV[0][2][1](2, 2) == Approx(0));
    REQUIRE(d3psiV[0][2][2](2, 2) == Approx(0));


    //Move electron 0 to some nontrivial position:
    ParticleSet::SingleParticlePos_t disp(0.02, -0.1, 0.05);
    elec.makeMove(0, disp);


    SPOSet::GradMatrix_t dionpsi(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::HessMatrix_t diongradpsi(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GradMatrix_t dionlaplpsi(elec.R.size(), sposet->getOrbitalSetSize());

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 0, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  0  Component  0 ===================
    REQUIRE(dionpsi[0][0][0] == Approx(0.0453112082));
    REQUIRE(diongradpsi[0][0](0, 0) == Approx(-0.2943513994));
    REQUIRE(diongradpsi[0][0](0, 1) == Approx(0.030468047));
    REQUIRE(diongradpsi[0][0](0, 2) == Approx(-0.0152340235));
    REQUIRE(dionlaplpsi[0][0][0] == Approx(1.333755581));
    REQUIRE(dionpsi[0][1][0] == Approx(-0.0006473819623));
    REQUIRE(diongradpsi[0][1](0, 0) == Approx(0.0004713407512));
    REQUIRE(diongradpsi[0][1](0, 1) == Approx(-0.0001254975603));
    REQUIRE(diongradpsi[0][1](0, 2) == Approx(6.274878013e-05));
    REQUIRE(dionlaplpsi[0][1][0] == Approx(0.001057940846));
    REQUIRE(dionpsi[0][2][0] == Approx(0.265578336));
    REQUIRE(diongradpsi[0][2](0, 0) == Approx(-0.08685804115));
    REQUIRE(diongradpsi[0][2](0, 1) == Approx(0.05438178417));
    REQUIRE(diongradpsi[0][2](0, 2) == Approx(-0.02719089209));
    REQUIRE(dionlaplpsi[0][2][0] == Approx(-1.882489819));
    REQUIRE(dionpsi[0][3][0] == Approx(-0.06444305979));
    REQUIRE(diongradpsi[0][3](0, 0) == Approx(-0.002013151923));
    REQUIRE(diongradpsi[0][3](0, 1) == Approx(-0.002535923431));
    REQUIRE(diongradpsi[0][3](0, 2) == Approx(0.001267961716));
    REQUIRE(dionlaplpsi[0][3][0] == Approx(0.4547401581));
    REQUIRE(dionpsi[0][4][0] == Approx(0.1454357726));
    REQUIRE(diongradpsi[0][4](0, 0) == Approx(-0.2330499431));
    REQUIRE(diongradpsi[0][4](0, 1) == Approx(0.09667641762));
    REQUIRE(diongradpsi[0][4](0, 2) == Approx(-0.04833820881));
    REQUIRE(dionlaplpsi[0][4][0] == Approx(-0.9197558839));
    REQUIRE(dionpsi[0][5][0] == Approx(-0.04329985085));
    REQUIRE(diongradpsi[0][5](0, 0) == Approx(0.09051993304));
    REQUIRE(diongradpsi[0][5](0, 1) == Approx(0.382375474));
    REQUIRE(diongradpsi[0][5](0, 2) == Approx(-0.07043361927));
    REQUIRE(dionlaplpsi[0][5][0] == Approx(0.2201672051));
    REQUIRE(dionpsi[0][6][0] == Approx(0.01207541177));
    REQUIRE(diongradpsi[0][6](0, 0) == Approx(-0.02524405435));
    REQUIRE(diongradpsi[0][6](0, 1) == Approx(0.0800332842));
    REQUIRE(diongradpsi[0][6](0, 2) == Approx(0.3929818664));
    REQUIRE(dionlaplpsi[0][6][0] == Approx(-0.0614000824));


    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 1, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  1  Component  1 ===================
    REQUIRE(dionpsi[0][0][1] == Approx(0.0001412373768));
    REQUIRE(diongradpsi[0][0](1, 0) == Approx(0.0001124265646));
    REQUIRE(diongradpsi[0][0](1, 1) == Approx(-0.001383378615));
    REQUIRE(diongradpsi[0][0](1, 2) == Approx(-1.449757545e-05));
    REQUIRE(dionlaplpsi[0][0][1] == Approx(-0.001252043663));
    REQUIRE(dionpsi[0][1][1] == Approx(-0.01029290716));
    REQUIRE(diongradpsi[0][1](1, 0) == Approx(-0.06178485148));
    REQUIRE(diongradpsi[0][1](1, 1) == Approx(0.0971577216));
    REQUIRE(diongradpsi[0][1](1, 2) == Approx(0.002885675005));
    REQUIRE(dionlaplpsi[0][1][1] == Approx(-0.1403103458));
    REQUIRE(dionpsi[0][2][1] == Approx(-0.0230872583));
    REQUIRE(diongradpsi[0][2](1, 0) == Approx(-0.02537847709));
    REQUIRE(diongradpsi[0][2](1, 1) == Approx(0.2268946564));
    REQUIRE(diongradpsi[0][2](1, 2) == Approx(0.001988963201));
    REQUIRE(dionlaplpsi[0][2][1] == Approx(0.2028851421));
    REQUIRE(dionpsi[0][3][1] == Approx(0.01850231814));
    REQUIRE(diongradpsi[0][3](1, 0) == Approx(0.05709948475));
    REQUIRE(diongradpsi[0][3](1, 1) == Approx(-0.1776515965));
    REQUIRE(diongradpsi[0][3](1, 2) == Approx(-0.003685792479));
    REQUIRE(dionlaplpsi[0][3][1] == Approx(-0.1280699725));
    REQUIRE(dionpsi[0][4][1] == Approx(-0.02136209962));
    REQUIRE(diongradpsi[0][4](1, 0) == Approx(-0.03836586276));
    REQUIRE(diongradpsi[0][4](1, 1) == Approx(0.2084578148));
    REQUIRE(diongradpsi[0][4](1, 2) == Approx(0.002581590766));
    REQUIRE(dionlaplpsi[0][4][1] == Approx(0.1792683544));
    REQUIRE(dionpsi[0][5][1] == Approx(-0.1942343714));
    REQUIRE(diongradpsi[0][5](1, 0) == Approx(-0.3037357197));
    REQUIRE(diongradpsi[0][5](1, 1) == Approx(-0.09561978734));
    REQUIRE(diongradpsi[0][5](1, 2) == Approx(0.02118492506));
    REQUIRE(dionlaplpsi[0][5][1] == Approx(0.6410434658));
    REQUIRE(dionpsi[0][6][1] == Approx(-0.03930992259));
    REQUIRE(diongradpsi[0][6](1, 0) == Approx(-0.06331544695));
    REQUIRE(diongradpsi[0][6](1, 1) == Approx(-0.002807368817));
    REQUIRE(diongradpsi[0][6](1, 2) == Approx(-0.02801340823));
    REQUIRE(dionlaplpsi[0][6][1] == Approx(0.1369061053));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 2, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  2  Component  2 ===================
    REQUIRE(dionpsi[0][0][2] == Approx(1.302648961e-06));
    REQUIRE(diongradpsi[0][0](2, 0) == Approx(1.865129579e-06));
    REQUIRE(diongradpsi[0][0](2, 1) == Approx(6.142092043e-08));
    REQUIRE(diongradpsi[0][0](2, 2) == Approx(2.602225618e-05));
    REQUIRE(dionlaplpsi[0][0][2] == Approx(1.234692903e-06));
    REQUIRE(dionpsi[0][1][2] == Approx(3.248738084e-07));
    REQUIRE(diongradpsi[0][1](2, 0) == Approx(-2.044420189e-06));
    REQUIRE(diongradpsi[0][1](2, 1) == Approx(-7.011145137e-08));
    REQUIRE(diongradpsi[0][1](2, 2) == Approx(6.532522353e-06));
    REQUIRE(dionlaplpsi[0][1][2] == Approx(-6.10958506e-06));
    REQUIRE(dionpsi[0][2][2] == Approx(3.264249981e-06));
    REQUIRE(diongradpsi[0][2](2, 0) == Approx(2.820971234e-05));
    REQUIRE(diongradpsi[0][2](2, 1) == Approx(9.405184964e-07));
    REQUIRE(diongradpsi[0][2](2, 2) == Approx(6.481420782e-05));
    REQUIRE(dionlaplpsi[0][2][2] == Approx(5.73961989e-05));
    REQUIRE(dionpsi[0][3][2] == Approx(0.0001288974413));
    REQUIRE(diongradpsi[0][3](2, 0) == Approx(0.0002840756879));
    REQUIRE(diongradpsi[0][3](2, 1) == Approx(9.281700408e-06));
    REQUIRE(diongradpsi[0][3](2, 2) == Approx(0.002573308008));
    REQUIRE(dionlaplpsi[0][3][2] == Approx(0.0003025314443));
    REQUIRE(dionpsi[0][4][2] == Approx(-7.300043903e-05));
    REQUIRE(diongradpsi[0][4](2, 0) == Approx(-0.0001000016834));
    REQUIRE(diongradpsi[0][4](2, 1) == Approx(-3.233243534e-06));
    REQUIRE(diongradpsi[0][4](2, 2) == Approx(-0.001458391774));
    REQUIRE(dionlaplpsi[0][4][2] == Approx(-3.546690719e-05));
    REQUIRE(dionpsi[0][5][2] == Approx(2.910525987e-06));
    REQUIRE(diongradpsi[0][5](2, 0) == Approx(1.307065133e-05));
    REQUIRE(diongradpsi[0][5](2, 1) == Approx(1.560390706e-06));
    REQUIRE(diongradpsi[0][5](2, 2) == Approx(-2.92731811e-06));
    REQUIRE(dionlaplpsi[0][5][2] == Approx(3.797816228e-05));
    REQUIRE(dionpsi[0][6][2] == Approx(-1.56074936e-05));
    REQUIRE(diongradpsi[0][6](2, 0) == Approx(-7.009049656e-05));
    REQUIRE(diongradpsi[0][6](2, 1) == Approx(-2.048666792e-06));
    REQUIRE(diongradpsi[0][6](2, 2) == Approx(2.967709412e-06));
    REQUIRE(dionlaplpsi[0][6][2] == Approx(-0.0002018111858));

    //Same tests as before, but for the gradient only.

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 0, dionpsi);
    //============== Ion  0  Component  0 ===================
    REQUIRE(dionpsi[0][0][0] == Approx(0.0453112082));
    REQUIRE(dionpsi[0][1][0] == Approx(-0.0006473819623));
    REQUIRE(dionpsi[0][2][0] == Approx(0.265578336));
    REQUIRE(dionpsi[0][3][0] == Approx(-0.06444305979));
    REQUIRE(dionpsi[0][4][0] == Approx(0.1454357726));
    REQUIRE(dionpsi[0][5][0] == Approx(-0.04329985085));
    REQUIRE(dionpsi[0][6][0] == Approx(0.01207541177));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 1, dionpsi);
    //============== Ion  1  Component  1 ===================
    REQUIRE(dionpsi[0][0][1] == Approx(0.0001412373768));
    REQUIRE(dionpsi[0][1][1] == Approx(-0.01029290716));
    REQUIRE(dionpsi[0][2][1] == Approx(-0.0230872583));
    REQUIRE(dionpsi[0][3][1] == Approx(0.01850231814));
    REQUIRE(dionpsi[0][4][1] == Approx(-0.02136209962));
    REQUIRE(dionpsi[0][5][1] == Approx(-0.1942343714));
    REQUIRE(dionpsi[0][6][1] == Approx(-0.03930992259));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 2, dionpsi);
    //============== Ion  2  Component  2 ===================
    REQUIRE(dionpsi[0][0][2] == Approx(1.302648961e-06));
    REQUIRE(dionpsi[0][1][2] == Approx(3.248738084e-07));
    REQUIRE(dionpsi[0][2][2] == Approx(3.264249981e-06));
    REQUIRE(dionpsi[0][3][2] == Approx(0.0001288974413));
    REQUIRE(dionpsi[0][4][2] == Approx(-7.300043903e-05));
    REQUIRE(dionpsi[0][5][2] == Approx(2.910525987e-06));
    REQUIRE(dionpsi[0][6][2] == Approx(-1.56074936e-05));
  }
}

TEST_CASE("ReadMolecularOrbital GTO HCN", "[wavefunction]") { test_HCN(false); }

TEST_CASE("ReadMolecularOrbital Numerical HCN", "[wavefunction]") { test_HCN(true); }

} // namespace qmcplusplus

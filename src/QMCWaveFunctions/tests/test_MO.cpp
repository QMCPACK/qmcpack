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
#include "Particle/DistanceTableData.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Numerics/GaussianBasisSet.h"
#ifdef ENABLE_SOA
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#else
#include "QMCWaveFunctions/MolecularOrbitals/LocalizedBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/LCOrbitalSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/SphericalBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#endif
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{

void test_He(bool transform)
{
  std::ostringstream section_name;
  section_name << "He, transform orbitals to grid: " << (transform ? "T":"F");

  SECTION(section_name.str()) {
    OHMMS::Controller->initialize(0, NULL);
    Communicate *c = OHMMS::Controller;

    ParticleSet elec;
    std::vector<int> agroup(2);
    agroup[0] = 1;
    agroup[1] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0] = 0.0;

    SpeciesSet &tspecies = elec.getSpeciesSet();
    int upIdx = tspecies.addSpecies("u");
    int downIdx = tspecies.addSpecies("d");
    int massIdx = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx) = 1.0;
    tspecies(massIdx, downIdx) = 1.0;

    ParticleSet ions;
    ions.setName("ion0");
    ions.create(1);
    ions.R[0] = 0.0;
    SpeciesSet &ispecies = ions.getSpeciesSet();
    int heIdx = ispecies.addSpecies("He");
    ions.update();


  #ifdef ENABLE_SOA
    elec.addTable(ions,DT_SOA);
  #else
    elec.addTable(ions,DT_AOS);
  #endif
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("he_sto3g.wfj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    TrialWaveFunction psi(c);

    OrbitalBuilderBase::PtclPoolType particle_set_map;
    particle_set_map["e"] = &elec;
    particle_set_map["ion0"] = &ions;


    SPOSetBuilderFactory bf(elec, psi, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);
    if (transform) {
      // input file is set to transform GTO's to numerical orbitals by default
    } else {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar *)"transform", (const xmlChar *)"no");
      xmlSetProp(MO_base[0], (const xmlChar *)"key", (const xmlChar *)"GTO");
    }

    SPOSetBuilder *bb = bf.createSPOSetBuilder(MO_base[0]);
    REQUIRE(bb != NULL);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    bb->loadBasisSetFromXML(MO_base[0]);
    SPOSetBase *sposet = bb->createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSetBase::ValueVector_t values;
    SPOSetBase::GradVector_t dpsi;
    SPOSetBase::ValueVector_t d2psi;
    values.resize(1);
    dpsi.resize(1);
    d2psi.resize(1);

    // Call makeMove to compute the distances
    ParticleSet::SingleParticlePos_t newpos(0.0001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluate(elec, 0, values);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    REQUIRE(values[0] == Approx(   0.9996037001));

    sposet->evaluate(elec, 0, values, dpsi, d2psi);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    REQUIRE(values[0] == Approx(   0.9996037001));
    REQUIRE(dpsi[0][0] == Approx(-0.0006678035459));
    REQUIRE(dpsi[0][1] == Approx(              0));
    REQUIRE(dpsi[0][2] == Approx(              0));
    REQUIRE(d2psi[0] == Approx(   -20.03410564));


    ParticleSet::SingleParticlePos_t disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);

    sposet->evaluate(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(   0.2315567641));
    REQUIRE(dpsi[0][0] == Approx(  -0.3805431885));
    REQUIRE(dpsi[0][1] == Approx(              0));
    REQUIRE(dpsi[0][2] == Approx(              0));
    REQUIRE(d2psi[0] == Approx(  -0.2618497452));

    SPOSetBuilderFactory::clear();
  }
}

TEST_CASE("ReadMolecularOrbital GTO He","[wavefunction]")
{
  test_He(false);
}

TEST_CASE("ReadMolecularOrbital Numerical He","[wavefunction]")
{
  test_He(true);
}


void test_Ne(bool transform)
{
  std::ostringstream section_name;
  section_name << "Neon, transform orbitals to grid: " << (transform ? "T":"F");

  SECTION(section_name.str()) {
    OHMMS::Controller->initialize(0, NULL);
    Communicate *c = OHMMS::Controller;

    ParticleSet elec;;
    std::vector<int> agroup(2);
    agroup[0] = 1;
    agroup[1] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0]= 0.0;

    SpeciesSet &tspecies = elec.getSpeciesSet();
    int upIdx = tspecies.addSpecies("u");
    int downIdx = tspecies.addSpecies("d");
    int massIdx = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx) = 1.0;
    tspecies(massIdx, downIdx) = 1.0;

    ParticleSet ions;
    ions.setName("ion0");
    ions.create(1);
    ions.R[0][0] = 0.0;
    ions.R[0][1] = 0.0;
    ions.R[0][2] = 0.0;
    SpeciesSet &ispecies = ions.getSpeciesSet();
    int heIdx = ispecies.addSpecies("Ne");
    ions.update();


  #ifdef ENABLE_SOA
    elec.addTable(ions,DT_SOA);
  #else
    elec.addTable(ions,DT_AOS);
  #endif
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("ne_def2_svp.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    TrialWaveFunction psi(c);

    OrbitalBuilderBase::PtclPoolType particle_set_map;
    particle_set_map["e"] = &elec;
    particle_set_map["ion0"] = &ions;


    SPOSetBuilderFactory bf(elec, psi, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (transform) {
      // input file is set to transform GTO's to numerical orbitals by default
    } else {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar *)"transform", (const xmlChar *)"no");
      xmlSetProp(MO_base[0], (const xmlChar *)"key", (const xmlChar *)"GTO");
    }

    SPOSetBuilder *bb = bf.createSPOSetBuilder(MO_base[0]);
    REQUIRE(bb != NULL);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    bb->loadBasisSetFromXML(MO_base[0]);
    SPOSetBase *sposet = bb->createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSetBase::ValueVector_t values;
    SPOSetBase::GradVector_t dpsi;
    SPOSetBase::ValueVector_t d2psi;
    values.resize(5);
    dpsi.resize(5);
    d2psi.resize(5);

    ParticleSet::SingleParticlePos_t newpos(0.00001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluate(elec, 0, values);

    std::cout << "values = " << values << std::endl;

    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    REQUIRE(values[0] == Approx(   -16.11819042));

    sposet->evaluate(elec, 0, values, dpsi, d2psi);
    std::cout << "values = " << values << std::endl;
    std::cout << "dpsi = " << dpsi << std::endl;
    std::cout << "d2psi = " << d2psi << std::endl;
    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    REQUIRE(values[0] == Approx(   -16.11819042));
    REQUIRE(dpsi[0][0] == Approx(   0.1747261458));
    REQUIRE(dpsi[0][1] == Approx(              0));
    REQUIRE(dpsi[0][2] == Approx(              0));


    ParticleSet::SingleParticlePos_t disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);
    sposet->evaluate(elec, 0, values);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-0.005041631374));


    sposet->evaluate(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx(-0.005041631374));
    REQUIRE(dpsi[0][0] == Approx(  0.01862216578));
    REQUIRE(dpsi[0][1] == Approx(              0));
    REQUIRE(dpsi[0][2] == Approx(              0));
    REQUIRE(d2psi[0] == Approx( -0.01551755818));


    SPOSetBuilderFactory::clear();
  }
}

TEST_CASE("ReadMolecularOrbital GTO Ne","[wavefunction]")
{
  test_Ne(false);
}

TEST_CASE("ReadMolecularOrbital Numerical Ne","[wavefunction]")
{
  test_Ne(true);
}

TEST_CASE("ReadMolecularOrbital HCN","[wavefunction]")
{
}

void test_HCN(bool transform)
{
  std::ostringstream section_name;
  section_name << "HCN, transform orbitals to grid: " << (transform ? "T":"F");

  SECTION(section_name.str()) {
    OHMMS::Controller->initialize(0, NULL);
    Communicate *c = OHMMS::Controller;

    Libxml2Document doc;
    bool okay = doc.parse("hcn.structure.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();
    Tensor<int, 3> tmat;
    tmat(0,0) = 1;
    tmat(1,1) = 1;
    tmat(2,2) = 1;

    ParticleSet ions;
    XMLParticleParser parse_ions(ions, tmat);
    OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
    REQUIRE(particleset_ion.size() == 1);
    parse_ions.put(particleset_ion[0]);

    REQUIRE(ions.groups() == 3);
    REQUIRE(ions.R.size() == 3);
    ions.update();

    ParticleSet elec;
    XMLParticleParser parse_elec(elec, tmat);
    OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
    REQUIRE(particleset_elec.size() == 1);
    parse_elec.put(particleset_elec[0]);

    REQUIRE(elec.groups() == 2);
    REQUIRE(elec.R.size() == 14);

    elec.R = 0.0;

  #ifdef ENABLE_SOA
    elec.addTable(ions,DT_SOA);
  #else
    elec.addTable(ions,DT_AOS);
  #endif
    elec.update();

    Libxml2Document doc2;
    okay = doc2.parse("hcn.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root2 = doc2.getRoot();

    TrialWaveFunction psi(c);

    OrbitalBuilderBase::PtclPoolType particle_set_map;
    particle_set_map["e"] = &elec;
    particle_set_map["ion0"] = &ions;


    SPOSetBuilderFactory bf(elec, psi, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (transform) {
      // input file is set to transform GTO's to numerical orbitals by default
    } else {
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], (const xmlChar *)"transform", (const xmlChar *)"no");
      xmlSetProp(MO_base[0], (const xmlChar *)"key", (const xmlChar *)"GTO");
    }

    SPOSetBuilder *bb = bf.createSPOSetBuilder(MO_base[0]);
    REQUIRE(bb != NULL);

    OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
    bb->loadBasisSetFromXML(MO_base[0]);
    SPOSetBase *sposet = bb->createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSetBase::ValueVector_t values;
    SPOSetBase::GradVector_t dpsi;
    SPOSetBase::ValueVector_t d2psi;
    values.resize(7);
    dpsi.resize(7);
    d2psi.resize(7);

    ParticleSet::SingleParticlePos_t newpos;
    elec.makeMove(0, newpos);

    sposet->evaluate(elec, 0, values);

    //typedef LCOrbitalSet<LocalizedBasisSet<SphericalBasisSet<GaussianCombo<double>>>, false> OrbType;
    //OrbType *lcob = dynamic_cast<OrbType *>(sposet);
    //REQUIRE(lcob != NULL);

    //std::cout << "atomic orbs = " << lcob->myBasisSet->Phi << std::endl;
    //std::cout << "values = " << values << std::endl;
    REQUIRE(values[0] == Approx( 0.009452265234));
    REQUIRE(values[1] == Approx(  0.02008357407));
    REQUIRE(values[2] == Approx(   0.4163749594));
    REQUIRE(values[3] == Approx( -0.08854428343));
    REQUIRE(values[4] == Approx(    0.273158705));
    REQUIRE(values[5] == Approx(              0));
    REQUIRE(values[6] == Approx(              0));

    sposet->evaluate(elec, 0, values, dpsi, d2psi);


    // Generated form gen_mo.py for position [0.0, 0.0, 0.0]
    REQUIRE(values[0] == Approx( 0.009452265234));
    REQUIRE(dpsi[0][0] == Approx( -0.05400764372));
    REQUIRE(dpsi[0][1] == Approx(              0));
    REQUIRE(dpsi[0][2] == Approx(              0));
    REQUIRE(d2psi[0] == Approx(   0.2532157143));

    REQUIRE(values[1] == Approx(  0.02008357407));
    REQUIRE(dpsi[1][0] == Approx(   0.1009262252));
    REQUIRE(dpsi[1][1] == Approx(              0));
    REQUIRE(dpsi[1][2] == Approx(              0));
    REQUIRE(d2psi[1] == Approx(   0.3423520138));

    REQUIRE(values[2] == Approx(   0.4163749594));
    REQUIRE(dpsi[2][0] == Approx(  -0.1202256419));
    REQUIRE(dpsi[2][1] == Approx(              0));
    REQUIRE(dpsi[2][2] == Approx(              0));
    REQUIRE(d2psi[2] == Approx(   -1.178149899));

    REQUIRE(values[3] == Approx( -0.08854428343));
    REQUIRE(dpsi[3][0] == Approx(-0.004505552544));
    REQUIRE(dpsi[3][1] == Approx(              0));
    REQUIRE(dpsi[3][2] == Approx(              0));
    REQUIRE(d2psi[3] == Approx(   0.2838238091));

    REQUIRE(values[4] == Approx(    0.273158705));
    REQUIRE(dpsi[4][0] == Approx( -0.01125044248));
    REQUIRE(dpsi[4][1] == Approx(              0));
    REQUIRE(dpsi[4][2] == Approx(              0));
    REQUIRE(d2psi[4] == Approx(  -0.9173261582));

    REQUIRE(values[5] == Approx(              0));
    REQUIRE(dpsi[5][0] == Approx(              0));
    REQUIRE(dpsi[5][1] == Approx(   0.4221165864));
    REQUIRE(dpsi[5][2] == Approx( -0.08191634629));
    REQUIRE(d2psi[5] == Approx(              0));

    REQUIRE(values[6] == Approx(              0));
    REQUIRE(dpsi[6][0] == Approx(              0));
    REQUIRE(dpsi[6][1] == Approx(  0.08191634629));
    REQUIRE(dpsi[6][2] == Approx(   0.4221165864));
    REQUIRE(d2psi[6] == Approx(              0));

    SPOSetBuilderFactory::clear();
  }
}

TEST_CASE("ReadMolecularOrbital GTO HCN","[wavefunction]")
{
  test_HCN(false);
}

TEST_CASE("ReadMolecularOrbital Numerical HCN","[wavefunction]")
{
  test_HCN(true);
}

}


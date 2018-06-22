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

}


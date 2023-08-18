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
#include <ResourceCollection.h>

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
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& elec(*elec_ptr);
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

    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& ions(*ions_ptr);
    ions.setName("ion0");
    ions.create({1});
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

    WaveFunctionComponentBuilder::PSetMap particle_set_map;
    particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
    particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);
    if (!transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], castCharToXMLChar("transform"), castCharToXMLChar("no"));
      xmlSetProp(MO_base[0], castCharToXMLChar("key"), castCharToXMLChar("GTO"));
    }

    const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
    auto& bb(*bb_ptr);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    auto sposet = bb.createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    SPOSet::ValueVector values;
    SPOSet::GradVector dpsi;
    SPOSet::ValueVector d2psi;
    values.resize(1);
    dpsi.resize(1);
    d2psi.resize(1);

    // Call makeMove to compute the distances
    ParticleSet::SingleParticlePos newpos(0.0001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    CHECK(values[0] == Approx(0.9996037001));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);

    // Generated from gen_mo.py for position [0.0001, 0.0, 0.0]
    CHECK(values[0] == Approx(0.9996037001));
    CHECK(dpsi[0][0] == Approx(-0.0006678035459));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));
    CHECK(d2psi[0] == Approx(-20.03410564));


    ParticleSet::SingleParticlePos disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    CHECK(values[0] == Approx(0.2315567641));
    CHECK(dpsi[0][0] == Approx(-0.3805431885));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));
    CHECK(d2psi[0] == Approx(-0.2618497452));
  }
}

TEST_CASE("ReadMolecularOrbital GTO He", "[wavefunction]") { test_He(false); }

TEST_CASE("ReadMolecularOrbital Numerical He", "[wavefunction]") { test_He(true); }

void test_He_mw(bool transform)
{
  // set up ion particle set as normal
  Communicate* c = OHMMS::Controller;

  const SimulationCell simulation_cell;
  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto& elec(*elec_ptr);
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

  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto& ions(*ions_ptr);
  ions.setName("ion0");
  ions.create({1});
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

  WaveFunctionComponentBuilder::PSetMap particle_set_map;
  particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
  particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

  SPOSetBuilderFactory bf(c, elec, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
  REQUIRE(MO_base.size() == 1);
  if (!transform)
  {
    // input file is set to transform GTO's to numerical orbitals by default
    // use direct evaluation of GTO's
    xmlSetProp(MO_base[0], castCharToXMLChar("transform"), castCharToXMLChar("no"));
    xmlSetProp(MO_base[0], castCharToXMLChar("key"), castCharToXMLChar("GTO"));
  }

  const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
  auto& bb(*bb_ptr);

  OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
  auto sposet = bb.createSPOSet(slater_base[0]);

  //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

  SPOSet::ValueVector psi;
  SPOSet::GradVector dpsi;
  SPOSet::ValueVector d2psi;
  psi.resize(1);
  dpsi.resize(1);
  d2psi.resize(1);

  // Call makeMove to compute the distances
  ParticleSet::SingleParticlePos newpos(0.0001, 0.0, 0.0);
  elec.makeMove(0, newpos);
  // set up second walkers
  // auto elec2 = elec.makeClone();

  sposet->evaluateVGL(elec, 0, psi, dpsi, d2psi);
  CHECK(std::real(psi[0]) == Approx(0.9996037001));
  CHECK(std::real(dpsi[0][0]) == Approx(-0.000667803579));
  CHECK(std::real(dpsi[0][1]) == Approx(0));
  CHECK(std::real(dpsi[0][2]) == Approx(0));
  CHECK(std::real(d2psi[0]) == Approx(-20.0342132));


  // vectors of SPOSets, ParticleSets, V/G/L (leading dim of each == nwalkers)
  RefVectorWithLeader<SPOSet> spo_list(*sposet);
  spo_list.push_back(*sposet);

  RefVectorWithLeader<ParticleSet> P_list(elec);
  P_list.push_back(elec);

  RefVector<SPOSet::ValueVector> psi_list;
  RefVector<SPOSet::GradVector> dpsi_list;
  RefVector<SPOSet::ValueVector> d2psi_list;

  // create V,G,L arrays for walker 1
  SPOSet::ValueVector psi_1(sposet->getOrbitalSetSize());
  SPOSet::GradVector dpsi_1(sposet->getOrbitalSetSize());
  SPOSet::ValueVector d2psi_1(sposet->getOrbitalSetSize());

  psi_list.push_back(psi_1);
  dpsi_list.push_back(dpsi_1);
  d2psi_list.push_back(d2psi_1);

  //second walker
  // interchange positions and shift y instead of x
  ParticleSet::SingleParticlePos dy(0.0, 0.0001, 0.0);
  ParticleSet elec_2(elec);
  elec_2.R[0] = elec.R[1];
  elec_2.R[1] = elec.R[0];
  elec_2.update();
  elec_2.makeMove(0, dy);

  std::unique_ptr<SPOSet> sposet_2(sposet->makeClone());
  SPOSet::ValueVector psi_2(sposet->getOrbitalSetSize());
  SPOSet::GradVector dpsi_2(sposet->getOrbitalSetSize());
  SPOSet::ValueVector d2psi_2(sposet->getOrbitalSetSize());
  spo_list.push_back(*sposet_2);
  P_list.push_back(elec_2);
  psi_list.push_back(psi_2);
  dpsi_list.push_back(dpsi_2);
  d2psi_list.push_back(d2psi_2);

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec.createResource(pset_res);
  sposet->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, P_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

  sposet->mw_evaluateVGL(spo_list, P_list, 0, psi_list, dpsi_list, d2psi_list);

  CHECK(std::real(psi_list[0].get()[0]) == Approx(psi[0]));
  CHECK(std::real(dpsi_list[0].get()[0][0]) == Approx(dpsi[0][0]));
  CHECK(std::real(dpsi_list[0].get()[0][1]) == Approx(dpsi[0][1]));
  CHECK(std::real(dpsi_list[0].get()[0][2]) == Approx(dpsi[0][2]));
  CHECK(std::real(d2psi_list[0].get()[0]) == Approx(d2psi[0]));

  CHECK(std::real(psi_list[1].get()[0]) == Approx(psi[0]));
  CHECK(std::real(dpsi_list[1].get()[0][0]) == Approx(dpsi[0][1])); // x, y switched here
  CHECK(std::real(dpsi_list[1].get()[0][1]) == Approx(dpsi[0][0]));
  CHECK(std::real(dpsi_list[1].get()[0][2]) == Approx(dpsi[0][2]));
  CHECK(std::real(d2psi_list[1].get()[0]) == Approx(d2psi[0]));
}

TEST_CASE("mw_evaluate Numerical He", "[wavefunction]") { test_He_mw(true); }

void test_EtOH_mw(bool transform)
{
  // set up ion particle set as normal
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("ethanol.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  const SimulationCell simulation_cell;
  auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto& ions(*ions_ptr);
  XMLParticleParser parse_ions(ions);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.readXML(particleset_ion[0]);

  REQUIRE(ions.groups() == 3);
  REQUIRE(ions.R.size() == 9);
  ions.update();

  auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
  auto& elec(*elec_ptr);
  XMLParticleParser parse_elec(elec);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.readXML(particleset_elec[0]);

  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 26);

  elec.R = 0.0;

  elec.addTable(ions);
  elec.update();

  Libxml2Document doc2;
  okay = doc2.parse("ethanol.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  WaveFunctionComponentBuilder::PSetMap particle_set_map;
  particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
  particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

  SPOSetBuilderFactory bf(c, elec, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);
  if (!transform)
  {
    // input file is set to transform GTO's to numerical orbitals by default
    // use direct evaluation of GTO's
    xmlSetProp(MO_base[0], castCharToXMLChar("transform"), castCharToXMLChar("no"));
    xmlSetProp(MO_base[0], castCharToXMLChar("key"), castCharToXMLChar("GTO"));
  }

  xmlSetProp(MO_base[0], castCharToXMLChar("cuspCorrection"), castCharToXMLChar("no"));

  const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
  auto& bb(*bb_ptr);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  auto sposet = bb.createSPOSet(slater_base[0]);


  //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;
  size_t n_mo = sposet->getOrbitalSetSize();
  SPOSet::ValueVector psiref_0(n_mo);
  SPOSet::GradVector dpsiref_0(n_mo);
  SPOSet::ValueVector d2psiref_0(n_mo);
  // Call makeMove to compute the distances
  //ParticleSet::SingleParticlePos newpos(0.0001, 0.0, 0.0);
  //std::cout << elec.R[0] << std::endl;
  //elec.makeMove(0, newpos);
  //elec.update();
  elec.R[0] = {0.0001, 0.0, 0.0};
  elec.update();
  std::cout << elec.R[0] << std::endl;
  // set up second walkers
  // auto elec2 = elec.makeClone();
  sposet->evaluateVGL(elec, 0, psiref_0, dpsiref_0, d2psiref_0);

  CHECK(std::real(psiref_0[0]) == Approx(-0.001664403313));
  CHECK(std::real(psiref_0[1]) == Approx(0.01579976715));
  CHECK(std::real(dpsiref_0[0][0]) == Approx(-0.0001961749098));
  CHECK(std::real(dpsiref_0[0][1]) == Approx(0.003340392074));
  CHECK(std::real(dpsiref_0[0][2]) == Approx(9.461877818e-05));
  CHECK(std::real(dpsiref_0[1][0]) == Approx(-0.005476152264));
  CHECK(std::real(dpsiref_0[1][1]) == Approx(-0.06648077046));
  CHECK(std::real(dpsiref_0[1][2]) == Approx(2.086541402e-05));
  CHECK(std::real(d2psiref_0[0]) == Approx(0.01071299243));
  CHECK(std::real(d2psiref_0[1]) == Approx(0.4121970776));

  //ParticleSet::SingleParticlePos newpos2(0.0, 0.04, 0.02);
  //elec.makeMove(1, newpos2);
  elec.R[1] = {0.0, 0.04, 0.02};
  elec.update();
  SPOSet::ValueVector psiref_1(n_mo);
  SPOSet::GradVector dpsiref_1(n_mo);
  SPOSet::ValueVector d2psiref_1(n_mo);
  sposet->evaluateVGL(elec, 1, psiref_1, dpsiref_1, d2psiref_1);

  CHECK(std::real(psiref_1[0]) == Approx(-0.001528135727));
  CHECK(std::real(psiref_1[0]) == Approx(-0.001528135727));
  CHECK(std::real(psiref_1[1]) == Approx(0.01351541907));
  CHECK(std::real(d2psiref_1[0]) == Approx(0.01001796854));
  CHECK(std::real(d2psiref_1[1]) == Approx(0.2912963205));
  CHECK(std::real(dpsiref_1[0][0]) == Approx(-0.0004235196101));
  CHECK(std::real(dpsiref_1[0][1]) == Approx(0.003351193375));
  CHECK(std::real(dpsiref_1[0][2]) == Approx(0.0001374796409));
  CHECK(std::real(dpsiref_1[1][0]) == Approx(-0.003873067027));
  CHECK(std::real(dpsiref_1[1][1]) == Approx(-0.0483167767));
  CHECK(std::real(dpsiref_1[1][2]) == Approx(-0.0008320732335));

  // vectors of SPOSets, ParticleSets, V/G/L (leading dim of each == nwalkers)
  RefVectorWithLeader<SPOSet> spo_list(*sposet);
  std::unique_ptr<SPOSet> sposet_2(sposet->makeClone());
  spo_list.push_back(*sposet_2);
  spo_list.push_back(*sposet);
  // second walker
  // exchange positions
  ParticleSet elec_2(elec);
  elec_2.R[0] = elec.R[1];
  elec_2.R[1] = elec.R[0];
  elec_2.update();
  RefVectorWithLeader<ParticleSet> P_list(elec);
  P_list.push_back(elec);
  P_list.push_back(elec_2);
  // create V,G,L arrays for walker 1
  SPOSet::ValueVector psi_1(n_mo);
  SPOSet::GradVector dpsi_1(n_mo);
  SPOSet::ValueVector d2psi_1(n_mo);
  SPOSet::ValueVector psi_2(n_mo);
  SPOSet::GradVector dpsi_2(n_mo);
  SPOSet::ValueVector d2psi_2(n_mo);
  RefVector<SPOSet::ValueVector> psi_list   = {psi_1, psi_2};
  RefVector<SPOSet::GradVector> dpsi_list   = {dpsi_1, dpsi_2};
  RefVector<SPOSet::ValueVector> d2psi_list = {d2psi_1, d2psi_2};

  size_t nw = psi_list.size();
  SPOSet::ValueVector psi_v_1(n_mo);
  SPOSet::ValueVector psi_v_2(n_mo);
  RefVector<SPOSet::ValueVector> psi_v_list{psi_v_1, psi_v_2};

  ResourceCollection pset_res("test_pset_res");
  ResourceCollection spo_res("test_spo_res");

  elec.createResource(pset_res);
  sposet->createResource(spo_res);

  ResourceCollectionTeamLock<ParticleSet> mw_pset_lock(pset_res, P_list);
  ResourceCollectionTeamLock<SPOSet> mw_sposet_lock(spo_res, spo_list);

  sposet->mw_evaluateVGL(spo_list, P_list, 0, psi_list, dpsi_list, d2psi_list);
  sposet->mw_evaluateValue(spo_list, P_list, 0, psi_v_list);

  for (size_t iorb = 0; iorb < n_mo; iorb++)
  {
    for (size_t iw = 0; iw < nw; iw++)
    {
      // test values from OffloadMWVArray impl.
      CHECK(std::real(psi_v_list[iw].get()[iorb]) == Approx(psi_list[iw].get()[iorb]));
    }
    CHECK(std::real(psi_list[0].get()[iorb]) == Approx(psiref_0[iorb]));
    CHECK(std::real(psi_list[1].get()[iorb]) == Approx(psiref_1[iorb]));
    CHECK(std::real(d2psi_list[0].get()[iorb]) == Approx(d2psiref_0[iorb]));
    CHECK(std::real(d2psi_list[1].get()[iorb]) == Approx(d2psiref_1[iorb]));
    for (size_t idim = 0; idim < SPOSet::DIM; idim++)
    {
      CHECK(std::real(dpsi_list[0].get()[iorb][idim]) == Approx(dpsiref_0[iorb][idim]));
      CHECK(std::real(dpsi_list[1].get()[iorb][idim]) == Approx(dpsiref_1[iorb][idim]));
    }
  }
}

TEST_CASE("mw_evaluate Numerical EtOH", "[wavefunction]") { test_EtOH_mw(true); }
TEST_CASE("mw_evaluate GTO EtOH", "[wavefunction]") { test_EtOH_mw(false); }

void test_Ne(bool transform)
{
  std::ostringstream section_name;
  section_name << "Neon, transform orbitals to grid: " << (transform ? "T" : "F");

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    const SimulationCell simulation_cell;
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& elec(*elec_ptr);

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

    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& ions(*ions_ptr);
    ions.setName("ion0");
    ions.create({1});
    ions.R[0]            = {0.0, 0.0, 0.0};
    SpeciesSet& ispecies = ions.getSpeciesSet();
    int heIdx            = ispecies.addSpecies("Ne");
    ions.update();


    elec.addTable(ions);
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("ne_def2_svp.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    WaveFunctionComponentBuilder::PSetMap particle_set_map;
    particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
    particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (!transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], castCharToXMLChar("transform"), castCharToXMLChar("no"));
      xmlSetProp(MO_base[0], castCharToXMLChar("key"), castCharToXMLChar("GTO"));
    }

    const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
    auto& bb(*bb_ptr);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    auto sposet = bb.createSPOSet(slater_base[0]);

    //std::cout << "basis set size = " << sposet->getBasisSetSize() << std::endl;

    const int norbs = 5;
    SPOSet::ValueVector values;
    SPOSet::GradVector dpsi;
    SPOSet::ValueVector d2psi;
    values.resize(norbs);
    dpsi.resize(norbs);
    d2psi.resize(norbs);

    ParticleSet::SingleParticlePos newpos(0.00001, 0.0, 0.0);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    std::cout << "values = " << values << std::endl;

    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    CHECK(values[0] == Approx(-16.11819042));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    std::cout << "values = " << values << std::endl;
    std::cout << "dpsi = " << dpsi << std::endl;
    std::cout << "d2psi = " << d2psi << std::endl;
    // Generated from gen_mo.py for position [1e-05, 0.0, 0.0]
    CHECK(values[0] == Approx(-16.11819042));
    CHECK(dpsi[0][0] == Approx(0.1747261458));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));


    ParticleSet::SingleParticlePos disp(1.0, 0.0, 0.0);
    elec.makeMove(0, disp);
    sposet->evaluateValue(elec, 0, values);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    CHECK(values[0] == Approx(-0.005041631374));


    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);
    // Generated from gen_mo.py for position [1.0, 0.0, 0.0]
    CHECK(values[0] == Approx(-0.005041631374));
    CHECK(dpsi[0][0] == Approx(0.01862216578));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));
    CHECK(d2psi[0] == Approx(-0.01551755818));

    // when a determinant only contains a single particle.
    SPOSet::ValueVector phi(1), phiinv(1);
    phiinv[0] = 100;
    VirtualParticleSet VP(elec, 2);
    std::vector<ParticleSet::SingleParticlePos> newpos2(2);
    std::vector<SPOSet::ValueType> ratios2(2);
    newpos2[0] = disp;
    newpos2[1] = -disp;
    VP.makeMoves(elec, 0, newpos2);
    sposet->evaluateDetRatios(VP, phi, phiinv, ratios2);
    CHECK(ratios2[0] == Approx(-0.504163137)); // values[0] * phiinv[0]
    CHECK(ratios2[1] == Approx(-0.504163137)); // symmetric move
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
    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& ions(*ions_ptr);
    XMLParticleParser parse_ions(ions);
    OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
    REQUIRE(particleset_ion.size() == 1);
    parse_ions.readXML(particleset_ion[0]);

    REQUIRE(ions.groups() == 3);
    REQUIRE(ions.R.size() == 3);
    ions.update();

    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& elec(*elec_ptr);
    XMLParticleParser parse_elec(elec);
    OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
    REQUIRE(particleset_elec.size() == 1);
    parse_elec.readXML(particleset_elec[0]);

    REQUIRE(elec.groups() == 2);
    REQUIRE(elec.R.size() == 14);

    elec.R = 0.0;

    elec.addTable(ions);
    elec.update();

    Libxml2Document doc2;
    okay = doc2.parse("hcn.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root2 = doc2.getRoot();

    WaveFunctionComponentBuilder::PSetMap particle_set_map;
    particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
    particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    if (!transform)
    {
      // input file is set to transform GTO's to numerical orbitals by default
      // use direct evaluation of GTO's
      xmlSetProp(MO_base[0], castCharToXMLChar("transform"), castCharToXMLChar("no"));
      xmlSetProp(MO_base[0], castCharToXMLChar("key"), castCharToXMLChar("GTO"));
    }

    const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
    auto& bb(*bb_ptr);

    OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
    auto sposet = bb.createSPOSet(slater_base[0]);


    SPOSet::ValueVector values;
    SPOSet::GradVector dpsi;
    SPOSet::ValueVector d2psi;
    values.resize(7);
    dpsi.resize(7);
    d2psi.resize(7);

    ParticleSet::SingleParticlePos newpos;
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    CHECK(values[0] == Approx(0.009452265234));
    CHECK(values[1] == Approx(0.02008357407));
    CHECK(values[2] == Approx(0.4163749594));
    CHECK(values[3] == Approx(-0.08854428343));
    CHECK(values[4] == Approx(0.273158705));
    CHECK(values[5] == Approx(0));
    CHECK(values[6] == Approx(0));

    sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);


    // Generated form gen_mo.py for position [0.0, 0.0, 0.0]
    CHECK(values[0] == Approx(0.009452265234));
    CHECK(dpsi[0][0] == Approx(-0.05400764372));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));
    CHECK(d2psi[0] == Approx(0.2532157143));

    CHECK(values[1] == Approx(0.02008357407));
    CHECK(dpsi[1][0] == Approx(0.1009262252));
    CHECK(dpsi[1][1] == Approx(0));
    CHECK(dpsi[1][2] == Approx(0));
    CHECK(d2psi[1] == Approx(0.3423520138));

    CHECK(values[2] == Approx(0.4163749594));
    CHECK(dpsi[2][0] == Approx(-0.1202256419));
    CHECK(dpsi[2][1] == Approx(0));
    CHECK(dpsi[2][2] == Approx(0));
    CHECK(d2psi[2] == Approx(-1.178149899));

    CHECK(values[3] == Approx(-0.08854428343));
    CHECK(dpsi[3][0] == Approx(-0.004505552544));
    CHECK(dpsi[3][1] == Approx(0));
    CHECK(dpsi[3][2] == Approx(0));
    CHECK(d2psi[3] == Approx(0.2838238091));

    CHECK(values[4] == Approx(0.273158705));
    CHECK(dpsi[4][0] == Approx(-0.01125044248));
    CHECK(dpsi[4][1] == Approx(0));
    CHECK(dpsi[4][2] == Approx(0));
    CHECK(d2psi[4] == Approx(-0.9173261582));

    CHECK(values[5] == Approx(0));
    CHECK(dpsi[5][0] == Approx(0));
    CHECK(dpsi[5][1] == Approx(0.4221165864));
    CHECK(dpsi[5][2] == Approx(-0.08191634629));
    CHECK(d2psi[5] == Approx(0));

    CHECK(values[6] == Approx(0));
    CHECK(dpsi[6][0] == Approx(0));
    CHECK(dpsi[6][1] == Approx(0.08191634629));
    CHECK(dpsi[6][2] == Approx(0.4221165864));
    CHECK(d2psi[6] == Approx(0));

    //==========Hessian and Grad Hessian Test==========
    SPOSet::HessVector dhpsi;
    dhpsi.resize(7);

    sposet->evaluateVGH(elec, 0, values, dpsi, dhpsi);

    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    CHECK(values[0] == Approx(0.009452265234));
    CHECK(dpsi[0][0] == Approx(-0.05400764372));
    CHECK(dpsi[0][1] == Approx(0));
    CHECK(dpsi[0][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    CHECK(dhpsi[0](0, 0) == Approx(0.3523924743));
    CHECK(dhpsi[0](0, 1) == Approx(0));
    CHECK(dhpsi[0](0, 2) == Approx(0));
    CHECK(dhpsi[0](1, 1) == Approx(-0.04958838002));
    CHECK(dhpsi[0](1, 2) == Approx(0));
    CHECK(dhpsi[0](2, 2) == Approx(-0.04958838002));

    /////////////////////////////////////////////////////////////////////////////
    //Now we're going to test the API's called by SPOSet for higher derivatives//
    /////////////////////////////////////////////////////////////////////////////
    SPOSet::ValueMatrix psiM(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GradMatrix dpsiM(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::HessMatrix hesspsiV(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GGGMatrix d3psiV(elec.R.size(), sposet->getOrbitalSetSize());

    sposet->evaluate_notranspose(elec, 0, elec.R.size(), psiM, dpsiM, hesspsiV, d3psiV);


    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    CHECK(psiM[0][0] == Approx(0.009452265234));
    CHECK(dpsiM[0][0][0] == Approx(-0.05400764372));
    CHECK(dpsiM[0][0][1] == Approx(0));
    CHECK(dpsiM[0][0][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    CHECK(hesspsiV[0][0](0, 0) == Approx(0.3523924743));
    CHECK(hesspsiV[0][0](0, 1) == Approx(0));
    CHECK(hesspsiV[0][0](0, 2) == Approx(0));
    CHECK(hesspsiV[0][0](1, 1) == Approx(-0.04958838002));
    CHECK(hesspsiV[0][0](1, 2) == Approx(0));
    CHECK(hesspsiV[0][0](2, 2) == Approx(-0.04958838002));

    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    CHECK(d3psiV[0][0][0](0, 0) == Approx(-2.241965465));
    CHECK(d3psiV[0][0][0](0, 1) == Approx(0));
    CHECK(d3psiV[0][0][0](0, 2) == Approx(0));
    CHECK(d3psiV[0][0][0](1, 1) == Approx(0.3714481861));
    CHECK(d3psiV[0][0][0](1, 2) == Approx(0));
    CHECK(d3psiV[0][0][0](2, 2) == Approx(0.3714481861));
    CHECK(d3psiV[0][0][1](1, 1) == Approx(0));
    CHECK(d3psiV[0][0][1](1, 2) == Approx(0));
    CHECK(d3psiV[0][0][1](2, 2) == Approx(0));
    CHECK(d3psiV[0][0][2](2, 2) == Approx(0));

    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    CHECK(psiM[0][1] == Approx(0.02008357407));
    CHECK(dpsiM[0][1][0] == Approx(0.1009262252));
    CHECK(dpsiM[0][1][1] == Approx(0));
    CHECK(dpsiM[0][1][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    CHECK(hesspsiV[0][1](0, 0) == Approx(0.5298289497));
    CHECK(hesspsiV[0][1](0, 1) == Approx(0));
    CHECK(hesspsiV[0][1](0, 2) == Approx(0));
    CHECK(hesspsiV[0][1](1, 1) == Approx(-0.09373846794));
    CHECK(hesspsiV[0][1](1, 2) == Approx(0));
    CHECK(hesspsiV[0][1](2, 2) == Approx(-0.09373846794));

    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    CHECK(d3psiV[0][1][0](0, 0) == Approx(2.594787656));
    CHECK(d3psiV[0][1][0](0, 1) == Approx(0));
    CHECK(d3psiV[0][1][0](0, 2) == Approx(0));
    CHECK(d3psiV[0][1][0](1, 1) == Approx(-0.5720485625));
    CHECK(d3psiV[0][1][0](1, 2) == Approx(0));
    CHECK(d3psiV[0][1][0](2, 2) == Approx(-0.5720485625));
    CHECK(d3psiV[0][1][1](1, 1) == Approx(0));
    CHECK(d3psiV[0][1][1](1, 2) == Approx(0));
    CHECK(d3psiV[0][1][1](2, 2) == Approx(0));
    CHECK(d3psiV[0][1][2](2, 2) == Approx(0));
    // Generated from gen_mo.py for position [0.0, 0.0, 0.0]
    CHECK(psiM[0][2] == Approx(0.4163749594));
    CHECK(dpsiM[0][2][0] == Approx(-0.1202256419));
    CHECK(dpsiM[0][2][1] == Approx(0));
    CHECK(dpsiM[0][2][2] == Approx(0));
    //Hessian (xx,xy,xz,yy,yz,zz)
    CHECK(hesspsiV[0][2](0, 0) == Approx(-0.02607695984));
    CHECK(hesspsiV[0][2](0, 1) == Approx(0));
    CHECK(hesspsiV[0][2](0, 2) == Approx(0));
    CHECK(hesspsiV[0][2](1, 1) == Approx(-0.5760364698));
    CHECK(hesspsiV[0][2](1, 2) == Approx(0));
    CHECK(hesspsiV[0][2](2, 2) == Approx(-0.5760364698));
    //GradHessian (xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz)
    CHECK(d3psiV[0][2][0](0, 0) == Approx(-0.227147312));
    CHECK(d3psiV[0][2][0](0, 1) == Approx(0));
    CHECK(d3psiV[0][2][0](0, 2) == Approx(0));
    CHECK(d3psiV[0][2][0](1, 1) == Approx(0.2992015499));
    CHECK(d3psiV[0][2][0](1, 2) == Approx(0));
    CHECK(d3psiV[0][2][0](2, 2) == Approx(0.2992015499));
    CHECK(d3psiV[0][2][1](1, 1) == Approx(0));
    CHECK(d3psiV[0][2][1](1, 2) == Approx(0));
    CHECK(d3psiV[0][2][1](2, 2) == Approx(0));
    CHECK(d3psiV[0][2][2](2, 2) == Approx(0));


    //Move electron 0 to some nontrivial position:
    ParticleSet::SingleParticlePos disp(0.02, -0.1, 0.05);
    elec.makeMove(0, disp);


    SPOSet::GradMatrix dionpsi(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::HessMatrix diongradpsi(elec.R.size(), sposet->getOrbitalSetSize());
    SPOSet::GradMatrix dionlaplpsi(elec.R.size(), sposet->getOrbitalSetSize());

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 0, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  0  Component  0 ===================
    CHECK(dionpsi[0][0][0] == Approx(0.0453112082));
    CHECK(diongradpsi[0][0](0, 0) == Approx(-0.2943513994));
    CHECK(diongradpsi[0][0](0, 1) == Approx(0.030468047));
    CHECK(diongradpsi[0][0](0, 2) == Approx(-0.0152340235));
    CHECK(dionlaplpsi[0][0][0] == Approx(1.333755581));
    CHECK(dionpsi[0][1][0] == Approx(-0.0006473819623));
    CHECK(diongradpsi[0][1](0, 0) == Approx(0.0004713407512));
    CHECK(diongradpsi[0][1](0, 1) == Approx(-0.0001254975603));
    CHECK(diongradpsi[0][1](0, 2) == Approx(6.274878013e-05));
    CHECK(dionlaplpsi[0][1][0] == Approx(0.001057940846));
    CHECK(dionpsi[0][2][0] == Approx(0.265578336));
    CHECK(diongradpsi[0][2](0, 0) == Approx(-0.08685804115));
    CHECK(diongradpsi[0][2](0, 1) == Approx(0.05438178417));
    CHECK(diongradpsi[0][2](0, 2) == Approx(-0.02719089209));
    CHECK(dionlaplpsi[0][2][0] == Approx(-1.882489819));
    CHECK(dionpsi[0][3][0] == Approx(-0.06444305979));
    CHECK(diongradpsi[0][3](0, 0) == Approx(-0.002013151923));
    CHECK(diongradpsi[0][3](0, 1) == Approx(-0.002535923431));
    CHECK(diongradpsi[0][3](0, 2) == Approx(0.001267961716));
    CHECK(dionlaplpsi[0][3][0] == Approx(0.4547401581));
    CHECK(dionpsi[0][4][0] == Approx(0.1454357726));
    CHECK(diongradpsi[0][4](0, 0) == Approx(-0.2330499431));
    CHECK(diongradpsi[0][4](0, 1) == Approx(0.09667641762));
    CHECK(diongradpsi[0][4](0, 2) == Approx(-0.04833820881));
    CHECK(dionlaplpsi[0][4][0] == Approx(-0.9197558839));
    CHECK(dionpsi[0][5][0] == Approx(-0.04329985085));
    CHECK(diongradpsi[0][5](0, 0) == Approx(0.09051993304));
    CHECK(diongradpsi[0][5](0, 1) == Approx(0.382375474));
    CHECK(diongradpsi[0][5](0, 2) == Approx(-0.07043361927));
    CHECK(dionlaplpsi[0][5][0] == Approx(0.2201672051));
    CHECK(dionpsi[0][6][0] == Approx(0.01207541177));
    CHECK(diongradpsi[0][6](0, 0) == Approx(-0.02524405435));
    CHECK(diongradpsi[0][6](0, 1) == Approx(0.0800332842));
    CHECK(diongradpsi[0][6](0, 2) == Approx(0.3929818664));
    CHECK(dionlaplpsi[0][6][0] == Approx(-0.0614000824));


    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 1, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  1  Component  1 ===================
    CHECK(dionpsi[0][0][1] == Approx(0.0001412373768));
    CHECK(diongradpsi[0][0](1, 0) == Approx(0.0001124265646));
    CHECK(diongradpsi[0][0](1, 1) == Approx(-0.001383378615));
    CHECK(diongradpsi[0][0](1, 2) == Approx(-1.449757545e-05));
    CHECK(dionlaplpsi[0][0][1] == Approx(-0.001252043663));
    CHECK(dionpsi[0][1][1] == Approx(-0.01029290716));
    CHECK(diongradpsi[0][1](1, 0) == Approx(-0.06178485148));
    CHECK(diongradpsi[0][1](1, 1) == Approx(0.0971577216));
    CHECK(diongradpsi[0][1](1, 2) == Approx(0.002885675005));
    CHECK(dionlaplpsi[0][1][1] == Approx(-0.1403103458));
    CHECK(dionpsi[0][2][1] == Approx(-0.0230872583));
    CHECK(diongradpsi[0][2](1, 0) == Approx(-0.02537847709));
    CHECK(diongradpsi[0][2](1, 1) == Approx(0.2268946564));
    CHECK(diongradpsi[0][2](1, 2) == Approx(0.001988963201));
    CHECK(dionlaplpsi[0][2][1] == Approx(0.2028851421));
    CHECK(dionpsi[0][3][1] == Approx(0.01850231814));
    CHECK(diongradpsi[0][3](1, 0) == Approx(0.05709948475));
    CHECK(diongradpsi[0][3](1, 1) == Approx(-0.1776515965));
    CHECK(diongradpsi[0][3](1, 2) == Approx(-0.003685792479));
    CHECK(dionlaplpsi[0][3][1] == Approx(-0.1280699725));
    CHECK(dionpsi[0][4][1] == Approx(-0.02136209962));
    CHECK(diongradpsi[0][4](1, 0) == Approx(-0.03836586276));
    CHECK(diongradpsi[0][4](1, 1) == Approx(0.2084578148));
    CHECK(diongradpsi[0][4](1, 2) == Approx(0.002581590766));
    CHECK(dionlaplpsi[0][4][1] == Approx(0.1792683544));
    CHECK(dionpsi[0][5][1] == Approx(-0.1942343714));
    CHECK(diongradpsi[0][5](1, 0) == Approx(-0.3037357197));
    CHECK(diongradpsi[0][5](1, 1) == Approx(-0.09561978734));
    CHECK(diongradpsi[0][5](1, 2) == Approx(0.02118492506));
    CHECK(dionlaplpsi[0][5][1] == Approx(0.6410434658));
    CHECK(dionpsi[0][6][1] == Approx(-0.03930992259));
    CHECK(diongradpsi[0][6](1, 0) == Approx(-0.06331544695));
    CHECK(diongradpsi[0][6](1, 1) == Approx(-0.002807368817));
    CHECK(diongradpsi[0][6](1, 2) == Approx(-0.02801340823));
    CHECK(dionlaplpsi[0][6][1] == Approx(0.1369061053));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 2, dionpsi, diongradpsi, dionlaplpsi);

    //============== Ion  2  Component  2 ===================
    CHECK(dionpsi[0][0][2] == Approx(1.302648961e-06));
    CHECK(diongradpsi[0][0](2, 0) == Approx(1.865129579e-06));
    CHECK(diongradpsi[0][0](2, 1) == Approx(6.142092043e-08));
    CHECK(diongradpsi[0][0](2, 2) == Approx(2.602225618e-05));
    CHECK(dionlaplpsi[0][0][2] == Approx(1.234692903e-06));
    CHECK(dionpsi[0][1][2] == Approx(3.248738084e-07));
    CHECK(diongradpsi[0][1](2, 0) == Approx(-2.044420189e-06));
    CHECK(diongradpsi[0][1](2, 1) == Approx(-7.011145137e-08));
    CHECK(diongradpsi[0][1](2, 2) == Approx(6.532522353e-06));
    CHECK(dionlaplpsi[0][1][2] == Approx(-6.10958506e-06));
    CHECK(dionpsi[0][2][2] == Approx(3.264249981e-06));
    CHECK(diongradpsi[0][2](2, 0) == Approx(2.820971234e-05));
    CHECK(diongradpsi[0][2](2, 1) == Approx(9.405184964e-07));
    CHECK(diongradpsi[0][2](2, 2) == Approx(6.481420782e-05));
    CHECK(dionlaplpsi[0][2][2] == Approx(5.73961989e-05));
    CHECK(dionpsi[0][3][2] == Approx(0.0001288974413));
    CHECK(diongradpsi[0][3](2, 0) == Approx(0.0002840756879));
    CHECK(diongradpsi[0][3](2, 1) == Approx(9.281700408e-06));
    CHECK(diongradpsi[0][3](2, 2) == Approx(0.002573308008));
    CHECK(dionlaplpsi[0][3][2] == Approx(0.0003025314443));
    CHECK(dionpsi[0][4][2] == Approx(-7.300043903e-05));
    CHECK(diongradpsi[0][4](2, 0) == Approx(-0.0001000016834));
    CHECK(diongradpsi[0][4](2, 1) == Approx(-3.233243534e-06));
    CHECK(diongradpsi[0][4](2, 2) == Approx(-0.001458391774));
    CHECK(dionlaplpsi[0][4][2] == Approx(-3.546690719e-05));
    CHECK(dionpsi[0][5][2] == Approx(2.910525987e-06));
    CHECK(diongradpsi[0][5](2, 0) == Approx(1.307065133e-05));
    CHECK(diongradpsi[0][5](2, 1) == Approx(1.560390706e-06));
    CHECK(diongradpsi[0][5](2, 2) == Approx(-2.92731811e-06));
    CHECK(dionlaplpsi[0][5][2] == Approx(3.797816228e-05));
    CHECK(dionpsi[0][6][2] == Approx(-1.56074936e-05));
    CHECK(diongradpsi[0][6](2, 0) == Approx(-7.009049656e-05));
    CHECK(diongradpsi[0][6](2, 1) == Approx(-2.048666792e-06));
    CHECK(diongradpsi[0][6](2, 2) == Approx(2.967709412e-06));
    CHECK(dionlaplpsi[0][6][2] == Approx(-0.0002018111858));

    //Same tests as before, but for the gradient only.

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 0, dionpsi);
    //============== Ion  0  Component  0 ===================
    CHECK(dionpsi[0][0][0] == Approx(0.0453112082));
    CHECK(dionpsi[0][1][0] == Approx(-0.0006473819623));
    CHECK(dionpsi[0][2][0] == Approx(0.265578336));
    CHECK(dionpsi[0][3][0] == Approx(-0.06444305979));
    CHECK(dionpsi[0][4][0] == Approx(0.1454357726));
    CHECK(dionpsi[0][5][0] == Approx(-0.04329985085));
    CHECK(dionpsi[0][6][0] == Approx(0.01207541177));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 1, dionpsi);
    //============== Ion  1  Component  1 ===================
    CHECK(dionpsi[0][0][1] == Approx(0.0001412373768));
    CHECK(dionpsi[0][1][1] == Approx(-0.01029290716));
    CHECK(dionpsi[0][2][1] == Approx(-0.0230872583));
    CHECK(dionpsi[0][3][1] == Approx(0.01850231814));
    CHECK(dionpsi[0][4][1] == Approx(-0.02136209962));
    CHECK(dionpsi[0][5][1] == Approx(-0.1942343714));
    CHECK(dionpsi[0][6][1] == Approx(-0.03930992259));

    sposet->evaluateGradSource(elec, 0, elec.R.size(), ions, 2, dionpsi);
    //============== Ion  2  Component  2 ===================
    CHECK(dionpsi[0][0][2] == Approx(1.302648961e-06));
    CHECK(dionpsi[0][1][2] == Approx(3.248738084e-07));
    CHECK(dionpsi[0][2][2] == Approx(3.264249981e-06));
    CHECK(dionpsi[0][3][2] == Approx(0.0001288974413));
    CHECK(dionpsi[0][4][2] == Approx(-7.300043903e-05));
    CHECK(dionpsi[0][5][2] == Approx(2.910525987e-06));
    CHECK(dionpsi[0][6][2] == Approx(-1.56074936e-05));


    //Finally, going to test evaluateGradSourceRow.  Same template and reference
    //values as above.
    SPOSet::GradVector dionpsivec;
    dionpsivec.resize(7);

    sposet->evaluateGradSourceRow(elec, 0, ions, 0, dionpsivec);
    //============== Ion  0  Component  0 ===================
    CHECK(dionpsivec[0][0] == Approx(0.0453112082));
    CHECK(dionpsivec[1][0] == Approx(-0.0006473819623));
    CHECK(dionpsivec[2][0] == Approx(0.265578336));
    CHECK(dionpsivec[3][0] == Approx(-0.06444305979));
    CHECK(dionpsivec[4][0] == Approx(0.1454357726));
    CHECK(dionpsivec[5][0] == Approx(-0.04329985085));
    CHECK(dionpsivec[6][0] == Approx(0.01207541177));

    sposet->evaluateGradSourceRow(elec, 0, ions, 1, dionpsivec);
    //============== Ion  1  Component  1 ===================
    CHECK(dionpsivec[0][1] == Approx(0.0001412373768));
    CHECK(dionpsivec[1][1] == Approx(-0.01029290716));
    CHECK(dionpsivec[2][1] == Approx(-0.0230872583));
    CHECK(dionpsivec[3][1] == Approx(0.01850231814));
    CHECK(dionpsivec[4][1] == Approx(-0.02136209962));
    CHECK(dionpsivec[5][1] == Approx(-0.1942343714));
    CHECK(dionpsivec[6][1] == Approx(-0.03930992259));

    sposet->evaluateGradSourceRow(elec, 0, ions, 2, dionpsivec);
    //============== Ion  2  Component  2 ===================
    CHECK(dionpsivec[0][2] == Approx(1.302648961e-06));
    CHECK(dionpsivec[1][2] == Approx(3.248738084e-07));
    CHECK(dionpsivec[2][2] == Approx(3.264249981e-06));
    CHECK(dionpsivec[3][2] == Approx(0.0001288974413));
    CHECK(dionpsivec[4][2] == Approx(-7.300043903e-05));
    CHECK(dionpsivec[5][2] == Approx(2.910525987e-06));
    CHECK(dionpsivec[6][2] == Approx(-1.56074936e-05));
  }
}

TEST_CASE("ReadMolecularOrbital GTO HCN", "[wavefunction]") { test_HCN(false); }

TEST_CASE("ReadMolecularOrbital Numerical HCN", "[wavefunction]") { test_HCN(true); }

} // namespace qmcplusplus

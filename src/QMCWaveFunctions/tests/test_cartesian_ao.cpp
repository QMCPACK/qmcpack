//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
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
void test_cartesian_ao()
{
  std::ostringstream section_name;
  section_name << "Cartesian AO ordering";

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    const SimulationCell simulation_cell;
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& elec(*elec_ptr);
    std::vector<int> agroup(2);
    agroup[0] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0] = 0.0;

    SpeciesSet& tspecies     = elec.getSpeciesSet();
    int upIdx                = tspecies.addSpecies("u");
    int massIdx              = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx) = 1.0;

    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& ions(*ions_ptr);
    ions.setName("ion0");
    ions.create({1});
    ions.R[0]            = 0.0;
    SpeciesSet& ispecies = ions.getSpeciesSet();
    int hIdx             = ispecies.addSpecies("H");
    ions.update();

    elec.addTable(ions);
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("cartesian_order.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    WaveFunctionComponentBuilder::PSetMap particle_set_map;
    particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
    particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));

    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
    auto& bb(*bb_ptr);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    auto sposet = bb.createSPOSet(slater_base[0]);

    SPOSet::ValueVector values;
    values.resize(1);

    // Call makeMove to compute the distances
    ParticleSet::SingleParticlePos newpos(0.1, -0.3, 0.2);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    //generated from ao_order_test.py
    REQUIRE(values[0] == Approx(0.48224527310155046).epsilon(1E-6));
  }
}

void test_dirac_ao()
{
  std::ostringstream section_name;
  section_name << "Dirac AO ordering";

  SECTION(section_name.str())
  {
    Communicate* c = OHMMS::Controller;

    const SimulationCell simulation_cell;
    auto elec_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& elec(*elec_ptr);
    std::vector<int> agroup(2);
    agroup[0] = 1;
    elec.setName("e");
    elec.create(agroup);
    elec.R[0] = 0.0;

    SpeciesSet& tspecies     = elec.getSpeciesSet();
    int upIdx                = tspecies.addSpecies("u");
    int massIdx              = tspecies.addAttribute("mass");
    tspecies(massIdx, upIdx) = 1.0;

    auto ions_ptr = std::make_unique<ParticleSet>(simulation_cell);
    auto& ions(*ions_ptr);
    ions.setName("ion0");
    ions.create({1});
    ions.R[0]            = 0.0;
    SpeciesSet& ispecies = ions.getSpeciesSet();
    int hIdx             = ispecies.addSpecies("H");
    ions.update();

    elec.addTable(ions);
    elec.update();

    Libxml2Document doc;
    bool okay = doc.parse("dirac_order.wfnoj.xml");
    REQUIRE(okay);
    xmlNodePtr root = doc.getRoot();

    WaveFunctionComponentBuilder::PSetMap particle_set_map;
    particle_set_map.emplace(elec_ptr->getName(), std::move(elec_ptr));
    particle_set_map.emplace(ions_ptr->getName(), std::move(ions_ptr));


    SPOSetBuilderFactory bf(c, elec, particle_set_map);

    OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
    REQUIRE(MO_base.size() == 1);

    const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
    auto& bb(*bb_ptr);

    OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
    auto sposet = bb.createSPOSet(slater_base[0]);

    SPOSet::ValueVector values;
    values.resize(1);

    // Call makeMove to compute the distances
    ParticleSet::SingleParticlePos newpos(0.1, -0.3, 0.2);
    elec.makeMove(0, newpos);

    sposet->evaluateValue(elec, 0, values);

    //from ao_order_test.py
    REQUIRE(values[0] == Approx(0.35953790416302006).epsilon(1E-6));
  }
}

TEST_CASE("Cartesian Gaussian Ordering", "[wavefunction]") { test_cartesian_ao(); }
TEST_CASE("Dirac Cartesian Gaussian Ordering", "[wavefunction]") { test_dirac_ao(); }

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "LCAO/LCAOrbitalBuilder.h"
#include <ResourceCollection.h>

#include <stdio.h>
#include <string>
#include <limits>

using std::string;

namespace qmcplusplus
{

TEST_CASE("LCAO DiamondC_2x1x1", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  ParticleSet::ParticleLayout lattice;
  lattice.R = {6.7463223, 6.7463223, 0.0, 0.0, 3.37316115, 3.37316115, 3.37316115, 0.0, 3.37316115};


  ParticleSetPool ptcl = ParticleSetPool(c);
  ptcl.setSimulationCell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  auto elec_uptr = std::make_unique<ParticleSet>(ptcl.getSimulationCell());
  ParticleSet& ions_(*ions_uptr);
  ParticleSet& elec_(*elec_uptr);

  ions_.setName("ion0");
  ptcl.addParticleSet(std::move(ions_uptr));
  ions_.create({4});
  ions_.R[0]           = {0.0, 0.0, 0.0};
  ions_.R[1]           = {1.686580575, 1.686580575, 1.686580575};
  ions_.R[2]           = {3.37316115, 3.37316115, 0.0};
  ions_.R[3]           = {5.059741726, 5.059741726, 1.686580575};
  SpeciesSet& ispecies = ions_.getSpeciesSet();
  int Cidx             = ispecies.addSpecies("C");

  ions_.print(app_log());

  elec_.setName("elec");
  ptcl.addParticleSet(std::move(elec_uptr));
  elec_.create({8, 8});
  elec_.R[0] = {0.0, 1.0, 0.0};
  elec_.R[1] = {0.0, 1.1, 0.0};
  elec_.R[2] = {0.0, 1.2, 0.0};
  elec_.R[3] = {0.0, 1.3, 0.0};

  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  // const int ei_table_index = elec_.addTable(ions_);

  // diamondC_2x1x1
  // from tests/solids/diamondC_2x1x1-Gaussian_pp/C_Diamond-Twist0.wfj.xml
  const std::string wf_xml_str = R"(
  <wavefunction name="psi0" target="e">
      <sposet_collection type="MolecularOrbital" name="LCAOBSet" source="ion0" twist="0  0  0" href="C_Diamond.h5" PBCimages="5  5  5">
      <basisset name="LCAOBSet" key="GTO" transform="no">
        <grid type="log" ri="1.e-6" rf="1.e2" npts="1001"/>
      </basisset>
      <sposet basisset="LCAOBSet" name="spo-up" size="8">
        <occupation mode="ground"/>
        <coefficient size="116" spindataset="0"/>
      </sposet>
      <sposet basisset="LCAOBSet" name="spo-dn" size="8">
        <occupation mode="ground"/>
        <coefficient size="116" spindataset="0"/>
      </sposet>
    </sposet_collection>
  </wavefunction>
  )";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml_str);
  REQUIRE(okay);

  xmlNodePtr root    = doc.getRoot();
  xmlNodePtr spo_xml = xmlFirstElementChild(root);

  LCAOrbitalBuilder lcaoSet(elec_, ions_, c, spo_xml);
  auto spo = lcaoSet.createSPOSetFromXML(spo_xml);
  REQUIRE(spo);
  auto& lcao_spos = dynamic_cast<const LCAOrbitalSet&>(*spo);
  REQUIRE(!lcao_spos.isIdentity());
}
} // namespace qmcplusplus

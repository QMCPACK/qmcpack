//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"

#include "LCAO/LCAOrbitalBuilder.h"
#include "LCAO/MultiQuinticSpline1D.h"
#include "Numerics/SoaCartesianTensor.h"
#include "Numerics/SoaSphericalTensor.h"
#include "LCAO/SoaAtomicBasisSet.h"
#include "LCAO/SoaLocalizedBasisSet.h"
#include "LCAO/AOBasisBuilder.h"
#include "LCAO/MultiFunctorAdapter.h"


namespace qmcplusplus
{


TEST_CASE("LCAOrbitalBuilder", "[wavefunction][LCAO]")
{
  Communicate* c = OHMMS::Controller;

  using Real      = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;


  const SimulationCell sim_cell;
  ParticleSet elec(sim_cell);
  elec.setName("e");
  ParticleSet ions(sim_cell);
  ions.setName("ion0");
  ions.create({1});
  SpeciesSet& source_species(ions.getSpeciesSet());
  source_species.addSpecies("C");
  ions.update();


  // Radial basis: Numerical (Transform to grid)    Angular part: Cartesian
  auto wf_xml_num_cart = R"(
    <tmp>
     <basisset keyword="STO" transform="yes">
        <atomicBasisSet angular="cartesian" elementType="C" normalized="no" type="STO">
          <basisGroup l="0" m="0" n="0" rid="1S1" type="Slater">
            <radfunc contraction="1.0" exponent="11.4" n="1"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  Libxml2Document doc;
  bool okay = doc.parseFromString(wf_xml_num_cart);
  REQUIRE(okay);

  using BasisSet_t = LCAOrbitalBuilder::BasisSet_t;
  LCAOrbitalBuilder lcaob_num_cart(elec, ions, c, doc.getRoot());
  const auto& bs = lcaob_num_cart.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<
            SoaLocalizedBasisSet<SoaAtomicBasisSet<MultiQuinticSpline1D<Real>, SoaCartesianTensor<Real>>, ValueType>*>(
            bs.get()) != nullptr);


  // Radial basis: Numerical (Transform to grid)    Angular part: Spherical
  auto wf_xml_num_sph = R"(
    <tmp>
     <basisset keyword="STO" transform="yes">
        <atomicBasisSet angular="spherical" elementType="C" normalized="no" type="STO">
          <basisGroup l="0" m="0" n="0" rid="1S1" type="Slater">
            <radfunc contraction="1.0" exponent="11.4" n="1"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  okay = doc.parseFromString(wf_xml_num_sph);
  REQUIRE(okay);

  LCAOrbitalBuilder lcaob_num_sph(elec, ions, c, doc.getRoot());
  const auto& bs2 = lcaob_num_sph.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<
            SoaLocalizedBasisSet<SoaAtomicBasisSet<MultiQuinticSpline1D<Real>, SoaSphericalTensor<Real>>, ValueType>*>(
            bs2.get()) != nullptr);


  // Radial basis: GTO    Angular part: Cartesian
  auto wf_xml_gto_cart = R"(
    <tmp>
      <basisset keyword="GTO" name="LCAOBSet">
        <atomicBasisSet angular="cartesian" elementType="C" name="Gaussian-G2" normalized="no" type="Gaussian">
          <grid npts="1001" rf="1.e2" ri="1.e-6" type="log" />
          <basisGroup l="0" n="0" rid="H00" type="Gaussian">
            <radfunc contraction="0.002006000801" exponent="82.64" />
            <radfunc contraction="0.01534300612" exponent="12.41" />
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  okay = doc.parseFromString(wf_xml_gto_cart);
  REQUIRE(okay);

  LCAOrbitalBuilder lcaob_gto_cart(elec, ions, c, doc.getRoot());
  const auto& bs3 = lcaob_gto_cart.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<SoaLocalizedBasisSet<
            SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<Real>>, SoaCartesianTensor<Real>>, ValueType>*>(
            bs3.get()) != nullptr);


  // Radial basis: GTO    Angular part: Spherical
  auto wf_xml_gto_sph = R"(
    <tmp>
      <basisset keyword="GTO" name="LCAOBSet">
        <atomicBasisSet angular="spherical" elementType="C" name="Gaussian-G2" normalized="no" type="Gaussian">
          <grid npts="1001" rf="1.e2" ri="1.e-6" type="log" />
          <basisGroup l="0" n="0" rid="H00" type="Gaussian">
            <radfunc contraction="0.002006000801" exponent="82.64" />
            <radfunc contraction="0.01534300612" exponent="12.41" />
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  okay = doc.parseFromString(wf_xml_gto_sph);
  REQUIRE(okay);

  LCAOrbitalBuilder lcaob_gto_sph(elec, ions, c, doc.getRoot());
  const auto& bs4 = lcaob_gto_sph.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<SoaLocalizedBasisSet<
            SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<Real>>, SoaSphericalTensor<Real>>, ValueType>*>(
            bs4.get()) != nullptr);


  // Radial basis: STO    Angular part: Cartesian
  auto wf_xml_sto_cart = R"(
    <tmp>
     <basisset keyword="STO" transform="no">
        <atomicBasisSet angular="cartesian" elementType="C" normalized="no" type="STO">
          <basisGroup l="0" m="0" n="0" rid="1S1" type="Slater">
            <radfunc contraction="1.0" exponent="11.4" n="1"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  okay = doc.parseFromString(wf_xml_sto_cart);
  REQUIRE(okay);

  LCAOrbitalBuilder lcaob_sto_cart(elec, ions, c, doc.getRoot());
  const auto& bs5 = lcaob_sto_cart.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<SoaLocalizedBasisSet<
            SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<Real>>, SoaCartesianTensor<Real>>, ValueType>*>(
            bs5.get()) != nullptr);


  // Radial basis: STO    Angular part: Spherical
  auto wf_xml_sto_sph = R"(
    <tmp>
     <basisset keyword="STO" transform="no">
        <atomicBasisSet angular="spherical" elementType="C" normalized="no" type="STO">
          <basisGroup l="0" m="0" n="0" rid="1S1" type="Slater">
            <radfunc contraction="1.0" exponent="11.4" n="1"/>
          </basisGroup>
        </atomicBasisSet>
      </basisset>
    </tmp>
  )";

  okay = doc.parseFromString(wf_xml_sto_sph);
  REQUIRE(okay);

  LCAOrbitalBuilder lcaob_sto_sph(elec, ions, c, doc.getRoot());
  const auto& bs6 = lcaob_sto_sph.getBasissetMap().at("LCAOBSet");
  CHECK(dynamic_cast<SoaLocalizedBasisSet<
            SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<Real>>, SoaSphericalTensor<Real>>, ValueType>*>(
            bs6.get()) != nullptr);
}

} // namespace qmcplusplus


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

#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/lcao/CuspCorrection.h"

#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
TEST_CASE("readCuspInfo", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  typedef OneDimGridBase<double> GridType;

  Matrix<CuspCorrectionParameters> info;
  int num_center       = 3;
  int orbital_set_size = 7;
  info.resize(num_center, orbital_set_size);

  bool okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);
  REQUIRE(okay);

  // N
  REQUIRE(info(0, 0).redo == Approx(0.0));                   // redo
  REQUIRE(info(0, 0).C == Approx(0.0));                      // C
  REQUIRE(info(0, 0).sg == Approx(1.0));                     // sg
  REQUIRE(info(0, 0).Rc == Approx(0.0769130700800000));      // rc
  REQUIRE(info(0, 0).alpha[0] == Approx(2.29508580995773));  // a1
  REQUIRE(info(0, 0).alpha[1] == Approx(-7.00028778782666)); // a2
  REQUIRE(info(0, 0).alpha[2] == Approx(0.834942828252775)); // a3
  REQUIRE(info(0, 0).alpha[3] == Approx(-4.61597420905980)); // a4
  REQUIRE(info(0, 0).alpha[4] == Approx(31.6558091872316));  // a5

  // Spot check a few values from these centers
  // C
  REQUIRE(info(0, 6).C == Approx(0.0));        // C
  REQUIRE(info(0, 6).alpha[4] == Approx(0.0)); // a5

  // H
  REQUIRE(info(2, 4).alpha[4] == Approx(-404.733151049101)); // a5
}


TEST_CASE("applyCuspInfo", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("hcn.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  Tensor<int, 3> tmat;
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

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

  elec.addTable(ions, DT_SOA);
  elec.update();

  Libxml2Document doc2;
  okay = doc2.parse("hcn.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet* sposet = bb->createSPOSet(slater_base[0]);

  LCAOrbitalSet* lcob = dynamic_cast<LCAOrbitalSet*>(sposet);
  REQUIRE(lcob != NULL);


  LCAOrbitalSet phi = LCAOrbitalSet(lcob->myBasisSet);
  phi.setOrbitalSetSize(lcob->OrbitalSetSize);
  phi.BasisSetSize = lcob->BasisSetSize;
  phi.setIdentity(false);

  LCAOrbitalSet eta = LCAOrbitalSet(lcob->myBasisSet);
  eta.setOrbitalSetSize(lcob->OrbitalSetSize);
  eta.BasisSetSize = lcob->BasisSetSize;
  eta.setIdentity(false);

  *(eta.C) = *(lcob->C);
  *(phi.C) = *(lcob->C);


  int num_center = 3;
  std::vector<bool> corrCenter(num_center, "true");

  // N is first atom
  int center_idx = 0;

  typedef QMCTraits::RealType RealType;

  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(1.00180500));
  CHECK((*eta.C)(0, 0) == Approx(0.0));

  int orbital_set_size = 7;
  Matrix<CuspCorrectionParameters> info;
  info.resize(num_center, orbital_set_size);
  okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);
  REQUIRE(okay);

  *(eta.C) = *(lcob->C);
  *(phi.C) = *(lcob->C);

  // C is second atom
  center_idx = 1;
  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(0.0));
  CHECK((*eta.C)(0, 0) == Approx(1.00180500));


  removeSTypeOrbitals(corrCenter, *lcob);

  CHECK((*lcob->C)(0, 0) == Approx(0.0));
  CHECK((*lcob->C)(0, 1) == Approx(0.0));
  CHECK((*lcob->C)(0, 2) == Approx(0.0));
  CHECK((*lcob->C)(0, 3) != 0.0);

  SPOSetBuilderFactory::clear();
}

} // namespace qmcplusplus

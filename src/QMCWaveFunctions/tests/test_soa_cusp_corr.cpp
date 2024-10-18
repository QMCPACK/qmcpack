
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

#include "QMCWaveFunctions/LCAO/LCAOrbitalSet.h"
#include "QMCWaveFunctions/LCAO/CuspCorrectionConstruction.h"

#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
TEST_CASE("readCuspInfo", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;

  using GridType = OneDimGridBase<double>;

  Matrix<CuspCorrectionParameters> info;
  int num_center       = 3;
  int orbital_set_size = 7;
  info.resize(num_center, orbital_set_size);

  bool okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);
  REQUIRE(okay);

  // N
  CHECK(info(0, 0).redo == Approx(0.0));                   // redo
  CHECK(info(0, 0).C == Approx(0.0));                      // C
  CHECK(info(0, 0).sg == Approx(1.0));                     // sg
  CHECK(info(0, 0).Rc == Approx(0.0769130700800000));      // rc
  CHECK(info(0, 0).alpha[0] == Approx(2.29508580995773));  // a1
  CHECK(info(0, 0).alpha[1] == Approx(-7.00028778782666)); // a2
  CHECK(info(0, 0).alpha[2] == Approx(0.834942828252775)); // a3
  CHECK(info(0, 0).alpha[3] == Approx(-4.61597420905980)); // a4
  CHECK(info(0, 0).alpha[4] == Approx(31.6558091872316));  // a5

  // Spot check a few values from these centers
  // C
  CHECK(info(0, 6).C == Approx(0.0));        // C
  CHECK(info(0, 6).alpha[4] == Approx(0.0)); // a5

  // H
  CHECK(info(2, 4).alpha[4] == Approx(-404.733151049101)); // a5
}


TEST_CASE("applyCuspInfo", "[wavefunction]")
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

  const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
  auto& bb(*bb_ptr);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  auto sposet = bb.createSPOSet(slater_base[0]);

  LCAOrbitalSet& lcob = dynamic_cast<LCAOrbitalSet&>(*sposet);

  LCAOrbitalSet phi("phi", std::unique_ptr<LCAOrbitalSet::basis_type>(lcob.myBasisSet->makeClone()),
                    lcob.getOrbitalSetSize(), lcob.isIdentity(), lcob.isOMPoffload());

  LCAOrbitalSet eta("eta", std::unique_ptr<LCAOrbitalSet::basis_type>(lcob.myBasisSet->makeClone()),
                    lcob.getOrbitalSetSize(), lcob.isIdentity(), lcob.isOMPoffload());

  *(eta.C) = *(lcob.C);
  *(phi.C) = *(lcob.C);


  int num_center = 3;
  std::vector<bool> corrCenter(num_center, "true");

  // N is first atom
  int center_idx = 0;

  using RealType = QMCTraits::RealType;

  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(1.00180500));
  CHECK((*eta.C)(0, 0) == Approx(0.0));

  int orbital_set_size = 7;
  Matrix<CuspCorrectionParameters> info;
  info.resize(num_center, orbital_set_size);
  okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);

  REQUIRE(okay);
  Vector<RealType> xgrid;
  Vector<RealType> rad_orb;
  int ngrid = 10;
  xgrid.resize(ngrid);
  for (int i = 0; i < ngrid; i++)
  {
    xgrid[i] = 0.012 * (i + 1);
  }

  rad_orb.resize(ngrid);

  int mo_idx = 0;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  // Comparisons generated from gen_cusp_corr.py
  //  Center  0  MO 0 rc =  0.07691307008
  CHECK(rad_orb[0] == Approx(9.1266186340)); // x = 0.012
  CHECK(rad_orb[1] == Approx(8.3939106599)); // x = 0.024
  CHECK(rad_orb[2] == Approx(7.7213972780)); // x = 0.036
  CHECK(rad_orb[3] == Approx(7.1039662640)); // x = 0.048
  CHECK(rad_orb[4] == Approx(6.5370601478)); // x = 0.06
  CHECK(rad_orb[5] == Approx(6.0165935481)); // x = 0.072
  CHECK(rad_orb[6] == Approx(5.5390213984)); // x = 0.084
  CHECK(rad_orb[7] == Approx(5.1023814795)); // x = 0.096
  CHECK(rad_orb[8] == Approx(4.7033287383)); // x = 0.108
  CHECK(rad_orb[9] == Approx(4.3370522377)); // x = 0.12


  mo_idx = 1;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  //  Center  0  MO 1 rc =  0.060909477888
  CHECK(rad_orb[0] == Approx(-0.0099816961)); // x = 0.012
  CHECK(rad_orb[1] == Approx(-0.0092950723)); // x = 0.024
  CHECK(rad_orb[2] == Approx(-0.0086498844)); // x = 0.036
  CHECK(rad_orb[3] == Approx(-0.0080440071)); // x = 0.048
  CHECK(rad_orb[4] == Approx(-0.0074778482)); // x = 0.06
  CHECK(rad_orb[5] == Approx(-0.0069529708)); // x = 0.072
  CHECK(rad_orb[6] == Approx(-0.0064707256)); // x = 0.084
  CHECK(rad_orb[7] == Approx(-0.0060313791)); // x = 0.096
  CHECK(rad_orb[8] == Approx(-0.0056312867)); // x = 0.108
  CHECK(rad_orb[9] == Approx(-0.0052652668)); // x = 0.12


  // Reset the MO matrices for another center

  *(eta.C) = *(lcob.C);
  *(phi.C) = *(lcob.C);


  // C is second atom
  center_idx = 1;
  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(0.0));
  CHECK((*eta.C)(0, 0) == Approx(1.00180500));

  mo_idx = 0;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  //  Center  1  MO 0 rc =  0.105
  CHECK(rad_orb[0] == Approx(0.0017535517)); // x = 0.012
  CHECK(rad_orb[1] == Approx(0.0016496533)); // x = 0.024
  CHECK(rad_orb[2] == Approx(0.0015544835)); // x = 0.036
  CHECK(rad_orb[3] == Approx(0.0014678130)); // x = 0.048
  CHECK(rad_orb[4] == Approx(0.0013891000)); // x = 0.06
  CHECK(rad_orb[5] == Approx(0.0013175785)); // x = 0.072
  CHECK(rad_orb[6] == Approx(0.0012523246)); // x = 0.084
  CHECK(rad_orb[7] == Approx(0.0011923038)); // x = 0.096
  CHECK(rad_orb[8] == Approx(0.0011364095)); // x = 0.108
  CHECK(rad_orb[9] == Approx(0.0010837868)); // x = 0.12


  removeSTypeOrbitals(corrCenter, lcob);

  CHECK((*lcob.C)(0, 0) == Approx(0.0));
  CHECK((*lcob.C)(0, 1) == Approx(0.0));
  CHECK((*lcob.C)(0, 2) == Approx(0.0));
  CHECK((*lcob.C)(0, 3) != 0.0);
}

TEST_CASE("HCN MO with cusp", "[wavefunction]")
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

  xmlSetProp(MO_base[0], castCharToXMLChar("cuspCorrection"), castCharToXMLChar("yes"));

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

  elec.R = 0.0;
  elec.update();
  ParticleSet::SingleParticlePos newpos;
  elec.makeMove(0, newpos);

  sposet->evaluateValue(elec, 0, values);

  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(0.00945227));
  CHECK(values[1] == Approx(0.0200836));
  CHECK(values[2] == Approx(0.416375));
  CHECK(values[3] == Approx(-0.0885443));
  CHECK(values[4] == Approx(0.273159));
  CHECK(values[5] == Approx(0));
  CHECK(values[6] == Approx(0));

  // Put electron near N atom
  elec.R[0][0] = -1.09;
  elec.update();
  elec.makeMove(0, newpos);

  values = 0.0;
  sposet->evaluateValue(elec, 0, values);
  //std::cout << "values = " << values << std::endl;
  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(9.5150713253));
  CHECK(values[1] == Approx(-0.0086731542));
  CHECK(values[2] == Approx(-1.6426151116));
  CHECK(values[3] == Approx(0.6569242017));
  CHECK(values[4] == Approx(0.9775522176));
  CHECK(values[5] == Approx(0.0000000000));
  CHECK(values[6] == Approx(0.0000000000));


  values = 0.0;
  sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);

  //std::cout << "values = " << values << std::endl;
  //std::cout << "dpsi = " << dpsi << std::endl;
  //std::cout << "d2psi = " << d2psi << std::endl;

  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(9.5150713253));
  CHECK(values[1] == Approx(-0.0086731542));
  CHECK(values[2] == Approx(-1.6426151116));
  CHECK(values[3] == Approx(0.6569242017));
  CHECK(values[4] == Approx(0.9775522176));
  CHECK(values[5] == Approx(0.0000000000));
  CHECK(values[6] == Approx(0.0000000000));

  CHECK(dpsi[0][0] == Approx(-66.5007223213));
  CHECK(dpsi[0][1] == Approx(0.0000000000));
  CHECK(dpsi[0][2] == Approx(0.0000000000));
  CHECK(d2psi[0] == Approx(-21540.9990552510));

  CHECK(values[1] == Approx(-0.0086731542));
  CHECK(dpsi[1][0] == Approx(0.0616909346));
  CHECK(dpsi[1][1] == Approx(0.0000000000));
  CHECK(dpsi[1][2] == Approx(0.0000000000));
  CHECK(d2psi[1] == Approx(19.8720529007));


  SPOSet::ValueMatrix all_values;
  SPOSet::GradMatrix all_grad;
  SPOSet::ValueMatrix all_lap;
  all_values.resize(7, 7);
  all_grad.resize(7, 7);
  all_lap.resize(7, 7);


  sposet->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(9.5150713253));

  CHECK(all_values[0][0] == Approx(9.5150713253));
  CHECK(all_grad[0][0][0] == Approx(-66.5007223213));
  CHECK(all_grad[0][0][1] == Approx(0.0000000000));
  CHECK(all_grad[0][0][2] == Approx(0.0000000000));
  CHECK(all_lap[0][0] == Approx(-21540.9990552510));

  CHECK(all_values[0][1] == Approx(-0.0086731542));
  CHECK(all_grad[0][1][0] == Approx(0.0616909346));
  CHECK(all_grad[0][1][1] == Approx(0.0000000000));
  CHECK(all_grad[0][1][2] == Approx(0.0000000000));
  CHECK(all_lap[0][1] == Approx(19.8720529007));


  // Test the makeClone method
  std::unique_ptr<SPOSet> sposet_clone(sposet->makeClone());

  sposet_clone->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(9.5150713253));

  CHECK(all_values[0][0] == Approx(9.5150713253));
  CHECK(all_grad[0][0][0] == Approx(-66.5007223213));
  CHECK(all_grad[0][0][1] == Approx(0.0000000000));
  CHECK(all_grad[0][0][2] == Approx(0.0000000000));
  CHECK(all_lap[0][0] == Approx(-21540.9990552510));

  CHECK(all_values[0][1] == Approx(-0.0086731542));
  CHECK(all_grad[0][1][0] == Approx(0.0616909346));
  CHECK(all_grad[0][1][1] == Approx(0.0000000000));
  CHECK(all_grad[0][1][2] == Approx(0.0000000000));
  CHECK(all_lap[0][1] == Approx(19.8720529007));
}

// Test case with multiple atoms of the same type
TEST_CASE("Ethanol MO with cusp", "[wavefunction]")
{
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

  xmlSetProp(MO_base[0], castCharToXMLChar("cuspCorrection"), castCharToXMLChar("yes"));

  const auto bb_ptr = bf.createSPOSetBuilder(MO_base[0]);
  auto& bb(*bb_ptr);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  auto sposet = bb.createSPOSet(slater_base[0]);

  SPOSet::ValueVector values;
  SPOSet::GradVector dpsi;
  SPOSet::ValueVector d2psi;
  values.resize(13);
  dpsi.resize(13);
  d2psi.resize(13);

  elec.R = 0.0;
  // Put electron near O atom
  elec.R[0][0] = -2.10;
  elec.R[0][1] = 0.50;

  elec.update();
  ParticleSet::SingleParticlePos newpos;
  elec.makeMove(0, newpos);

  sposet->evaluateValue(elec, 0, values);

  // Values from gen_cusp_corr.py
  CHECK(values[0] == Approx(4.3617329704));
  CHECK(values[1] == Approx(0.0014119853));
  CHECK(values[2] == Approx(0.0001156461));
  CHECK(values[3] == Approx(-0.6722670611));
  CHECK(values[4] == Approx(0.2762949842));
  CHECK(values[5] == Approx(0.2198735778));
  CHECK(values[6] == Approx(0.0659454461));
  CHECK(values[7] == Approx(0.2952071056));
  CHECK(values[8] == Approx(0.0322071389));
  CHECK(values[9] == Approx(0.0877981239));
  CHECK(values[10] == Approx(-0.2151873873));
  CHECK(values[11] == Approx(0.4250074750));
  CHECK(values[12] == Approx(0.0767950823));

  sposet->evaluateVGL(elec, 0, values, dpsi, d2psi);

  CHECK(values[0] == Approx(4.3617329704));
  CHECK(values[1] == Approx(0.0014119853));
  CHECK(values[2] == Approx(0.0001156461));
  CHECK(values[3] == Approx(-0.6722670611));
  CHECK(values[4] == Approx(0.2762949842));
  CHECK(values[5] == Approx(0.2198735778));
  CHECK(values[6] == Approx(0.0659454461));
  CHECK(values[7] == Approx(0.2952071056));
  CHECK(values[8] == Approx(0.0322071389));
  CHECK(values[9] == Approx(0.0877981239));
  CHECK(values[10] == Approx(-0.2151873873));
  CHECK(values[11] == Approx(0.4250074750));
  CHECK(values[12] == Approx(0.0767950823));

  CHECK(dpsi[0][0] == Approx(-27.2844138432));
  CHECK(dpsi[0][1] == Approx(15.9958208598));
  CHECK(dpsi[0][2] == Approx(0.0195317131));
  CHECK(d2psi[0] == Approx(-293.2869628790));

  CHECK(dpsi[12][0] == Approx(1.7548511775));
  CHECK(dpsi[12][1] == Approx(2.2759333828));
  CHECK(dpsi[12][2] == Approx(-1.4878277937));
  CHECK(d2psi[12] == Approx(-4.3399821309));


  SPOSet::ValueMatrix all_values;
  SPOSet::GradMatrix all_grad;
  SPOSet::ValueMatrix all_lap;
  all_values.resize(10, 13);
  all_grad.resize(10, 13);
  all_lap.resize(10, 13);

  sposet->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  CHECK(all_values[0][0] == Approx(4.3617329704));
  CHECK(all_grad[0][0][0] == Approx(-27.2844138432));
  CHECK(all_grad[0][0][1] == Approx(15.9958208598));
  CHECK(all_grad[0][0][2] == Approx(0.0195317131));
  CHECK(all_lap[0][0] == Approx(-293.2869628790));

  CHECK(all_values[0][11] == Approx(0.4250074750));
  CHECK(all_grad[0][11][0] == Approx(-0.3947036210));
  CHECK(all_grad[0][11][1] == Approx(0.9883840215));
  CHECK(all_grad[0][11][2] == Approx(1.7863218842));
  CHECK(all_lap[0][11] == Approx(-33.5202249813));
}


TEST_CASE("broadcastCuspInfo", "[wavefunction]")
{
  Communicate* c = OHMMS::Controller;
  CuspCorrectionParameters cp;
  int root = 0;
  if (c->rank() == root)
  {
    cp.Rc       = 2.0;
    cp.C        = 3.0;
    cp.sg       = -1.0;
    cp.alpha[0] = 1.1;
    cp.alpha[1] = 1.2;
    cp.alpha[2] = 1.3;
    cp.alpha[3] = 1.4;
    cp.alpha[4] = 1.5;
    cp.redo     = 1;
  }

  broadcastCuspInfo(cp, *c, root);

  CHECK(cp.Rc == Approx(2.0));
  CHECK(cp.C == Approx(3.0));
  CHECK(cp.sg == Approx(-1.0));
  CHECK(cp.alpha[0] == Approx(1.1));
  CHECK(cp.alpha[1] == Approx(1.2));
  CHECK(cp.alpha[2] == Approx(1.3));
  CHECK(cp.alpha[3] == Approx(1.4));
  CHECK(cp.alpha[4] == Approx(1.5));
  REQUIRE(cp.redo == 1);
}


} // namespace qmcplusplus

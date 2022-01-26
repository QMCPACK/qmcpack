//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Numerics/Quadrature.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "QMCHamiltonians/SOECPComponent.h"

//for wavefunction
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "QMCWaveFunctions/SpinorSet.h"
//for nonlocal moves
#include "QMCHamiltonians/NonLocalTOperator.h"


//for Hamiltonian manipulations.
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"

#ifdef QMC_COMPLEX //This is for the spinor test.
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#endif

namespace qmcplusplus
{
QMCTraits::RealType getSplinedSOPot(SOECPComponent* so_comp, int l, double r) { return so_comp->sopp_m[l]->splint(r); }

TEST_CASE("CheckSphericalIntegration", "[hamiltonian]")
{
  // Use the built-in quadrature rule check
  for (int quadrature_index = 1; quadrature_index < 8; quadrature_index++)
  {
    Quadrature3D<QMCTraits::RealType> myRule(quadrature_index, false);
    REQUIRE(myRule.quad_ok);
  }
}

TEST_CASE("ReadFileBuffer_no_file", "[hamiltonian]")
{
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("does_not_exist");
  REQUIRE(open_okay == false);
}

TEST_CASE("ReadFileBuffer_simple_serial", "[hamiltonian]")
{
  // Initializing with no Communicate pointer under MPI,
  //   this will read the file on every node.  Should be okay
  //   for testing purposes.
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);
  REQUIRE(std::string("File contents\n") == buf.contents());
}

TEST_CASE("ReadFileBuffer_simple_mpi", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  ReadFileBuffer buf(c);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);
  REQUIRE(std::string("File contents\n") == buf.contents());
}

TEST_CASE("ReadFileBuffer_ecp", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  ECPComponentBuilder ecp("test_read_ecp", c);

  bool okay = ecp.read_pp_file("C.BFD.xml");
  REQUIRE(okay);

  REQUIRE(ecp.Zeff == 4);

  // TODO: add more checks that pseudopotential file was read correctly
}

TEST_CASE("ReadFileBuffer_sorep", "[hamiltonian]")
{
  Communicate* c = OHMMS::Controller;

  ECPComponentBuilder ecp("test_read_sorep", c);

  bool okay = ecp.read_pp_file("so_ecp_test.xml");
  REQUIRE(okay);

  REQUIRE(ecp.Zeff == 13);

  double rlist[5] = {0.001, 0.500, 1.000, 2.000, 10.000};
  double so_p[5]  = {0.999999000005, 0.778800783071, 0.3678794411714, 0.01831563888873418, 0.000};
  double so_d[5]  = {9.99998000e-01, 6.06530660e-01, 1.35335283e-01, 3.35462628e-04, 0.000};
  double so_f[5]  = {9.99997000e-01, 4.72366553e-01, 4.97870684e-02, 6.14421235e-06, 0.000};

  for (int i = 0; i < 5; i++)
  {
    double r        = rlist[i];
    double so_p_ref = so_p[i];
    double so_d_ref = so_d[i];
    double so_f_ref = so_f[i];

    double so_p_val = getSplinedSOPot(ecp.pp_so.get(), 0, r);
    double so_d_val = getSplinedSOPot(ecp.pp_so.get(), 1, r);
    double so_f_val = getSplinedSOPot(ecp.pp_so.get(), 2, r);

    REQUIRE(so_p_val == Approx(so_p_ref));
    REQUIRE(so_d_val == Approx(so_d_ref));
    REQUIRE(so_f_val == Approx(so_f_ref));
  }

  // TODO: add more checks that pseudopotential file was read correctly
}


TEST_CASE("ReadFileBuffer_reopen", "[hamiltonian]")
{
  // Initializing with no Communicate pointer under MPI,
  //   this will read the file on every node.  Should be okay
  //   for testing purposes.
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);

  open_okay = buf.open_file("C.BFD.xml");
  REQUIRE(open_okay == true);

  read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length > 14);
}

void copyGridUnrotatedForTest(NonLocalECPComponent& nlpp) { nlpp.rrotsgrid_m = nlpp.sgridxyz_m; }
void copyGridUnrotatedForTest(SOECPComponent& sopp) { sopp.rrotsgrid_m = sopp.sgridxyz_m; }

TEST_CASE("Evaluate_ecp", "[hamiltonian]")
{
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R.diagonal(20);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion0");
  ions.create(2);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;
  ions.R[1][0] = 6.0;
  ions.R[1][1] = 0.0;
  ions.R[1][2] = 0.0;

  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("Na");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 1;
  ion_species(iatnumber, pIdx)  = 11;
  ions.createSK();

  elec.setName("e");
  std::vector<int> agroup(2, 1);
  elec.create(agroup);
  elec.R[0][0] = 2.0;
  elec.R[0][1] = 0.0;
  elec.R[0][2] = 0.0;
  elec.R[1][0] = 3.0;
  elec.R[1][1] = 0.0;
  elec.R[1][2] = 0.0;


  SpeciesSet& tspecies         = elec.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  int massIdx                  = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  tspecies(massIdx, upIdx)     = 1.0;
  tspecies(massIdx, downIdx)   = 1.0;

  elec.createSK();

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  //Cool.  Now to construct a wavefunction with 1 and 2 body jastrow (no determinant)
  TrialWaveFunction psi;

  //Add the two body jastrow
  const char* particles = "<tmp> \
  <jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\" gpu=\"no\">  \
      <correlation speciesA=\"u\" speciesB=\"d\" rcut=\"10\" size=\"8\"> \
          <coefficients id=\"ud\" type=\"Array\"> 2.015599059 1.548994099 1.17959447 0.8769687661 0.6245736507 0.4133517767 0.2333851935 0.1035636904</coefficients> \
        </correlation> \
  </jastrow> \
  </tmp> \
  ";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  xmlNodePtr jas2 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow(c, elec);
  psi.addComponent(jastrow.buildComponent(jas2));
  // Done with two body jastrow.

  //Add the one body jastrow.
  const char* particles2 = "<tmp> \
  <jastrow name=\"J1\" type=\"One-Body\" function=\"Bspline\" source=\"ion0\" print=\"yes\"> \
        <correlation elementType=\"Na\" rcut=\"10\" size=\"10\" cusp=\"0\"> \
          <coefficients id=\"eNa\" type=\"Array\"> 1.244201343 -1.188935609 -1.840397253 -1.803849126 -1.612058635 -1.35993202 -1.083353212 -0.8066295188 -0.5319252448 -0.3158819772</coefficients> \
        </correlation> \
      </jastrow> \
  </tmp> \
  ";
  bool okay3             = doc.parseFromString(particles2);
  REQUIRE(okay3);

  root = doc.getRoot();

  xmlNodePtr jas1 = xmlFirstElementChild(root);

  RadialJastrowBuilder jastrow1bdy(c, elec, ions);
  psi.addComponent(jastrow1bdy.buildComponent(jas1));

  //Now we set up the nonlocal ECP component.
  ECPComponentBuilder ecp("test_read_ecp", c);

  bool okay2 = ecp.read_pp_file("Na.BFD.xml");

  NonLocalECPComponent* nlpp = ecp.pp_nonloc.get();

  REQUIRE(nlpp != nullptr);

  //This line is required because the randomized quadrature Lattice is set by
  //random number generator in NonLocalECPotential.  We take the unrotated
  //quadrature Lattice instead...
  copyGridUnrotatedForTest(*nlpp);

  const int myTableIndex = elec.addTable(ions);

  const auto& myTable = elec.getDistTableAB(myTableIndex);

  // update all distance tables
  ions.update();
  elec.update();

  //Need to set up temporary data for this configuration in trial wavefunction.  Needed for ratios.
  double logpsi = psi.evaluateLog(elec);
  REQUIRE(logpsi == Approx(5.1497823982));

  auto test_evaluateOne = [&]()
  {
    double Value1(0.0);
    //Using SoA distance tables, hence the guard.
    for (int jel = 0; jel < elec.getTotalNum(); jel++)
    {
      const auto& dist  = myTable.getDistRow(jel);
      const auto& displ = myTable.getDisplRow(jel);
      for (int iat = 0; iat < ions.getTotalNum(); iat++)
        if (nlpp != nullptr && dist[iat] < nlpp->getRmax())
          Value1 += nlpp->evaluateOne(elec, iat, psi, jel, dist[iat], -displ[iat], false);
    }
    //These numbers are validated against an alternate code path via wavefunction tester.
    CHECK(Value1 == Approx(6.9015710211e-02));
  };

  {
    test_evaluateOne();
    nlpp->initVirtualParticle(elec);
    test_evaluateOne();
    nlpp->deleteVirtualParticle();
  }

  opt_variables_type optvars;
  std::vector<ValueType> dlogpsi;
  std::vector<ValueType> dhpsioverpsi;

  psi.checkInVariables(optvars);
  optvars.resetIndex();
  const int NumOptimizables(optvars.size());
  psi.checkOutVariables(optvars);
  auto test_evaluateValueAndDerivatives = [&]()
  {
    dlogpsi.resize(NumOptimizables, ValueType(0));
    dhpsioverpsi.resize(NumOptimizables, ValueType(0));
    psi.evaluateDerivatives(elec, optvars, dlogpsi, dhpsioverpsi);
    REQUIRE(std::real(dlogpsi[0]) == Approx(-0.2211666667));
    REQUIRE(std::real(dlogpsi[2]) == Approx(-0.1215));
    REQUIRE(std::real(dlogpsi[3]) == Approx(0.0));
    REQUIRE(std::real(dlogpsi[9]) == Approx(-0.0853333333));
    REQUIRE(std::real(dlogpsi[10]) == Approx(-0.745));

    REQUIRE(std::real(dhpsioverpsi[0]) == Approx(-0.6463306581));
    REQUIRE(std::real(dhpsioverpsi[2]) == Approx(1.5689981479));
    REQUIRE(std::real(dhpsioverpsi[3]) == Approx(0.0));
    REQUIRE(std::real(dhpsioverpsi[9]) == Approx(0.279561213));
    REQUIRE(std::real(dhpsioverpsi[10]) == Approx(-0.3968828778));

    double Value1 = 0.0;
    //Using SoA distance tables, hence the guard.
    for (int jel = 0; jel < elec.getTotalNum(); jel++)
    {
      const auto& dist  = myTable.getDistRow(jel);
      const auto& displ = myTable.getDisplRow(jel);
      for (int iat = 0; iat < ions.getTotalNum(); iat++)
        if (nlpp != nullptr && dist[iat] < nlpp->getRmax())
          Value1 += nlpp->evaluateValueAndDerivatives(elec, iat, psi, jel, dist[iat], -displ[iat], optvars, dlogpsi,
                                                      dhpsioverpsi);
    }
    CHECK(Value1 == Approx(6.9015710211e-02));

    CHECK(std::real(dhpsioverpsi[0]) == Approx(-0.6379341942));
    CHECK(std::real(dhpsioverpsi[2]) == Approx(1.5269279991));
    CHECK(std::real(dhpsioverpsi[3]) == Approx(-0.0355730676));
    CHECK(std::real(dhpsioverpsi[9]) == Approx(0.279561213));
    CHECK(std::real(dhpsioverpsi[10]) == Approx(-0.3968763604));
  };

  {
    test_evaluateValueAndDerivatives();
    nlpp->initVirtualParticle(elec);
    test_evaluateValueAndDerivatives();
    nlpp->deleteVirtualParticle();
  }

  double Value2(0.0);
  double Value3(0.0);
  ParticleSet::ParticlePos PulayTerm, HFTerm, HFTerm2;
  HFTerm.resize(ions.getTotalNum());
  HFTerm2.resize(ions.getTotalNum());
  PulayTerm.resize(ions.getTotalNum());
  HFTerm    = 0;
  HFTerm2   = 0;
  PulayTerm = 0;

  for (int jel = 0; jel < elec.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.getDistRow(jel);
    const auto& displ = myTable.getDisplRow(jel);
    for (int iat = 0; iat < ions.getTotalNum(); iat++)
      if (nlpp != nullptr && dist[iat] < nlpp->getRmax())
      {
        Value2 += nlpp->evaluateOneWithForces(elec, iat, psi, jel, dist[iat], -displ[iat], HFTerm[iat]);
        Value3 +=
            nlpp->evaluateOneWithForces(elec, ions, iat, psi, jel, dist[iat], -displ[iat], HFTerm2[iat], PulayTerm);
      }
  }
  //These values are validated against print statements.
  //Two-body jastrow-only wave functions agree with finite difference of NLPP to machine precision.
  //  These numbers assume the Hellman Feynmann implementation is correct, and dump the values
  //  when a one body term is added in.

  REQUIRE(Value2 == Approx(6.9015710211e-02));
  REQUIRE(Value3 == Approx(6.9015710211e-02));

  //The total force (HFTerm+PulayTerm) is validated against finite difference of nonlocal PP w.r.t
  //ion coordinates. delta=1e-6.  Should be good up to 7th or 8th sig fig. These are:
  // F[0][0]= 0.3474359
  // F[0][1]= 0
  // F[0][2]= 0
  // F[1][0]=-0.002734064
  // F[1][1]= 0
  // F[1][2]= 0

  REQUIRE(HFTerm[0][0] == Approx(-0.3557369031));
  REQUIRE(HFTerm[0][1] == Approx(0.0));
  REQUIRE(HFTerm[0][2] == Approx(0.0));
  REQUIRE(HFTerm[1][0] == Approx(0.001068673105));
  REQUIRE(HFTerm[1][1] == Approx(0.0));
  REQUIRE(HFTerm[1][2] == Approx(0.0));

  REQUIRE(HFTerm2[0][0] == Approx(-0.3557369031));
  REQUIRE(HFTerm2[0][1] == Approx(0.0));
  REQUIRE(HFTerm2[0][2] == Approx(0.0));
  REQUIRE(HFTerm2[1][0] == Approx(0.001068673105));
  REQUIRE(HFTerm2[1][1] == Approx(0.0));
  REQUIRE(HFTerm2[1][2] == Approx(0.0));

  REQUIRE(PulayTerm[0][0] == Approx(0.008300993315));
  REQUIRE(PulayTerm[0][1] == Approx(0.0));
  REQUIRE(PulayTerm[0][2] == Approx(0.0));
  REQUIRE(PulayTerm[1][0] == Approx(0.001665391103));
  REQUIRE(PulayTerm[1][1] == Approx(0.0));
  REQUIRE(PulayTerm[1][2] == Approx(0.0));

  //Comparing against finite difference results above, here's what we get.
  //HFTerm[0][0]+PulayTerm[0][0] = −0.34743591
  //HFTerm[0][1]+PulayTerm[0][1] =  0.0
  //HFTerm[0][2]+PulayTerm[0][2] =  0.0
  //HFTerm[1][0]+PulayTerm[1][0] =  0.002734064
  //HFTerm[1][1]+PulayTerm[1][1] =  0.0
  //HFTerm[1][2]+PulayTerm[1][2] =  0.0
}

#ifdef QMC_COMPLEX
TEST_CASE("Evaluate_soecp", "[hamiltonian]")
{
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!! Evaluate SOECPComponent !!!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  using RealType  = QMCTraits::RealType;
  using ValueType = QMCTraits::ValueType;
  using PosType   = QMCTraits::PosType;

  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = false; // periodic
  lattice.R.diagonal(20);
  lattice.LR_dim_cutoff = 15;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  auto ions_uptr = std::make_unique<ParticleSet>(simulation_cell);
  auto elec_uptr = std::make_unique<ParticleSet>(simulation_cell);
  ParticleSet& ions(*ions_uptr);
  ParticleSet& elec(*elec_uptr);

  ions.setName("ion0");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;


  SpeciesSet& ion_species       = ions.getSpeciesSet();
  int pIdx                      = ion_species.addSpecies("H");
  int pChargeIdx                = ion_species.addAttribute("charge");
  int iatnumber                 = ion_species.addAttribute("atomic_number");
  ion_species(pChargeIdx, pIdx) = 0;
  ion_species(iatnumber, pIdx)  = 1;
  ions.createSK();


  elec.setName("e");
  elec.create(1);
  elec.R[0][0]  = 0.138;
  elec.R[0][1]  = -0.24;
  elec.R[0][2]  = 0.216;
  elec.spins[0] = 0.0;

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();

  ions.resetGroups();
  elec.resetGroups();

  TrialWaveFunction psi;

  std::vector<PosType> kup, kdn;
  std::vector<RealType> k2up, k2dn;
  QMCTraits::IndexType nelec = elec.getTotalNum();
  REQUIRE(nelec == 1);

  kup.resize(nelec);
  kup[0] = PosType(1, 1, 1);

  k2up.resize(nelec);
  //For some goofy reason, EGOSet needs to be initialized with:
  //1.) A k-vector list (fine).
  //2.) A list of -|k|^2.  To save on expensive - sign multiplication apparently.
  k2up[0] = -dot(kup[0], kup[0]);

  kdn.resize(nelec);
  kdn[0] = PosType(2, 2, 2);

  k2dn.resize(nelec);
  k2dn[0] = -dot(kdn[0], kdn[0]);

  auto spo_up = std::make_unique<EGOSet>(kup, k2up);
  auto spo_dn = std::make_unique<EGOSet>(kdn, k2dn);

  auto spinor_set = std::make_unique<SpinorSet>();
  spinor_set->set_spos(std::move(spo_up), std::move(spo_dn));
  QMCTraits::IndexType norb = spinor_set->getOrbitalSetSize();
  REQUIRE(norb == 1);

  auto dd = std::make_unique<DiracDeterminant<>>(std::move(spinor_set), 0, nelec);

  psi.addComponent(std::move(dd));

  //Now we set up the SO ECP component.
  ECPComponentBuilder ecp("test_read_soecp", c);

  bool okay3 = ecp.read_pp_file("so_ecp_test.xml");
  REQUIRE(okay3);

  SOECPComponent* sopp = ecp.pp_so.get();
  REQUIRE(sopp != nullptr);
  copyGridUnrotatedForTest(*sopp);

  const int myTableIndex = elec.addTable(ions);

  const auto& myTable = elec.getDistTableAB(myTableIndex);

  // update all distance tables
  ions.update();
  elec.update();

  //Need to set up temporary data for this configuration in trial wavefunction.  Needed for ratios.
  double logpsi = psi.evaluateLog(elec);

  RealType Value1(0.0);

  for (int jel = 0; jel < elec.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.getDistRow(jel);
    const auto& displ = myTable.getDisplRow(jel);
    for (int iat = 0; iat < ions.getTotalNum(); iat++)
    {
      if (sopp != nullptr && dist[iat] < sopp->getRmax())
      {
        Value1 += sopp->evaluateOne(elec, iat, psi, jel, dist[iat], RealType(-1) * displ[iat]);
      }
    }
  }
  REQUIRE(Value1 == Approx(-0.3214176962));
}
#endif


} // namespace qmcplusplus

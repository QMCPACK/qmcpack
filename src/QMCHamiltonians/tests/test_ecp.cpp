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

//for wavefunction
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
//for nonlocal moves
#include "QMCHamiltonians/NonLocalTOperator.h"

//for Hamiltonian manipulations.
#include "Lattice/ParticleBConds.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "QMCApp/ParticleSetPool.h"
#include "LongRange/EwaldHandler3D.h"

namespace qmcplusplus
{
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
  OHMMS::Controller->initialize(0, NULL);
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
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  ECPComponentBuilder ecp("test_read_ecp", c);

  bool okay = ecp.read_pp_file("C.BFD.xml");
  REQUIRE(okay);

  REQUIRE(ecp.Zeff == 4);

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

void copyGridUnrotatedForTest(NonLocalECPComponent& nlpp)
{
  nlpp.rrotsgrid_m = nlpp.sgridxyz_m;
}

TEST_CASE("Evaluate_ecp", "[hamiltonian]")
{
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::PosType PosType;

  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  //Cell definition:

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R.diagonal(20);
  Lattice.LR_dim_cutoff = 15;
  Lattice.reset();


  ParticleSet ions;
  ParticleSet elec;

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
  ions.Lattice = Lattice;
  ions.createSK();


  elec.Lattice = Lattice;
  elec.setName("e");
  elec.create(2);
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

  ParticleSetPool ptcl = ParticleSetPool(c);

  ions.resetGroups();

  // The call to resetGroups is needed transfer the SpeciesSet
  // settings to the ParticleSet
  elec.resetGroups();

  //Cool.  Now to construct a wavefunction with 1 and 2 body jastrow (no determinant)
  TrialWaveFunction psi = TrialWaveFunction(c);

  //Add the two body jastrow
  const char* particles = "<tmp> \
  <jastrow name=\"J2\" type=\"Two-Body\" function=\"Bspline\" print=\"yes\">  \
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

  RadialJastrowBuilder jastrow(elec, psi);
  bool build_okay = jastrow.put(jas2);
  REQUIRE(build_okay);
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

  RadialJastrowBuilder jastrow1bdy(elec, psi, ions);
  bool build_okay2 = jastrow1bdy.put(jas1);
  REQUIRE(build_okay2);

  //Now we set up the nonlocal ECP component.
  ECPComponentBuilder ecp("test_read_ecp", c);

  bool okay2 = ecp.read_pp_file("Na.BFD.xml");

  NonLocalECPComponent* nlpp = ecp.pp_nonloc;

  //This line is required because the randomized quadrature Lattice is set by
  //random number generator in NonLocalECPotential.  We take the unrotated
  //quadrature Lattice instead...
  copyGridUnrotatedForTest(*nlpp);


  //Not testing nonlocal moves here, but the PP functions take this as an argument
  NonLocalTOperator nonLocalOps(elec.getTotalNum());
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);

  const int myTableIndex = elec.addTable(ions, DT_SOA_PREFERRED);

  const auto& myTable = elec.getDistTable(myTableIndex);

  // update all distance tables
  ions.update();
  elec.update();

  //Need to set up temporary data for this configuration in trial wavefunction.  Needed for ratios.
  double logpsi = psi.evaluateLog(elec);

  double Value1(0.0);
#ifdef ENABLE_SOA
  //Forces are only implemented in SOA version, hence the guard.
  double Value2(0.0);
  double Value3(0.0);
  ParticleSet::ParticlePos_t PulayTerm, HFTerm, HFTerm2;
  HFTerm.resize(ions.getTotalNum());
  HFTerm2.resize(ions.getTotalNum());
  PulayTerm.resize(ions.getTotalNum());
  HFTerm    = 0;
  HFTerm2   = 0;
  PulayTerm = 0;

  for (int jel = 0; jel < elec.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.Distances[jel];
    const auto& displ = myTable.Displacements[jel];
    for (int iat = 0; iat < ions.getTotalNum(); iat++)
      if (nlpp != nullptr && dist[iat] < nlpp->getRmax())
      {
        Value1 += nlpp->evaluateOne(elec, iat, psi, jel, dist[iat], RealType(-1) * displ[iat], 0, Txy);

        Value2 +=
            nlpp->evaluateOneWithForces(elec, iat, psi, jel, dist[iat], RealType(-1) * displ[iat], HFTerm[iat], 0, Txy);
        Value3 += nlpp->evaluateOneWithForces(elec, ions, iat, psi, jel, dist[iat], RealType(-1) * displ[iat],
                                              HFTerm2[iat], PulayTerm, 0, Txy);
      }
  }
  //These numbers are validated against an alternate code path via wavefunction tester.
  REQUIRE(Value1 == Approx(6.9015710211e-02));

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


#else

  for (int iat = 0; iat < ions.getTotalNum(); iat++)
  {
    if (nlpp == nullptr)
      continue;
    for (int nn = myTable.M[iat], iel = 0; nn < myTable.M[iat + 1]; nn++, iel++)
    {
      const RealType r(myTable.r(nn));
      if (r > nlpp->getRmax())
        continue;
      Value1 += nlpp->evaluateOne(elec, iat, psi, iel, r, myTable.dr(nn), 0, Txy);
    }
  }
  REQUIRE(Value1 == Approx(6.9015710211e-02));

#endif
}
} // namespace qmcplusplus

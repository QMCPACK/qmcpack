#include "catch.hpp"

#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/CoulombPBCAA.h"
#include "LongRange/EwaldHandler2D.h"

namespace qmcplusplus
{

TEST_CASE("Coulomb PBC A-A Ewald2D square", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  const double vmad_sq = -1.95013246;
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R.diagonal(1.0);
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.update();

  bool is_active = true;
  CoulombPBCAA caa = CoulombPBCAA(elec, is_active);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_sq));
}

TEST_CASE("Coulomb PBC A-A Ewald2D triangular", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0; // this is crutial!
  const double vmad_tri = -1.1061025865191676;
  const double alat = std::sqrt(2.0*M_PI/std::sqrt(3));

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true; // periodic
  lattice.ndim = 2;
  lattice.R = 0.0;
  lattice.R(0, 0) = alat;
  lattice.R(1, 0) = -1.0/2*alat;
  lattice.R(1, 1) = std::sqrt(3)/2*alat;
  lattice.R(2, 2) = 2*alat;
  lattice.LR_dim_cutoff = 30.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1});
  elec.R[0] = {0.0, 0.0, 0.0};

  SpeciesSet& tspecies       = elec.getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int chargeIdx              = tspecies.addAttribute("charge");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(chargeIdx, upIdx) = -1;
  tspecies(massIdx, upIdx)   = 1.0;

  elec.createSK();
  elec.update();

  bool is_active = true;
  CoulombPBCAA caa = CoulombPBCAA(elec, is_active);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_tri));
}
} // qmcplusplus

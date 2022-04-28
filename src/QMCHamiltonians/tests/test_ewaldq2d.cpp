#include "catch.hpp"

#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/CoulombPBCAA.h"

namespace qmcplusplus
{

TEST_CASE("Coulomb PBC A-A Ewald Quasi2D square", "[hamiltonian]")
{
  LRCoulombSingleton::CoulombHandler = 0;
  LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
  const double vmad_sq = -1.95013246;
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = true;
  lattice.BoxBConds[2] = false; // ppn
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

  CoulombPBCAA caa = CoulombPBCAA(elec, true, false, false);

  double val = caa.evaluate(elec);
  CHECK(val == Approx(vmad_sq));
}

} // qmcplusplus

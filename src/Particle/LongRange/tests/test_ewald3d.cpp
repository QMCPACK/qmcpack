//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 and QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, University of Illinois Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Configuration.h"
#include "Lattice/CrystalLattice.h"
#include "Particle/ParticleSet.h"
#include "LongRange/EwaldHandler3D.h"

namespace qmcplusplus
{
using mRealType = EwaldHandler3D::mRealType;

/** evalaute bare Coulomb using EwaldHandler3D
 */
TEST_CASE("ewald3d", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Lattice.LR_rc == Approx(2.5));
  REQUIRE(Lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  EwaldHandler3D handler(ref, Lattice.LR_kc);

  // make sure initBreakup changes the default sigma
  REQUIRE(handler.Sigma == Approx(Lattice.LR_kc));
  handler.initBreakup(ref);
  REQUIRE(handler.Sigma == Approx(std::sqrt(Lattice.LR_kc / (2.0 * Lattice.LR_rc))));

  REQUIRE(handler.MaxKshell == 78);
  REQUIRE(handler.LR_rc == Approx(2.5));
  REQUIRE(handler.LR_kc == Approx(12));

  mRealType r, dr, rinv;
  mRealType vsr, vlr;
  int nr = 101;
  dr     = 5.0 / nr; // L/[# of grid points]
  for (int ir = 1; ir < nr; ir++)
  {
    r    = ir * dr;
    rinv = 1. / r;
    vsr  = handler.evaluate(r, rinv);
    vlr  = handler.evaluateLR(r);
    // short-range part must vanish after rcut
    if (r > 2.5)
      REQUIRE(vsr == Approx(0.0));
    // sum must recover the Coulomb potential
    REQUIRE(vsr + vlr == Approx(rinv));
  }
}

/** evalaute bare Coulomb derivative using EwaldHandler3D
 */
TEST_CASE("ewald3d df", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Lattice.LR_rc == Approx(2.5));
  REQUIRE(Lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  EwaldHandler3D handler(ref, Lattice.LR_kc);

  // make sure initBreakup changes the default sigma
  REQUIRE(handler.Sigma == Approx(Lattice.LR_kc));
  handler.initBreakup(ref);
  REQUIRE(handler.Sigma == Approx(std::sqrt(Lattice.LR_kc / (2.0 * Lattice.LR_rc))));

  REQUIRE(handler.MaxKshell == 78);
  REQUIRE(handler.LR_rc == Approx(2.5));
  REQUIRE(handler.LR_kc == Approx(12));

  mRealType r, dr, rinv;
  mRealType rm, rp; // minus (m), plus (p)
  mRealType vsrm, vsrp, dvsr, vlrm, vlrp, dvlr;
  dr                           = 0.00001; // finite difference step
  std::vector<mRealType> rlist = {0.1, 0.5, 1.0, 2.0, 2.5};
  for (auto it = rlist.begin(); it != rlist.end(); ++it)
  {
    r = *it;
    // test short-range piece
    rm   = r - dr;
    rinv = 1. / rm;
    vsrm = handler.evaluate(rm, rinv);
    vlrm = handler.evaluateLR(rm);
    rp   = r + dr;
    rinv = 1. / rp;
    vsrp = handler.evaluate(rp, rinv);
    vlrp = handler.evaluateLR(rp);
    dvsr = (vsrp - vsrm) / (2 * dr);
    rinv = 1. / r;
    REQUIRE(handler.srDf(r, rinv) == Approx(dvsr));
    // test long-range piece
    dvlr = (vlrp - vlrm) / (2 * dr);
    REQUIRE(handler.lrDf(r) == Approx(dvlr));
  }
}

} // namespace qmcplusplus

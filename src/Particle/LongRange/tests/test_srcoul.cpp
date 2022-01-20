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
#include "LongRange/LRHandlerSRCoulomb.h"

namespace qmcplusplus
{
using mRealType = LRHandlerBase::mRealType;

struct EslerCoulomb3D_ForSRCOUL
{ // stripped down version of LRCoulombSingleton::CoulombFunctor for 3D
  double norm;
  inline double operator()(double r, double rinv) const { return rinv; }
  inline double Vk(double k) const { return 1. / (k * k); }
  inline double dVk_dk(double k) const { return -2 * norm / (k * k * k); }
  void reset(ParticleSet& ref) { norm = 4.0 * M_PI / ref.getLRBox().Volume; }
  void reset(ParticleSet& ref, double rs) { reset(ref); } // ignore rs
  inline double df(double r) const { return -1. / (r * r); }
};

/** evalaute bare Coulomb in 3D
 */
TEST_CASE("srcoul", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Approx(Lattice.LR_rc) == 2.5);
  REQUIRE(Approx(Lattice.LR_kc) == 12);

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.SK.getKLists()
  ref.createSK();
  LRHandlerSRCoulomb<EslerCoulomb3D_ForSRCOUL, LPQHISRCoulombBasis> handler(ref);

  handler.initBreakup(ref);
  REQUIRE(handler.MaxKshell == 78);
  REQUIRE(Approx(handler.LR_rc) == 2.5);
  REQUIRE(Approx(handler.LR_kc) == 12);

  mRealType r, dr, rinv;
  mRealType vsr;
  int nr = 101;
  dr     = 5.0 / nr; // L/[# of grid points]
  for (int ir = 1; ir < nr; ir++)
  {
    r    = ir * dr;
    rinv = 1. / r;
    vsr  = handler.evaluate(r, rinv);
    // short-range part must vanish after rcut
    if (r > 2.5)
      REQUIRE(vsr == Approx(0.0));
    //// !!!! SR values not validated, see "srcoul df" test
    //REQUIRE(vsr == Approx(rinv));
  }
}

/** evalaute bare Coulomb derivative in 3D
 */
TEST_CASE("srcoul df", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Approx(Lattice.LR_rc) == 2.5);
  REQUIRE(Approx(Lattice.LR_kc) == 12);

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.SK.getKLists()
  ref.createSK();
  LRHandlerSRCoulomb<EslerCoulomb3D_ForSRCOUL, LPQHISRCoulombBasis> handler(ref);

  handler.initBreakup(ref);
  REQUIRE(handler.MaxKshell == 78);
  REQUIRE(Approx(handler.LR_rc) == 2.5);
  REQUIRE(Approx(handler.LR_kc) == 12);

  EslerCoulomb3D_ForSRCOUL fref;
  fref.reset(ref);
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
    // test total derivative
    REQUIRE(dvsr + dvlr == Approx(fref.df(r)));
  }
}

} // namespace qmcplusplus

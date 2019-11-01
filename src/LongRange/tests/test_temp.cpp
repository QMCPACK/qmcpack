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
#include "LongRange/LRHandlerTemp.h"

namespace qmcplusplus
{

using mRealType = LRHandlerBase::mRealType;

struct EslerCoulomb3D
{ // stripped down version of LRCoulombSingleton::CoulombFunctor for 3D
  double norm;
  inline double operator()(double r, double rinv) {return rinv;}
  void reset(ParticleSet& ref) {norm=4.0*M_PI/ref.LRBox.Volume;}
  inline double Xk(double k, double rc) {return -norm/(k*k)*std::cos(k*rc);}
  inline double Fk(double k, double rc) {return -Xk(k, rc);}
  inline double integrate_r2(double r) const {return 0.5*r*r;}
  inline double df(double r) {return 0;} // ignore derivatives for now
  void reset(ParticleSet& ref, double rs) {reset(ref);} // ignore rs
};

/** evalaute bare Coulomb in 3D using LRHandlerTemp
 */
TEST_CASE("temp3d", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  REQUIRE(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  REQUIRE(Lattice.LR_rc == Approx(2.5));
  REQUIRE(Lattice.LR_kc == Approx(12));

  ParticleSet ref; // handler needs ref.SK.KLists
  ref.Lattice = Lattice;  // !!!! crucial for access to Volume
  ref.LRBox = Lattice;  // !!!! crucial for S(k) update
  StructFact *SK = new StructFact(ref, Lattice.LR_kc);
  ref.SK = SK;
  LRHandlerTemp<EslerCoulomb3D, LPQHIBasis> handler(ref);

  handler.initBreakup(ref);
  REQUIRE(handler.MaxKshell == 78);
  REQUIRE(handler.LR_rc == Approx(2.5));
  REQUIRE(handler.LR_kc == 12);

  mRealType r, dr, rinv;
  mRealType vsr, vlr;
  int nr = 101;
  dr = 5.0/nr;  // L/[# of grid points]
  for (int ir=1;ir<nr;ir++)
  {
    r = ir*dr;
    rinv = 1./r;
    vsr = handler.evaluate(r, rinv);
    vlr = handler.evaluateLR(r);
    // short-range part must vanish after rcut
    if (r>2.5) REQUIRE(vsr == Approx(0.0));
    // sum must recover the Coulomb potential
    REQUIRE(vsr+vlr == Approx(rinv));
  }
}

} // namespace qmcplusplus

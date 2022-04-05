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
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
using pRealType = qmcplusplus::LRHandlerBase::pRealType;

struct CoulombF2
{
  inline double operator()(double k2) { return 4 * M_PI / k2; }
};

/** evalaute bare Coulomb using DummyLRHandler
 */
TEST_CASE("dummy", "[lrhandler]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds     = true;
  Lattice.LR_dim_cutoff = 30.;
  Lattice.R.diagonal(5.0);
  Lattice.reset();
  CHECK(Lattice.Volume == Approx(125));
  Lattice.SetLRCutoffs(Lattice.Rv);
  //Lattice.printCutoffs(app_log());
  CHECK(Lattice.LR_rc == Approx(2.5));
  CHECK(Lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(Lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  DummyLRHandler<CoulombF2> handler(Lattice.LR_kc);

  handler.initBreakup(ref);

  std::cout << "handler.MaxKshell is " << handler.MaxKshell << std::endl;
  CHECK( (std::is_same<OHMMS_PRECISION, OHMMS_PRECISION_FULL>::value ?
     handler.MaxKshell == 78 : handler.MaxKshell >= 117 && handler.MaxKshell <= 128 ));
  CHECK(handler.LR_kc == Approx(12));
  CHECK(handler.LR_rc == Approx(0));

  std::vector<pRealType> rhok1(handler.MaxKshell);
  std::vector<pRealType> rhok2(handler.MaxKshell);
  CoulombF2 fk;
  double norm = 4 * M_PI / Lattice.Volume;
  // no actual LR breakup happened in DummyLRHandler,
  //  the full Coulomb potential should be retained in kspace
  for (int ish = 0; ish < handler.MaxKshell; ish++)
  {
    int ik           = ref.getSimulationCell().getKLists().kshell[ish];
    double k2        = ref.getSimulationCell().getKLists().ksq[ik];
    double fk_expect = fk(k2);
    CHECK(handler.Fk_symm[ish] == Approx(norm * fk_expect));
  }
  // ?? cannot access base class method, too many overloads?
  // handler.evaluate(SK->getKLists().kshell, rhok1.data(), rhok2.data());
}

} // namespace qmcplusplus

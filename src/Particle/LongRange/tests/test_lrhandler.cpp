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
  Lattice lattice;
  lattice.BoxBConds     = true;
  lattice.LR_dim_cutoff = 30.;
  lattice.R.diagonal(5.0);
  lattice.reset();
  CHECK(lattice.Volume == Approx(125));
  lattice.SetLRCutoffs(lattice.Rv);
  //lattice.printCutoffs(app_log());
  CHECK(lattice.LR_rc == Approx(2.5));
  CHECK(lattice.LR_kc == Approx(12));

  const SimulationCell simulation_cell(lattice);
  ParticleSet ref(simulation_cell);       // handler needs ref.getSimulationCell().getKLists()
  ref.createSK();
  DummyLRHandler<CoulombF2> handler(lattice.LR_kc);

  handler.initBreakup(ref);

  std::cout << "handler.MaxKshell is " << handler.MaxKshell << std::endl;
  CHECK( handler.MaxKshell == 78);
  CHECK(handler.LR_kc == Approx(12));
  CHECK(handler.LR_rc == Approx(0));

  std::vector<pRealType> rhok1(handler.MaxKshell);
  std::vector<pRealType> rhok2(handler.MaxKshell);
  CoulombF2 fk;
  double norm = 4 * M_PI / lattice.Volume;
  // no actual LR breakup happened in DummyLRHandler,
  //  the full Coulomb potential should be retained in kspace
  for (int ish = 0; ish < handler.MaxKshell; ish++)
  {
    int ik           = ref.getSimulationCell().getKLists().getKShell()[ish];
    double k2        = ref.getSimulationCell().getKLists().getKSQWorking()[ik];
    double fk_expect = fk(k2);
    CHECK(handler.Fk_symm[ish] == Approx(norm * fk_expect));
  }
  // ?? cannot access base class method, too many overloads?
  // handler.evaluate(SK->getKLists().kshell, rhok1.data(), rhok2.data());
}

} // namespace qmcplusplus

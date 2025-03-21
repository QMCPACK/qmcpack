//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////
#include "catch.hpp"

#include "Particle/LongRange/KContainer.h"
#include "OhmmsPETE/TinyVector.h"
#include "Particle/Lattice/CrystalLattice.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "Utilities/for_testing/checkVector.hpp"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"

namespace qmcplusplus
{
using PosType = TinyVector<double, 3>;

TEST_CASE("kcontainer at gamma in 3D", "[longrange]")
{
  const int ndim    = 3;
  const double alat = 1.0;
  const double blat = 2 * M_PI / alat;

  // check first 3 shells of kvectors
  const std::vector<double> kcs = {blat, std::sqrt(2) * blat, std::sqrt(3) * blat};
  const std::vector<int> nks    = {6, 18, 26};

  Lattice lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik = 0; ik < kcs.size(); ik++)
  {
    const double kc = kcs[ik] + 1e-6;
    klists.updateKLists(lattice, kc,
                                                                                                         ndim);
    CHECK(klists.getKpts().size() == nks[ik]);
  }
  const int mxk    = klists.getKpts().size();
  int gvecs[26][3] = {{-1, 0, 0},  {0, -1, 0},  {0, 0, -1}, {0, 0, 1},   {0, 1, 0},    {1, 0, 0},   {-1, -1, 0},
                      {-1, 0, -1}, {-1, 0, 1},  {-1, 1, 0}, {0, -1, -1}, {0, -1, 1},   {0, 1, -1},  {0, 1, 1},
                      {1, -1, 0},  {1, 0, -1},  {1, 0, 1},  {1, 1, 0},   {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1},
                      {-1, 1, 1},  {1, -1, -1}, {1, -1, 1}, {1, 1, -1},  {1, 1, 1}};

  for (int ik = 0; ik < mxk; ik++)
  {
    for (int ldim = 0; ldim < ndim; ldim++)
    {
      CHECK(klists.getKpts()[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.getKpts()[ik][ldim] * blat == Approx(klists.getKptsCartWorking()[ik][ldim]));
    }
  }
}

TEST_CASE("kcontainer at twist in 3D", "[longrange]")
{
  const int ndim    = 3;
  const double alat = 1.0;
  const double blat = 2 * M_PI / alat;

  // twist one shell of kvectors
  const double kc = blat + 1e-6;

  Lattice lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;

  PosType twist;
  twist[0] = 0.1;
  klists.updateKLists(lattice, kc,
                                                                                                       ndim, twist);
  CHECK(klists.getKpts().size() == 1);
  CHECK(klists.getKptsCartWorking()[0][0] == Approx(blat * (twist[0] - 1)));

  twist = {-0.5, 0, 0.5};
  klists.updateKLists(lattice, kc,
                                                                                                       ndim, twist);
  int gvecs[3][3] = {{0, 0, -1}, {1, 0, -1}, {1, 0, 0}};
  CHECK(klists.getKpts().size() == 3);
  for (int ik = 0; ik < klists.getKpts().size(); ik++)
    for (int ldim = 0; ldim < ndim; ldim++)
      CHECK(klists.getKptsCartWorking()[ik][ldim] == Approx(blat * (twist[ldim] + gvecs[ik][ldim])));
}

TEST_CASE("kcontainer at gamma in 2D", "[longrange]")
{
  const int ndim    = 2;
  const double alat = 1.0;
  const double blat = 2 * M_PI / alat;

  // check first 3 shells of kvectors
  const std::vector<double> kcs = {blat, std::sqrt(2) * blat, 2 * blat};
  const std::vector<int> nks    = {4, 8, 12};

  Lattice lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik = 0; ik < kcs.size(); ik++)
  {
    const double kc = kcs[ik] + 1e-6;
    klists.updateKLists(lattice, kc,
                                                                                                         ndim);
    CHECK(klists.getKpts().size() == nks[ik]);
  }
  const int mxk    = klists.getKpts().size();
  int gvecs[12][3] = {
      {-1, 0, 0}, {0, -1, 0}, {0, 1, 0},  {1, 0, 0},  {-1, -1, 0}, {-1, 1, 0},
      {1, -1, 0}, {1, 1, 0},  {-2, 0, 0}, {0, -2, 0}, {0, 2, 0},   {2, 0, 0},
  };

  for (int ik = 0; ik < mxk; ik++)
  {
    for (int ldim = 0; ldim < ndim; ldim++)
    {
      CHECK(klists.getKpts()[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.getKpts()[ik][ldim] * blat == Approx(klists.getKptsCartWorking()[ik][ldim]));
    }
  }
}

TEST_CASE("kcontainer at twist in 2D", "[longrange]")
{
  const int ndim    = 2;
  const double alat = 1.0;
  const double blat = 2 * M_PI / alat;

  // twist one shell of kvectors
  const double kc = blat + 1e-6;

  Lattice lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;

  PosType twist;
  twist[0] = 0.1;
  klists.updateKLists(lattice, kc,
                                                                                                       ndim, twist);
  CHECK(klists.getKpts().size() == 1);
  CHECK(klists.getKptsCartWorking()[0][0] == Approx(blat * (twist[0] - 1)));

  twist[1] = 0.1;
  klists.updateKLists(lattice, kc,
                                                                                                       ndim, twist);
  CHECK(klists.getKpts().size() == 2);
  CHECK(klists.getKptsCartWorking()[0][0] == Approx(blat * (twist[0] - 1)));
  CHECK(klists.getKptsCartWorking()[0][1] == Approx(blat * twist[1]));
  CHECK(klists.getKptsCartWorking()[1][0] == Approx(blat * (twist[0])));
  CHECK(klists.getKptsCartWorking()[1][1] == Approx(blat * (twist[1] - 1)));

  twist = {-0.5, 0.5, 0};
  klists.updateKLists(lattice, kc,
                                                                                                       ndim, twist);
  CHECK(klists.getKpts().size() == 3);
  //for (int ik=0;ik<3;ik++)
  //  app_log() << klists.getKptsCartWorking()[ik] << std::endl;
  CHECK(klists.getKptsCartWorking()[0][0] == Approx(blat * (twist[0] - 0)));
  CHECK(klists.getKptsCartWorking()[0][1] == Approx(blat * (twist[1] - 1)));
  CHECK(klists.getKptsCartWorking()[1][0] == Approx(blat * (twist[0] + 1)));
  CHECK(klists.getKptsCartWorking()[1][1] == Approx(blat * (twist[1] - 1)));
  CHECK(klists.getKptsCartWorking()[2][0] == Approx(blat * (twist[0] + 1)));
  CHECK(klists.getKptsCartWorking()[2][1] == Approx(blat * twist[1]));
}

TEST_CASE("kcontainer for diamond", "[longrange]")
{
  int ndim = 3;

  using Real         = QMCTraits::RealType;
  using FullPrecReal = QMCTraits::FullPrecRealType;

  Lattice lattice;
  lattice.BoxBConds = {true, true, true};
  Tensor<double, 3> lattice_mat{3.37316115, 3.37316115, 0.00000000, 0.00000000, 3.37316115,
                                3.37316115, 3.37316115, 0.00000000, 3.37316115};
  lattice.set(lattice_mat);
  lattice.LR_dim_cutoff = 15;
  SimulationCell cell(lattice);

  const KContainer& klists = cell.getKLists();
  std::remove_cv_t<std::remove_reference_t<decltype(KContainer().getKpts())>> kpoint_lists = {
      {-1, -1, -1}, {-1, 0, 0},   {0, -1, 0},   {0, 0, -1},   {0, 0, 1},    {0, 1, 0},    {1, 0, 0},    {1, 1, 1},
      {-1, -1, 0},  {-1, 0, -1},  {0, -1, -1},  {0, 1, 1},    {1, 0, 1},    {1, 1, 0},    {-2, -1, -1}, {-1, -2, -1},
      {-1, -1, -2}, {-1, 0, 1},   {-1, 1, 0},   {0, -1, 1},   {0, 1, -1},   {1, -1, 0},   {1, 0, -1},   {1, 1, 2},
      {1, 2, 1},    {2, 1, 1},    {-2, -2, -1}, {-2, -1, -2}, {-2, -1, 0},  {-2, 0, -1},  {-1, -2, -2}, {-1, -2, 0},
      {-1, -1, 1},  {-1, 0, -2},  {-1, 1, -1},  {-1, 1, 1},   {0, -2, -1},  {0, -1, -2},  {0, 1, 2},    {0, 2, 1},
      {1, -1, -1},  {1, -1, 1},   {1, 0, 2},    {1, 1, -1},   {1, 2, 0},    {1, 2, 2},    {2, 0, 1},    {2, 1, 0},
      {2, 1, 2},    {2, 2, 1},    {-2, -2, -2}, {-2, 0, 0},   {0, -2, 0},   {0, 0, -2},   {0, 0, 2},    {0, 2, 0},
      {2, 0, 0},    {2, 2, 2},    {-2, -2, 0},  {-2, 0, -2},  {0, -2, -2},  {0, 2, 2},    {2, 0, 2},    {2, 2, 0},
      {-3, -2, -2}, {-3, -1, -1}, {-2, -3, -2}, {-2, -2, -3}, {-2, 0, 1},   {-2, 1, 0},   {-1, -3, -1}, {-1, -1, -3},
      {-1, 0, 2},   {-1, 2, 0},   {0, -2, 1},   {0, -1, 2},   {0, 1, -2},   {0, 2, -1},   {1, -2, 0},   {1, 0, -2},
      {1, 1, 3},    {1, 3, 1},    {2, -1, 0},   {2, 0, -1},   {2, 2, 3},    {2, 3, 2},    {3, 1, 1},    {3, 2, 2},
      {-3, -2, -1}, {-3, -1, -2}, {-2, -3, -1}, {-2, -1, -3}, {-2, -1, 1},  {-2, 1, -1},  {-1, -3, -2}, {-1, -2, -3},
      {-1, -2, 1},  {-1, 1, -2},  {-1, 1, 2},   {-1, 2, 1},   {1, -2, -1},  {1, -1, -2},  {1, -1, 2},   {1, 2, -1},
      {1, 2, 3},    {1, 3, 2},    {2, -1, 1},   {2, 1, -1},   {2, 1, 3},    {2, 3, 1},    {3, 1, 2},    {3, 2, 1},
      {-3, -3, -2}, {-3, -2, -3}, {-3, -1, 0},  {-3, 0, -1},  {-2, -3, -3}, {-2, 1, 1},   {-1, -3, 0},  {-1, -1, 2},
      {-1, 0, -3},  {-1, 2, -1},  {0, -3, -1},  {0, -1, -3},  {0, 1, 3},    {0, 3, 1},    {1, -2, 1},   {1, 0, 3},
      {1, 1, -2},   {1, 3, 0},    {2, -1, -1},  {2, 3, 3},    {3, 0, 1},    {3, 1, 0},    {3, 2, 3},    {3, 3, 2},
      {-3, -3, -3}, {-3, -3, -1}, {-3, -2, 0},  {-3, -1, -3}, {-3, 0, -2},  {-3, 0, 0},   {-2, -3, 0},  {-2, -2, 1},
      {-2, 0, -3},  {-2, 1, -2},  {-1, -3, -3}, {-1, 2, 2},   {0, -3, -2},  {0, -3, 0},   {0, -2, -3},  {0, 0, -3},
      {0, 0, 3},    {0, 2, 3},    {0, 3, 0},    {0, 3, 2},    {1, -2, -2},  {1, 3, 3},    {2, -1, 2},   {2, 0, 3},
      {2, 2, -1},   {2, 3, 0},    {3, 0, 0},    {3, 0, 2},    {3, 1, 3},    {3, 2, 0},    {3, 3, 1},    {3, 3, 3},
      {-4, -2, -2}, {-2, -4, -2}, {-2, -2, -4}, {-2, 0, 2},   {-2, 2, 0},   {0, -2, 2},   {0, 2, -2},   {2, -2, 0},
      {2, 0, -2},   {2, 2, 4},    {2, 4, 2},    {4, 2, 2},    {-4, -3, -2}, {-4, -2, -3}, {-4, -2, -1}, {-4, -1, -2},
      {-3, -4, -2}, {-3, -2, -4}, {-3, -1, 1},  {-3, 1, -1},  {-2, -4, -3}, {-2, -4, -1}, {-2, -3, -4}, {-2, -1, -4},
      {-2, -1, 2},  {-2, 1, 2},   {-2, 2, -1},  {-2, 2, 1},   {-1, -4, -2}, {-1, -3, 1},  {-1, -2, -4}, {-1, -2, 2},
      {-1, 1, -3},  {-1, 1, 3},   {-1, 2, -2},  {-1, 3, 1},   {1, -3, -1},  {1, -2, 2},   {1, -1, -3},  {1, -1, 3},
      {1, 2, -2},   {1, 2, 4},    {1, 3, -1},   {1, 4, 2},    {2, -2, -1},  {2, -2, 1},   {2, -1, -2},  {2, 1, -2},
      {2, 1, 4},    {2, 3, 4},    {2, 4, 1},    {2, 4, 3},    {3, -1, 1},   {3, 1, -1},   {3, 2, 4},    {3, 4, 2},
      {4, 1, 2},    {4, 2, 1},    {4, 2, 3},    {4, 3, 2},    {-4, -3, -3}, {-4, -1, -1}, {-3, -4, -3}, {-3, -3, -4},
      {-3, -3, 0},  {-3, 0, -3},  {-3, 0, 1},   {-3, 1, 0},   {-1, -4, -1}, {-1, -1, -4}, {-1, 0, 3},   {-1, 3, 0},
      {0, -3, -3},  {0, -3, 1},   {0, -1, 3},   {0, 1, -3},   {0, 3, -1},   {0, 3, 3},    {1, -3, 0},   {1, 0, -3},
      {1, 1, 4},    {1, 4, 1},    {3, -1, 0},   {3, 0, -1},   {3, 0, 3},    {3, 3, 0},    {3, 3, 4},    {3, 4, 3},
      {4, 1, 1},    {4, 3, 3},    {-4, -3, -1}, {-4, -1, -3}, {-3, -4, -1}, {-3, -2, 1},  {-3, -1, -4}, {-3, 1, -2},
      {-2, -3, 1},  {-2, 1, -3},  {-1, -4, -3}, {-1, -3, -4}, {-1, 2, 3},   {-1, 3, 2},   {1, -3, -2},  {1, -2, -3},
      {1, 3, 4},    {1, 4, 3},    {2, -1, 3},   {2, 3, -1},   {3, -1, 2},   {3, 1, 4},    {3, 2, -1},   {3, 4, 1},
      {4, 1, 3},    {4, 3, 1},    {-4, -4, -3}, {-4, -3, -4}, {-4, -1, 0},  {-4, 0, -1},  {-3, -4, -4}, {-3, 1, 1},
      {-1, -4, 0},  {-1, -1, 3},  {-1, 0, -4},  {-1, 3, -1},  {0, -4, -1},  {0, -1, -4},  {0, 1, 4},    {0, 4, 1},
      {1, -3, 1},   {1, 0, 4},    {1, 1, -3},   {1, 4, 0},    {3, -1, -1},  {3, 4, 4},    {4, 0, 1},    {4, 1, 0},
      {4, 3, 4},    {4, 4, 3},    {-4, -4, -2}, {-4, -2, -4}, {-4, -2, 0},  {-4, 0, -2},  {-2, -4, -4}, {-2, -4, 0},
      {-2, -2, 2},  {-2, 0, -4},  {-2, 2, -2},  {-2, 2, 2},   {0, -4, -2},  {0, -2, -4},  {0, 2, 4},    {0, 4, 2},
      {2, -2, -2},  {2, -2, 2},   {2, 0, 4},    {2, 2, -2},   {2, 4, 0},    {2, 4, 4},    {4, 0, 2},    {4, 2, 0},
      {4, 2, 4},    {4, 4, 2},    {-4, -4, -4}, {-4, 0, 0},   {0, -4, 0},   {0, 0, -4},   {0, 0, 4},    {0, 4, 0},
      {4, 0, 0},    {4, 4, 4},    {-5, -3, -3}, {-5, -2, -2}, {-4, -4, -1}, {-4, -3, 0},  {-4, -1, -4}, {-4, 0, -3},
      {-3, -5, -3}, {-3, -4, 0},  {-3, -3, -5}, {-3, -3, 1},  {-3, 0, -4},  {-3, 0, 2},   {-3, 1, -3},  {-3, 2, 0},
      {-2, -5, -2}, {-2, -2, -5}, {-2, 0, 3},   {-2, 3, 0},   {-1, -4, -4}, {-1, 3, 3},   {0, -4, -3},  {0, -3, -4},
      {0, -3, 2},   {0, -2, 3},   {0, 2, -3},   {0, 3, -2},   {0, 3, 4},    {0, 4, 3},    {1, -3, -3},  {1, 4, 4},
      {2, -3, 0},   {2, 0, -3},   {2, 2, 5},    {2, 5, 2},    {3, -2, 0},   {3, -1, 3},   {3, 0, -2},   {3, 0, 4},
      {3, 3, -1},   {3, 3, 5},    {3, 4, 0},    {3, 5, 3},    {4, 0, 3},    {4, 1, 4},    {4, 3, 0},    {4, 4, 1},
      {5, 2, 2},    {5, 3, 3},    {-5, -3, -2}, {-5, -2, -3}, {-3, -5, -2}, {-3, -2, -5}, {-3, -1, 2},  {-3, 2, -1},
      {-2, -5, -3}, {-2, -3, -5}, {-2, 1, 3},   {-2, 3, 1},   {-1, -3, 2},  {-1, 2, -3},  {1, -2, 3},   {1, 3, -2},
      {2, -3, -1},  {2, -1, -3},  {2, 3, 5},    {2, 5, 3},    {3, -2, 1},   {3, 1, -2},   {3, 2, 5},    {3, 5, 2},
      {5, 2, 3},    {5, 3, 2},    {-5, -4, -3}, {-5, -3, -4}, {-5, -2, -1}, {-5, -1, -2}, {-4, -5, -3}, {-4, -3, -5},
      {-4, -1, 1},  {-4, 1, -1},  {-3, -5, -4}, {-3, -4, -5}, {-3, 1, 2},   {-3, 2, 1},   {-2, -5, -1}, {-2, -1, -5},
      {-2, -1, 3},  {-2, 3, -1},  {-1, -5, -2}, {-1, -4, 1},  {-1, -2, -5}, {-1, -2, 3},  {-1, 1, -4},  {-1, 1, 4},
      {-1, 3, -2},  {-1, 4, 1},   {1, -4, -1},  {1, -3, 2},   {1, -1, -4},  {1, -1, 4},   {1, 2, -3},   {1, 2, 5},
      {1, 4, -1},   {1, 5, 2},    {2, -3, 1},   {2, 1, -3},   {2, 1, 5},    {2, 5, 1},    {3, -2, -1},  {3, -1, -2},
      {3, 4, 5},    {3, 5, 4},    {4, -1, 1},   {4, 1, -1},   {4, 3, 5},    {4, 5, 3},    {5, 1, 2},    {5, 2, 1},
      {5, 3, 4},    {5, 4, 3},    {-5, -4, -4}, {-5, -4, -2}, {-5, -3, -1}, {-5, -2, -4}, {-5, -1, -3}, {-5, -1, -1},
      {-4, -5, -4}, {-4, -5, -2}, {-4, -4, -5}, {-4, -2, -5}, {-4, -2, 1},  {-4, 0, 1},   {-4, 1, -2},  {-4, 1, 0},
      {-3, -5, -1}, {-3, -2, 2},  {-3, -1, -5}, {-3, 2, -2},  {-2, -5, -4}, {-2, -4, -5}, {-2, -4, 1},  {-2, -3, 2},
      {-2, 1, -4},  {-2, 2, -3},  {-2, 2, 3},   {-2, 3, 2},   {-1, -5, -3}, {-1, -5, -1}, {-1, -3, -5}, {-1, -1, -5},
      {-1, 0, 4},   {-1, 2, 4},   {-1, 4, 0},   {-1, 4, 2},   {0, -4, 1},   {0, -1, 4},   {0, 1, -4},   {0, 4, -1},
      {1, -4, -2},  {1, -4, 0},   {1, -2, -4},  {1, 0, -4},   {1, 1, 5},    {1, 3, 5},    {1, 5, 1},    {1, 5, 3},
      {2, -3, -2},  {2, -2, -3},  {2, -2, 3},   {2, -1, 4},   {2, 3, -2},   {2, 4, -1},   {2, 4, 5},    {2, 5, 4},
      {3, -2, 2},   {3, 1, 5},    {3, 2, -2},   {3, 5, 1},    {4, -1, 0},   {4, -1, 2},   {4, 0, -1},   {4, 2, -1},
      {4, 2, 5},    {4, 4, 5},    {4, 5, 2},    {4, 5, 4},    {5, 1, 1},    {5, 1, 3},    {5, 2, 4},    {5, 3, 1},
      {5, 4, 2},    {5, 4, 4},    {-4, -4, 0},  {-4, 0, -4},  {0, -4, -4},  {0, 4, 4},    {4, 0, 4},    {4, 4, 0},
      {-5, -5, -3}, {-5, -3, -5}, {-5, -2, 0},  {-5, 0, -2},  {-3, -5, -5}, {-3, 2, 2},   {-2, -5, 0},  {-2, -2, 3},
      {-2, 0, -5},  {-2, 3, -2},  {0, -5, -2},  {0, -2, -5},  {0, 2, 5},    {0, 5, 2},    {2, -3, 2},   {2, 0, 5},
      {2, 2, -3},   {2, 5, 0},    {3, -2, -2},  {3, 5, 5},    {5, 0, 2},    {5, 2, 0},    {5, 3, 5},    {5, 5, 3},
      {-5, -5, -4}, {-5, -4, -5}, {-5, -4, -1}, {-5, -1, -4}, {-5, -1, 0},  {-5, 0, -1},  {-4, -5, -5}, {-4, -5, -1},
      {-4, -3, 1},  {-4, -1, -5}, {-4, 1, -3},  {-4, 1, 1},   {-3, -4, 1},  {-3, 1, -4},  {-1, -5, -4}, {-1, -5, 0},
      {-1, -4, -5}, {-1, -1, 4},  {-1, 0, -5},  {-1, 3, 4},   {-1, 4, -1},  {-1, 4, 3},   {0, -5, -1},  {0, -1, -5},
      {0, 1, 5},    {0, 5, 1},    {1, -4, -3},  {1, -4, 1},   {1, -3, -4},  {1, 0, 5},    {1, 1, -4},   {1, 4, 5},
      {1, 5, 0},    {1, 5, 4},    {3, -1, 4},   {3, 4, -1},   {4, -1, -1},  {4, -1, 3},   {4, 1, 5},    {4, 3, -1},
      {4, 5, 1},    {4, 5, 5},    {5, 0, 1},    {5, 1, 0},    {5, 1, 4},    {5, 4, 1},    {5, 4, 5},    {5, 5, 4},
  };
  {
    INFO("Checking kpoint_lists");
    auto check = checkVector(klists.getKpts(), kpoint_lists, true);
    CHECKED_ELSE(check.result) { FAIL(check.result_message); }
  }
}


} // namespace qmcplusplus

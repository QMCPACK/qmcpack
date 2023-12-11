//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////
#include "catch.hpp"

#include "Particle/LongRange/KContainer.h"
#include "OhmmsPETE/TinyVector.h"

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

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik = 0; ik < kcs.size(); ik++)
  {
    const double kc = kcs[ik] + 1e-6;
    klists.updateKLists(lattice, kc, ndim);
    CHECK(klists.kpts.size() == nks[ik]);
  }
  const int mxk    = klists.kpts.size();
  int gvecs[26][3] = {{-1, 0, 0},  {0, -1, 0},  {0, 0, -1}, {0, 0, 1},   {0, 1, 0},    {1, 0, 0},   {-1, -1, 0},
                      {-1, 0, -1}, {-1, 0, 1},  {-1, 1, 0}, {0, -1, -1}, {0, -1, 1},   {0, 1, -1},  {0, 1, 1},
                      {1, -1, 0},  {1, 0, -1},  {1, 0, 1},  {1, 1, 0},   {-1, -1, -1}, {-1, -1, 1}, {-1, 1, -1},
                      {-1, 1, 1},  {1, -1, -1}, {1, -1, 1}, {1, 1, -1},  {1, 1, 1}};

  for (int ik = 0; ik < mxk; ik++)
  {
    for (int ldim = 0; ldim < ndim; ldim++)
    {
      CHECK(klists.kpts[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.kpts[ik][ldim] * blat == Approx(klists.kpts_cart[ik][ldim]));
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

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;

  PosType twist;
  twist[0] = 0.1;
  klists.updateKLists(lattice, kc, ndim, twist);
  CHECK(klists.kpts.size() == 1);
  CHECK(klists.kpts_cart[0][0] == Approx(blat * (twist[0] - 1)));

  twist = {-0.5, 0, 0.5};
  klists.updateKLists(lattice, kc, ndim, twist);
  int gvecs[3][3] = {{0, 0, -1}, {1, 0, -1}, {1, 0, 0}};
  CHECK(klists.kpts.size() == 3);
  for (int ik = 0; ik < klists.kpts.size(); ik++)
    for (int ldim = 0; ldim < ndim; ldim++)
      CHECK(klists.kpts_cart[ik][ldim] == Approx(blat * (twist[ldim] + gvecs[ik][ldim])));
}

TEST_CASE("kcontainer at gamma in 2D", "[longrange]")
{
  const int ndim    = 2;
  const double alat = 1.0;
  const double blat = 2 * M_PI / alat;

  // check first 3 shells of kvectors
  const std::vector<double> kcs = {blat, std::sqrt(2) * blat, 2 * blat};
  const std::vector<int> nks    = {4, 8, 12};

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik = 0; ik < kcs.size(); ik++)
  {
    const double kc = kcs[ik] + 1e-6;
    klists.updateKLists(lattice, kc, ndim);
    CHECK(klists.kpts.size() == nks[ik]);
  }
  const int mxk    = klists.kpts.size();
  int gvecs[12][3] = {
      {-1, 0, 0}, {0, -1, 0}, {0, 1, 0},  {1, 0, 0},  {-1, -1, 0}, {-1, 1, 0},
      {1, -1, 0}, {1, 1, 0},  {-2, 0, 0}, {0, -2, 0}, {0, 2, 0},   {2, 0, 0},
  };

  for (int ik = 0; ik < mxk; ik++)
  {
    for (int ldim = 0; ldim < ndim; ldim++)
    {
      CHECK(klists.kpts[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.kpts[ik][ldim] * blat == Approx(klists.kpts_cart[ik][ldim]));
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

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;

  PosType twist;
  twist[0] = 0.1;
  klists.updateKLists(lattice, kc, ndim, twist);
  CHECK(klists.kpts.size() == 1);
  CHECK(klists.kpts_cart[0][0] == Approx(blat * (twist[0] - 1)));

  twist[1] = 0.1;
  klists.updateKLists(lattice, kc, ndim, twist);
  CHECK(klists.kpts.size() == 2);
  CHECK(klists.kpts_cart[0][0] == Approx(blat * (twist[0] - 1)));
  CHECK(klists.kpts_cart[0][1] == Approx(blat * twist[1]));
  CHECK(klists.kpts_cart[1][0] == Approx(blat * (twist[0])));
  CHECK(klists.kpts_cart[1][1] == Approx(blat * (twist[1] - 1)));

  twist = {-0.5, 0.5, 0};
  klists.updateKLists(lattice, kc, ndim, twist);
  CHECK(klists.kpts.size() == 3);
  //for (int ik=0;ik<3;ik++)
  //  app_log() << klists.kpts_cart[ik] << std::endl;
  CHECK(klists.kpts_cart[0][0] == Approx(blat * (twist[0] - 0)));
  CHECK(klists.kpts_cart[0][1] == Approx(blat * (twist[1] - 1)));
  CHECK(klists.kpts_cart[1][0] == Approx(blat * (twist[0] + 1)));
  CHECK(klists.kpts_cart[1][1] == Approx(blat * (twist[1] - 1)));
  CHECK(klists.kpts_cart[2][0] == Approx(blat * (twist[0] + 1)));
  CHECK(klists.kpts_cart[2][1] == Approx(blat * twist[1]));
}

} // namespace qmcplusplus

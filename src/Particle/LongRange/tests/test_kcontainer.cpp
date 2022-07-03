#include "catch.hpp"

#include "Particle/LongRange/KContainer.h"

namespace qmcplusplus
{

TEST_CASE("kcontainer at gamma in 3D", "[longrange]")
{
  const int ndim = 3;
  const double alat = 1.0;
  const double blat = 2*M_PI/alat;

  // check first 3 shells of kvectors
  const std::vector<double> kcs = {blat, std::sqrt(2)*blat, std::sqrt(3)*blat};
  const std::vector<int> nks = {6, 18, 26};

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik=0;ik<kcs.size();ik++)
  {
    const double kc = kcs[ik]+1e-6;
    klists.updateKLists(lattice, kc, ndim);
    CHECK(klists.kpts.size() == nks[ik]);
  }
  const int mxk = klists.kpts.size();
  int gvecs[mxk][OHMMS_DIM] = {
    {-1,  0,  0},
    { 0, -1,  0},
    { 0,  0, -1},
    { 0,  0,  1},
    { 0,  1,  0},
    { 1,  0,  0},
    {-1, -1,  0},
    {-1,  0, -1},
    {-1,  0,  1},
    {-1,  1,  0},
    { 0, -1, -1},
    { 0, -1,  1},
    { 0,  1, -1},
    { 0,  1,  1},
    { 1, -1,  0},
    { 1,  0, -1},
    { 1,  0,  1},
    { 1,  1,  0},
    {-1, -1, -1},
    {-1, -1,  1},
    {-1,  1, -1},
    {-1,  1,  1},
    { 1, -1, -1},
    { 1, -1,  1},
    { 1,  1, -1},
    { 1,  1,  1}
  };

  for (int ik=0;ik<mxk;ik++)
  {
    for (int ldim=0;ldim<ndim;ldim++)
    {
      CHECK(klists.kpts[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.kpts[ik][ldim]*blat == klists.kpts_cart[ik][ldim]);
    }
  }
}

TEST_CASE("kcontainer at gamma in 2D", "[longrange]")
{
  const int ndim = 2;
  const double alat = 1.0;
  const double blat = 2*M_PI/alat;

  // check first 3 shells of kvectors
  const std::vector<double> kcs = {blat, std::sqrt(2)*blat, 2*blat};
  const std::vector<int> nks = {4, 8, 12};

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.R.diagonal(1.0);
  lattice.set(lattice.R); // compute Rv and Gv from R

  KContainer klists;
  for (int ik=0;ik<kcs.size();ik++)
  {
    const double kc = kcs[ik]+1e-6;
    klists.updateKLists(lattice, kc, ndim);
    CHECK(klists.kpts.size() == nks[ik]);
  }
  const int mxk = klists.kpts.size();
  int gvecs[mxk][OHMMS_DIM] = {
    {-1,  0,  0},
    { 0, -1,  0},
    { 0,  1,  0},
    { 1,  0,  0},
    {-1, -1,  0},
    {-1,  1,  0},
    { 1, -1,  0},
    { 1,  1,  0},
    {-2,  0,  0},
    { 0, -2,  0},
    { 0,  2,  0},
    { 2,  0,  0},
  };

  for (int ik=0;ik<mxk;ik++)
  {
    for (int ldim=0;ldim<ndim;ldim++)
    {
      CHECK(klists.kpts[ik][ldim] == gvecs[ik][ldim]);
      CHECK(klists.kpts[ik][ldim]*blat == klists.kpts_cart[ik][ldim]);
    }
  }
}

} // qmcplusplus

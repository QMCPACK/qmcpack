//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file distancetables_soa.cpp
 *
 * Test code for the accuracy of AoS to SoA transformation of distance tables.
 */
#include <Configuration.h>
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTable.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "random.hpp"
#include "mpi/collectives.h"
#include "Sandbox/input.hpp"
#include "Sandbox/pseudo.hpp"
#include "Utilities/Timer.h"
#include "Sandbox/common.hpp"
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller->initialize(env);
#endif
  Communicate* myComm = OHMMS::Controller;
  if (OHMMS::Controller->rank() != 0)
  {
    outputManager.shutOff();
  }

  bool ionode = (myComm->rank() == 0);
  int na      = 4;
  int nb      = 4;
  int nc      = 1;
  int nsteps  = 1;
  int iseed   = 11;
  int nx = 48, ny = 48, nz = 60;
  //thread blocking
  //int ncrews=1; //default is 1
  int tileSize = -1;
  int ncrews   = 1;

  char* g_opt_arg;
  int opt;
  while ((opt = getopt(argc, argv, "hs:g:i:")) != -1)
  {
    switch (opt)
    {
    case 'h':
      printf("[-g \"n0 n1 n2\"]\n");
      return 1;
    case 'g': //tiling1 tiling2 tiling3
      sscanf(optarg, "%d %d %d", &na, &nb, &nc);
      break;
    case 's': //random seed
      iseed = atoi(optarg);
      break;
    case 'i': //number of iterations
      nsteps = atoi(optarg);
      break;
    }
  }

  Random.init(0, 1, iseed);

  typedef QMCTraits::RealType RealType;
  typedef ParticleSet::ParticlePos_t ParticlePos_t;
  typedef ParticleSet::ParticleLayout_t LatticeType;
  typedef ParticleSet::TensorType TensorType;
  typedef ParticleSet::PosType PosType;

  Tensor<int, 3> tmat(na, 0, 0, 0, nb, 0, 0, 0, nc);
  double t0 = 0.0, t1 = 0.0;
  constexpr OHMMS_PRECISION scale(1);

  RandomGenerator<RealType> random_th(11);

  ParticleSet ions, els;

  ions.setName("ion0");
  ions.Lattice.BoxBConds = 1;
  ions.Lattice.LR_rc     = 5;
  tile_cell(ions, tmat, scale);
  ions.setCoordinates(ions.R); //this needs to be handled internally

  const int nions = ions.getTotalNum();
  const int nels  = count_electrons(ions);
  const int nels3 = 3 * nels;

  { //create up/down electrons
    els.Lattice.BoxBConds = 1;
    els.Lattice.LR_rc     = 5;
    els.Lattice           = ions.Lattice;
    vector<int> ud(2);
    ud[0] = nels / 2;
    ud[1] = nels - ud[0];
    els.create(ud);
    els.R.InUnit = PosUnit::Lattice;
    random_th.generate_uniform(&els.R[0][0], nels3);
    els.convert2Cart(els.R);   // convert to Cartiesian
    els.setCoordinates(els.R); //this needs to be handled internally
    els.setName("e");
  }

  constexpr RealType eps = numeric_limits<float>::epsilon();

  //copy of ParticleSet for validations
  ParticleSet::ParticlePos_t Rcopy(els.R);

  const auto& d_ee = els.getDistTableAA(els.addTable(els));
  const auto& d_ie = els.getDistTableAB(els.addTable(ions));

  RealType Rsim = els.Lattice.WignerSeitzRadius;

  //SoA version does not need update if PbyP
  els.update();

  //SoA check symmetry
  double sym_err = 0.0;
  int nn         = 0;
  for (int iel = 0; iel < nels; ++iel)
    for (int jel = iel + 1; jel < nels; ++jel)
    {
      RealType dsym = std::abs(d_ee.getDistRow(jel)[iel] - d_ee.getDistRow(iel)[jel]);
      sym_err += dsym;
      ++nn;
    }
  cout << "---------------------------------" << endl;
  cout << "AA SoA(upper) - SoA(lower) distances     = " << sym_err / nn << endl;

  ParticlePos_t delta(nels);

  //main particle-by-particle update
  RealType sqrttau = 0.2;
  for (int s = 0; s < nsteps; ++s)
  {
    random_th.generate_normal(&delta[0][0], nels3);
    for (int iel = 0; iel < nels; ++iel)
    {
      PosType dr      = sqrttau * delta[iel];
      bool valid_move = els.makeMoveAndCheck(iel, dr);

      if (valid_move && Random() < 0.5)
        els.rejectMove(iel);
      else
        els.acceptMove(iel);
    }
  }

  els.donePbyP();

  {
    double r_err = 0.0;
    for (int iat = 0; iat < nels; ++iat)
    {
      PosType d  = els.R[iat] - Rcopy[iat];
      RealType r = sqrt(dot(d, d));
      r_err += r;
    }
    cout << "Done with the sweep. Diffusion |els.R-R0|^2/nels = " << r_err / nels << endl;
  }

  OHMMS::Controller->finalize();

  return 0;
}

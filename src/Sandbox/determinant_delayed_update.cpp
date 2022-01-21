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
/** @file moveonsphere.cpp
 * @brief Miniapp to capture the computation in NonLocalPP.
 *
 * Using only einspine SPO part + Jastrow as a wavefunction.
 */
#include <Configuration.h>
#include "Utilities/PrimeNumberSet.h"
#include "Utilities/Timer.h"
#include "random.hpp"
#include "mpi/collectives.h"
#include <getopt.h>
using namespace std;
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/simd.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/DelayedUpdate.h"
using namespace qmcplusplus;

template<typename RNG, typename T>
inline void generate(RNG& rng, T* restrict data, size_t n)
{
  constexpr T shift(0.5);
  std::generate(data, data + n, rng);
  for (int i = 0; i < n; ++i)
    data[i] -= shift;
}

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller->initialize(env);
#endif
  Communicate* myComm = OHMMS::Controller;

  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::ValueType ValueType;
#if defined(QMC_COMPLEX)
  typedef std::complex<OHMMS_PRECISION_FULL> mValueType;
#else
  typedef OHMMS_PRECISION_FULL mValueType;
#endif
  //use the global generator

  bool ionode   = (myComm->rank() == 0);
  int nels      = 128;
  int iseed     = 11;
  int nsteps    = 100;
  int ncrews    = 1;
  int nsubsteps = 1;
  int delay     = 4;

  PrimeNumberSet<uint32_t> myPrimes;

  bool debug = false;
  char* g_opt_arg;
  int opt;
  while ((opt = getopt(argc, argv, "hdn:i:c:k:s:")) != -1)
  {
    switch (opt)
    {
    case 'h':
      printf("[-n int=64]\n");
      return 1;
    case 'd': //debug
      debug = true;
      break;
    case 'n': //number of MC steps
      nels = atoi(optarg);
      break;
    case 'i': //number of MC steps
      nsteps = atoi(optarg);
      break;
    case 'k': //number of MC steps
      delay = atoi(optarg);
      break;
    case 's': //the number of sub steps for drift/diffusion
      nsubsteps = atoi(optarg);
      break;
    case 'c': //number of crews per team
      ncrews = atoi(optarg);
      break;
    }
  }

  Random.init(iseed);

  //turn off output
  if (omp_get_max_threads() > 1)
  {
    outputManager.pause();
  }

  Timer bigClock;
  bigClock.restart();
  double t_compute = 0.0, t_ratio = 0.0, t_accept = 0.0, error = 0.0;
  int naccepted = 0;
#pragma omp parallel reduction(+ : error, naccepted, t_ratio, t_accept, t_compute)
  {
    Timer clock, clock_mc;
    clock.restart();
    double t_compute_loc = 0.0, t_ratio_loc = 0.0, t_accept_loc = 0.0;

    const int np = omp_get_num_threads();
    const int ip = omp_get_thread_num();

    const int teamID = ip / ncrews;
    const int crewID = ip % ncrews;

    RandomGenerator random_th(myPrimes[ip]);

    Matrix<ValueType> psiM(nels, nels), psiM_inv(nels, nels);
    Vector<ValueType> psiV(nels), invRow(nels);

    DiracMatrix<ValueType> detEng;
    DelayedUpdate<ValueType, QMCTraits::QTFull::ValueType> FahyEng;
    DelayedUpdate<ValueType, QMCTraits::QTFull::ValueType> delayedEng;

    FahyEng.resize(nels, 1);
    delayedEng.resize(nels, delay);

    generate(random_th, psiM.data(), nels * nels);
    std::complex<RealType> logdet;
    detEng.invert_transpose(psiM, psiM_inv, logdet);

    if (debug)
    {
      Matrix<ValueType> psiM0(nels, nels);
      psiM0 = psiM_inv;

      ValueType ratio_0, ratio_1;
      double err = 0.0;
      for (int iel = 0; iel < nels; ++iel)
      {
        clock_mc.restart();
        generate(random_th, psiV.data(), nels);
        FahyEng.getInvRow(psiM0, iel, invRow);
        ratio_0 = simd::dot(invRow.data(), psiV.data(), invRow.size());
        delayedEng.getInvRow(psiM_inv, iel, invRow);
        ratio_1 = simd::dot(invRow.data(), psiV.data(), invRow.size());

        err += std::abs(ratio_1 - ratio_0);
        if (std::abs(ratio_0) > 0.5 * random_th())
        {
          FahyEng.acceptRow(psiM0, iel, psiV, ratio_0);
          delayedEng.acceptRow(psiM_inv, iel, psiV, ratio_1);
        }
      }
      delayedEng.updateInvMat(psiM_inv);
      error += err;
    }

    int naccepted_loc = 0;
    if (delay > 1)
    { //use delayed update
      ValueType ratio;
      for (int mc = 0; mc < nsteps; ++mc)
      {
        for (int iel = 0; iel < nels; ++iel)
        {
          generate(random_th, psiV.data(), nels);
          clock_mc.restart();
          delayedEng.getInvRow(psiM_inv, iel, invRow);
          ratio = simd::dot(invRow.data(), psiV.data(), invRow.size());
          t_ratio_loc += clock_mc.elapsed();

          if (std::abs(ratio) > 0.5 * random_th())
          {
            naccepted_loc++;
            clock_mc.restart();
            delayedEng.acceptRow(psiM_inv, iel, psiV, ratio);
            t_accept_loc += clock_mc.elapsed();
          }
        }
        if (delayedEng.getDelayCount() > 0)
        {
          clock_mc.restart();
          delayedEng.updateInvMat(psiM_inv);
          t_accept_loc += clock_mc.elapsed();
        }
      }
    }
    else
    {
      ValueType ratio;
      for (int mc = 0; mc < nsteps; ++mc)
      {
        for (int iel = 0; iel < nels; ++iel)
        {
          generate(random_th, psiV.data(), nels);
          clock_mc.restart();
          FahyEng.getInvRow(psiM_inv, iel, invRow);
          ratio = simd::dot(invRow.data(), psiV.data(), invRow.size());
          t_ratio_loc += clock_mc.elapsed();
          if (std::abs(ratio) > 0.5 * random_th())
          {
            naccepted_loc++;
            clock_mc.restart();
            FahyEng.acceptRow(psiM_inv, iel, psiV, ratio);
            t_accept_loc += clock_mc.elapsed();
          }
        }
      }
    }

    naccepted += naccepted_loc;
    t_compute += t_compute_loc;
    t_ratio += t_ratio_loc;
    t_accept += t_accept_loc;

    //using both float/float
    //DiracDet<OHMMS_PRECISION,OHMMS_PRECISION> det_lp(nels);
    //det_lp.initialize(random_th);
    //if(debug) det_lp.debug(random_th);

  } //end of omp parallel

  int nthreads   = omp_get_max_threads();
  double omp_fac = 1.0 / nthreads;

  t_compute *= omp_fac;
  t_ratio *= omp_fac;
  t_accept *= omp_fac;

  cout.setf(std::ios::scientific, std::ios::floatfield);
  cout.precision(4);

  int nthreads_nested = 1;
#pragma omp parallel
  {
#pragma omp master
    nthreads_nested = omp_get_max_threads();
  }

  if (myComm->rank() == 0)
  {
    cout << "# determinant " << nels << " rank " << delay << " Total accepted " << naccepted << " /" << nels * nsteps
         << " " << naccepted / static_cast<double>(nels * nsteps) << " error " << error * omp_fac << endl;
    cout << "# N K MPI OMP-walker OMP-det T_accept T_ratio T_total T_accept/call T_ratio/call T_total/step " << endl;
    cout << "Det " << nels << " " << delay << " " << myComm->size() << " " << nthreads << " " << nthreads_nested << " "
         << t_accept << " " << t_ratio << " " << (t_ratio + t_accept) << " " << t_accept / naccepted << " "
         << t_ratio / (nsteps * nels) << " " << (t_ratio + t_accept) / (nsteps / nsubsteps) << endl;
  }
  //t_diffusion*=1.0/static_cast<double>(nsteps*nsubsteps*nthreads);
  //t_pseudo   *=1.0/static_cast<double>(nsteps*nthreads);
  //cout << "#per MC step steps " << nsteps << " substeps " << nsubsteps << endl;
  //cout << "diffusion_mc " << t_diffusion << " pseudo_mc  " << t_pseudo << endl;

  OHMMS::Controller->finalize();

  return 0;
}

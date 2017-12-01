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
#include <Utilities/PrimeNumberSet.h>
#include <Utilities/Timer.h>
#include <random/random.hpp>
#include <mpi/collectives.h>
#include <getopt.h>
using namespace std;
#include <simd/allocator.hpp>
#include <simd/simd.hpp>
#include <simd/algorithm.hpp>
#include <QMCWaveFunctions/Fermion/DiracMatrix.h>
using namespace qmcplusplus;

template<typename RNG, typename T>
inline void generate(RNG& rng, T* restrict data, size_t n)
{
  constexpr T shift(0.5);
  rng.generate_uniform(data,n);
  for(int i=0; i<n; ++i) data[i]-=shift;
}

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
  Communicate* mycomm=OHMMS::Controller;

  typedef OHMMS_PRECISION REAL_T;

  //use the global generator

  bool ionode=(mycomm->rank() == 0);
  int nels=8;
  int iseed=11;
  int nsteps=100;
  int ncrews=1;
  int nsubsteps=1;
  int delay=4;

  PrimeNumberSet<uint32_t> myPrimes;

  bool debug=false;
  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hdn:i:c:k:s:")) != -1)
  {
    switch(opt)
    {
      case 'h':
        printf("[-n int=64]\n");
        return 1;
      case 'd': //debug
        debug=true;
        break;
      case 'n': //number of MC steps
        nels=atoi(optarg);
        break;
      case 'i': //number of MC steps
        nsteps=atoi(optarg);
        break;
      case 'k': //number of MC steps
        delay=atoi(optarg);
        break;
      case 's'://the number of sub steps for drift/diffusion
        nsubsteps=atoi(optarg);
        break;
      case 'c'://number of crews per team
        ncrews=atoi(optarg);
        break;
    }
  }

  Random.init(0,1,iseed);

  //turn off output
  if(omp_get_max_threads()>1)
  {
    OhmmsInfo::Log->turnoff();
    OhmmsInfo::Warn->turnoff();
  }

  Timer bigClock;
  bigClock.restart();
  double t_compute=0.0, t_ratio=0.0, t_accept=0.0, error=0.0;
  int naccepted=0;
#pragma omp parallel reduction(+:t_compute, t_ratio, t_accept, error, naccepted)
  {
    Timer clock, clock_mc;
    clock.restart();
    double t_compute_loc=0.0, t_ratio_loc=0.0, t_accept_loc=0.0;

    const int np=omp_get_num_threads();
    const int ip=omp_get_thread_num();

    const int teamID=ip/ncrews;
    const int crewID=ip%ncrews;

    RandomGenerator<REAL_T> random_th(myPrimes[ip]);

    Matrix<REAL_T> psiM(nels,nels),psiM_inv(nels,nels);
    Vector<REAL_T> psiV(nels), Ainv_row(nels);

    DiracMatrix<REAL_T> detEng;
    DelayedUpdate<REAL_T,double> delayedEng;

    delayedEng.resize(nels,delay);

    generate(random_th,psiM.data(),nels*nels);
    simd::transpose(psiM.data(),nels,nels,psiM_inv.data(),nels,nels);
    detEng.invert(psiM_inv,true);


    if(debug)
    {
      REAL_T ratio_0, ratio_1;
      double err=0.0;
      for(int iel=0; iel<nels; ++iel)
      {
        clock_mc.restart();
        generate(random_th, psiV.data(), nels);
        ratio_0=simd::dot(psiM_inv[iel],psiV.data(),nels);

        delayedEng.getInvRow(psiM_inv,iel,Ainv_row.data());
        ratio_1=simd::dot(Ainv_row.data(),psiV.data(),nels);
        err += std::abs(ratio_1-ratio_0);
        if(ratio_0>0 && ratio_0>0.5*random_th())
        {
          detEng.updateRow(psiM_inv,psiV.data(),iel,ratio_0);
          delayedEng.acceptRow(psiM_inv,psiV.data(),iel);
        }
      }
      error += err;
    }

    int naccepted_loc=0;
    if(delay>1)
    {//use delayed update
      REAL_T ratio;
      for(int mc=0; mc<nsteps; ++mc)
      {
        for(int iel=0; iel<nels; ++iel)
        {
          generate(random_th, psiV.data(), nels);
          clock_mc.restart();
          delayedEng.getInvRow(psiM_inv,iel,Ainv_row.data());
          ratio=simd::dot(Ainv_row.data(),psiV.data(),nels);
          t_ratio_loc+=clock_mc.elapsed();

          if(ratio>0 && ratio>0.5*random_th())
          {
            naccepted_loc++;
            clock_mc.restart();
            delayedEng.acceptRow(psiM_inv,psiV.data(),iel);
            t_accept_loc+=clock_mc.elapsed();
          }
        }
        delayedEng.updateInvMat(psiM_inv);
      }
    }
    else
    {
      REAL_T ratio;
      for(int mc=0; mc<nsteps; ++mc)
      {
        for(int iel=0; iel<nels; ++iel)
        {
          generate(random_th, psiV.data(), nels);
          clock_mc.restart();
          ratio=simd::dot(psiM_inv[iel],psiV.data(),nels);
          t_ratio_loc+=clock_mc.elapsed();
          if(ratio>0 && ratio>0.5*random_th())
          {
            naccepted_loc++;
            clock_mc.restart();
            detEng.updateRow(psiM_inv,psiV.data(),iel,ratio);
            t_accept_loc+=clock_mc.elapsed();
          }
        }
      }
    }

    naccepted+= naccepted_loc;
    t_compute+= t_compute_loc;
    t_ratio  += t_ratio_loc;
    t_accept += t_accept_loc;

    //using both float/float
    //DiracDet<OHMMS_PRECISION,OHMMS_PRECISION> det_lp(nels);
    //det_lp.initialize(random_th);
    //if(debug) det_lp.debug(random_th);

  } //end of omp parallel

  int nthreads=omp_get_max_threads();
  double omp_fac=1.0/nthreads;

  t_compute*=omp_fac;
  t_ratio  *=omp_fac;
  t_accept *=omp_fac;

  cout.setf(std::ios::scientific, std::ios::floatfield);
  cout.precision(4);

  cout << "determinant " << nels << " rank " << delay << " Total accepted " << naccepted << " /" << nels*nsteps << " " 
    << naccepted/static_cast<double>(nels*nsteps) << " error " << error*omp_fac << endl;
  cout << "Total " << (t_ratio+t_accept) << " ratio " << t_ratio  << " accept " << t_accept << endl;
  cout << "Per   " << (t_ratio+t_accept)/(nsteps/nsubsteps) << " ratio " << t_ratio/(nsteps*nels)  << " accept " << t_accept/naccepted << endl;
  //t_diffusion*=1.0/static_cast<double>(nsteps*nsubsteps*nthreads);
  //t_pseudo   *=1.0/static_cast<double>(nsteps*nthreads);
  //cout << "#per MC step steps " << nsteps << " substeps " << nsubsteps << endl;
  //cout << "diffusion_mc " << t_diffusion << " pseudo_mc  " << t_pseudo << endl;

  OHMMS::Controller->finalize();

  return 0;
}

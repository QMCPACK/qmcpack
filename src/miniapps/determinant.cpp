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
#include <miniapps/determinant.hpp>
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  if (OHMMS::Controller->rank()!=0) {
    outputManager.shutOff();
  }
  Communicate* mycomm=OHMMS::Controller;

  //use the global generator

  bool ionode=(mycomm->rank() == 0);
  int nels=8;
  int iseed=11;
  int nsteps=100;
  int ncrews=1;
  int nsubsteps=1;

  PrimeNumberSet<uint32_t> myPrimes;

  bool debug=false;
  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hdn:i:c:s:")) != -1)
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
    outputManager.pause();
  }

  Timer bigClock;
  bigClock.restart();
  double t_compute=0.0, t_ratio=0.0, t_accept=0.0;
  int naccepted=0;
#pragma omp parallel reduction(+:t_compute, t_ratio, t_accept, naccepted)
  {
    Timer clock, clock_mc;
    clock.restart();
    double t_compute_loc=0.0, t_ratio_loc=0.0, t_accept_loc=0.0;

    const int np=omp_get_num_threads();
    const int ip=omp_get_thread_num();

    const int teamID=ip/ncrews;
    const int crewID=ip%ncrews;

    RandomGenerator<OHMMS_PRECISION> random_th(myPrimes[ip]);

    DiracDet<OHMMS_PRECISION> det(nels);
    det.initialize(random_th);
    if(debug) det.debug();

    int naccepted_loc=0;
    for(int mc=0; mc<nsteps; ++mc)
    {
      if(mc%nsubsteps ==0)
      {
        clock_mc.restart();
        det.recompute();
        t_compute_loc+=clock_mc.elapsed();
      }
      for(int iel=0; iel<nels; ++iel)
      {
        clock_mc.restart();
        auto ratio=det.ratio(iel);
        t_ratio_loc+=clock_mc.elapsed();
        if(ratio>0 && ratio>0.5*random_th())
        {
          naccepted_loc++;
          clock_mc.restart();
          det.accept(iel);
          t_accept_loc+=clock_mc.elapsed();
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

  cout << "determinant " << nels << " Total accepted " << naccepted << " /" << nels*nsteps << " " 
    << naccepted/static_cast<double>(nels*nsteps) << endl;
  cout << "Total recompute " << t_compute << " ratio " << t_ratio  << " accept " << t_accept << endl;
  cout << "Per   recompute " << t_compute/(nsteps/nsubsteps) << " ratio " << t_ratio/(nsteps*nels)  << " accept " << t_accept/naccepted << endl;
  //t_diffusion*=1.0/static_cast<double>(nsteps*nsubsteps*nthreads);
  //t_pseudo   *=1.0/static_cast<double>(nsteps*nthreads);
  //cout << "#per MC step steps " << nsteps << " substeps " << nsubsteps << endl;
  //cout << "diffusion_mc " << t_diffusion << " pseudo_mc  " << t_pseudo << endl;

  OHMMS::Controller->finalize();

  return 0;
}

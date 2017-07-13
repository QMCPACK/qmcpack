//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file einspline_smp.cpp
 * @brief Benchmark einspline. 
 */
// Affinity for this run is to be set, KMP_AFFINITY=compact,granularity=thread,verbose

#include <Configuration.h>
#include <Utilities/OhmmsInfo.h>
#include <Message/Communicate.h>
#include <mpi/collectives.h>
#include <Utilities/Timer.h>
#include <spline2/tests/einspline_shared_woNested.h>
#include <getopt.h>
#include <stdio.h>
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("bspline_tiledThreading",mycomm->rank());

  int nx=48,ny=48,nz=48;
  int num_splines=128;
  int nsamples=512;
  int niters=10;
  int nThreadsLevel2 =1;

  int tileSize = -1; 

  int opt;
  while((opt = getopt(argc, argv, "hg:x:y:z:i:s:p:t:a:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations\n");
      return 1;
    case 'g':
      nx=ny=nz=atoi(optarg);
      break;
    case 'x':
      nx=atoi(optarg);
      break;
    case 'y':
      ny=atoi(optarg);
      break;
    case 'z':
      nz=atoi(optarg);
      break;
    case 's':
      num_splines=atoi(optarg);
      break;
    case 'p':
      nsamples=atoi(optarg);
      break;
    case 'i':
      niters=atoi(optarg);
      break;
    case 't':
      nThreadsLevel2=atoi(optarg);
      break;
    case 'a':
      tileSize=atoi(optarg);
      break;
    }
  }

  if ( tileSize == -1 ) {
    tileSize = num_splines;
  }
  if ( tileSize > num_splines || num_splines%tileSize != 0 ) {
    cout << "\nTile size should be less than and multiple of num_splines. " << tileSize << " " << num_splines << endl;
    exit(0);
  } 
  int numTiles = num_splines/tileSize;
  cout << "Number of tiles " << numTiles << ", tile size " << tileSize << endl;

  typedef TinyVector<double,3> timer_type;
  timer_type s_timer_t;
  Timer big_clock;

  if ( num_splines%numTiles != 0 ) {
    cout << "\nNumber of splines should be divisible by numTiles " << numTiles << "\n" ;
    exit(0);
  }

  omp_set_nested(false);
  int nThreadsLevel1 = omp_get_max_threads()/nThreadsLevel2;
  if ( omp_get_max_threads() % nThreadsLevel2 != 0 ) {
    cout << "\nNumber of level 2 threads should be multiple. \n";
    exit(0);
  }
  std::cout << "\nTotal Threads = " << omp_get_max_threads() << ", At level 1: " << nThreadsLevel1 << ", At level 2: " << nThreadsLevel2 << std::endl;

  int factor = numTiles/nThreadsLevel2;
  if ( factor*nThreadsLevel2 != numTiles ) {
    std::cout << "\nNumber of tiles should be multiple of number of threads at level 2. " << std::endl;
    exit(0);
  }
  if ( numTiles < nThreadsLevel2) {
    std::cout << "\nNumber of tiles should be greater than number of threads at level 2." << std::endl;
    exit(0);
  }

  // create copies of benchamark to run in different teams.
  einspline3d_benchmark<float> s_bench[numTiles];
  for ( int index = 0; index < numTiles; index++ ) {
    s_bench[index].set(nx, ny, nz, num_splines/numTiles, 1);
  }
  einspline3d_benchmark<float> s_benchTeam[nThreadsLevel1][numTiles];
  #pragma omp parallel for collapse(2)
  for ( int i = 0; i < nThreadsLevel1; i++ ) {
    for ( int index = 0; index < numTiles; index++ ) {
      s_benchTeam[i][index] = s_bench[index];
    }
  } 

  double t_init=big_clock.elapsed();
  cout << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
    << " num_splines = " << num_splines << " num_samples = " << nsamples
    << " iterations = " << niters
    << " number of operations in millions " << endl;
  cout << "#   mpi   openmp    datatype     "
    << "     value       vgl         vgh_op  " << endl;

  std::flush(std::cout);

  double am_total_time = 0.0;
  big_clock.restart();
  #pragma omp parallel reduction (+: am_total_time )
  {
    int iThreadID = omp_get_thread_num();
    //For KNC with 240 threads - The objective is to make 60 teams each with 4 threads (called creews)
    int iTeamID = iThreadID/nThreadsLevel2;   // 5/4=1 will go to team 1, 120/4=30 will go to team number 30
    int iCrewID = iThreadID%nThreadsLevel2;   // 6%2=2 will be crew number 2, 239%2=3 will be crew number 3.

    random_position_generator<float> s_pos(nsamples,iTeamID);
    timer_type s_timer(0.0);

    Timer am_big_clock; 
    am_big_clock.restart();    
    for(int i=0; i<niters; ++i)
    {
      s_pos.randomize();
      int iBlockingFactor = numTiles/nThreadsLevel2;
      for ( int index = 0; index < numTiles; index++ ) {

        int iExecutionID = index/iBlockingFactor;
        if ( iExecutionID == iCrewID ) {
#ifdef BSPLINE_EVAL_VALUE
          s_timer+=s_benchTeam[iTeamID][index].test_onlyValue(s_pos.Vpos, s_pos.VGLpos, s_pos.VGHpos, 0 ); 
#else
          s_timer+=s_benchTeam[iTeamID][index].test_all(s_pos.Vpos, s_pos.VGLpos, s_pos.VGHpos, 0 ); 
#endif
        }
      }
    }
    am_total_time += am_big_clock.elapsed();
#pragma omp critical
    s_timer_t += s_timer;
  }
  cout << endl;
  double t_comp=big_clock.elapsed();

  mpi::reduce(*mycomm,s_timer_t);

  if(mycomm->rank()==0)
  {
    double nops=double(nThreadsLevel1*mycomm->size())*double(num_splines)*double(nsamples)*1.e-6;
    double tfac=1.0/static_cast<double>(mycomm->size()*omp_get_max_threads()*niters);
    s_timer_t*=tfac;
    cout.setf(std::ios::scientific, std::ios::floatfield);
    cout.precision(6);
    cout << argv[0] << setw(4) << mycomm->size() << setw(4) << nThreadsLevel1 <<  " single Ops  "<< std::fixed << nops/double(s_timer_t[0]) << "   " << nops/double(s_timer_t[1]) << "   " << nops/double(s_timer_t[2])<< endl;
    cout << argv[0] << setw(4) << mycomm->size() <<  setw(4) << nThreadsLevel1 <<  " single sec  "<< std::fixed  << s_timer_t << endl;
    cout << " Initialization = " << std::fixed  << t_init << endl;
    cout << " Total time = " << std::fixed << am_total_time/omp_get_max_threads() << endl;
  }

  OHMMS::Controller->finalize();
  return 0;
}


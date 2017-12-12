//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file einspline_spo.cpp
 * @brief Miniapp to capture the computation in einspline + NonLocalPP.
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <random/random.hpp>
#include <mpi/collectives.h>
#include <miniapps/input.hpp>
#include <miniapps/pseudo.hpp>
#include <Utilities/Timer.h>
#include <miniapps/common.hpp>
#include <miniapps/einspline_spo.hpp>
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  if (OHMMS::Controller->rank() != 0)
  {
    outputManager.shutOff();
  }
  Communicate* mycomm=OHMMS::Controller;

  typedef QMCTraits::RealType           RealType;
  typedef ParticleSet::ParticlePos_t    ParticlePos_t;
  typedef ParticleSet::ParticleLayout_t LatticeType;
  typedef ParticleSet::TensorType       TensorType;
  typedef ParticleSet::PosType          PosType;

  //use the global generator

  bool ionode=(mycomm->rank() == 0);
#if USE_NIO
  int na=1; int nb=1; int nc=1;
  int nx=37,ny=37,nz=37;
#else
  int na=4; int nb=4; int nc=1;
  int nx=48,ny=48,nz=60;
#endif
  int nsteps=100;
  int iseed=11;
  //thread blocking
  //int ncrews=1; //default is 1
  int tileSize=-1;
  int ncrews=1;

  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hsg:i:b:c:a:")) != -1)
  {
    switch(opt)
    {
      case 'h':
        printf("[-g \"n0 n1 n2\"]\n");
        return 1;
      case 'g': //tiling1 tiling2 tiling3
        sscanf(optarg,"%d %d %d",&na,&nb,&nc);
        break;
      case 'i': //number of MC steps
        nsteps=atoi(optarg);
        break;
      case 's'://random seed
        iseed=atoi(optarg);
        break;
      case 'c'://number of crews per team
        ncrews=atoi(optarg);
        break;
      case 'a':
        tileSize=atoi(optarg);
        break;
    }
  }

  //Random.init(0,1,iseed);
  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);

  //turn off output
//  if(omp_get_max_threads()>1)
  {
    outputManager.pause();
  }

  int nptcl=0;
  int nknots_copy = 0;
  double t0=0.0;
  OHMMS_PRECISION ratio=0.0;

  using spo_type=einspline_spo<OHMMS_PRECISION>;
  spo_type spo_main;
  int nTiles=1;

  {
    Tensor<OHMMS_PRECISION,3> lattice_b;
    ParticleSet ions;
    OHMMS_PRECISION scale=1.0;
    lattice_b=tile_cell(ions,tmat,scale);
    const int nions=ions.getTotalNum();
    const int nels=count_electrons(ions)/2;
    tileSize=(tileSize>0)?tileSize:nels;
    nTiles=nels/tileSize;
    if(ionode)
      cout << "\nNumber of orbitals/splines = " << nels << " and Tile size = " << tileSize << " and Number of tiles = " << nTiles << " and Iterations = " << nsteps << endl;
    spo_main.set(nx,ny,nz,nels,nTiles);
    spo_main.Lattice.set(lattice_b);
  }

  double tInit = 0.0;
  double vgh_t=0.0, val_t=0.0;
  double nspheremoves=0;
  double dNumVGHCalls = 0;

//  const int time_array_size=17;
//  int timeIterArray[time_array_size] = {1, 10, 20, 30, 40, 50, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000}; // c++11 array.

  const int time_array_size=6;
  int timeIterArray[time_array_size] = {1, 10, 20, 30, 40, 50 }; // c++11 array.


  double vgh_t_iter[time_array_size] = {0}; // c++11 array.
  double val_t_iter[time_array_size] = {0}; // c++11 array.

  Timer bigClock;
  bigClock.restart();
#pragma omp parallel reduction(+:t0,tInit,ratio,vgh_t,val_t,nspheremoves,dNumVGHCalls)
  {
    Timer initClock;
    initClock.restart();
    const int np=omp_get_num_threads();
    const int ip=omp_get_thread_num();
    const int teamID=ip/ncrews;
    const int crewID=ip%ncrews;

    //create generator within the thread
    RandomGenerator<RealType> random_th(MakeSeed(teamID,np));

    ParticleSet ions, els;
    const OHMMS_PRECISION scale=1.0;
    ions.Lattice.BoxBConds=1;  
    tile_cell(ions,tmat,scale);

    const int nions=ions.getTotalNum();
    const int nels=count_electrons(ions);
    const int nels3=3*nels;

#pragma omp master
    nptcl=nels;

    {//create up/down electrons
      els.Lattice.BoxBConds=1;   els.Lattice.set(ions.Lattice);
      vector<int> ud(2); ud[0]=nels/2; ud[1]=nels-ud[0];
      els.create(ud);
      els.R.InUnit=1;
      random_th.generate_uniform(&els.R[0][0],nels3);
      els.convert2Cart(els.R); // convert to Cartiesian
    }

    //update content: compute distance tables and structure factor
    els.update();
    ions.update(); 

    //create pseudopp
    NonLocalPP<OHMMS_PRECISION> ecp(random_th);
    //create spo per thread
    spo_type spo(spo_main,ncrews,crewID); 

    //use teams
    //if(ncrews>1 && ncrews>=nTiles ) spo.set_range(ncrews,ip%ncrews);

    //this is the cutoff from the non-local PP
    const RealType Rmax(1.7);
    const int nknots(ecp.size());
    const RealType tau=2.0;

    ParticlePos_t delta(nels);
    ParticlePos_t rOnSphere(nknots);

#pragma omp master
    nknots_copy = nknots;

    RealType sqrttau=2.0;
    RealType accept=0.5;

    vector<RealType> ur(nels);
    random_th.generate_uniform(ur.data(),nels);
    const double zval=1.0*static_cast<double>(nels)/static_cast<double>(nions);

//test random numbers
//    Stat<RealType> measure;
//    auto res=measure.apply(ur.data(),nels);
//#pragma omp critical
//    printf("%d average= %.4f +- %.4f\n",omp_get_thread_num(),res.first,std::sqrt(res.second/nels));

    Timer clock, clock_a;
    double tt=0, vgh_t_loc=0, v_t_loc=0;      // local to the thread.
    int my_accepted=0, my_vals=0;

    double vgh_t_loc_iter[time_array_size] = {0};
    double v_t_loc_iter[time_array_size] = {0};

    tInit += initClock.elapsed();
    clock_a.restart();
    for(int mc=0; mc<nsteps; ++mc)
    {
#pragma omp barrier
      random_th.generate_normal(&delta[0][0],nels3);
      random_th.generate_uniform(ur.data(),nels);


      //VMC
      for(int iel=0; iel<nels; ++iel)
      {
        PosType pos=els.R[iel]+sqrttau*delta[iel];;
        clock.restart();
        spo.evaluate_vgh(els.R[iel]);
        vgh_t_loc+=clock.elapsed();
        if(ur[iel]>accept) 
        {
          els.R[iel]=pos;
          my_accepted++;
        }
       }
#pragma omp barrier

      random_th.generate_uniform(ur.data(),nels);
      ecp.randomize(rOnSphere); // pick random sphere
      for(int iat=0,kat=0; iat<nions; ++iat)
      {
        const int nnF=static_cast<int>(ur[kat++]*zval);
        RealType r=Rmax*ur[kat++];
        auto centerP=ions.R[iat];
        my_vals += (nnF*nknots);

        for(int nn=0; nn<nnF; ++nn)
        {
          for (int k=0; k < nknots ; k++)
          {
            PosType pos=centerP+r*rOnSphere[k];
            clock.restart();
            spo.evaluate_v(pos);
            v_t_loc+=clock.elapsed();
          }
        } // els
      } //ions

      for ( int p = 0; p < time_array_size; p++ ) {
           if ( mc == timeIterArray[p]-1 ){ 
             vgh_t_loc_iter[p] = vgh_t_loc;
             v_t_loc_iter[p] = v_t_loc;
           }
      }

    } // steps.
    tt+=clock_a.elapsed();

    ratio += RealType(my_accepted)/RealType(nels*nsteps);
    t0   +=tt/nsteps;
    vgh_t+=vgh_t_loc/nsteps; // time accumulated for all the threads for each time step.
    for ( int p = 0; p < time_array_size; p++ ) {
        vgh_t_loc_iter[p] /= timeIterArray[p];
        v_t_loc_iter[p] /= timeIterArray[p];
    }
    val_t+=v_t_loc/nsteps;
    nspheremoves+=RealType(my_vals)/RealType(nsteps);
    dNumVGHCalls += nels;

#pragma omp critical 
    {
      for ( int p = 0; p < time_array_size; p++ ) {
        vgh_t_iter[p] += vgh_t_loc_iter[p]; 
        val_t_iter[p] += v_t_loc_iter[p]; 
      }
    }
  } //end of omp parallel
  double tBigClock = bigClock.elapsed();

  double dTotalThreads = omp_get_max_threads();
  // Find out time, per walker.
  ///////////////////////
  double fac=1.0/dTotalThreads;
  t0 *= fac;
  vgh_t *= fac;
  val_t *= fac;
  tInit *= fac;
  ///////////////////////

  //collect timing and normalized by the number of ranks
  typedef TinyVector<double,4> timer_type;
  timer_type global_t(t0,vgh_t,val_t,0.0);
  timer_type global_t_1(tInit,tBigClock,0.0,0.0);

  mpi::reduce(*mycomm,global_t);
  mpi::reduce(*mycomm,global_t_1);

  const int nmpi=mycomm->size();
  t0=global_t[0]/nmpi;
  vgh_t=global_t[1]/nmpi;
  val_t=global_t[2]/nmpi;

  tInit=global_t_1[0]/nmpi;
  tBigClock=global_t_1[1]/nmpi;

  dNumVGHCalls = dNumVGHCalls/dTotalThreads;
  nspheremoves = nspheremoves/dTotalThreads;
  // Time per call:
  ////////////
  vgh_t = vgh_t/nptcl;
  val_t = val_t/nspheremoves;  

  for ( int p = 0; p < time_array_size; p++ ) {
    vgh_t_iter[p] *= fac;
    vgh_t_iter[p] /= nptcl;

    val_t_iter[p] *= fac;
    val_t_iter[p] /= nspheremoves;
  }

  /////////////

  // Number of operations per call.
  ///////////////
  double nMajorThreads=omp_get_max_threads()/ncrews;
  double numSplines = 0.5*nptcl;   // internally, these many splines are worked upon.
  double throughput_vgh = (nmpi*nMajorThreads*numSplines*1.e-06)/(vgh_t);    // Per vgh_t rate. 
  double throughput_v =  (nmpi*nMajorThreads*numSplines*1e-06)/(val_t);    // Per vgh_t rate. 
  ///////////////

  if(ionode)
  {
    cout << "Grid size: " << nx << "x" << ny << "x" << nz << " and Graphite N = " << nptcl << std::fixed << " and Ratio = " << ratio*fac << " and SphereMoves/ion = " << double(nspheremoves*4)/double(nknots_copy*nptcl); 

//    cout << "\nIterations: " << nsteps;
    cout << "\nMPI: " << nmpi << " and Threads: " << omp_get_max_threads(); 
    cout << "\nLevel 1 threads: " << int(nMajorThreads) << " and Nested threads = " << ncrews; 
    cout << "\nNumber of calls to Value = " << setprecision(1) << std::fixed << nspheremoves << " and VGH = " << dNumVGHCalls << endl;
    cout << "\nTime per call: Value and VGH: "<< setprecision(4) << std::scientific << val_t << "    " << vgh_t ; 
    cout << "\nThroughput: Value and VGH: " << setprecision(2) << std::fixed << throughput_v << "    " << throughput_vgh << endl;

    cout << "\nStats: MPI Threads_L1 Threads_L2 throughput_v throughput_vgh: " << nmpi << "    " << int(nMajorThreads) << "    " << ncrews << "    " << throughput_v << "    " << throughput_vgh << endl;
    cout << "\nTime Stats: VAL VGH  " << val_t*nsteps*nspheremoves << "   " << vgh_t*nsteps*nptcl << endl;

    cout << "\nTotal initialization time = " << setprecision(5) << tInit;
    cout << "\nTotal time of interest = " << setprecision(2) << t0*nsteps;
    cout << "\nTotal time of full app = " << setprecision(2) << tBigClock;
    cout << endl; 
    
    double numerator_th = (nMajorThreads*numSplines*1.e-06);
    for ( int p = 0; p < time_array_size; p++ ) {
       cout << timeIterArray[p] << " " << numerator_th/val_t_iter[p]  << " " << numerator_th/vgh_t_iter[p] << endl;
    }
  }


  OHMMS::Controller->finalize();

  return 0;
}

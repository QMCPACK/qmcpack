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
/** @file einspline_spo_nested.cpp
 * @brief Derived einspline_spo with nested parallelism
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
  if (OHMMS::Controller->rank() != 0) {
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
  int na=4;
  int nb=4;
  int nc=1;
  int nsteps=100;
  int iseed=11;
  int nx=48,ny=48,nz=60;
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
      case 'a':
        tileSize=atoi(optarg);
        break;
    }
  }

  //Random.init(0,1,iseed);
  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);

  //turn off output
  if(omp_get_max_threads()>1)
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

  Timer bigClock;
  bigClock.restart();
#pragma omp parallel reduction(+:t0,tInit,ratio,vgh_t,val_t,nspheremoves,dNumVGHCalls)
  {
    Timer initClock;
    initClock.restart();
    const int np=omp_get_num_threads();
    const int ip=omp_get_thread_num();

    //create generator within the thread
    RandomGenerator<RealType> random_th(MakeSeed(ip,np));

    ParticleSet ions, els;
    const OHMMS_PRECISION scale=1.0;
    tile_cell(ions,tmat,scale);

    const int nions=ions.getTotalNum();
    const int nels=count_electrons(ions);
    const int nels3=3*nels;

#pragma omp master
    {
      nptcl=nels;
      ncrews=omp_get_max_threads();
    }

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

    spo_type spo(spo_main,1,0);

    Timer clock_a;
    double tt=0, vgh_t_loc=0, v_t_loc=0;
    int my_accepted=0, my_vals=0;

    tInit += initClock.elapsed();
    clock_a.restart();

    //starting a nested region which corresponds to advanceWalker
    //random number generators are initialized so that all the threads generate
    //the same random moves and all the data are private except for
    //the position of the particle that is moving
#pragma omp parallel reduction(+:vgh_t_loc,v_t_loc,my_vals)
    {
      //this is the cutoff from the non-local PP
      const RealType Rmax(1.7);
      const RealType tau=2.0;

      RandomGenerator<RealType> my_random(random_th);
      NonLocalPP<OHMMS_PRECISION> ecp(random_th);

      const int nknots(ecp.size());
      ParticlePos_t delta(nels);
      ParticlePos_t rOnSphere(nknots);

      RealType sqrttau=2.0;
      RealType accept=0.5;

      vector<RealType> ur(nels);
      random_th.generate_uniform(ur.data(),nels);
      const double zval=2.0*static_cast<double>(nels)/static_cast<double>(nions);

      Timer clock;

      double vgh_t_loc2=0.0;
      double v_t_loc2=0.0;
      int my_vals2=0;

      for(int mc=0; mc<nsteps; ++mc)
      {
        my_random.generate_normal(&delta[0][0],nels3);
        my_random.generate_uniform(ur.data(),nels);

        //drift-diffusion stage
        for(int iel=0; iel<nels; ++iel)
        {
          PosType pos=els.R[iel]+sqrttau*delta[iel];;

          clock.restart();
          spo.evaluate_vgh_pfor(els.R[iel]); //internally using omp for over the blocks
          vgh_t_loc2+=clock.elapsed();

#pragma omp master
          if(ur[iel]>accept) 
          {
            els.R[iel]=pos;
            my_accepted++;
          }
        }

#pragma omp barrier

        my_random.generate_uniform(ur.data(),nels);
        ecp.randomize(rOnSphere); // pick random sphere
        for(int iat=0,kat=0; iat<nions; ++iat)
        {
          const int nnF=static_cast<int>(ur[kat++]*zval);
          RealType r=Rmax*ur[kat++];
          auto centerP=ions.R[iat];
          my_vals2 += (nnF*nknots);

          for(int nn=0; nn<nnF; ++nn)
          {
            for (int k=0; k < nknots ; k++)
            {
              PosType pos=centerP+r*rOnSphere[k];
              clock.restart();
              spo.evaluate_v_pfor(pos);
              v_t_loc2+=clock.elapsed();
            }
          } // els
        } //ions
      } // steps.

      vgh_t_loc+=vgh_t_loc2;
      v_t_loc+=v_t_loc2;
      my_vals+= my_vals2;

#pragma omp master
      nknots_copy = nknots;

    }//parallel region
    tt+=clock_a.elapsed();

    ratio += RealType(my_accepted)/RealType(nels*nsteps);
    t0   +=tt/nsteps;

    vgh_t+=vgh_t_loc/(nsteps*ncrews); // time accumulated for all the threads for each time step.
    val_t+=v_t_loc/(nsteps*ncrews);
    nspheremoves+=RealType(my_vals)/RealType(nsteps*ncrews);
    dNumVGHCalls += nels;

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
  /////////////

  // Number of operations per call.
  ///////////////
  double nMajorThreads=omp_get_max_threads();
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

    cout << "\nTotal initialization time = " << setprecision(5) << tInit;
    cout << "\nTotal time of interest = " << setprecision(2) << t0*nsteps;
    cout << "\nTotal time of full app = " << setprecision(2) << tBigClock;
    cout << endl; 
  }


  OHMMS::Controller->finalize();

  return 0;
}

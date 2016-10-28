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
/** @file moveonsphere.cpp
 * @brief Miniapp to capture the computation in NonLocalPP.
 *
 * Using only einspine SPO part + Jastrow as a wavefunction.
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/DistanceTable.h>
#include <Utilities/PrimeNumberSet.h>
#include <Utilities/Timer.h>
#include <random/random.hpp>
#include <mpi/collectives.h>
#include <miniapps/graphite.hpp>
#include <miniapps/pseudo.hpp>
#include <miniapps/common.hpp>
#include <miniapps/einspline_spo.hpp>
#include <miniapps/FakeWaveFunction.h>
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo welcome(argc,argv,OHMMS::Controller->rank());
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
  RealType Rmax(1.7);
  bool useSoA=false;

  PrimeNumberSet<uint32_t> myPrimes;

  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hos:g:i:b:c:a:r:")) != -1)
  {
    switch(opt)
    {
      case 'h':
        printf("[-g \"n0 n1 n2\"]\n");
        return 1;
      case 'o'://optimized method using SoA distance table and J2
        useSoA=true;
        break;
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
      case 'r'://rmax
        Rmax=atof(optarg);
        break;
      case 'a':
        tileSize=atoi(optarg);
        break;
    }
  }

  Random.init(0,1,iseed);
  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);

  //turn off output
  if(omp_get_max_threads()>1)
  {
    OhmmsInfo::Log->turnoff();
    OhmmsInfo::Warn->turnoff();
  }

  int nptcl=0;
  int nknots_copy=0;
  double t0=0.0,t1=0.0;
  OHMMS_PRECISION ratio=0.0;

  using spo_type=einspline_spo<OHMMS_PRECISION>;
  spo_type spo_main;
  int nTiles=1;

  {
    Tensor<OHMMS_PRECISION,3> lattice_b;
    ParticleSet ions;
    OHMMS_PRECISION scale=1.0;
    lattice_b=tile_graphite(ions,tmat,scale);
    const int nions=ions.getTotalNum();
    const int nels=2*nions;
    tileSize=(tileSize>0)?tileSize:nels;
    nTiles=nels/tileSize;
    if(ionode) {
      cout << "\nNumber of orbitals/splines = " << nels << " and Tile size = " << tileSize << " and Number of tiles = " << nTiles << " and Iterations = " << nsteps << endl;
    cout << "Rmax " << Rmax << endl;
    }
    spo_main.set(nx,ny,nz,nels,nTiles);
    spo_main.Lattice.set(lattice_b);
  }

  if(ionode) 
  {
    if(useSoA)
      cout << "Using SoA distance table and Jastrow + einspilne " << endl;
    else
      cout << "Using AoS distance table and Jastrow + einspilne of the original implementation " << endl;
  }

  double tInit = 0.0;
  double vgh_t=0.0, val_t=0.0, d_t0=0, d_t1=0, j_t0=0, j_t1=0, j_t2, j_t3=0, j_t4=0, j_t5=0;
  double nspheremoves=0;
  double dNumVGHCalls = 0;

  const int time_array_size=7;
  int timeIterArray[time_array_size] = {1, 5, 10, 20, 30, 40, 50 }; // c++11 array.
  double vgh_t_iter[time_array_size] = {0}; // c++11 array.
  double val_t_iter[time_array_size] = {0}; // c++11 array.

  Timer bigClock;
  bigClock.restart();
#pragma omp parallel reduction(+:t0,tInit,ratio,vgh_t,val_t,nspheremoves,dNumVGHCalls,d_t0,d_t1,j_t0,j_t1, j_t2, j_t3, j_t4, j_t5) 
  {
    Timer initClock;
    initClock.restart();

    ParticleSet ions, els;
    const OHMMS_PRECISION scale=1.0;

    const int np=omp_get_num_threads();
    const int ip=omp_get_thread_num();

    const int teamID=ip/ncrews;
    const int crewID=ip%ncrews;

    //create spo per thread
    spo_type spo(spo_main,ncrews,crewID); 

    //create generator within the thread
    RandomGenerator<RealType> random_th(myPrimes[ip]);

    ions.Lattice.BoxBConds=1;
    tile_graphite(ions,tmat,scale);

    const int nions=ions.getTotalNum();
    const int nels=4*nions;
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

    FakeWaveFunctionBase* Jastrow;

    if(useSoA) 
      Jastrow=new SoAWaveFunction(ions,els);
    else
      Jastrow=new AoSWaveFunction(ions,els);

    //create pseudopp
    NonLocalPP<OHMMS_PRECISION> ecp(random_th);

    //this is the cutoff from the non-local PP
    const int nknots(ecp.size());
    const RealType tau=2.0;

    ParticlePos_t delta(nels);
    ParticlePos_t rOnSphere(nknots);

#pragma omp master
    nknots_copy = nknots;

    RealType sqrttau=2.0;
    RealType accept=0.5;

    aligned_vector<RealType> ur(nels);
    random_th.generate_uniform(ur.data(),nels);

    double vgh_t_loc_iter[time_array_size] = {0};
    double v_t_loc_iter[time_array_size] = {0};

    constexpr RealType czero(0);

    els.update();
    Jastrow->evaluateLog(els);

    int my_accepted=0;
    for(int mc=0; mc<nsteps; ++mc)
    {
      random_th.generate_normal(&delta[0][0],nels3);
      random_th.generate_uniform(ur.data(),nels);

      Jastrow->evaluateLog(els);

      //drift-and-diffusion
      for(int l=0; l<5; ++l) 
      {
        for(int iel=0; iel<nels; ++iel)
        {
          //compute G[iel] with the current position to make a move
          els.setActive(iel);
          PosType grad_now=Jastrow->evalGrad(els,iel);

          //move iel el and compute the ratio
          PosType dr=sqrttau*delta[iel];
          bool isValid=els.makeMoveAndCheck(iel,dr); 

          if(!isValid) continue;

          PosType grad_new;
          RealType j2_ratio=Jastrow->ratioGrad(els,iel,grad_new);

          spo.evaluate_vgh(els.R[iel]);

          if(ur[iel]>accept) //MC
          { 
            els.acceptMove(iel);           
            Jastrow->acceptMove(els,iel);
            my_accepted++; 
          }
          else 
          { 
            els.rejectMove(iel); 
            Jastrow->restore(iel);
          }
        } // iel
      } //sub branch
      els.donePbyP();
      Jastrow->evaluateGL(els);
    }

    //cleanup
    delete Jastrow;
  } //end of omp parallel

  OHMMS::Controller->finalize();

  return 0;
}

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
/** @file j2debug.cpp
 * @brief Debugging J2OribitalSoA 
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/DistanceTable.h>
#include <OhmmsSoA/VectorSoaContainer.h>
#include <Utilities/PrimeNumberSet.h>
#include <random/random.hpp>
#include <Sandbox/input.hpp>
#include <Sandbox/pseudo.hpp>
#include <Utilities/Timer.h>
#include <Sandbox/common.hpp>
#include <Sandbox/FakeWaveFunction.h>
#include <getopt.h>

using namespace std;
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(0, NULL);
  Communicate* mycomm=OHMMS::Controller;

  typedef QMCTraits::RealType           RealType;
  typedef ParticleSet::ParticlePos_t    ParticlePos_t;
  typedef ParticleSet::ParticleLayout_t LatticeType;
  typedef ParticleSet::TensorType       TensorType;
  typedef ParticleSet::PosType          PosType;

  //use the global generator

  int na=4;
  int nb=4;
  int nc=1;
  int nsteps=100;
  int iseed=11;
  RealType Rmax(2.7);

  char *g_opt_arg;
  int opt;
  while((opt = getopt(argc, argv, "hg:i:r:")) != -1)
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
      case 'r'://rmax
        Rmax=atof(optarg);
        break;
    }
  }

//  Random.init(0,1,iseed);
  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);

  //turn off output
  if(omp_get_max_threads()>1)
  {
    outputManager.pause();
  }

  int nptcl=0;
  double t0=0.0,t1=0.0;
  OHMMS_PRECISION ratio=0.0;

  PrimeNumberSet<uint32_t> myPrimes;

  {
    ParticleSet ions, els;
    ions.setName("ion");
    els.setName("e");
    OHMMS_PRECISION scale=1.0;

    int np=omp_get_num_threads();
    int ip=omp_get_thread_num();

    //create generator within the thread
    RandomGenerator<RealType> random_th(myPrimes[ip]);

    tile_cell(ions,tmat,scale);
    ions.RSoA=ions.R; //fill the SoA

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
      els.RSoA=els.R;
    }

    ParticleSet els_aos(els);

    //create tables
    DistanceTableData* d_ee=DistanceTable::add(els,DT_SOA);

    DistanceTableData* d_ee_aos=DistanceTable::add(els_aos,DT_AOS);

    ParticlePos_t delta(nels);

    RealType sqrttau=2.0;
    RealType accept=0.5;

    vector<RealType> ur(nels);
    random_th.generate_uniform(ur.data(),nels);

    SoAWaveFunction J(ions,els);
    AoSWaveFunction J_aos(ions,els_aos);

    DistanceTableData* d_ie=DistanceTable::add(ions,els,DT_SOA);

    constexpr RealType czero(0);

    constexpr RealType small=std::numeric_limits<RealType>::epsilon();

    //compute distance tables
    els.update();
    els_aos.update();

    //for(int mc=0; mc<nsteps; ++mc)
    {
      J.evaluateLog(els);
      J_aos.evaluateLog(els_aos);

      cout << "Check values " << J.LogValue << " " << els.G[0] << " " << els.L[0] << endl;
      cout << "evaluateLog::V Error = " << (J.LogValue-J_aos.LogValue)/nels << endl;
      {
        double g_err=0.0;
        for(int iel=0; iel<nels; ++iel)
        {
          PosType dr= (els.G[iel]-els_aos.G[iel]);
          RealType d=sqrt(dot(dr,dr));
          g_err += d;
        }
        cout << "evaluateLog::G Error = " << g_err/nels << endl;
      }
      {
        double l_err=0.0;
        for(int iel=0; iel<nels; ++iel)
        {
          l_err += abs(els.L[iel]-els_aos.L[iel]);
        }
        cout << "evaluateLog::L Error = " << l_err/nels << endl;
      }

      random_th.generate_normal(&delta[0][0],nels3);
      double g_eval=0.0;
      double r_ratio=0.0;
      double g_ratio=0.0;

      int naccepted=0;
      for(int iel=0; iel<nels; ++iel)
      {
        els.setActive(iel);
        PosType grad_soa=J.evalGrad(els,iel);

        els_aos.setActive(iel);
        PosType grad_aos=J_aos.evalGrad(els_aos,iel)-grad_soa;

        g_eval+=sqrt(dot(grad_aos,grad_aos));

        PosType dr=sqrttau*delta[iel];

        grad_soa=0;
        els.makeMoveAndCheck(iel,dr); 
        RealType r_soa=J.ratioGrad(els,iel,grad_soa);

        grad_aos=0;
        els_aos.makeMoveAndCheck(iel,dr); 
        RealType r_aos=J_aos.ratioGrad(els_aos,iel,grad_aos);

        grad_aos-=grad_soa;
        g_ratio+=sqrt(dot(grad_aos,grad_aos));
        r_ratio += abs(r_soa/r_aos-1);

        if(Random() < r_aos)
        {
          J.acceptMove(els,iel);
          els.acceptMove(iel);

          els_aos.acceptMove(iel);
          J_aos.acceptMove(els_aos,iel);
          naccepted++;
        }
        else
        {
          els.rejectMove(iel);
          els_aos.rejectMove(iel);
        }
      }
      cout << "Accepted " << naccepted << "/" << nels << endl;
      cout << "evalGrad::G      Error = " << g_eval/nels << endl;
      cout << "ratioGrad::G     Error = " << g_ratio/nels << endl;
      cout << "ratioGrad::Ratio Error = " << r_ratio/nels << endl;

      els.donePbyP();
      els_aos.donePbyP();

      J.evaluateGL(els);
      J_aos.evaluateGL(els_aos);

      {
        double g_err=0.0;
        for(int iel=0; iel<nels; ++iel)
        {
          PosType dr= (els.G[iel]-els_aos.G[iel]);
          RealType d=sqrt(dot(dr,dr));
          g_err += d;
        }
        cout << "evaluteGL::G Error = " << g_err/nels << endl;
      }
      {
        double l_err=0.0;
        for(int iel=0; iel<nels; ++iel)
        {
          l_err += abs(els.L[iel]-els_aos.L[iel]);
        }
        cout << "evaluteGL::L Error = " << l_err/nels << endl;
      }

      //now ratio only
      r_ratio=0.0;
      constexpr int nknots=12;
      int nsphere=0;
      for(int jel=0; jel<nels; ++jel)
      {
        const auto &dist = d_ie->Distances[jel];
        for(int iat=0; iat<nions; ++iat)
          if(dist[iat]<Rmax)
          {
            nsphere++;
            random_th.generate_uniform(&delta[0][0],nknots*3);
            for(int k=0; k<nknots;++k)
            {
              els.makeMoveOnSphere(jel,delta[k]);
              RealType r_soa=J.ratio(els,jel);
              els.rejectMove(jel);

              els_aos.makeMoveOnSphere(jel,delta[k]);
              RealType r_aos=J_aos.ratio(els_aos,jel);
              els_aos.rejectMove(jel);
              r_ratio += abs(r_soa/r_aos-1);
            }
          }
      }
      cout << "ratio with SphereMove  Error = " << r_ratio/nsphere << " # of moves =" << nsphere << endl;
    }
  } //end of omp parallel

  OHMMS::Controller->finalize();

  return 0;
}

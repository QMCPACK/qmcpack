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
/** @file exchange_walker.cpp
 * @brief Test boost::mpi 
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/MCWalkerConfiguration.h>
#include <Utilities/PrimeNumberSet.h>
#include <random/random.hpp>
#include <miniapps/graphite.hpp>
#include <Utilities/Timer.h>
#include <miniapps/common.hpp>
#include <boost/mpi.hpp>
#include <getopt.h>

using namespace std;

namespace qmcplusplus
{

  namespace mpi=boost::mpi;

  /** dummy walker to show how boost::mpi and boost::serialization work
   */
  template<typename T, unsigned D>
    struct DummyWalker
    {
      using PosType=TinyVector<T,D>;
      long ID;
      ParticleAttrib<PosType> R;
      Matrix<double> Properties;

      DummyWalker() 
      {
        Properties.resize(2,8);
      }

      DummyWalker(const DummyWalker&)=default;

      inline void resize(size_t n)
      {
        R.resize(n);
      }

      template<class Archive>
        inline void serialize(Archive & ar, const unsigned int version)
        {
          ar & ID;
          ar & boost::serialization::make_array(&R[0][0],D*R.size());
          ar & boost::serialization::make_array(Properties.data(),Properties.size());
        }
    };

}


int main(int argc, char** argv)
{
  using namespace qmcplusplus;

  mpi::environment env(mpi::threading::funneled);
  mpi::communicator world;
  //create two MPI groups
  mpi::communicator rows=world.split(world.rank()/2);


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

  size_t nptcl=0;
  double t0=0.0,t1=0.0;
  OHMMS_PRECISION ratio=0.0;

  PrimeNumberSet<uint32_t> myPrimes;

  {
    ParticleSet ions; 
    MCWalkerConfiguration els;
    OHMMS_PRECISION scale=1.0;

    int np=omp_get_num_threads();
    int ip=omp_get_thread_num();

    //create generator within the thread
    RandomGenerator<RealType> random_th(myPrimes[world.rank()]);

    tile_graphite(ions,tmat,scale);
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

    size_t nw=4;
    using walker_t=DummyWalker<OHMMS_PRECISION,3>;
    vector<walker_t> walkers(nw);
    for(size_t i=0; i<nw; ++i)
    {
      walkers[i].resize(nptcl);
      walkers[i].ID=world.rank();
      walkers[i].R=0;
      walkers[i].Properties=i*world.size()+world.rank();
    }

    walkers[0].R=els.R;

    char fname[128];
    sprintf(fname,"debug.p%d",world.rank());
    ofstream fout(fname);

    //bcast a string
    string message("dummy");
    if(world.rank()==0) 
    {
      message="exchange_walker";
      walkers[0].R=els.R;
    }
    broadcast(world, message, 0);

    fout << message << endl << endl;
    fout << world.rank() << " " << rows.rank() << " " << rows.size() << endl;

    //send the skeleton
    broadcast(world, mpi::skeleton(walkers[0]), 0);

    mpi::request reqs[2];
    if(rows.rank()==1)
    {
      mpi::content c = mpi::get_content(walkers[1]);
      reqs[0]=rows.irecv(0,911,c);
    }
    else if(rows.rank()==0)
    {
      mpi::content c = mpi::get_content(walkers[0]);
      reqs[0]=rows.isend(1,911,c);
    }
    mpi::wait_any(reqs,reqs+1);

    fout << "Properties " << endl;
    fout << walkers[0].Properties << endl;
    fout << endl;
    fout << walkers[1].Properties << endl;
    fout << endl;
    fout << "ID " << walkers[1].ID << endl;
    for(size_t e=0; e<4; ++e)
      fout << walkers[0].R[e] << walkers[1].R[e] << endl;
  } //end of omp parallel

  return 0;
}

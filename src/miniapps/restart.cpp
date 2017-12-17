//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file restart.cpp
 * @brief developing restart IO
 */

#include <Configuration.h>
#include <Message/CommOperators.h>
#include <Particle/MCWalkerConfiguration.h>
#include <Particle/HDFWalkerOutput.h>
#include <HDFVersion.h>
#include <Particle/HDFWalkerInput_0_4.h>
#include <OhmmsApp/RandomNumberControl.h>
#include <random/random.hpp>
#include <miniapps/input.hpp>
#include <miniapps/pseudo.hpp>
#include <Utilities/Timer.h>
#include <miniapps/common.hpp>
#include <getopt.h>
#include <mpi/collectives.h>

using namespace std;
using namespace qmcplusplus;

void setWalkerOffsets(MCWalkerConfiguration& W, Communicate* myComm)
{
  std::vector<int> nw(myComm->size(),0),nwoff(myComm->size()+1,0);
  nw[myComm->rank()]=W.getActiveWalkers();
  myComm->allreduce(nw);
  for(int ip=0; ip<myComm->size(); ip++)
    nwoff[ip+1]=nwoff[ip]+nw[ip];
  W.setGlobalNumWalkers(nwoff[myComm->size()]);
  W.setWalkerOffsets(nwoff);
}

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(0, NULL);
  Communicate* myComm=OHMMS::Controller;
  myComm->setName("restart");

  typedef QMCTraits::RealType              RealType;
  typedef ParticleSet::ParticlePos_t       ParticlePos_t;
  typedef ParticleSet::ParticleLayout_t    LatticeType;
  typedef ParticleSet::TensorType          TensorType;
  typedef ParticleSet::PosType             PosType;
  typedef RandomGenerator_t::uint_type     uint_type;
  typedef MCWalkerConfiguration::Walker_t  Walker_t;

  //use the global generator

  int na=4;
  int nb=4;
  int nc=1;
  int nsteps=100;
  int iseed=11;
  int AverageWalkersPerNode=0;
  int nwtot;
  std::vector<int> wPerNode;
  RealType Rmax(1.7);

  const int NumThreads=omp_get_max_threads();

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
      case 'w'://the number of walkers
        AverageWalkersPerNode=atoi(optarg);
        break;
      case 'r'://rmax
        Rmax=atof(optarg);
        break;
    }
  }

  // set the number of walkers equal to the threads.
  if(!AverageWalkersPerNode) AverageWalkersPerNode=NumThreads;
  //set nwtot, to be random
  nwtot=std::abs(AverageWalkersPerNode+myComm->rank()%5-2);
  FairDivideLow(nwtot,NumThreads,wPerNode);

  //Random.init(0,1,iseed);
  Tensor<int,3> tmat(na,0,0,0,nb,0,0,0,nc);

  //turn off output
  if(myComm->rank())
  {
    outputManager.shutOff();
  }

  int nptcl=0;
  double t0=0.0,t1=0.0;

  RandomNumberControl::make_seeds();
  std::vector<RandomGenerator_t> myRNG(NumThreads);
  std::vector<uint_type> mt(Random.state_size(),0);
  std::vector<MCWalkerConfiguration> elecs(NumThreads);

  ParticleSet ions;
  OHMMS_PRECISION scale=1.0;
  tile_cell(ions,tmat,scale);

  #pragma omp parallel reduction(+:t0)
  {
    int ip=omp_get_thread_num();

    MCWalkerConfiguration& els=elecs[ip];

    //create generator within the thread
    myRNG[ip]=*RandomNumberControl::Children[ip];
    RandomGenerator_t& random_th=myRNG[ip];

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

    if(!ip) elecs[0].createWalkers(nwtot);
    #pragma omp barrier

    for(MCWalkerConfiguration::iterator wi=elecs[0].begin()+wPerNode[ip]; wi!=elecs[0].begin()+wPerNode[ip+1]; wi++)
      els.saveWalker(**wi);

    // save random seeds and electron configurations.
    *RandomNumberControl::Children[ip]=myRNG[ip];
    //MCWalkerConfiguration els_save(els);

  } //end of omp parallel
  Random.save(mt);

  setWalkerOffsets(elecs[0], myComm);

  //storage variables for timers
  double h5write = 0.0, h5read = 0.0; //random seed R/W speeds
  double walkerWrite = 0.0, walkerRead =0.0; //walker R/W speeds
  Timer h5clock; //timer for the program

  // dump random seeds
  myComm->barrier();
  h5clock.restart(); //start timer
  RandomNumberControl::write("restart",myComm);
  myComm->barrier();
  h5write += h5clock.elapsed(); //store timer

  // flush random seeds to zero
  #pragma omp parallel
  {
    int ip=omp_get_thread_num();
    RandomGenerator_t& random_th=*RandomNumberControl::Children[ip];
    std::vector<uint_type> vt(random_th.state_size(),0);
    random_th.load(vt);
  }
  std::vector<uint_type> mt_temp(Random.state_size(),0);
  Random.load(mt_temp);

  // load random seeds
  myComm->barrier();
  h5clock.restart(); //start timer
  RandomNumberControl::read("restart",myComm);
  myComm->barrier();
  h5read += h5clock.elapsed(); //store timer

  // validate random seeds
  int mismatch_count=0;
  #pragma omp parallel reduction(+:mismatch_count)
  {
    int ip=omp_get_thread_num();
    RandomGenerator_t& random_th=myRNG[ip];
    std::vector<uint_type> vt_orig(random_th.state_size());
    std::vector<uint_type> vt_load(random_th.state_size());
    random_th.save(vt_orig);
    RandomNumberControl::Children[ip]->save(vt_load);
    for(int i=0; i<random_th.state_size(); i++)
      if(vt_orig[i]!=vt_load[i]) mismatch_count++;
  }
  Random.save(mt_temp);
  for(int i=0; i<Random.state_size(); i++)
    if(mt_temp[i]!=mt[i]) mismatch_count++;

  myComm->allreduce(mismatch_count);

  if(!myComm->rank())
  {
    if(mismatch_count!=0)
      std::cout << "Fail: random seeds mismatch between write and read!\n"
                << "  state_size= " << myRNG[0].state_size() << " mismatch_cout=" << mismatch_count << std::endl;
    else
      std::cout << "Pass: random seeds match exactly between write and read!\n";
  }

  // dump electron coordinates.
  HDFWalkerOutput wOut(elecs[0],"restart",myComm);
  myComm->barrier();
  h5clock.restart(); //start timer
  wOut.dump(elecs[0],1);
  myComm->barrier();
  walkerWrite += h5clock.elapsed(); //store timer
  if(!myComm->rank()) std::cout << "Walkers are dumped!\n";

  // save walkers before destroying them
  std::vector<Walker_t> saved_walkers;
  for(int wi=0; wi<elecs[0].getActiveWalkers(); wi++)
    saved_walkers.push_back(*elecs[0][wi]);
  elecs[0].destroyWalkers(elecs[0].begin(),elecs[0].end());

  // load walkers
  const char *restart_input = \
"<tmp> \
  <mcwalkerset fileroot=\"restart\" node=\"-1\" version=\"3 0\" collected=\"yes\"/> \
</tmp> \
";

  Libxml2Document doc;
  bool okay = doc.parseFromString(restart_input);
  xmlNodePtr root = doc.getRoot();
  xmlNodePtr restart_leaf = xmlFirstElementChild(root);

  HDFVersion in_version(0,4);
  HDFWalkerInput_0_4 wIn(elecs[0],myComm,in_version);
  myComm->barrier();
  h5clock.restart(); //start timer
  wIn.put(restart_leaf);
  myComm->barrier();
  walkerRead += h5clock.elapsed(); //store time spent

  if(saved_walkers.size()!=elecs[0].getActiveWalkers())
    std::cout << "Fail: Rank " << myComm->rank() << " had " << saved_walkers.size()
              << "  walkers but loaded only " << elecs[0].getActiveWalkers() << " walkers from file!" << std::endl;

  mismatch_count=0;
  for(int wi=0; wi<saved_walkers.size(); wi++)
  {
    saved_walkers[wi].R=saved_walkers[wi].R-elecs[0][wi]->R;
    if(Dot(saved_walkers[wi].R,saved_walkers[wi].R)>std::numeric_limits<RealType>::epsilon()) mismatch_count++;
  }
  myComm->allreduce(mismatch_count);

  if(!myComm->rank())
  {
    if(mismatch_count!=0)
      std::cout << "Fail: electron coordinates mismatch between write and read!\n"
                << "  mismatch_cout=" << mismatch_count << std::endl;
    else
      std::cout << "Pass: electron coordinates match exactly between write and read!\n";
  }

  //print out hdf5 R/W times
  TinyVector<double,4> timers(h5read, h5write, walkerRead, walkerWrite);
  mpi::reduce(*myComm, timers);
  h5read = timers[0]/myComm->size();
  h5write = timers[1]/myComm->size();
  walkerRead = timers[2]/myComm->size();
  walkerWrite = timers[3]/myComm->size();
  if(myComm->rank() == 0)
  {
    cout << "\nTotal time of writing random seeds to HDF5 file: " << setprecision(6) << h5write << "\n";
    cout << "\nTotal time of reading random seeds in HDF5 file: " << setprecision(6) << h5read << "\n";
    cout << "\nTotal time of writing walkers to HDF5 file: " << setprecision(6) << walkerWrite << "\n";
    cout << "\nTotal time of reading walkers in HDF5 file: " << setprecision(6) << walkerRead << "\n";
  }

  if(myComm->size()>1)
  {
    Communicate* subComm = new Communicate(*myComm, 2);
    subComm->setName("restart2");

    if(subComm->getGroupID() == 0)
    {
      elecs[0].destroyWalkers(elecs[0].begin(),elecs[0].end());
      HDFWalkerInput_0_4 subwIn(elecs[0],subComm,in_version);
      subwIn.put(restart_leaf);
      subComm->barrier();
      if(!subComm->rank()) std::cout << "Walkers are loaded again by the subgroup!\n";
      setWalkerOffsets(elecs[0], subComm);
      HDFWalkerOutput subwOut(elecs[0],"XXXX",subComm);
      subwOut.dump(elecs[0],1);
      if(!subComm->rank()) std::cout << "Walkers are dumped again by the subgroup!\n";
    }
    delete subComm;
  }

  OHMMS::Controller->finalize();

  return 0;
}

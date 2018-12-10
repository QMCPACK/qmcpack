#include<tuple>
#include<map>
#include<string>
#include<iomanip>

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include<Message/MPIObjectBase.h>
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"
#include <qmc_common.h>

#include "AFQMC/config.h"
#include "AFQMC/Drivers/BenchmarkDriver.h"

namespace qmcplusplus {

bool BenchmarkDriver::run()
{

  app_log()<<"***********************************************************\n"
           <<"************ Starting Benchmark/Tests/Timings *************\n"
           <<"***********************************************************\n";
 
  prop0->benchmark(benchmark_list,maxnW,delnW,nrepeat,wlkBucket);

  wlkBucket->benchmark(benchmark_list,maxnW,delnW,nrepeat);
  
  return true;
}

bool BenchmarkDriver::parse(xmlNodePtr cur)
{
  if(cur==NULL) return false;

  std::string str,str1;

  ncores_per_TG=1; 
  ParameterSet m_param;
  m_param.add(benchmark_list,"list","string");
  m_param.add(maxnW,"maxnw","int");
  m_param.add(delnW,"delnw","int");
  m_param.add(nrepeat,"repeat","int");
  m_param.add(ncores_per_TG,"ncores_per_TG","int");
  m_param.add(ncores_per_TG,"ncores","int");
  m_param.add(ncores_per_TG,"cores","int");
  m_param.add(dt,"timestep","double");
  m_param.put(cur);

  return true;
}


bool BenchmarkDriver::setup(HamPtr h0, WSetPtr w0, PropPtr p0, WfnPtr wf0)
{

  ham0=h0;
  wlkBucket=w0;
  prop0=p0;
  wfn0=wf0;
  restarted=false;

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Beginning Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  app_log()<<" Using " <<ncores_per_TG <<" cores per node in a TaskGroup. \n";
  // right now this TG is not used. It is needed for setup purposes and to  
  // get a unique TG number for every group of cores on a node (used in the WalkerSet)
  TG.setup(ncores_per_TG,1,false);  
  std::vector<int> TGdata(5);
  TG.getSetupInfo(TGdata); 
   
  // setup local-to-node MPI Comm
  // TGdata[0]: node_number
  myComm->split_comm(TGdata[0],MPI_COMM_NODE_LOCAL);
  TG.setNodeCommLocal(MPI_COMM_NODE_LOCAL);
  int key = TG.getTGNumber(); // This works because the TG used has nnodes_per_TG=1 
  myComm->split_comm(key,MPI_COMM_TG_LOCAL);
  TG.setTGCommLocal(MPI_COMM_TG_LOCAL);
  key = TG.getCoreRank();  
  myComm->split_comm(key,MPI_COMM_TG_LOCAL_HEADS);

  key = TG.getCoreID();
  myComm->split_comm(key,MPI_COMM_HEAD_OF_NODES);
  TG.setHeadOfNodesComm(MPI_COMM_HEAD_OF_NODES);

  CommBuffer.setup(TG.getCoreRank()==0,std::string("COMMBuffer_")+std::to_string(TG.getTGNumber()),MPI_COMM_TG_LOCAL);
  TG.setBuffer(&CommBuffer);

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Hamiltonian \n"
           <<"****************************************************\n"
           <<std::endl;

  // hamiltonian
  if(!ham0->init(TGdata,&CommBuffer,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,MPI_COMM_HEAD_OF_NODES)) {
    app_error()<<"Error initializing Hamiltonian in BenchmarkDriver::setup" <<std::endl; 
    return false; 
  }   

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Wavefunction \n"
           <<"****************************************************\n"
           <<std::endl;

  hdf_archive read(myComm);
  if(!wfn0->init(TGdata,&CommBuffer,read,std::string(""),MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,MPI_COMM_HEAD_OF_NODES)) {
    app_error()<<"Error initializing Wavefunction in BenchmarkDriver::setup" <<std::endl; 
    return false; 
  }   
  if(!wfn0->setup(ham0)) {
    app_error()<<"Error in WavefunctionHandler::setup in BenchmarkDriver::setup" <<std::endl; 
    return false; 
  }   

  app_log()<<"\n****************************************************\n"
           <<"              Initializating Propagator \n"
           <<"****************************************************\n"
           <<std::endl;

  // propagator
  if(!prop0->setup(TGdata,&CommBuffer,ham0,wfn0,dt,read,std::string(""),MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,MPI_COMM_HEAD_OF_NODES)) { 
    app_error()<<"Error in PropagatorBase::setup in BenchmarkDriver::setup" <<std::endl; 
    return false; 
  }   

  // you will also need the Propagator's TG to handle the distibuted case 
  wfn0->setupFactorizedHamiltonian(prop0->is_vn_sparse(),prop0->getSpvn(),prop0->getDvn(),dt,prop0->getTG());

  app_log()<<"\n****************************************************\n"
           <<"             Initializating Walker Handler \n"
           <<"****************************************************\n"
           <<std::endl;

  // walker set
  wlkBucket->setup(TG.getCoreRank(),ncores_per_TG,TG.getTGNumber(),MPI_COMM_TG_LOCAL_HEADS,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,&LocalTimer);
  wlkBucket->setHF(wfn0->getHF());
  wlkBucket->initWalkers(maxnW);
  wfn0->evaluateLocalEnergyAndOverlap("ImportanceSampling",-1,wlkBucket);

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Finished Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  myComm->barrier(); 

  return true;
}

}

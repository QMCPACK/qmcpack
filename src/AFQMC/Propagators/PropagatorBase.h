
#ifndef QMCPLUSPLUS_AFQMC_PROPAGATORBASE_H
#define QMCPLUSPLUS_AFQMC_PROPAGATORBASE_H

#include<fstream> 

#include "OhmmsData/libxmldefs.h"
#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h" 
//#include "AFQMC/Walkers/SlaterDetWalker.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "AFQMC/Estimators/EstimatorHandler.h"
#include "Utilities/RandomGenerator.h"
#include "AFQMC/Estimators/SlaterDetOperations.h"

namespace qmcplusplus
{


class PropagatorBase: public MPIObjectBase, public AFQMCInfo
{
  public:
       
  PropagatorBase(Communicate *c, RandomGenerator_t* r): MPIObjectBase(c),TG(c,"PropagatorTG"),rng(r),Order_Taylor_Expansion(6),name(""),hdf_write_tag(""),hdf_write_file(""),hdf_read_tag(""),hdf_read_file(""),parallel_factorization(true),ncores_per_TG(1),nnodes_per_TG(1),parallelPropagation(true),distributeSpvn(false),core_rank(0),sparsePropagator(true) 
  {
  }

  ~PropagatorBase() {}

//  virtual void Propagate(int n, SlaterDetWalker&, RealType& E1, const RealType E2=0)=0;

  virtual void Propagate(int steps, int& steps_total, WalkerHandlerBase*, RealType& E1)=0;

  virtual bool parse(xmlNodePtr)=0;

  virtual bool setup(std::vector<int>&,SPComplexSMVector*,HamiltonianBase*,WavefunctionHandler*, RealType dt, hdf_archive&, const std::string&, MPI_Comm, MPI_Comm, MPI_Comm)=0;

  virtual bool hdf_write(hdf_archive&, const std::string&)=0;

  virtual bool hdf_read(hdf_archive&,const std::string&)=0;

  bool is_vn_sparse() { return sparsePropagator; }

  afqmc::TaskGroup* getTG() { return &TG; } 

  virtual SPValueSMVector* getDvn()=0;

  virtual SPValueSMSpMat* getSpvn()=0;

  virtual void benchmark(std::string&,int,int,int,WalkerHandlerBase*)=0;

  SlaterDetOperations* SDetOps;

  // timestep
  RealType dt; 

  afqmc::TaskGroup TG;
  int ncores_per_TG,nnodes_per_TG;
  int core_rank;

  bool parallelPropagation;
  bool distributeSpvn; 

  bool sparsePropagator;

  bool parallel_factorization;
  bool head_of_nodes;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;

  // used to sort snD values using only indexes 
  _mySort_snD_ mySort;

  RandomGenerator_t* rng;

  int Order_Taylor_Expansion;

  // name of object
  std::string name; 

  // id of HDF group of this object
  // The actual datasets will be stored on:
  //   /Propagators/ACTUAL_PROPAGATOR/hdf_tag 
  std::string hdf_write_tag;
  // hdf file where data will be stored. if a filename is found
  // on the xml section of the propagator, it will be used.
  // otherwise the one from the driver will be used.  
  std::string hdf_write_file; 

  std::string hdf_read_tag; 
  std::string hdf_read_file; 

};


}

#endif


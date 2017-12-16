#ifndef QMCPLUSPLUS_AFQMC_ESTIMATORBASE_H
#define QMCPLUSPLUS_AFQMC_ESTIMATORBASE_H

#include<Message/MPIObjectBase.h>
#include"AFQMC/config.h"
#include<vector>
#include<iostream>
#include<fstream>

#include "io/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Walkers/WalkerHandlerBase.h"

namespace qmcplusplus
{

class EstimatorBase: public MPIObjectBase, public AFQMCInfo
{

  public:

  typedef HamiltonianBase* HamPtr;
  typedef WavefunctionHandler* WfnPtr;
  typedef WalkerHandlerBase* WSetPtr;

  EstimatorBase(Communicate *c):MPIObjectBase(c) {}

  ~EstimatorBase() {}

  virtual void accumulate_block(WSetPtr wlks)=0;

  virtual void accumulate_step(WSetPtr wlks, std::vector<ComplexType>& curData)=0;

  //virtual void accumulate_substep(WSetPtr wlks){};

  virtual void print(std::ofstream& out,WalkerHandlerBase* wlks)=0; 

  virtual void tags(std::ofstream& out)=0; 

  virtual void average(WSetPtr wlks)=0;

  virtual bool parse(xmlNodePtr)=0;

  virtual bool setup(std::vector<int>& TGdata, SPComplexSMVector *v,HamiltonianBase*,WavefunctionHandler*,myTimer* LocalTimer, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm cm)=0;

  virtual double getEloc() {return 0;}

  virtual double getEloc_step() {return 0;}

  virtual void setTargetWeight(RealType w0) {} 

  protected:

  bool writer;

  std::string filename, filetype, ID;

  HamPtr ham0;

  WfnPtr wfn0;

};
}

#endif

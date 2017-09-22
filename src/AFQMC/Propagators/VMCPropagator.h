
#ifndef QMCPLUSPLUS_AFQMC_VMCPROPAGATOR_H
#define QMCPLUSPLUS_AFQMC_VMCPROPAGATOR_H

#include "OhmmsData/libxmldefs.h"
#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Propagators/PropagatorBase.h"
#include "AFQMC/Wavefunctions/WavefunctionHandler.h" 
#include "AFQMC/Walkers/WalkerHandlerBase.h"
#include "Utilities/RandomGenerator.h"
#include "AFQMC/Estimators/SlaterDetOperations.h"
#include "AFQMC/Hamiltonians/ProjectorBase.h"

namespace qmcplusplus
{


class VMCPropagator: public PropagatorBase
{
  public:
       
  VMCPropagator(Communicate *c, RandomGenerator_t* r):PropagatorBase(c,r),use_eig(false),cutoff(1e-6) 
  {
  }

  ~VMCPropagator() {}

  void Propagate(int steps, int& steps_total, WalkerHandlerBase*, RealType& E1);

  bool parse(xmlNodePtr);

  bool setup(std::vector<int>& TGdata, SPComplexSMVector *v,HamiltonianBase*,WavefunctionHandler*, RealType dt, hdf_archive&, const std::string&,MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm);

  bool hdf_write(hdf_archive&, const std::string&);

  bool hdf_read(hdf_archive&,const std::string&);

  void benchmark(std::string&,int,int,int,WalkerHandlerBase*) {}

  SPValueSMVector* getDvn() { return NULL; } 

  SPValueSMSpMat* getSpvn() { return NULL; }

  private:

  ProjectorBase* proj0;

  bool use_eig;
 
  double cutoff; 

  ComplexMatrix vHS;

  ComplexSMSpMat Spvn;

  SlaterDetOperations* Sdet;

  // local storage
  ComplexMatrix T1;
  ComplexMatrix T2;
  
  ComplexMatrix SL, SR;

  std::vector<ComplexType> sigmaL;
  std::vector<ComplexType> sigmaR;

  void applyHSPropagator(ComplexType*, std::vector<ComplexType>& , int order=-1);

  void sampleGaussianFields(std::vector<ComplexType>&);

};


}

#endif


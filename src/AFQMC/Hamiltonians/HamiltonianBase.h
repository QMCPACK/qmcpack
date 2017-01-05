
#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANBASE_H
#define QMCPLUSPLUS_AFQMC_HAMILTONIANBASE_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include<Message/MPIObjectBase.h>
#include "OhmmsData/libxmldefs.h"
#include "AFQMC/Utilities/taskgroup.h"

#include"AFQMC/config.h"

namespace qmcplusplus
{

class HamiltonianBase: public MPIObjectBase, public AFQMCInfo
{

  public:
 
  HamiltonianBase(Communicate *c):MPIObjectBase(c),TG(c,"HamiltonianTG"),name(""),filetype("undefined"),filename("undefined"),test_breakup(false),head_of_nodes(false),distribute_Ham(false),min_i(0),max_i(0),nnodes_per_TG(1)
  {
    FrozenCoreEnergy = NuclearCoulombEnergy = ValueType(0.0);
  }

  ~HamiltonianBase() {}

  inline int getNMO_FULL() {return NMO_FULL;}
  inline int getNAEA() { return NAEA;}
  inline int getNAEB() { return NAEB;}
  inline int getNCA() { return NCA;}
  inline int getNCB() { return NCB;}
  inline bool RHF() {return spinRestricted;}
  
  bool init(std::vector<int>& TGdata, ComplexSMVector* TGbuff, MPI_Comm tg_comm, MPI_Comm node_comm )
  {
    if(nnodes_per_TG > 1) distribute_Ham=true;
    if(!TG.quick_setup(TGdata[4],nnodes_per_TG,TGdata[0],TGdata[1],TGdata[2],TGdata[3]))
      return false;
    TG.setBuffer(TGbuff);
    TG.setNodeCommLocal(node_comm);
    TG.setTGCommLocal(tg_comm); 

    if(filetype == "fcidump" || filetype == "ascii")
      return initFromASCII(filename); 
    else if(filetype == "xml")
      return initFromXML(filename);
    else if(filetype == "hdf5")
      return initFromHDF5(filename);
    else {
      app_error()<<"Unknown filetype in HamiltonianBase::init(): " <<filetype <<std::endl; 
      return false;
    }
  } 

  virtual void calculateHSPotentials(const RealType cut, const RealType dt, ComplexSMSpMat&, TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool paral )=0; 

  virtual void calculateHSPotentials_Diagonalization(const RealType cut, const RealType dt, ComplexSMSpMat&, TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool paral)=0; 

  virtual void calculateOneBodyPropagator(const RealType cut, const RealType dt, ComplexMatrix& Hadd, std::vector<s2D<ComplexType> >& Pkin)=0;

  virtual bool generateFullHamiltonianForME()=0; 

  virtual bool getFullHam(std::vector<s1D<ValueType> >*& h, ComplexSMSpMat*& v)=0;
  
  // parse xml input node
  virtual bool parse(xmlNodePtr cur)=0; 

  // check object
  virtual bool checkObject()=0;

  // keep a copy of Reference State in FCIDUMP
  ComplexMatrix RefWFn; 

  // name of the object
  std::string name;

  // nuclear coulomb term
  ValueType NuclearCoulombEnergy; 
  ValueType FrozenCoreEnergy; 

  // timestep
  RealType dt; 

  void setHeadComm(bool hd, MPI_Comm comm) {
    head_of_nodes=hd;
    MPI_COMM_HEAD_OF_NODES = comm;
  }

//  void setTGComm(bool hd, MPI_Comm comm) {
//    head_of_local_tg=hd;
//    MPI_COMM_LOCAL_TG = comm;
//  }

  protected:

  // for hamiltonian distribution 
  
  TaskGroup TG; 
  int nnodes_per_TG;

  std::string filetype;
  std::string filename;
  bool test_breakup;
  bool head_of_nodes;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;
  bool distribute_Ham;  // implement assuming factorized Ham first  
  int min_i, max_i;
  

  virtual bool initFromASCII(const std::string& fileName)=0; 

  virtual bool initFromXML(const std::string& fileName)=0;

  virtual bool initFromHDF5(const std::string& fileName)=0; 

};
}

#endif


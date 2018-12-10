
#ifndef QMCPLUSPLUS_AFQMC_PROJECTORBASE_H
#define QMCPLUSPLUS_AFQMC_PROJECTORBASE_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include<Message/MPIObjectBase.h>
#include "OhmmsData/libxmldefs.h"

#include"AFQMC/config.h"
#include"AFQMC/Hamiltonians/HamiltonianBase.h"

namespace qmcplusplus
{

class ProjectorBase: public MPIObjectBase, public AFQMCInfo
{

  typedef HamiltonianBase* HamPtr;

  public:
 
  ProjectorBase(Communicate *c):MPIObjectBase(c),name(""),filetype("undefined"),filename("undefined"),test_breakup(false),head_of_nodes(false),cutoff_sparse(1e-5),eigcut(1e-5) 
  {
    
  }

  ~ProjectorBase() {}

  bool init(HamPtr h)
  {
    ham0 = h;
    if(filetype == "ascii")
      return initFromASCII(filename); 
    else if(filetype == "xml")
      return initFromXML(filename);
    else if(filetype == "hdf5")
      return initFromHDF5(filename);
    else
      return initFromGuess();
  } 

  virtual void calculateHSPotentials(ComplexSMSpMat&)=0; 

  virtual void calculateHSPotentials_Diagonalization(ComplexSMSpMat&)=0; 

  // parse xml input node
  virtual bool parse(xmlNodePtr cur)=0; 

  // check object
  virtual bool checkObject()=0;

  // name of the object
  std::string name;

  void setHeadComm(bool hd, MPI_Comm comm) {
    head_of_nodes=hd;
    MPI_COMM_HEAD_OF_NODES = comm;
  }

  virtual void hdf_write()=0;

  protected:
  
  HamPtr ham0;

  std::string hdf_write_file;
  RealType cutoff_sparse;
  RealType eigcut;
  std::string filetype;
  std::string filename;
  bool test_breakup;
  bool head_of_nodes;
  MPI_Comm MPI_COMM_HEAD_OF_NODES;

  virtual bool initFromASCII(const std::string& fileName)=0; 

  virtual bool initFromXML(const std::string& fileName)=0;

  virtual bool initFromHDF5(const std::string& fileName)=0; 

  virtual bool initFromGuess()=0; 

};
}

#endif


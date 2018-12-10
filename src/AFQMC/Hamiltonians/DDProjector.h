
#ifndef QMCPLUSPLUS_AFQMC_DDPROJECTOR_H
#define QMCPLUSPLUS_AFQMC_DDPROJECTOR_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include<Message/MPIObjectBase.h>
#include "OhmmsData/libxmldefs.h"

#include"AFQMC/config.h"

namespace qmcplusplus
{

class DDProjector: public ProjectorBase 
{

  typedef HamiltonianBase* HamPtr;

  public:
 
  DDProjector(Communicate *c):ProjectorBase(c)
  {
    
  }

  ~DDProjector() {}

  void calculateHSPotentials(ComplexSMSpMat&);

  void calculateHSPotentials_Diagonalization(ComplexSMSpMat&);

  // parse xml input node
  bool parse(xmlNodePtr cur); 

  // check object
  bool checkObject();

  void setHeadComm(bool hd, MPI_Comm comm) {
    head_of_nodes=hd;
    MPI_COMM_HEAD_OF_NODES = comm;
  }

  void hdf_write();

  protected:

  bool parallel_factorization=false;
  bool use_eig=true;  

  // sum_ij P(i,j) n_i n_j , where n_k = c+_k c_k  
  // P(i,j) == 0, P(i,j)==P(j,i)
  ValueMatrix Pmat;
  
  bool initFromASCII(const std::string& fileName) { return false;}; 

  bool initFromXML(const std::string& fileName) { return false; };

  bool initFromHDF5(const std::string& fileName); 
  
  bool initFromGuess();

};
}

#endif


#ifndef QMCPLUSPLUS_AFQMC_SELECTEDCI_H
#define QMCPLUSPLUS_AFQMC_SELECTEDCI_H

#include<Message/MPIObjectBase.h>
#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Drivers/Driver.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"

namespace qmcplusplus
{

class selectedCI: public Driver 
{

  typedef HamiltonianBase* HamPtr;
  typedef AFQMCInfo* InfoPtr;

  public:

    selectedCI(Communicate *c):Driver(c),build_full_hamiltonian(true),maxit(5),
                               diag_in_steps(-1),cutoff_list(0.1),cutoff_diag(0.1),
                               output_filename("selectedCI")
    {
      name = "selectedCI";
      project_title = "selectedCI";
    }

    ~selectedCI() {}

    bool run();

    bool parse(xmlNodePtr); 

    bool setup(HamPtr,WSetPtr,PropPtr,WfnPtr);

    bool checkpoint(int,int);

    bool restart(hdf_archive&); 

    bool clear();

    bool diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ValueMatrix& eigVec, std::vector<IndexType>& occ1, int nci1, std::vector<IndexType>& occ2, int nci2, bool eigV=true);

  protected:  

    SparseGeneralHamiltonian* sHam;

    int maxit;
    int diag_in_steps;
    double cutoff_list;
    double cutoff_diag;
    std::string output_filename;
    bool build_full_hamiltonian; 
    ValueType NuclearCoulombEnergy; 

    std::vector<IndexType> occ_orbs;
    std::vector<ValueType> ci;
    std::vector<IndexType> val_at_pivot;

    void sort_multiple(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right, std::vector<ValueType>::iterator v1, std::vector<ValueType>::iterator v2);
    void sort_list(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right);
    void remove_repeated(std::vector<IndexType>& vnew, std::vector<IndexType>& vold);
    inline bool list_equal(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right)
    { 
      for(int i=0; i<NAEA+NAEB; i++)
        if(!(*(left++) == *(right++))) return false;
      return true;
    }
    inline bool list_order(std::vector<IndexType>::iterator left, std::vector<IndexType>::iterator right)
    { 
      for(int i=0; i<NAEA+NAEB; i++,left++,right++)
        if(*(left) == *(right)) 
          continue; 
        else  
          return *left < *right;
      return false; // they are equal 
    }
    inline bool mysearch(std::vector<IndexType>::iterator& first, std::vector<IndexType>::iterator last, std::vector<IndexType>::iterator val)
    {

      std::vector<IndexType>::iterator left = first;
      std::vector<IndexType>::iterator middle, right = last;
      IndexType ne=NAEA+NAEB;
      while (left <= right) {
        middle = left;
        std::advance(middle,(std::distance(left,right)/ne)/2*ne);
        if (list_equal(middle,val)) {
          first=middle;
          return true;
        } else if(list_order(middle,val))
          left = middle + (NAEA+NAEB);
        else 
          right = middle - (NAEA+NAEB);
      }
      first = left;
      if(right < left) first = right; 
      return false;
    } 

};
}

#endif

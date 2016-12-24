#ifndef QMCPLUSPLUS_AFQMC_GENERALSINGLEDETERMINANT_H
#define QMCPLUSPLUS_AFQMC_GENERALSINGLEDETERMINANT_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "AFQMC/config.h"
#include <Message/MPIObjectBase.h>
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Hamiltonians/SparseGeneralHamiltonian.h"

namespace qmcplusplus
{

/*
 * Class that implements a pure single determinant trial wave-function.
 * General means that the Slater matrix is composed of 0,1 only, 
 * meaning that it represents an eigenstate of the mean-field 
 * solution used to construct the hamiltonian. 
 *
 */  
class GeneralSingleDeterminant: public WavefunctionBase
{

  typedef WavefunctionBase* WfnPtr;
  typedef GeneralSingleDeterminant ThisWfn;
  typedef GeneralSingleDeterminant* ThisWfnPtr;
  typedef HamiltonianBase* HamPtr;
  typedef std::vector<IndexType>::iterator        VIndexit; 
  typedef std::vector<s1D<ValueType> >::iterator  s1Dit; 
  typedef std::vector<s2D<ValueType> >::iterator  s2Dit; 
  typedef std::vector<s4D<ValueType> >::iterator  s4Dit; 

  public:

    GeneralSingleDeterminant(Communicate *c):WavefunctionBase(c),trialDensityMatrix_needsupdate(true),cutoff(1e-5),wfntype(0),uhf_walker(true)
    {}

    ~GeneralSingleDeterminant() {}

    bool setup(HamPtr cur); 

    bool parse(xmlNodePtr ); 

    bool init(hdf_archive& read, const std::string& tag)
    {
      if(filetype == "none" || filetype == "" || init_type=="diagH1" || init_type=="ground")
        return setup_local();
      else if(filetype == "fcidump" || filetype == "ascii" || filetype == "sqc_ascii")
        return initFromAscii(filename);
      else if(filetype == "xml")
        return initFromXML(filename);
      else if(filetype == "hdf5") {
        hdf_archive readF(myComm);
        if(head_of_nodes)
          if(!readF.open(filename,H5F_ACC_RDONLY,false))
            APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
        if(!initFromHDF5(readF,hdf_read_tag)) {
          app_error()<<" Problems reading hdf5 file in WavefunctionBase::init()";
          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
          return false;
        }
        readHamFromFile=true;
        if(head_of_nodes) readF.close();
        return true;
      } else {
        if(!initFromHDF5(read,tag)) {
          app_error()<<" Problems reading restart file in WavefunctionBase::init()";
          APP_ABORT(" Problems reading hdf5 file in WavefunctionBase::init()");
          return false;
        }
        readHamFromFile=true;
        return true;
      }
      app_error()<<" Could not find a wavefunction initialization type. \n";
      return false;
    }

    //bool hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors=true);
    bool hdf_write();

    void evaluateMeanFields()  {}

    void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G) {}

    void evaluateTwoBodyMixedDensityMatrix() {}

    void evaluateLocalEnergy(const ComplexType* , ComplexType& , ComplexType&, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void evaluateOverlap(const ComplexType* , ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1);
    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1);

    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1 );
    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1 );

  //private:
 
    bool setup_local();  
  
    bool initFromAscii(std::string fileName); 
    bool initFromHDF5(hdf_archive&,const std::string&); 
    bool initFromXML(std::string fileName) {} 

    bool getHamiltonian(HamPtr ); 

    // evaluates and stores mixed density matrix in mixed_density_matrix
    // this evaluates the mixed density matrix in reduced format (only NAEX*NMO non-zero sector)
    void local_evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta);

    void local_evaluateOneBodyTrialDensityMatrix();

    ValueType NuclearCoulombEnergy; 

    RealType cutoff;

    // vector that contains the list of orbitals occupied in the given Slater Matrix 
    // Indexes are 0-based.
//    std::vector<IndexType> occup_alpha; 
//    std::vector<IndexType> occup_beta; 
//    std::vector<IndexType> virtual_alpha; 
//    std::vector<IndexType> virtual_beta; 

    // alternative storage for occupied states for easy access/lookup
    // each MO (from 0...NMO-1) is mapped to either true/false based on occupation
//    std::map<IndexType,bool> isOcc_alpha; 
//    std::map<IndexType,bool> isOcc_beta; 

    // 1RDM of the trial density matrix.
    // Used to calculate mean-field energy and mean-field potentials 
    // This is a permutation of the identity matrix, do we need to store it??? 
    bool trialDensityMatrix_needsupdate; 
    ComplexMatrix trial_density_matrix;

    int wfntype;
    bool uhf_walker;

    // Local storage  
    // Notice that all these matrices are NMOxNAEX, since the rest of the NMO-NAEX columns are zero.
    // Careful must be taken when returning results to other objects in the code.
    ComplexMatrix overlap_inv;
    ComplexMatrix mixed_density_matrix;

    // temporary storage
    ComplexMatrix S0,S1,SS0,T0, SM; 
    ComplexVector V0; 

    std::vector<ComplexType> Cwork;
    std::vector<int> pivot;

    // Slater matrices for the trial wavefunction
    ComplexMatrix OrbMat;   

    // One-Body Hamiltonian. Stored in sparse form. 
    std::vector<s1D<ValueType> > hij;

    // Storage for two body hamiltonian 
    // Tensor is stored in sparse form. 
    ComplexSpMat SpHijkl;

    ComplexSMSpMat SMSpHijkl;

    // 
    SparseGeneralHamiltonian* sHam;

    /*
    // This is only for debugging and setup.
    // Should never be used inside execution loop
    ValueType H(IndexType I, IndexType J) {
      if( I < NMO && J < NMO ) {
        IndexType indx = I*NMO+J;  
        // use binary search later
        for(s1Dit it = hij.begin(); it<hij.end(); it++)
          if(std::get<0>(*it) == indx) return std::get<1>(*it);
      } else if( I >= NMO && J >= NMO ) {
        IndexType indx = I*NMO+J-NMO;  
        // use binary search later
        for(s1Dit it = hij.begin(); it<hij.end(); it++)
          if(std::get<0>(*it) == indx) return std::get<1>(*it);
      } 
      return static_cast<ValueType>(0.0);
    }

    // This is only for debugging and setup.
    // Should never be used inside execution loop
    ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) {

#ifdef AFQMC_DEBUG
// probably a good idea to check that (I,K) and (J,L) belong to same spin
#endif

      ValueType scl = static_cast<ValueType>(1.0);
      IndexType indx1 = (I<NMO)?(I*NMO+K):(I*NMO+K-NMO);
      IndexType indx2 = (J<NMO)?(J*NMO+L):(J*NMO+L-NMO);
      // only do this is you eliminate redundant pairs from list by *2.0
      //if( !(I==J && K==L) && ((I<NMO && J<NMO) || (I>=NMO && J>=NMO))  ) scl *= static_cast<ValueType>(0.5);  
      if( !(I==J && K==L) ) scl *= static_cast<ValueType>(0.5);  
      for(s2Dit it = Vijkl.begin(); it<Vijkl.end(); it++)
        if( (std::get<0>(*it) == indx1 && std::get<1>(*it) == indx2) ||
            (std::get<0>(*it) == indx2 && std::get<1>(*it) == indx1)   ) return std::get<2>(*it)*scl;
      return static_cast<ValueType>(0.0);
    }
    */

    inline IndexType Index2Mat(IndexType I, IndexType J) {
      return (J<NMO)?(I*NMO+J):(I*NMO+J-NMO);
    }

};
}

#endif


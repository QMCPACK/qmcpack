#ifndef QMCPLUSPLUS_AFQMC_MULTIPURESINGLEDETERMINANT_H
#define QMCPLUSPLUS_AFQMC_MULTIPURESINGLEDETERMINANT_H

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
 * Class that implements a linear combination of
 * pure single determinant trial wave-function.
 * Pure means that the Slater matrix is composed of 0,1 only, 
 * meaning that it represents an eigenstate of the mean-field 
 * solution used to construct the hamiltonian. 
 * Many tricks are used to minimize computation and storage
 * including sparse evalutation and low-rank updates
 */  
class MultiPureSingleDeterminant: public WavefunctionBase
{

  typedef WavefunctionBase* WfnPtr;
  typedef MultiPureSingleDeterminant ThisWfn;
  typedef MultiPureSingleDeterminant* ThisWfnPtr;
  typedef HamiltonianBase* HamPtr;
  typedef std::vector<IndexType>::iterator        VIndexit; 
  typedef std::vector<s1D<ValueType> >::iterator  s1Dit; 
  typedef std::vector<s2D<ValueType> >::iterator  s2Dit; 
  typedef std::vector<s4D<ValueType> >::iterator  s4Dit; 

  public:

    MultiPureSingleDeterminant(Communicate *c):WavefunctionBase(c),trialDensityMatrix_needsupdate(true),ref(0),max_excitation(0),cutoff(1e-5),runtype(0),rotated_hamiltonian(false),wfntype(0),diagHam(true),diag_in_steps(0),iterCI(false)
    {}

    ~MultiPureSingleDeterminant() {}

    bool parse(xmlNodePtr);

    bool setup (HamPtr cur) {
      return getHamiltonian(cur);
    } 

    bool hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors=true); 
    bool hdf_write();

    void evaluateMeanFields()  {}

    void evaluateTrialEnergy(ComplexType& ke, ComplexType& pe); 

    void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G) {}

    void evaluateTwoBodyMixedDensityMatrix() {}

    void evaluateLocalEnergy(const ComplexType* , ComplexType& , ComplexType&, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void evaluateOverlap(const ComplexType*  , ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1 );
    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1 );

    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG,  const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const ComplexType* buff, int ik0, int ikN, int pik0, ComplexSMSpMat&, std::vector<ComplexType>& v, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSpMat&, std::vector<ComplexType>& v, const int n=-1 );
    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<ComplexType> >& vn, const std::vector<IndexType>& vn_indx, ComplexSMSpMat&, std::vector<ComplexType>& v, const int n=-1 );

    ComplexType local_evaluateOverlapSlaterDet(int detn, const ComplexType* SlaterMat);

    int cntExcitations(std::vector<IndexType>&,std::vector<IndexType>&,IndexType&,IndexType&,IndexType&,IndexType&,std::vector<IndexType>&,RealType&);

    // very simple algorithm to generate a ci expansion in orthogonal space iteratively
    void iterativeCI(double cutoff, int nmax, int nmax_int, int maxit);

    bool check_occ_orbs() {
      if(rotated_hamiltonian) return true;
      int ne=NAEA+NAEB,nci = ci.size();
      if(occ_orbs.size()/ne != nci) {
        app_error()<<" Error in check_occ_orbs: nci != orbs.size(): " <<nci <<" " <<occ_orbs.size() <<std::endl;
        return false;
      }
      std::vector<IndexType>::iterator it = occ_orbs.begin();
      for(int i=0; i<nci; i++) {
        for(int k=0; k<NAEA; k++,it++) 
          if( *it < 0 || *it >= NMO ) {
            app_log()<<" Error in check_occ_orbs: det, occ: " <<i <<" ";
            for(int j=0; j<ne; j++) app_log()<<occ_orbs[i*ne+j] <<" ";  
            app_log()<<std::endl;
            return false;
          }
        for(int k=0; k<NAEB; k++,it++)
          if( *it < NMO || *it >= 2*NMO ) {
            app_error()<<" Error in check_occ_orbs: det, occ: " <<i <<" ";
            for(int j=0; j<ne; j++) app_error()<<occ_orbs[i*ne+j] <<" ";
            app_error()<<std::endl;
            return false;
          }
      }
      return true;
    }

  private:

    RealType cutoff;

    bool iterCI;
    int IterCI_maxit;
    double IterCI_cut; 

    bool initFromAscii(std::string fileName);
    bool initFromHDF5(hdf_archive&,const std::string&) {return true;}
    bool initFromXML(std::string fileName) {return true;}


    bool getHamiltonian(HamPtr ); 

    // evaluates and stores mixed density matrix in mixed_density_matrix
    // this evaluates the mixed density matrix in reduced format (only NAEX*NMO non-zero sector)
    void local_evaluateOneBodyMixedDensityMatrix(int, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, ComplexMatrix&, bool full=false);
  
    void local_evaluateOneBodyMixedDensityMatrixFull(const ComplexType* SlaterMat, ComplexType& ovl, ComplexMatrix&, bool full=false);

    void local_rankUpdateOneBodyMixedDensityMatrix(const int ndet, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool full=false);

    void local_evaluateOneBodyTrialDensityMatrix(bool full=false);

    void local_rankUpdateOneBodyTrialDensityMatrix(int n, bool full=false);

    bool diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ComplexMatrix& eigVec, std::vector<IndexType>& occ, int nci, bool eigV=true);

    bool diagHam;

    ValueType NuclearCoulombEnergy; 

    // for every determinant, 
    // vector that contains the list of orbitals occupied in the given Slater Matrix 
    // Indexes are 0-based.
    // Only used for setup/initialization, NOT FOR EXECUTION!!!! 
    //std::vector<vector<IndexType> > occup_alpha; 
    //std::vector<vector<IndexType> > occup_beta; 
    //std::vector<vector<IndexType> > virtual_alpha; 
    //std::vector<vector<IndexType> > virtual_beta; 

    // alternative storage for occupied states for easy access/lookup
    // each MO (from 0...NMO-1) is mapped to either true/false based on occupation
    //std::vector<map<IndexType,bool> > isOcc_alpha; 
    //std::vector<map<IndexType,bool> > isOcc_beta; 

    // for testing purposes only!!!
    ComplexMatrix StoreSM;


    bool rotated_hamiltonian;
    int wfntype;
    std::vector<ComplexType> OrbMat;
    int orbsize;

    int diag_in_steps;

    bool runtype;
    // 0: single copy of hamiltonian with extra computational cost  
    // 1: full storage of all hamiltonians (only option for rotated_hamiltonian) 

    int max_excitation;
    IndexType ref;
    // determinant coefficients
    std::vector<ComplexType> ci; 

    // stores differences wrt to reference determinant
    // stored continuously, tricky to access so careful 
    std::vector<IndexType> ijab; 
    // used to determine the number of excitations and their location in the list 
    std::vector<s2D<IndexType> > excitation_bounds;

    std::vector<IndexType> occ_orbs;
    std::vector<IndexType> occ_pairs;

    // storage for results 
    std::vector<ComplexType> overlaps, pe, ke;

    // 1RDM of the trial density matrix.
    // Used to calculate mean-field energy and mean-field potentials 
    // This is a permutation of the identity matrix, do we need to store it??? 
    bool trialDensityMatrix_needsupdate; 
    // this will always mean the reference determinant trial_density_matrix
    ComplexMatrix trial_density_matrix;
    ComplexMatrix rank_updated_trial_density_matrix;

    // Local storage  
    // Notice that all these matrices are NMOxNAEX, since the rest of the NMO-NAEX columns are zero.
    // Careful must be taken when returning results to other objects in the code.
    ComplexMatrix overlap_inv;
    // this will always mean the reference determinant mixed_density_matrix
    ComplexMatrix full_mixed_density_matrix;
    ComplexMatrix mixed_density_matrix;
    ComplexMatrix rank_updated_mixed_density_matrix;

    // temporary storage
    ComplexMatrix S0,S1,SS0; 
    ComplexVector V0;

    std::vector<IndexType> Iwork;

    std::vector<ComplexType> Cwork;
    std::vector<int> pivot;  

    // One-Body Hamiltonian. Stored in sparse form. 
    // Contains all possible terms used by all determinants
    // Terms are sorted to allow systematic access 
    std::vector<std::vector<s1D<ValueType> > > hij;

    std::vector<std::vector<s1D<ComplexType> > > haj;

    // list of pointers defining the beginning and end 
    // of the elements with first index i 
    std::vector<s1Dit > hij_indx;

    // Storage for two body hamiltonian 
    // Tensor is stored in sparse form. 
    // Contains all possible terms used by all determinants
    // Terms are sorted to allow systematic access 
    std::vector<s2D<ValueType> > Vijkl;

    std::vector<ComplexSMSpMat> SMSpHijkl;

    // list of pointers defining the beginning and end 
    // of the elements with first index i 
    std::vector<s2Dit > Vijkl_indx;
    std::vector<int> Vijkl_nterms_per_det; 

    // pointer to hamiltonian object 
    SparseGeneralHamiltonian* sHam;

    inline IndexType Index2Mat(IndexType I, IndexType J) {
      return (J<NMO)?(I*NMO+J):(I*NMO+J-NMO);
    }

    inline void Mat2Index(const IndexType IJ, IndexType& I, IndexType& J) {
      if( IJ < NMO*NMO) {
        I = IJ/NMO;
        J = IJ%NMO;
      } else {
        I = (IJ-NMO*NMO)/NMO+NMO;
        J = (IJ-NMO*NMO)%NMO+NMO;
      }
    }


};
}

#endif


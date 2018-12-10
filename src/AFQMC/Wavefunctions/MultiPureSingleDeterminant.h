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

    MultiPureSingleDeterminant(Communicate *c):WavefunctionBase(c),trialDensityMatrix_needsupdate(true),ref(0),cutoff(1e-6),runtype(0),rotated_hamiltonian(false),diagHam(true),diag_in_steps(0),iterCI(false),fast_alg(false),test_cnter(0),first_pass(true)
    {}

    ~MultiPureSingleDeterminant() {}

    bool parse(xmlNodePtr);

    bool setup (HamPtr cur) {
      return getHamiltonian(cur);
    } 

    bool hdf_write(hdf_archive& read, const std::string& tag, bool include_tensors=true); 
    bool hdf_write();

    int sizeOfInfoForDistributedPropagation()
    {
      if(closed_shell)
          return 2+NMO*NMO;  // green function becomes dense in evaluation of vbias with rotated_hamiltonian
      else
          return 2+2*NMO*NMO;
    }

    void evaluateMeanFields()  {}

    void evaluateTrialEnergy(ComplexType& ke, ComplexType& pe); 

    void evaluateOneBodyMixedDensityMatrix(const ComplexType* SlaterMat, ComplexMatrix& G) {}

    void evaluateTwoBodyMixedDensityMatrix() {}

    void evaluateLocalEnergy(const ComplexType* , ComplexType& , ComplexType&, ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void evaluateOverlap(const ComplexType*  , ComplexType& ovl_alpha, ComplexType& ovl_beta, const int n=-1 );

    void evaluateOneBodyMixedDensityMatrix(WalkerHandlerBase* wset, SPComplexSMVector* buf, int wlksz, int gfoffset, bool transposed, bool full=true);

    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMSpMat&, std::vector<SPComplexType>& v, const int n=-1 );
    void calculateMeanFieldMatrixElementOfOneBodyOperators(bool addBetaBeta, SPValueSMVector&, std::vector<SPComplexType>& v, const int n=-1 );

    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMSpMat&, SPComplexSMSpMat&, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const SPComplexType* GG, SPValueSMVector&, SPComplexSMVector&, std::vector<SPComplexType>& v, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMSpMat&, SPComplexSMSpMat&, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n=-1);
    void calculateMixedMatrixElementOfOneBodyOperatorsFromBuffer(bool addBetaBeta, const SPComplexType* buff, int ik0, int ikN, SPValueSMVector&, SPComplexSMVector&, std::vector<SPComplexType>& v, int walkerBlock, int nW, bool transposed, bool needsG, const int n=-1);

    void calculateMixedMatrixElementOfTwoBodyOperators(bool addBetaBeta, const ComplexType* SlaterMat, const std::vector<s4D<SPComplexType> >& vn, const std::vector<IndexType>& vn_indx, SPValueSMSpMat&, std::vector<SPComplexType>& v, const int n=-1 );

    ComplexType local_evaluateOverlapSlaterDet(int detn, const ComplexType* SlaterMat);

    int cntExcitations(std::vector<IndexType>&,std::vector<IndexType>&,IndexType&,IndexType&,IndexType&,IndexType&,std::vector<IndexType>&,RealType&);

    int countExct(int N, IndexType* D1, IndexType* D2, bool getIndx, IndexType* loc, IndexType* ik, IndexType* ak, RealType& psign);

    void prepare_excitations();

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
    bool initFromHDF5(hdf_archive&,const std::string&) { return false;}
    bool initFromXML(std::string fileName) { return false;}


    bool getHamiltonian(HamPtr ); 

    // evaluates and stores mixed density matrix in mixed_density_matrix
    // this evaluates the mixed density matrix in reduced format (only NAEX*NMO non-zero sector)
    void local_evaluateOneBodyMixedDensityMatrix(int, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, SPComplexMatrix&, bool full=false);
  
    void local_evaluateOneBodyMixedDensityMatrixFull(const ComplexType* SlaterMat, ComplexType& ovl, SPComplexMatrix&, bool full=false);

    void local_rankUpdateOneBodyMixedDensityMatrix(const int ndet, const ComplexType* SlaterMat, ComplexType& ovl_alpha, ComplexType& ovl_beta, bool full=false);

    void local_evaluateOneBodyTrialDensityMatrix();

    void local_rankUpdateOneBodyTrialDensityMatrix(int n, bool full=false);

    bool diagonalizeTrialWavefunction(std::vector<RealType>& eigVal, ComplexMatrix& eigVec, std::vector<IndexType>& occ, int nci, bool eigV=true);

    void calculate_unique_overlaps(ComplexType& ovl_ref, ComplexType& ovl, const ComplexType* SM);

    void calculate_QFull(ComplexType ovl_ref, bool getQFull, bool getGa, bool getGb);

    void prepare_hamiltonian_evaluation();

    void allocate_hamiltonian_evaluation();

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
    std::vector<ComplexType> OrbMat;
    int orbsize;

    int diag_in_steps;

    bool runtype;
    // 0: single copy of hamiltonian with extra computational cost  
    // 1: full storage of all hamiltonians (only option for rotated_hamiltonian) 

    bool fast_alg;

    IndexType ref;
    // determinant coefficients
    std::vector<ComplexType> ci; 

    // one for each excitation level
    std::vector<ComplexType> Kn;
    std::vector<SPComplexType> SPKn;
    std::vector<ComplexType> VB0;
    ComplexMatrix QFull; 
    ComplexMatrix G0; 

    // list of occupation numbers of all determinants
    std::vector<IndexType> occ_orbs;

    std::vector<RealType> det_sign;  // sign of determinant in iajb convention
    std::vector<ComplexType> ci_with_psign;   // ci coefficient including appropriate sign in iajb convention 

    int maxEx;
    int maxExa;
    int maxExb;
    int nunique_alpha;
    int nunique_beta;

    // unique determinants
    std::vector<IndexType> iajb_unique_alpha;    // compact structure with information on unique dets 
    std::vector<IndexType> iajb_unique_beta; 

    std::vector<std::pair<int,int>> map2unique;   // map from full determinant to unique lists

    std::vector<ComplexType> ovlp_unique_alpha;   // vector of unique overlaps
    std::vector<ComplexType> ovlp_unique_beta;

    SPComplexMatrix BB0inv_alpha;                 // SM*(SM_0)^-1
    SPComplexMatrix BB0inv_beta;
    SPComplexMatrix tBB0inv_alpha;                 
    SPComplexMatrix tBB0inv_beta;

    // storage for reference sector of green functions
    SPComplexMatrix refG;
    SPComplexMatrix Gia;
    SPComplexMatrix Gib;

    // storage for local energy evaluation
    std::vector<SPComplexType> vb;
    std::vector<SPComplexType> vb_helper;

    // storage for M*G
    SPComplexMatrix refP;
    SPComplexMatrix Pia;
    SPComplexMatrix Pib;

    // matrix of boundaries of SMSpHijkl
    Matrix<int> Hijkl_bounds; 
    Matrix<int> local_bounds; 

    SPValueMatrix TVn;

    int test_cnter;
    
    // debug
    bool first_pass; 
    SPValueSMSpMat H2;
    std::vector<s1D<ValueType> > thij;
    SPValueSMSpMat H0;
    SPValueSMSpMat H1;
    std::vector<s1D<ValueType> > tthij;
 

    // storage for results 
    std::vector<ComplexType> overlaps, pe, ke;

    // 1RDM of the trial density matrix.
    // Used to calculate mean-field energy and mean-field potentials 
    // This is a permutation of the identity matrix, do we need to store it??? 
    bool trialDensityMatrix_needsupdate; 
    // this will always mean the reference determinant trial_density_matrix
    ComplexMatrix trial_density_matrix;
    SPComplexMatrix SPtrial_density_matrix;
    ComplexMatrix rank_updated_trial_density_matrix;
    SPComplexMatrix SPrank_updated_trial_density_matrix;

    // Local storage  
    // Notice that all these matrices are NMOxNAEX, since the rest of the NMO-NAEX columns are zero.
    // Careful must be taken when returning results to other objects in the code.
    ComplexMatrix overlap_inv;
    // this will always mean the reference determinant mixed_density_matrix
    ComplexMatrix temp_density_matrix;
    ComplexMatrix temp_density_matrix_full;
    SPComplexMatrix full_mixed_density_matrix;
    SPComplexMatrix mixed_density_matrix;
    SPComplexMatrix rank_updated_mixed_density_matrix;

    std::vector<SPComplexType> cGF;

    std::vector<SPComplexType> local_buff;

    // ik breakup of Spvn
    IndexType ik0, ikN;   //  minimum and maximum values of ik index in Spvn
    IndexType pik0;  // locations of bounds of ik0 sector in Spvn 

    // temporary storage
    ComplexMatrix S0,S1,SS0; 
    SPComplexVector V0;

    std::vector<IndexType> Iwork;

    std::vector<ComplexType> Cwork;
    std::vector<int> pivot;  

    // One-Body Hamiltonian. Stored in sparse form. 
    // Contains all possible terms used by all determinants
    // Terms are sorted to allow systematic access 
    std::vector<std::vector<s1D<ValueType> > > hij;

    std::vector<std::vector<s1D<ComplexType> > > haj;

    std::vector<SPValueSMSpMat> SMSpHijkl;
    std::vector<SPComplexSMSpMat> SMSpHabkl;

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


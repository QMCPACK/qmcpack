
#ifndef QMCPLUSPLUS_AFQMC_SPARSEGENERALHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_SPARSEGENERALHAMILTONIAN_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include<Message/MPIObjectBase.h>
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/config.h"
#include "AFQMC/Hamiltonians/HamiltonianBase.h"
#include "AFQMC/Numerics/SparseMatrixOperations.h"

namespace qmcplusplus
{

/*
 * Note: 
 * 1. New policy, ValueType==real implies a purely complex Dvn, which is stored as real.
 *    For hamiltonians with real integrals but complex propagators (non positive definite), 
 *      you must perform the calculation with ValueType=complex if you want to use Dvn.
 */      


class SparseGeneralHamiltonian: public HamiltonianBase
{

//  template<typename spT> using ShmemAllocator = boost::interprocess::allocator<spT, boost::interprocess::managed_shared_memory::segment_manager>;
//  template<typename spT> using SMVector = boost::interprocess::vector<spT, ShmemAllocator<spT>>;

  typedef std::vector<s1D<ValueType> >::iterator  s1Dit;
  typedef std::vector<s2D<ValueType> >::iterator  s2Dit;
  typedef SMDenseVector<s4D<ValueType> >::iterator  s4Dit;

  public:
 
  SparseGeneralHamiltonian(Communicate *c):HamiltonianBase(c),orderStates(false),cutoff1bar(1e-8),cutoff_cholesky(1e-6),has_full_hamiltonian_for_matrix_elements(false),NMAX(-1),ascii_write_file(""),hdf_write_file(""),printEig(false),factorizedHamiltonian(false),v2full_transposed(false),test_2eint(false),zero_bad_diag_2eints(false),test_algo(true),has_hamiltonian_for_selCI(false),rotation(""),hdf_write_type("default"),inplace(true),skip_V2(false),useSpMSpM(false)
  {
  }

  ~SparseGeneralHamiltonian() {}

  // if ValueType=RealType and sparse=false, Dvn becomes the complex component of the propagator 
  // eventually this will also be true for the sparse format
  void calculateHSPotentials(const RealType cut, const RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool paral); 

  void calculateHSPotentials_Diagonalization(const RealType cut, const RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool paral);

  void calculateHSPotentials_FactorizedHam(const RealType cut, const RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, afqmc::TaskGroup& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool paral);

  void calculateOneBodyPropagator(const RealType cut, const RealType dt, ComplexMatrix& Hadd, std::vector<s2D<ComplexType> >& Pkin); 
  
  bool parse(xmlNodePtr cur); 

  // do a general check of parameters for consistency
  // make sure object was build consistently and correctly 
  bool checkObject(); 

  bool SparseHamiltonianFromFactorization( int indx, std::vector<OrbitalType>& jkl, std::vector<ValueType>& intgs, const RealType cut=1e-6);

  bool createHamiltonianForPureDeterminant(int type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b , std::vector<s1D<ValueType> >& , SPValueSMSpMat&, const RealType cut=1e-6);  

  bool createHamiltonianForGeneralDeterminant(int type, const ComplexMatrix& A,std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vabkl, const RealType cut=1e-6);

  inline bool generateFullHamiltonianForME() {
    
    std::map<IndexType,bool> all_alpha,all_beta;
    all_alpha.clear();
    all_beta.clear();
    for(IndexType i=0; i<NMO; i++) all_alpha[i]=true;
    for(IndexType i=NMO; i<2*NMO; i++) all_alpha[i]=false;
    for(IndexType i=0; i<NMO; i++) all_beta[i]=false;
    for(IndexType i=NMO; i<2*NMO; i++) all_beta[i]=true;

    SpH2_full_forME.setDims(2*NMO*NMO,2*NMO*NMO);
    SpH2_full_forME.setup(head_of_nodes,name+std::string("SpH2_full_forME"),TG.getNodeCommLocal());
    if(!createHamiltonianForPureDeterminant(1,false,all_alpha,all_beta,H1_full_forME,SpH2_full_forME,1e-5)) {
      app_error()<<"Error in createHamiltonianForPureDeterminant during call to generateFullHamiltonianForME(). \n";
      return false;
    }
    has_full_hamiltonian_for_matrix_elements = true;
    return true; 

  }

  inline bool getFullHam(std::vector<s1D<ValueType> > *& h, SPValueSMSpMat *& v) {
    if(!has_full_hamiltonian_for_matrix_elements) 
      generateFullHamiltonianForME();
    v = &SpH2_full_forME;
    h = &H1_full_forME;  
    return true;
  } 

  // should only be used with CIPSI like methods
//  ValueType H(IndexType I, IndexType J) {

 // }

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J) {
    if( (I>=NMO && J<NMO) || (I<NMO && J>=NMO) ) return ValueType(0);
    if(spinRestricted) {
      I = (I>=NMO)?(I-NMO):(I); 
      J = (J>=NMO)?(J-NMO):(J); 
    }
    if(I <= J) {
      s2Dit it = std::lower_bound( H1.begin(), H1.end(), std::forward_as_tuple(I,J,static_cast<ValueType>(0.0)),mySort);

      if (it != H1.end() && std::get<0>(*it) == I && std::get<1>(*it) == J) 
        return std::get<2>(*it);  
      else
        return static_cast<ValueType>(0.0);
    } else {
      s2Dit it = std::lower_bound( H1.begin(), H1.end(), std::forward_as_tuple(J,I,static_cast<ValueType>(0.0)),mySort);

      if (it != H1.end() && std::get<0>(*it) == J && std::get<1>(*it) == I)
        return myconj(std::get<2>(*it));
      else
        return static_cast<ValueType>(0.0);
    }
  }

  void generateIJ()
  {
      long N = spinRestricted?NMO:2*NMO;
      IJ.resize(N*(N+1)/2+1);
      long nt=0;
      OrbitalType i,j,k,l;
      OrbitalType i0=std::get<0>(V2[0]);
      OrbitalType j0=std::get<1>(V2[0]);
      long p2, p1 = mapUT(i0,j0,N);
      for(long ii=0; ii<=p1; ii++) IJ[ii]=0;
      ValueType V;
      for(s4Dit it = V2.begin(); it != V2.end(); it++, nt++) {
        std::tie (i,j,k,l,V) = *it;
        if(i!=i0 || j!=j0) {
          p2 = mapUT(i,j,N);
          for(long ii=p1+1; ii<=p2; ii++) IJ[ii]=nt;
          i0=i;
          j0=j;
          p1=p2;
        }
      }
      for(long ii=p1+1; ii<IJ.size(); ii++) IJ[ii]=nt;
  }

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) 
  {

    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);
    if(factorizedHamiltonian) { 

      if(spinRestricted) {
        I = (I>=NMO)?(I-NMO):(I);
        J = (J>=NMO)?(J-NMO):(J);
        K = (K>=NMO)?(K-NMO):(K);
        L = (L>=NMO)?(L-NMO):(L);
      }

      if(!V2_fact.isCompressed()) {
        app_error()<<" Error: Using uncompressed V2_fact in: SparseGeneralHamiltonian::H(I,J,K,L). " <<std::endl;
        APP_ABORT(" Error: Using uncompressed V2_fact in: SparseGeneralHamiltonian::H(I,J,K,L).");
      }
      SPValueSMSpMat::intPtr cols = V2_fact.column_data();
      SPValueSMSpMat::intPtr rows = V2_fact.row_data();
      SPValueSMSpMat::intPtr indx = V2_fact.row_index();
      SPValueSMSpMat::pointer vals = V2_fact.values();

      int ik = I*NMO+Index2Col(K); 
      int lj = L*NMO+Index2Col(J); 
      ValueType val;
// FIX FIX FIX: lj expansion has to be conjugated for complex orbitals !!!
#if defined(QMC_COMPLEX)
      val = SparseMatrixOperators::product_SpVSpVc<ValueType>(indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[lj+1]-indx[lj],cols+indx[lj],vals+indx[lj]);
#else
      val = SparseMatrixOperators::product_SpVSpV<ValueType>(indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[lj+1]-indx[lj],cols+indx[lj],vals+indx[lj]);
#endif
      return (std::abs(val)>cutoff1bar)?(val):(0);

    } else {

      if(spinRestricted) {
        I = (I>=NMO)?(I-NMO):(I); 
        J = (J>=NMO)?(J-NMO):(J); 
        K = (K>=NMO)?(K-NMO):(K); 
        L = (L>=NMO)?(L-NMO):(L); 
      }
      s4D<ValueType> s = std::make_tuple(I,J,K,L,ValueType(0.0));
      bool cjgt=find_smallest_permutation(s);
      long pp0=0,pp1=V2.size();
      if(IJ.size() != 0) {
        long N = spinRestricted?NMO:2*NMO;
        long p0 = mapUT(std::get<0>(s),std::get<1>(s),N);
        pp0 = IJ[p0];
        pp1 = IJ[p0+1];
      }
      s4Dit it = std::lower_bound( V2.begin()+pp0, V2.begin()+pp1, s, mySort);
  
      if (it != V2.end() &&  std::get<0>(*it)==std::get<0>(s) &&  std::get<1>(*it)==std::get<1>(s) && std::get<2>(*it)==std::get<2>(s) && std::get<3>(*it)==std::get<3>(s) ) {
#if defined(QMC_COMPLEX) 
        if(cjgt) 
          return std::conj(std::get<4>(*it));
        else
#endif
          return std::get<4>(*it);
      } else {
        return static_cast<ValueType>(0.0);
      }
    }
  } 

  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L, std::vector<s4D<ValueType> >& V, int NT ) 
  {

    if( (I>=NT && K<NT) || (I<NT && K>=NT) ) return ValueType(0);
    if( (J>=NT && L<NT) || (J<NT && L>=NT) ) return ValueType(0);
    if(factorizedHamiltonian) {                

      APP_ABORT(" Error: ValueType H(I,J,K,L,V,NT): not implemented with factorized hamiltonian. \n");

      if(spinRestricted) {
        I = (I>=NT)?(I-NT):(I);
        J = (J>=NT)?(J-NT):(J);
        K = (K>=NT)?(K-NT):(K);
        L = (L>=NT)?(L-NT):(L);
      }

      if(!V2_fact.isCompressed()) {
        app_error()<<" Error: Using uncompressed V2_fact in: SparseGeneralHamiltonian::H(I,J,K,L,NT). " <<std::endl;
        APP_ABORT(" Error: Using uncompressed V2_fact in: SparseGeneralHamiltonian::H(I,J,K,L,NT).");
      }
      SPValueSMSpMat::intPtr cols = V2_fact.column_data();
      SPValueSMSpMat::intPtr rows = V2_fact.row_data();
      SPValueSMSpMat::intPtr indx = V2_fact.row_index();
      SPValueSMSpMat::pointer vals = V2_fact.values();

      int ik = I*NT+K; 
      int lj = L*NT+J; 
      if(I>=NT) ik = (I-NT)*NT+(K-NT); 
      if(J>=NT) lj= (L-NT)*NT+(J-NT); 
      ValueType val;

#if defined(QMC_COMPLEX)
      val = SparseMatrixOperators::product_SpVSpVc<ValueType>(indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[lj+1]-indx[lj],cols+indx[lj],vals+indx[lj]);
#else
      val = SparseMatrixOperators::product_SpVSpV<ValueType>(indx[ik+1]-indx[ik],cols+indx[ik],vals+indx[ik],indx[lj+1]-indx[lj],cols+indx[lj],vals+indx[lj]);
#endif
      return (std::abs(val)>cutoff1bar)?(val):(0);

    } else {

      if(spinRestricted) {
        I = (I>=NT)?(I-NT):(I); 
        J = (J>=NT)?(J-NT):(J); 
        K = (K>=NT)?(K-NT):(K); 
        L = (L>=NT)?(L-NT):(L); 
      }
      //s4D<ValueType> s = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(I,J,K,L,static_cast<ValueType>(0.0)));
      s4D<ValueType> s = std::make_tuple(I,J,K,L,ValueType(0.0));
      bool cjgt = find_smallest_permutation(s);
      std::vector<s4D<ValueType> >::iterator it = std::lower_bound( V.begin(), V.end(), s, mySort);
  
      if (it != V.end() &&  std::get<0>(*it)==std::get<0>(s) &&  std::get<1>(*it)==std::get<1>(s) && std::get<2>(*it)==std::get<2>(s) && std::get<3>(*it)==std::get<3>(s) ) {
#if defined(QMC_COMPLEX) 
if(cjgt) 
        if(cjgt)
          return std::conj(std::get<4>(*it));
        else
//        s = find_smaller_equivalent_OneBar_for_integral_list(std::forward_as_tuple(I,J,K,L,std::get<4>(*it)));
//        return std::get<4>(s);
#endif
          return std::get<4>(*it);
      } else {
        return static_cast<ValueType>(0.0);
      }
    }
  } 


  ValueType H_2bar(IndexType I, IndexType J, IndexType K, IndexType L)
  {
    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);
    if( I==J || K==L ) return ValueType(0);
    return H(I,J,K,L) - H(I,J,L,K);
  }

  bool initializeCCProjector(ComplexMatrix& Pmat, RealType cut=1e-6);

  void generate_selCI_Ham(double cutoff); 

  void get_selCI_excitations(OrbitalType I, OrbitalType J, int spinSector, RealType cutoff, OrbitalType* occs, std::vector<OrbitalType>& KLs); 

  protected:

  // name of restart file
  std::string hdf_write_file; 

  std::string hdf_write_type; 

  std::string ascii_write_file; 

  // maximum number of MO in integral file
  int NMAX; 

  bool test_2eint;
  bool zero_bad_diag_2eints;
  bool test_algo; 
  bool inplace;  
  bool skip_V2;

  bool useSpMSpM;

  // stores one body integrals in s2D format 
  std::vector<s2D<ValueType> > H1;

  // shared memory vectors 
  SMDenseVector<s4D<ValueType> > V2;
//  SMDenseVector<s4D<ValueType> > V2_2bar; 

  bool v2full_transposed;
  std::vector<long> KL;
  std::vector<long> IJ;
  SMDenseVector<s4D<ValueType> > V2_full; 

  std::string rotation;;
  ValueMatrix rotationMatrix; 

  // 
  bool factorizedHamiltonian;  

  // factorized Ham : V2(ik,lj) = sum_n V2_fact(ik,n)*conj(V2_fact(lj,n)) 
  // similar to Spvn, without -dt/2 factor and without L+L* / L-L* rotation 
// NOTE: Make this identical to Spvn and return pointer to this object to propagator
// this avoids having 2 copied when reading in 3Index form 
  SPValueSMSpMat V2_fact; 

  std::vector<int> cholesky_residuals;

  bool has_hamiltonian_for_selCI;
// This is going to be a problem with enough orbitals, e.g. NMO ~> 1000
  int nmax_KL_selCI;
  IndexVector IJ_aa, IJ_bb, IJ_ab;  
  SMDenseVector<s2D<ValueType> > V2_selCI_aa; 
  SMDenseVector<s2D<ValueType> > V2_selCI_ab; 
  SMDenseVector<s2D<ValueType> > V2_selCI_bb; 

  bool has_full_hamiltonian_for_matrix_elements;
  SPValueSMSpMat SpH2_full_forME;
  std::vector<s1D<ValueType> > H1_full_forME;

  // do we need to sort the states according to eigenvalue?
  bool orderStates;

  bool printEig;

  bool close_shell;

  // cutoff to read 1bar terms in hamiltonian matrices 
  double cutoff1bar; 

  // cutoff for cholesky
  double cutoff_cholesky;  

  // needed only if we are ordering states (Molpro issue)
  // occup might be different from wfn, so do not mix
  std::vector<IndexType> occup_alpha;
  std::vector<IndexType> virtual_alpha;
  std::map<IndexType,bool> isOcc_alpha;
  std::vector<IndexType> occup_beta;
  std::vector<IndexType> virtual_beta;
  std::map<IndexType,bool> isOcc_beta;

  // right now assuming core states are the bottom NCX states in the list, (possibly after reordering).
  // Generalize this to partition states into (core,occ,virtual) based on arbitrary lists. 
  std::vector<IndexType> core_alpha;
  std::map<IndexType,bool> isCore_alpha;
  std::vector<IndexType> core_beta;
  std::map<IndexType,bool> isCore_beta;

  std::vector<ValueType> eig;

  std::vector<IndexType> orbSymm;
  std::vector<IndexType> occupPerSymm_alpha;
  std::vector<IndexType> occupPerSymm_beta;

  // store diagonal part of hamiltonian: Diag(i,k) = H(i,k,i,k);
  ValueMatrix DiagHam;

  bool initFromASCII(const std::string& fileName) 
  {
     // allow for different formats later on.
     // right now only FCIDUMP is allowed
     if(filetype == "fcidump")   
       return readFCIDUMP(fileName, false); //, H1, V2);
     else
       return false;
  } 

  bool initFromXML(const std::string& fileName) { return false;} 

  bool initFromHDF5(const std::string& fileName); 

  void hdf_write();

  void ascii_write();

//  bool generate_V2_2bar();

  // read integrals from ascii file, assuming FCIDUMP format
  // returns them in sparse form, using extended indexing to
  // encode alpha/beta, so index goes from 0-2*NMO-1 
  // initFromFCIDUMP takes care of storing them in the correct format of the 
  // derived class
  bool readFCIDUMP(const std::string& fileName, bool minimizeIO);//,
//      std::vector<s2D<ValueType> >& , std::vector<s4D<ValueType> >& );  

  // count number of elements in file
  bool countElementsFromFCIDUMP(std::ifstream&,int&,int&,int&,int&,int&,int&,int&,int&,std::map<IndexType,IndexType>&, std::map<IndexType,IndexType>&,int& n); 

  // read elements in FCIDUMP 
  bool readElementsFromFCIDUMP(std::ifstream&,std::vector<s2D<ValueType> >&, std::vector<s2D<ValueType> >&, SMDenseVector<s4D<ValueType> >&, std::vector<s4D<ValueType> >&, std::vector<s4D<ValueType> >&, ValueSMSpMat&, ValueSMSpMat&, ValueSMSpMat&, std::map<IndexType,IndexType>&, std::map<IndexType,IndexType>&); 

  // find all permutation of indexes among symmetry equivalent terms
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  void find_equivalent_OneBar_for_integral_list(s4D<ValueType> ijkl, std::vector<s4D<ValueType> >& v); 

  // find all permutation of indexes among symmetry equivalent terms
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  void find_equivalent_TwoBar_for_integral_list(s4D<ValueType> ijkl, std::vector<s4D<ValueType> >& v); 
  // find all permutation of indexes among symmetry equivalent terms
  // eliminates redundant terms and adjusts weight
  // redundant terms are those that are repeated when contracting against G*G
  // example: ik/jl and jl/ik. Instead of keeping 2, keep one and multiply V by 2.
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  void find_equivalent_OneBar_for_hamiltonian_generation(s4D<ValueType> ijkl, std::vector<s4D<ValueType> >& v); 

  // find smaller permutation of indexes among symmetry equivalent terms
  // This is the one stored in V2. 
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  void find_equivalent_TwoBar_for_hamiltonian_generation(s4D<ValueType> ijkl, std::vector<s4D<ValueType> >& v); 

  // more efficient version of find_smaller_equivalent_OneBar_for_integral_list
  bool find_smallest_permutation(s4D<ValueType>& ijkl);

  // find smaller permutation of indexes among symmetry equivalent terms
  // This is the one stored in V2. 
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  s4D<ValueType> find_smaller_equivalent_OneBar_for_integral_list(s4D<ValueType> ijkl); 

  // find smaller permutation of indexes among symmetry equivalent terms
  // This is the one stored in V2_2bar. 
  // NOT TU BE USED OUTSIDE INITIALIZATION! SLOW!
  s4D<ValueType> find_smaller_equivalent_TwoBar_for_integral_list(s4D<ValueType> ijkl); 

  //bool getFCIDUMPline(std::ifstream& in, ValueType& val, IndexType& ap, IndexType& bp, IndexType& cp, IndexType& dp);

  inline int getSpinSector(const IndexType& i, const IndexType& j, const IndexType& k, const IndexType& l) {
    if(i < NMO) {
      if(j < NMO) return 0;  // <alpha,alpha | alpha,alpha>
      else        return 1;  // <alpha,beta  | alpha,beta >
    } else {
      if(j < NMO) return 2;  // <beta,alpha | beta,alpha>
      else        return 3;  // <beta,beta  | beta,beta >
    }
  } 

  inline bool goodSpinSector(const IndexType& i, const IndexType& j, const IndexType& k, const IndexType& l, int NT) {
    if(i < NT) {
      if(j < NT) // <alpha,alpha | alpha,alpha> 
        return (k<NT&&l<NT); 
      else        // <alpha,beta  | alpha,beta >
        return (k<NT&&l>=NT); 
    } else {
      if(j < NT) // <beta,alpha | beta,alpha>
        return (k>=NT&&l<NT); 
      else        // <beta,beta  | beta,beta >
        return (k>=NT&&l>=NT); 
    }
  }

  inline int getSpinSector(const IndexType& i, const IndexType& j) {
    if(i < NMO) return 0; 
    return 1; 
  } 

  inline IndexType Index2Col(IndexType i) {
#if AFQMC_DEBUG
// assert( ((i<NMO)&&(j<NMO)) || ((i>NMO)&&(j>NMO))   )
#endif
   return (i<NMO)?(i):(i-NMO);
  }

  inline IndexType Index2Mat(IndexType i, IndexType j, bool GHF=false) {
#if AFQMC_DEBUG
// assert( ((i<NMO)&&(j<NMO)) || ((i>NMO)&&(j>NMO))   )
#endif
   if(GHF) 
     return i*2*NMO+j; 
   else
     return (i<NMO)?(i*NMO+j):(NMO*NMO+(i-NMO)*NMO+j-NMO);
  }

  // used to sort snD values using only indexes 
  _mySort_snD_ mySort;

  // used to identify equal index sets (value is not compared)  
  _myEqv_snD_ myEqv;

  void find_all_contributions_to_hamiltonian_closed_shell(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v); 
  void find_all_contributions_to_hamiltonian_spinRestricted(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v); 
  void find_all_contributions_to_hamiltonian_general(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v); 
  void find_all_contributions_to_hamiltonian_ghf(bool aa_only, OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType J1, ValueType J2, ValueType J3, ValueType J1a, ValueType J2a, ValueType J3a, double cut, std::vector<s4D<ValueType> >& v); 

  int count_allowed_terms(std::vector<s4D<ValueType> >& vs4D, std::map<IndexType,bool>&  occ_a, std::map<IndexType,bool>& occ_b)  
  {
    int cnt=0;
    for(std::vector<s4D<ValueType> >::iterator it = vs4D.begin(); it!=vs4D.end(); it++) 
      if( (occ_a[std::get<0>(*it)]||occ_b[std::get<0>(*it)]) && (occ_a[std::get<1>(*it)]||occ_b[std::get<1>(*it)]) ) cnt++;
    return cnt;
  }

  inline void push_ijkl(OrbitalType i, OrbitalType j, OrbitalType k, OrbitalType l, ValueType V, std::vector<s4D<ValueType> >& v, bool GHF=false) { 
    long ik = Index2Mat(i,k,GHF);  
    long jl = Index2Mat(j,l,GHF);
    if( ik == jl )
      v.push_back(std::make_tuple(i,j,k,l,V));
    else if(ik < jl)
      v.push_back(std::make_tuple(i,j,k,l,ValueType(2.0)*V));
    else
      v.push_back(std::make_tuple(j,i,l,k,ValueType(2.0)*V));
  }  

  int add_allowed_terms(std::vector<s4D<ValueType> >& vs4D, std::map<IndexType,bool>&  occ_a, std::map<IndexType,bool>& occ_b, SPValueSMSpMat& V , bool needs_locks=false, bool GHF=false)
  {
    int cnt=0;
    for(std::vector<s4D<ValueType> >::iterator it = vs4D.begin(); it!=vs4D.end(); it++) 
      if( (occ_a[std::get<0>(*it)]||occ_b[std::get<0>(*it)]) && (occ_a[std::get<1>(*it)]||occ_b[std::get<1>(*it)]) ) { 
        cnt++;
        V.add( Index2Mat(std::get<0>(*it),std::get<2>(*it),GHF) , Index2Mat(std::get<1>(*it),std::get<3>(*it),GHF) , static_cast<SPValueType>(std::get<4>(*it)), needs_locks);
      }
    return cnt;
  }

  void print_tuple(s1D<ValueType>& t) {
    std::cout<<"  -  " <<std::get<0>(t) <<" " <<std::get<1>(t) <<std::endl;
  }

  void print_tuple(s4D<ValueType>& v) {
    std::cout<<std::get<4>(v) <<" " 
             <<std::get<0>(v) <<" "
             <<std::get<1>(v) <<" "
             <<std::get<2>(v) <<" "
             <<std::get<3>(v) <<std::endl;
  }

  void print_Vs4D(std::vector<s4D<ValueType> >& v) {
    for(int i=0; i<v.size(); i++)
      std::cout<<std::get<4>(v[i]) <<" " 
               <<std::get<0>(v[i]) <<" "
               <<std::get<1>(v[i]) <<" "
               <<std::get<2>(v[i]) <<" "
               <<std::get<3>(v[i]) <<std::endl;
    std::cout<<std::endl;
  }

  inline long mapUT(long i, long j, long N) {
    if(j >= i)
      return N*i + j - (i*(i+1))/2;
    else
      return N*j + i - (j*(j+1))/2;
  }

  inline long mapUT_woD(long i, long j, long N) {
    if(j == i) {
      APP_ABORT(" Error in mapUT_woD: This should not happen. \n");
    } else if(j > i)
      return N*i + j - (i*(i+1))/2 - i-1;
    return N*j + i - (j*(j+1))/2 - j-1;
  }

  inline int mapUT(int i, int j, int N) {
    if(j >= i)
      return N*i + j - (i*(i+1))/2;
    return N*j + i - (j*(j+1))/2;
  }

  inline int mapUT_woD(int i, int j, int N) {
    if(j == i) {
      APP_ABORT(" Error in mapUT_woD: This should not happen. \n");
    } else if(j > i)
      return N*i + j - (i*(i+1))/2 - i-1;
    return N*j + i - (j*(j+1))/2 - j-1;
  }

  bool communicate_Vijkl(SPValueSMSpMat&);

#ifndef QMC_COMPLEX
  bool communicate_Vijkl(SPComplexSMSpMat&);
#endif

};
}

#endif


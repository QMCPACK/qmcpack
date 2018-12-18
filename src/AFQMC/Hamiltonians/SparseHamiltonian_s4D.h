
#ifndef QMCPLUSPLUS_AFQMC_SPARSEGENERALHAMILTONIAN_S4D_H
#define QMCPLUSPLUS_AFQMC_SPARSEGENERALHAMILTONIAN_S4D_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>
#include "OhmmsData/libxmldefs.h"

#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Hamiltonians/OneBodyHamiltonian.hpp"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

#include "AFQMC/Hamiltonians/Hamiltonian_Utilities.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"

namespace qmcplusplus
{
namespace afqmc
{

class SparseHamiltonian_s4D: public OneBodyHamiltonian
{

  typedef std::vector<s1D<ValueType> >::iterator  s1Dit;
  typedef std::vector<s2D<ValueType> >::iterator  s2Dit;
  typedef SMDenseVector<s4D<ValueType> >::iterator  s4Dit;

  public:

// NOTE: implement this by using a csr_matrix and map: {I,J,K,L,val} -> (mapUT(I,J), mapUT(K,L), val) 
 
  SparseHamiltonian_s4D(AFQMCInfo const& info, xmlNodePtr cur, std::vector<s2D<ValueType> >&& h, 
                        SpVType_shm_csr_matrix&& v2_, TaskGroup_& tg_, 
                        ValueType nucE=0, ValueType fzcE=0):
                            OneBodyHamiltonian(info,std::move(h),nucE,fzcE),
                            TG(tg_),
                            V2(std::move(v2_)),
                            cutoff1bar(1e-8),cutoff_cholesky(1e-6),
                            skip_V2(false),printEig(false),test_breakup(false),
                            v2_full_init(false),
                            V2_full({0,0},{0,0},0,
                                    boost::mpi3::intranode::allocator<SPValueType>(tg_.Node()))
  {
    distribute_Ham = (TG.getNumberOfTGs() > 1 );

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    std::string str("no");
    std::string str1("no");
    std::string str2("no");
    std::string str3("no");
    std::string str4("no");
    std::string str5("yes");
    ParameterSet m_param;
    m_param.add(cutoff1bar,"cutoff_1bar","double");
    m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
    m_param.add(cutoff_cholesky,"cutoff_factorization","double");
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(str,"skip_V2","std::string");
    m_param.add(str1,"printEig","std::string");
    m_param.add(str2,"test_2eint","std::string");
    m_param.add(str3,"fix_2eint","std::string");
    m_param.add(str4,"test_breakup","std::string");
    m_param.add(str5,"inplace","std::string");
    m_param.put(cur);

    std::transform(str.begin(),str.end(),str.begin(),(int (*)(int))tolower);
    std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int))tolower);
    std::transform(str2.begin(),str2.end(),str2.begin(),(int (*)(int))tolower);
    std::transform(str3.begin(),str3.end(),str3.begin(),(int (*)(int))tolower);
    std::transform(str4.begin(),str4.end(),str4.begin(),(int (*)(int))tolower);
    std::transform(str5.begin(),str5.end(),str5.begin(),(int (*)(int))tolower);
    if(str1 == "yes" || str1 == "true") printEig = true;
    if(str2 == "yes" || str2 == "true") test_2eint = true;
    if(str3 == "yes" || str3 == "true") zero_bad_diag_2eints = true;
    if(str4 == "yes" || str4 == "true") test_breakup = true; 
    if(str5 == "no" || str5 == "false") inplace = false;
    if(str == "yes" || str == "true") skip_V2 = true;

    if(skip_V2)
      app_log()<<" Skipping 2 electron integrals. Only correct if other objects are initialized from hdf5 files. \n";

  }

  ~SparseHamiltonian_s4D() {}

  SparseHamiltonian_s4D(SparseHamiltonian_s4D const& other) = delete;
  SparseHamiltonian_s4D(SparseHamiltonian_s4D && other) = default;
  SparseHamiltonian_s4D& operator=(SparseHamiltonian_s4D const& other) = delete;
  SparseHamiltonian_s4D& operator=(SparseHamiltonian_s4D && other) = default;

  boost::multi_array<ComplexType,2> getH1() const{ return OneBodyHamiltonian::getH1(); }

  // dummy routines for testing
  boost::multi_array<SPComplexType,1> halfRotatedHij(WALKER_TYPES type, PsiT_Matrix *Alpha, PsiT_Matrix *Beta) {
    APP_ABORT(" Finish Hij. \n");
    return boost::multi_array<SPComplexType,1>(extents[1]);
  }

  SpCType_shm_csr_matrix halfRotatedHijkl(WALKER_TYPES type, TaskGroup_& TGHam, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, const RealType cut=1e-6){
    APP_ABORT(" Finish Hijkl. \n");
    using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
    return SpCType_shm_csr_matrix({0,0},{0,0},0,Alloc(TG.Node()));
  }

  SpVType_shm_csr_matrix calculateHSPotentials(double cut, TaskGroup_& TGprop,
        boost::multi_array<ComplexType,2>& vn0) {
    APP_ABORT(" Finish HS. \n");
    using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
    return SpCType_shm_csr_matrix({0,0},{0,0},0,Alloc(TG.Node()));
  }

  HamiltonianOperations getHamiltonianOperations(bool pureSD, WALKER_TYPES type, 
            std::vector<PsiT_Matrix>& PsiT, double cutvn, double cutv2,
            TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& dump) {
    APP_ABORT("Calling getHamiltonianOperations S4D\n");
    return HamiltonianOperations{};
  }

  ValueType H(IndexType I, IndexType J) const  
  {  return OneBodyHamiltonian::H(I,J); }

  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) const 
  {

/*
    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);

    I = (I>=NMO)?(I-NMO):(I); 
    J = (J>=NMO)?(J-NMO):(J); 
    K = (K>=NMO)?(K-NMO):(K); 
    L = (L>=NMO)?(L-NMO):(L); 

    s4D<ValueType> s = std::make_tuple(I,J,K,L,ValueType(0.0));
    bool cjgt=find_smallest_permutation(s);
    long pp0=0,pp1=V2->size();
    if(IJ.size() != 0) {
      long N = NMO; 
      long p0 = mapUT(std::get<0>(s),std::get<1>(s),N);
      pp0 = IJ[p0];
      pp1 = IJ[p0+1];
    }
    SMDenseVector<s4D<ValueType> >::const_iterator it = std::lower_bound( V2.begin()+pp0, V2.begin()+pp1, s,                  
      [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs)
      {
        return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
      }
    );
  
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
*/
    return ValueType(0.0);
  } 

  ValueType H_2bar(IndexType I, IndexType J, IndexType K, IndexType L) const
  {
    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);
    if( I==J || K==L ) return ValueType(0);
    return H(I,J,K,L) - H(I,J,L,K);
  }

// For testing only, erased when classes are developed
  void calculateHSPotentials(RealType cut, const RealType dt, ComplexMatrix& vn0, SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, TaskGroup_& TGprop, std::vector<int>& nvec_per_node, bool sparse, bool paral)
  {  APP_ABORT(" Error: Calling old Hamiltonian routine from new classes. \n"); }

  bool createHamiltonianForPureDeterminant(int type, bool aa_only, std::map<IndexType,bool>& occ_a, std::map<IndexType,bool>& occ_b , std::vector<s1D<ValueType> >& , SPValueSMSpMat&, const RealType cut=1e-6)
  {  APP_ABORT(" Error: Calling old Hamiltonian routine from new classes. \n"); return false; }

  bool createHamiltonianForGeneralDeterminant(int type, const ComplexMatrix& A,std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vabkl, const RealType cut=1e-6)
  {  APP_ABORT(" Error: Calling old Hamiltonian routine from new classes. \n"); return false; }

  protected:

  // for hamiltonian distribution 
  TaskGroup_& TG;

  bool distribute_Ham;

  // cutoff to read 1bar terms in hamiltonian matrices 
  double cutoff1bar;

  // cutoff for cholesky
  double cutoff_cholesky;

  bool test_2eint;
  bool zero_bad_diag_2eints;
  bool inplace;  
  bool printEig;
  bool skip_V2;
  bool test_breakup;

  // shared memory vectors 
  SpVType_shm_csr_matrix V2;

  SpVType_shm_csr_matrix V2_full; 
  bool v2_full_init;

  // store diagonal part of hamiltonian: Diag(i,k) = H(i,k,i,k);
  boost::multi_array<ValueType,2> DiagHam;

};
}
}

#endif


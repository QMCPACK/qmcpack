
#ifndef QMCPLUSPLUS_AFQMC_SYMMETRICFACTORIZEDSPARSEHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_SYMMETRICFACTORIZEDSPARSEHAMILTONIAN_H

#include<iostream>
#include<vector> 
#include<map> 
#include<fstream>

#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"

#include "AFQMC/Hamiltonians/OneBodyHamiltonian.hpp"

#include "AFQMC/Hamiltonians/Hamiltonian_Utilities.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/Hamiltonians/rotateHamiltonian.hpp"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class SymmetricFactorizedSparseHamiltonian: public OneBodyHamiltonian 
{

  typedef std::vector<s1D<ValueType> >::iterator  s1Dit;
  typedef std::vector<s2D<ValueType> >::iterator  s2Dit;

  public:

  using shm_csr_matrix = SpVType_shm_csr_matrix;
  using csr_matrix_view = shm_csr_matrix::template matrix_view<int>; 

 
  SymmetricFactorizedSparseHamiltonian(AFQMCInfo const& info, xmlNodePtr cur, std::vector<s2D<ValueType> >&& h, 
                              shm_csr_matrix&& v2_, TaskGroup_& tg_, ValueType nucE=0, ValueType fzcE=0):
                                    OneBodyHamiltonian(info,std::move(h),nucE,fzcE),
                                    TG(tg_),V2_fact(std::move(v2_)),
                                    cutoff1bar(1e-8),cutoff_cholesky(1e-6),skip_V2(false),
                                    maximum_buffer_size(1024)
  {

    distribute_Ham = (TG.getNumberOfTGs() > 1 ); 

    if( distribute_Ham ) 
        APP_ABORT(" Error: Distributed SymmetricFactorizedHamiltonian not yet implemented.\n");
    
    if( !distribute_Ham ) 
      assert(V2_fact.global_origin()[0]==0 && V2_fact.global_origin()[1]==0); 

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    std::string str("no");
    ParameterSet m_param;
    m_param.add(cutoff1bar,"cutoff_1bar","double");
    m_param.add(cutoff_cholesky,"cutoff_decomposition","double");
    m_param.add(cutoff_cholesky,"cutoff_factorization","double");
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(str,"skip_V2","std::string");
    m_param.add(maximum_buffer_size,"buffer_size","int");
    m_param.put(cur);

    std::transform(str.begin(),str.end(),str.begin(),(int (*)(int))tolower);
    if(str == "yes" || str == "true") skip_V2 = true;

    if(skip_V2)
      app_log()<<" Skipping 2 electron integrals. Only correct if other objects are initialized from hdf5 files. \n";

  }

  ~SymmetricFactorizedSparseHamiltonian() {}

  SymmetricFactorizedSparseHamiltonian(SymmetricFactorizedSparseHamiltonian const& other) = delete;
  SymmetricFactorizedSparseHamiltonian(SymmetricFactorizedSparseHamiltonian && other) = default;
  SymmetricFactorizedSparseHamiltonian& operator=(SymmetricFactorizedSparseHamiltonian const& other) = delete;
  SymmetricFactorizedSparseHamiltonian& operator=(SymmetricFactorizedSparseHamiltonian && other) = default;

  boost::multi_array<ComplexType,2> getH1() const{ return OneBodyHamiltonian::getH1(); }

  // Haj[0:NAEA*NMO] for type CLOSED,
  // Haj[0:(NAEA+NAEB)*NMO] for type COLLINEAR,
  // Haj[0:(NAEA+NAEB)*2*NMO] for type NONCOLLINEAR
  boost::multi_array<SPComplexType,1> halfRotatedHij(WALKER_TYPES type, PsiT_Matrix *Alpha, PsiT_Matrix *Beta) {
    check_wavefunction_consistency(type,Alpha,Beta,NMO,NAEA,NAEB);
    return rotateHij(type,NMO,NAEA,NAEB,Alpha,Beta,OneBodyHamiltonian::H1);
  }

  SpCType_shm_csr_matrix halfRotatedHijkl(WALKER_TYPES type, TaskGroup_& TGHam, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, RealType const cut=1e-6);

  SpVType_shm_csr_matrix calculateHSPotentials(double cut, TaskGroup_& TGprop,
        boost::multi_array<ComplexType,2>& vn0); 

  HamiltonianOperations getHamiltonianOperations(bool pureSD, WALKER_TYPES type, 
            std::vector<PsiT_Matrix>& PsiT, double cutvn, double cutv2,
            TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& dump); 

  ValueType H(IndexType I, IndexType J) const 
  {  return OneBodyHamiltonian::H(I,J); }
 
  // this should never be used outside initialization routines.
  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) const 
  {

    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);

    I = (I>=NMO)?(I-NMO):(I);
    J = (J>=NMO)?(J-NMO):(J);
    K = (K>=NMO)?(K-NMO):(K);
    L = (L>=NMO)?(L-NMO):(L);

    int ik = I*NMO+Index2Col(NMO,K); 
    int jl = J*NMO+Index2Col(NMO,L); 
    ValueType val = csr::csrvv<ValueType>('N','N',V2_fact.sparse_row(ik),V2_fact.sparse_row(jl));
    return (std::abs(val)>cutoff1bar)?(val):(0);

  } 
  ValueType H_2bar(IndexType I, IndexType J, IndexType K, IndexType L) const
  {
    if( (I>=NMO && K<NMO) || (I<NMO && K>=NMO) ) return ValueType(0);
    if( (J>=NMO && L<NMO) || (J<NMO && L>=NMO) ) return ValueType(0);
    if( I==J || K==L ) return ValueType(0);
    return H(I,J,K,L) - H(I,J,L,K);
  }


  protected:

  // for hamiltonian distribution 
  TaskGroup_& TG;

  bool distribute_Ham;

  // factorized Ham : V2(ik,lj) = sum_n V2_fact(ik,n)*conj(V2_fact(lj,n)) 
  shm_csr_matrix V2_fact; 

  // store diagonal part of hamiltonian: Diag(i,k) = H(i,k,i,k);
  ValueMatrix DiagHam;

  // options read from xml
  double cutoff1bar;
  double cutoff_cholesky;
  bool skip_V2; 

  int maximum_buffer_size;
 
};

}
}

#endif


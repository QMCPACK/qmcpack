//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIAN_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIAN_HPP

#include<fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Hamiltonians/SymmetricFactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D.h"
//#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian_old.h"
//#include "AFQMC/Hamiltonians/SparseHamiltonian_s4D_old.h"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

namespace dummy
{
/*
 * Empty class to avoid need for default constructed Hamiltonians.
 * Throws is any visitor is called. 
 */
class dummy_Hamiltonian 
{
  public:
  dummy_Hamiltonian() {};

  ValueType getNuclearCoulombEnergy() const { 
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  } 

  ValueType H(IndexType i, IndexType j) {
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  } 

  ValueType H(IndexType i, IndexType j, IndexType k, IndexType l) {
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  } 

  boost::multi_array<ComplexType,2> getH1() const{ return boost::multi_array<ComplexType,2>{}; }

  void createHamiltonianForGeneralDeterminant(int type, const ComplexMatrix& A,
                    std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vabkl,
                    const RealType cut=1e-6) {
    throw std::runtime_error("calling visitor on dummy object");
  }

  boost::multi_array<SPComplexType,1> halfRotatedHij(WALKER_TYPES type, PsiT_Matrix *Alpha, PsiT_Matrix *Beta)
  {
    throw std::runtime_error("calling visitor on dummy object");
    return boost::multi_array<ComplexType,1>(extents[1]);
  }

  SpCType_shm_csr_matrix halfRotatedHijkl(WALKER_TYPES type, TaskGroup_& TGWfn, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, const RealType cut=1e-6) 
  {
    throw std::runtime_error("calling visitor on dummy object");
    using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
    return SpCType_shm_csr_matrix({0,0},{0,0},0,Alloc(TGWfn.Node()));
  }

  void calculateHSPotentials(RealType cut, const RealType dt, ComplexMatrix& vn0, 
        SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, TaskGroup_& TGprop, 
        std::vector<int>& nvec_per_node, bool sparse, bool paral)
  {
    throw std::runtime_error("calling visitor on dummy object");  
  }  

  SpVType_shm_csr_matrix calculateHSPotentials(double cut, TaskGroup_& TGprop, 
        boost::multi_array<ComplexType,2>& vn0)
  {
    throw std::runtime_error("calling visitor on dummy object");
    using Alloc = boost::mpi3::intranode::allocator<SPComplexType>;
    return SpCType_shm_csr_matrix({0,0},{0,0},0,Alloc(TGprop.Node()));
  }

  template<class... Args>
  HamiltonianOperations getHamiltonianOperations(Args&&... args) 
  //HamiltonianOperations getHamiltonianOperations(bool pureSD, WALKER_TYPES type, 
  //          std::vector<PsiT_Matrix>& PsiT, double cutvn, double cutv2,
  //          TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive *dump=nullptr)
  {
    throw std::runtime_error("calling visitor on dummy object");
    return HamiltonianOperations{}; 
  }

};
}

class Hamiltonian: public boost::variant<dummy::dummy_Hamiltonian,FactorizedSparseHamiltonian,SymmetricFactorizedSparseHamiltonian,SparseHamiltonian_s4D,THCHamiltonian> //,FactorizedSparseHamiltonian_old,SparseHamiltonian_s4D_old>
{
    using shm_csr_matrix = FactorizedSparseHamiltonian::shm_csr_matrix; 
    using csr_matrix_vew = FactorizedSparseHamiltonian::csr_matrix_view; 

    public: 

    Hamiltonian() { 
      APP_ABORT(" Error: Reached default constructor of Hamiltonian. \n");  
    } 
    explicit Hamiltonian(THCHamiltonian&& other) : variant(std::move(other)) {}
    explicit Hamiltonian(FactorizedSparseHamiltonian&& other) : variant(std::move(other)) {}
    explicit Hamiltonian(SymmetricFactorizedSparseHamiltonian&& other) : variant(std::move(other)) {}
    explicit Hamiltonian(SparseHamiltonian_s4D&& other) : variant(std::move(other)) {}
//    explicit Hamiltonian(FactorizedSparseHamiltonian_old&& other) : variant(std::move(other)) {}
//    explicit Hamiltonian(SparseHamiltonian_s4D_old&& other) : variant(std::move(other)) {}

    explicit Hamiltonian(THCHamiltonian const& other) = delete;
    explicit Hamiltonian(FactorizedSparseHamiltonian const& other) = delete;
    explicit Hamiltonian(SymmetricFactorizedSparseHamiltonian const& other) = delete;
    explicit Hamiltonian(SparseHamiltonian_s4D const& other) = delete; 
//    explicit Hamiltonian(FactorizedSparseHamiltonian_old const& other) = delete;
//    explicit Hamiltonian(SparseHamiltonian_s4D_old const& other) = delete; 

    Hamiltonian(Hamiltonian const& other) = delete; 
    Hamiltonian(Hamiltonian&& other) = default; 

    Hamiltonian& operator=(Hamiltonian const& other) = delete; 
    Hamiltonian& operator=(Hamiltonian&& other) = default; 

    ValueType getNuclearCoulombEnergy() {
        return boost::apply_visitor(
            [&](auto&& a){return a.getNuclearCoulombEnergy();},
            *this
        );
    }
    
    boost::multi_array<ComplexType,2> getH1() const{ 
        return boost::apply_visitor(
            [&](auto&& a){return a.getH1();},
            *this
        );
    }
    
    ValueType H(IndexType i, IndexType j) {
        return boost::apply_visitor(
            [&](auto&& a){return a.H(i,j);},
            *this
        );
    }

    ValueType H(IndexType i, IndexType j, IndexType k, IndexType l) {
        return boost::apply_visitor(
            [&](auto&& a){return a.H(i,j,k,l);},
            *this
        );
    }

    template<class... Args>
    HamiltonianOperations getHamiltonianOperations(Args&&... args) {
        return boost::apply_visitor(
            [&](auto&& a){return a.getHamiltonianOperations(std::forward<Args>(args)...);},
            *this
        );
    }

    // Everything below this point goes away when the code is complete. 
    // Right now it exists for testing mainly
    /*
    void createHamiltonianForGeneralDeterminant(int type, const ComplexMatrix& A, 
                    std::vector<s1D<ComplexType> >& hij, SPComplexSMSpMat& Vabkl, 
                    const RealType cut=1e-6) {
        boost::apply_visitor(
            [&](auto&& a){a.createHamiltonianForGeneralDeterminant(type,A,hij,Vabkl,cut);},
            *this
        );
    }

    boost::multi_array<SPComplexType,1> halfRotatedHij(WALKER_TYPES type, PsiT_Matrix *Alpha, PsiT_Matrix *Beta)
    {
        return boost::apply_visitor(
            [&](auto&& a){return a.halfRotatedHij(type,Alpha,Beta);},
            *this
        );
    }

    SpCType_shm_csr_matrix halfRotatedHijkl(WALKER_TYPES type, TaskGroup_& TGHam, PsiT_Matrix *Alpha, PsiT_Matrix *Beta, const RealType cut=1e-6)
    {
        return boost::apply_visitor(
            [&](auto&& a){return a.halfRotatedHijkl(type,TGHam,Alpha,Beta,cut);},
            *this
        );
    }

    void calculateHSPotentials(RealType cut, const RealType dt, ComplexMatrix& vn0, 
        SPValueSMSpMat& Spvn, SPValueSMVector& Dvn, TaskGroup_& TGprop, 
        std::vector<int>& nvec_per_node, bool sparse, bool paral) {
        boost::apply_visitor(
            [&](auto&& a){a.calculateHSPotentials(cut,dt,vn0,Spvn,Dvn,TGprop,
                    nvec_per_node,sparse,paral);},
            *this
        );
    }

    SpVType_shm_csr_matrix calculateHSPotentials(double cut, TaskGroup_& TGprop, 
        boost::multi_array<ComplexType,2>& vn0) {
        return boost::apply_visitor(
            [&](auto&& a){return a.calculateHSPotentials(cut,TGprop, vn0);},
            *this
        );
    }
    */
}; 

}

}

#endif

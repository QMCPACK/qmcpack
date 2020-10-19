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

#include <fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Hamiltonians/FactorizedSparseHamiltonian.h"
#include "AFQMC/Hamiltonians/THCHamiltonian.h"
#ifdef QMC_COMPLEX
#include "AFQMC/Hamiltonians/KPFactorizedHamiltonian.h"
//#include "AFQMC/Hamiltonians/KPTHCHamiltonian.h"
#else
#include "AFQMC/Hamiltonians/RealDenseHamiltonian.h"
#include "AFQMC/Hamiltonians/RealDenseHamiltonian_v2.h"
#endif
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
  dummy_Hamiltonian(){};

  ValueType getNuclearCoulombEnergy() const
  {
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  }

  ValueType H(IndexType i, IndexType j)
  {
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  }

  ValueType H(IndexType i, IndexType j, IndexType k, IndexType l)
  {
    throw std::runtime_error("calling visitor on dummy object");
    return 0;
  }

  boost::multi::array<ValueType, 2> getH1() const { return boost::multi::array<ValueType, 2>{}; }

  boost::multi::array<SPComplexType, 1> halfRotatedHij(WALKER_TYPES type, PsiT_Matrix* Alpha, PsiT_Matrix* Beta)
  {
    throw std::runtime_error("calling visitor on dummy object");
    return boost::multi::array<ComplexType, 1>(iextensions<1u>{1});
  }

  SpCType_shm_csr_matrix halfRotatedHijkl(WALKER_TYPES type,
                                          TaskGroup_& TGWfn,
                                          PsiT_Matrix* Alpha,
                                          PsiT_Matrix* Beta,
                                          const RealType cut = 1e-6)
  {
    throw std::runtime_error("calling visitor on dummy object");
    using Alloc = shared_allocator<SPComplexType>;
    return SpCType_shm_csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, Alloc(TGWfn.Node()));
  }

  SpVType_shm_csr_matrix calculateHSPotentials(double cut, TaskGroup_& TGprop, boost::multi::array<ComplexType, 2>& vn0)
  {
    throw std::runtime_error("calling visitor on dummy object");
    using Alloc = shared_allocator<SPComplexType>;
    return SpVType_shm_csr_matrix(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, Alloc(TGprop.Node()));
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
} // namespace dummy

#ifdef QMC_COMPLEX
class Hamiltonian
    : public boost::
          variant<dummy::dummy_Hamiltonian, FactorizedSparseHamiltonian, THCHamiltonian, KPFactorizedHamiltonian
                  //                                       ,KPTHCHamiltonian
                  >
#else
class Hamiltonian : public boost::variant<dummy::dummy_Hamiltonian,
                                          FactorizedSparseHamiltonian,
                                          THCHamiltonian,
                                          RealDenseHamiltonian,
                                          RealDenseHamiltonian_v2>
#endif
{
  using shm_csr_matrix = FactorizedSparseHamiltonian::shm_csr_matrix;
  using csr_matrix_vew = FactorizedSparseHamiltonian::csr_matrix_view;

public:
  Hamiltonian() { APP_ABORT(" Error: Reached default constructor of Hamiltonian. \n"); }
  explicit Hamiltonian(THCHamiltonian&& other) : variant(std::move(other)) {}
  explicit Hamiltonian(FactorizedSparseHamiltonian&& other) : variant(std::move(other)) {}
#ifdef QMC_COMPLEX
  explicit Hamiltonian(KPFactorizedHamiltonian&& other) : variant(std::move(other)) {}
//    explicit Hamiltonian(KPTHCHamiltonian&& other) : variant(std::move(other)) {}
#else
  explicit Hamiltonian(RealDenseHamiltonian&& other) : variant(std::move(other)) {}
  explicit Hamiltonian(RealDenseHamiltonian_v2&& other) : variant(std::move(other)) {}
#endif

  explicit Hamiltonian(THCHamiltonian const& other)              = delete;
  explicit Hamiltonian(FactorizedSparseHamiltonian const& other) = delete;
#ifdef QMC_COMPLEX
  explicit Hamiltonian(KPFactorizedHamiltonian const& other) = delete;
//    explicit Hamiltonian(KPTHCHamiltonian const& other) = delete;
#else
  explicit Hamiltonian(RealDenseHamiltonian const& other)    = delete;
  explicit Hamiltonian(RealDenseHamiltonian_v2 const& other) = delete;
#endif

  Hamiltonian(Hamiltonian const& other) = delete;
  Hamiltonian(Hamiltonian&& other)      = default;

  Hamiltonian& operator=(Hamiltonian const& other) = delete;
  Hamiltonian& operator=(Hamiltonian&& other) = default;

  ValueType getNuclearCoulombEnergy()
  {
    return boost::apply_visitor([&](auto&& a) { return a.getNuclearCoulombEnergy(); }, *this);
  }

  boost::multi::array<ValueType, 2> getH1() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.getH1(); }, *this);
  }

  ValueType H(IndexType i, IndexType j)
  {
    return boost::apply_visitor([&](auto&& a) { return a.H(i, j); }, *this);
  }

  ValueType H(IndexType i, IndexType j, IndexType k, IndexType l)
  {
    return boost::apply_visitor([&](auto&& a) { return a.H(i, j, k, l); }, *this);
  }

  template<class... Args>
  HamiltonianOperations getHamiltonianOperations(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.getHamiltonianOperations(std::forward<Args>(args)...); },
                                *this);
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif

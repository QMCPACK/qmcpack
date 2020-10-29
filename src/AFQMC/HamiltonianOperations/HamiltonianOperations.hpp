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

#ifndef QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_HPP
#define QMCPLUSPLUS_AFQMC_HAMILTONIANOPERATIONS_HPP

#include <fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/HamiltonianOperations/SparseTensor.hpp"
#include "AFQMC/HamiltonianOperations/THCOps.hpp"
#ifdef QMC_COMPLEX
#include "AFQMC/HamiltonianOperations/KP3IndexFactorization.hpp"
#include "AFQMC/HamiltonianOperations/KP3IndexFactorization_batched.hpp"
//#include "AFQMC/HamiltonianOperations/KPTHCOps.hpp"
#else
#include "AFQMC/HamiltonianOperations/Real3IndexFactorization.hpp"
//#include "AFQMC/HamiltonianOperations/Real3IndexFactorization_batched.hpp"
#include "AFQMC/HamiltonianOperations/Real3IndexFactorization_batched_v2.hpp"
#endif

namespace qmcplusplus
{
namespace afqmc
{
namespace dummy
{
/*
 * Empty class to avoid need for default constructed HamiltonianOperations.
 * Throws is any visitor is called.
 */
class dummy_HOps
{
public:
  dummy_HOps(){};

  template<class... Args>
  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return boost::multi::array<ComplexType, 2>{};
  }

  template<class... Args>
  void energy(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class... Args>
  void generalizedFockMatrix(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class... Args>
  void fast_energy(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  bool fast_ph_energy() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  template<class... Args>
  void vHS(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class... Args>
  void vbias(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  template<class... Args>
  void write2hdf(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
  }

  int number_of_ke_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }

  int local_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }

  int global_origin_cholesky_vector() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }

  int global_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return 0;
  }

  bool transposed_G_for_vbias() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool transposed_G_for_E() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool transposed_vHS() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  bool distribution_over_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return false;
  }

  boost::multi::array<ComplexType, 2> getHSPotentials()
  {
    throw std::runtime_error("calling visitor on dummy_HOps object");
    return boost::multi::array<ComplexType, 2>{};
  }

  HamiltonianTypes getHamType() const { return UNKNOWN; }
};

} // namespace dummy

using devMatrix = SPComplexMatrix<device_allocator<SPComplexType>>;
using shmMatrix = SPComplexMatrix<shared_allocator<SPComplexType>>;

#ifdef QMC_COMPLEX
class HamiltonianOperations : public boost::variant<dummy::dummy_HOps,
                                                    THCOps,
                                                    SparseTensor<ComplexType, ComplexType>,
                                                    KP3IndexFactorization,
                                                    KP3IndexFactorization_batched<devMatrix>,
                                                    KP3IndexFactorization_batched<shmMatrix>
                                                    //                              ,KPTHCOps
                                                    >
#else
class HamiltonianOperations : public boost::variant<dummy::dummy_HOps,
                                                    THCOps,
                                                    SparseTensor<RealType, RealType>,
                                                    SparseTensor<RealType, ComplexType>,
                                                    SparseTensor<ComplexType, RealType>,
                                                    SparseTensor<ComplexType, ComplexType>,
                                                    Real3IndexFactorization,
                                                    //                                  Real3IndexFactorization_batched,
                                                    Real3IndexFactorization_batched_v2>
#endif
{
#ifndef QMC_COMPLEX
  using STRR = SparseTensor<RealType, RealType>;
  using STRC = SparseTensor<RealType, ComplexType>;
  using STCR = SparseTensor<ComplexType, RealType>;
#endif
  using STCC = SparseTensor<ComplexType, ComplexType>;

public:
  HamiltonianOperations() : variant() {}
#ifndef QMC_COMPLEX
  explicit HamiltonianOperations(STRR&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(STRC&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(STCR&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(Real3IndexFactorization&& other) : variant(std::move(other)) {}
  //    explicit HamiltonianOperations(Real3IndexFactorization_batched&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(Real3IndexFactorization_batched_v2&& other) : variant(std::move(other)) {}
#else
  explicit HamiltonianOperations(KP3IndexFactorization&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(KP3IndexFactorization_batched<devMatrix>&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(KP3IndexFactorization_batched<shmMatrix>&& other) : variant(std::move(other)) {}
//    explicit HamiltonianOperations(KPTHCOps&& other) : variant(std::move(other)) {}
#endif
  explicit HamiltonianOperations(STCC&& other) : variant(std::move(other)) {}
  explicit HamiltonianOperations(THCOps&& other) : variant(std::move(other)) {}

#ifndef QMC_COMPLEX
  explicit HamiltonianOperations(STRR const& other)                    = delete;
  explicit HamiltonianOperations(STRC const& other)                    = delete;
  explicit HamiltonianOperations(STCR const& other)                    = delete;
  explicit HamiltonianOperations(Real3IndexFactorization const& other) = delete;
  //    explicit HamiltonianOperations(Real3IndexFactorization_batched const& other) = delete;
  explicit HamiltonianOperations(Real3IndexFactorization_batched_v2 const& other) = delete;
#else
  explicit HamiltonianOperations(KP3IndexFactorization const& other)                    = delete;
  explicit HamiltonianOperations(KP3IndexFactorization_batched<devMatrix> const& other) = delete;
  explicit HamiltonianOperations(KP3IndexFactorization_batched<shmMatrix> const& other) = delete;
//    explicit HamiltonianOperations(KPTHCOps const& other) = delete;
#endif
  explicit HamiltonianOperations(STCC const& other)   = delete;
  explicit HamiltonianOperations(THCOps const& other) = delete;

  HamiltonianOperations(HamiltonianOperations const& other) = delete;
  HamiltonianOperations(HamiltonianOperations&& other)      = default;

  HamiltonianOperations& operator=(HamiltonianOperations const& other) = delete;
  HamiltonianOperations& operator=(HamiltonianOperations&& other) = default;

  template<class... Args>
  boost::multi::array<ComplexType, 2> getOneBodyPropagatorMatrix(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.getOneBodyPropagatorMatrix(std::forward<Args>(args)...); },
                                *this);
  }

  template<class... Args>
  void write2hdf(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.write2hdf(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void energy(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.energy(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void fast_energy(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.fast_energy(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void generalizedFockMatrix(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.generalizedFockMatrix(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void vHS(Args&&... args)
  {
    boost::apply_visitor([&](auto&& s) { s.vHS(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void vbias(Args&&... args)
  {
    boost::apply_visitor([&](auto&& s) { s.vbias(std::forward<Args>(args)...); }, *this);
  }

  /*
    template<class MatA, class MatB>
    void vbias(const MatA& G, MatB&& v, double a=1., double c=0., int nd=0) {
        boost::apply_visitor(
            [&](auto&& s){
                if constexpr(typename std::is_convertible<std::>::value)
                    s.vbias(G,std::forward<MatB>(v),a,c,nd);
                else
                    throw 0;
            },
            *this
        );
    }
*/

  int local_number_of_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.local_number_of_cholesky_vectors(); }, *this);
  }

  int global_origin_cholesky_vector() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.global_origin_cholesky_vector(); }, *this);
  }

  int number_of_ke_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.number_of_ke_vectors(); }, *this);
  }

  int global_number_of_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.global_number_of_cholesky_vectors(); }, *this);
  }

  bool distribution_over_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.distribution_over_cholesky_vectors(); }, *this);
  }


  bool transposed_G_for_vbias() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.transposed_G_for_vbias(); }, *this);
  }

  bool transposed_G_for_E() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.transposed_G_for_E(); }, *this);
  }

  bool transposed_vHS() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.transposed_vHS(); }, *this);
  }

  bool fast_ph_energy() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.fast_ph_energy(); }, *this);
  }

  HamiltonianTypes getHamType() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.getHamType(); }, *this);
  }

  template<class... Args>
  boost::multi::array<ComplexType, 2> getHSPotentials()
  {
    return boost::apply_visitor([&](auto&& a) { return a.getHSPotentials(); }, *this);
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif

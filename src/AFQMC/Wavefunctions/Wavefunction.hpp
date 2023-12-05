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

#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTION_HPP
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTION_HPP

#include <fstream>

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"
#include "AFQMC/Wavefunctions/NOMSD.hpp"
#include "AFQMC/Wavefunctions/PHMSD.hpp"

namespace qmcplusplus
{
namespace afqmc
{
namespace dummy
{
/*
 * Empty class to avoid need for default constructed Wavefunctions.
 * Throws is any visitor is called. 
 */
class dummy_wavefunction
{
private:
  std::vector<ComplexType> ci;
  std::vector<PsiT_Matrix> orbs;

public:
  dummy_wavefunction(){};

  int size_of_G_for_vbias() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return 0;
  }
  int local_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return 0;
  }
  int global_number_of_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return 0;
  }
  int global_origin_cholesky_vector() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return 0;
  }
  int number_of_references_for_back_propagation() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return 0;
  }
  bool distribution_over_cholesky_vectors() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false;
  }

  WALKER_TYPES getWalkerType() const { return UNDEFINED_WALKER_TYPE; }

  bool transposed_G_for_vbias() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false;
  }

  bool transposed_G_for_E() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false;
  }

  bool transposed_vHS() const
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return false;
  }

  template<class Vec>
  void vMF(Vec&& v)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class MatG, class MatA>
  void vbias(const MatG& G, MatA&& v, double a = 1.0)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class MatX, class MatA>
  void vHS(MatX&& X, MatA&& v, double a = 1.0)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WSet>
  void Energy(WSet& wset)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class Mat, class TVec>
  void Energy(const WlkSet& wset, Mat&& E, TVec&& Ov)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class MatG>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, bool compact = true, bool transpose = false)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class MatG, class TVec>
  void MixedDensityMatrix(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact = true, bool transpose = false)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class MatG, class CVec1, class CVec2, class Mat1, class Mat2>
  void WalkerAveragedDensityMatrix(const WlkSet& wset,
                                   CVec1& wgt,
                                   MatG& G,
                                   CVec2& denom,
                                   Mat1&& Ovlp,
                                   Mat2&& DMsum,
                                   bool free_projection                          = false,
                                   boost::multi::array_ref<ComplexType, 3>* Refs = nullptr,
                                   boost::multi::array<ComplexType, 2>* detR     = nullptr)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class MatG>
  void MixedDensityMatrix_for_vbias(const WlkSet& wset, MatG&& G)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class... Args>
  void DensityMatrix(Args&&... args)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet, class TVec>
  void Overlap(const WlkSet& wset, TVec&& Ov)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet>
  void Overlap(WlkSet& wset)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  template<class WlkSet>
  void Orthogonalize(WlkSet& wset, bool impSamp)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  ComplexType getReferenceWeight(int i)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
    return ComplexType(0.0, 0.0);
  }

  template<class Mat>
  void getReferencesForBackPropagation(Mat&& A)
  {
    throw std::runtime_error("calling visitor on dummy_wavefunction object");
  }

  SlaterDetOperations SDet;
  SlaterDetOperations* getSlaterDetOperations() { return std::addressof(SDet); }

  HamiltonianOperations HOps;
  HamiltonianOperations* getHamiltonianOperations() { return std::addressof(HOps); }
};
} // namespace dummy

class Wavefunction : public boost::variant<dummy::dummy_wavefunction,
                                           NOMSD<devcsr_Matrix>,
                                           NOMSD<ComplexMatrix<node_allocator<ComplexType>>>,
                                           PHMSD>
{
public:
  Wavefunction() { APP_ABORT(" Error: Reached default constructor of Wavefunction. \n"); }
  explicit Wavefunction(NOMSD<devcsr_Matrix>&& other) : variant(std::move(other)) {}
  explicit Wavefunction(NOMSD<devcsr_Matrix> const& other) = delete;

  explicit Wavefunction(NOMSD<ComplexMatrix<node_allocator<ComplexType>>>&& other) : variant(std::move(other)) {}
  explicit Wavefunction(NOMSD<ComplexMatrix<node_allocator<ComplexType>>> const& other) = delete;

  explicit Wavefunction(PHMSD&& other) : variant(std::move(other)) {}
  explicit Wavefunction(PHMSD const& other) = delete;

  Wavefunction(Wavefunction const& other) = delete;
  Wavefunction(Wavefunction&& other)      = default;

  Wavefunction& operator=(Wavefunction const& other) = delete;
  Wavefunction& operator=(Wavefunction&& other) = default;

  int size_of_G_for_vbias() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.size_of_G_for_vbias(); }, *this);
  }

  int local_number_of_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.local_number_of_cholesky_vectors(); }, *this);
  }

  int global_number_of_cholesky_vectors() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.global_number_of_cholesky_vectors(); }, *this);
  }

  int global_origin_cholesky_vector() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.global_origin_cholesky_vector(); }, *this);
  }

  int number_of_references_for_back_propagation() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.number_of_references_for_back_propagation(); }, *this);
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

  WALKER_TYPES getWalkerType() const
  {
    return boost::apply_visitor([&](auto&& a) { return a.getWalkerType(); }, *this);
  }

  template<class... Args>
  void vMF(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.vMF(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void vbias(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.vbias(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void vHS(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.vHS(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void Energy(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.Energy(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void DensityMatrix(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.DensityMatrix(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void MixedDensityMatrix(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.MixedDensityMatrix(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void MixedDensityMatrix_for_vbias(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.MixedDensityMatrix_for_vbias(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void Overlap(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.Overlap(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void Orthogonalize(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.Orthogonalize(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  ComplexType getReferenceWeight(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.getReferenceWeight(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void getReferencesForBackPropagation(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.getReferencesForBackPropagation(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void WalkerAveragedDensityMatrix(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.WalkerAveragedDensityMatrix(std::forward<Args>(args)...); }, *this);
  }

  SlaterDetOperations* getSlaterDetOperations()
  {
    return boost::apply_visitor([&](auto&& a) { return a.getSlaterDetOperations(); }, *this);
  }

  HamiltonianOperations* getHamiltonianOperations()
  {
    return boost::apply_visitor([&](auto&& a) { return a.getHamiltonianOperations(); }, *this);
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif

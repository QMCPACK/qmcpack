//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP
#define QMCPLUSPLUS_AFQMC_SLATERDETOPERATIONS_HPP

#include "AFQMC/config.h"
#include "boost/variant.hpp"

#include "AFQMC/Memory/buffer_managers.h"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations_shared.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations_serial.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class SlaterDetOperations : public boost::variant<SlaterDetOperations_shared<ComplexType>,
                                                  SlaterDetOperations_serial<ComplexType, DeviceBufferManager>>
{
public:
  SlaterDetOperations() : variant()
  {
    app_warning() << (" WARNING: Building SlaterDetOperations with default constructor. \n");
  }

  explicit SlaterDetOperations(SlaterDetOperations_shared<ComplexType>&& other) : variant(std::move(other)) {}

  explicit SlaterDetOperations(SlaterDetOperations_shared<ComplexType> const& other) = delete;

  explicit SlaterDetOperations(SlaterDetOperations_serial<ComplexType, DeviceBufferManager> const& other) = delete;
  explicit SlaterDetOperations(SlaterDetOperations_serial<ComplexType, DeviceBufferManager>&& other)
      : variant(std::move(other))
  {}

  SlaterDetOperations(SlaterDetOperations const& other) = delete;
  SlaterDetOperations(SlaterDetOperations&& other)      = default;

  SlaterDetOperations& operator=(SlaterDetOperations const& other) = delete;
  SlaterDetOperations& operator=(SlaterDetOperations&& other) = default;

  // member functions visible outside the variant
  template<class... Args>
  ComplexType MixedDensityMatrix(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.MixedDensityMatrix(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BatchedMixedDensityMatrix(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BatchedMixedDensityMatrix(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BatchedDensityMatrices(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BatchedDensityMatrices(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BatchedOverlap(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BatchedOverlap(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BatchedPropagate(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BatchedPropagate(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  ComplexType MixedDensityMatrix_noHerm(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.MixedDensityMatrix_noHerm(std::forward<Args>(args)...); },
                                *this);
  }

  template<class... Args>
  ComplexType MixedDensityMatrixForWoodbury(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.MixedDensityMatrixForWoodbury(std::forward<Args>(args)...); },
                                *this);
  }

  template<class... Args>
  ComplexType MixedDensityMatrixFromConfiguration(Args&&... args)
  {
    return boost::
        apply_visitor([&](auto&& a) { return a.MixedDensityMatrixFromConfiguration(std::forward<Args>(args)...); },
                      *this);
  }

  template<class... Args>
  ComplexType Overlap(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.Overlap(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  ComplexType Overlap_noHerm(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.Overlap_noHerm(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  ComplexType OverlapForWoodbury(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.OverlapForWoodbury(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void Propagate(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.Propagate(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  ComplexType Orthogonalize(Args&&... args)
  {
    return boost::apply_visitor([&](auto&& a) { return a.Orthogonalize(std::forward<Args>(args)...); }, *this);
  }

  template<class... Args>
  void BatchedOrthogonalize(Args&&... args)
  {
    boost::apply_visitor([&](auto&& a) { a.BatchedOrthogonalize(std::forward<Args>(args)...); }, *this);
  }
};

} // namespace afqmc

} // namespace qmcplusplus

#endif

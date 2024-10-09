//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineOffload.hpp
 *
 * define classes MultiBsplineOffload
 * The evaluation functions are defined in MultiBsplineOffloadEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEOFFLOAD_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEOFFLOAD_HPP

#include "MultiBsplineBase.hpp"
#include "spline2/BsplineAllocator.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer and BsplineAllocator
 * @tparam T the precision of splines
 */
template<typename T>
class MultiBsplineOffload : public MultiBsplineBase<T>
{
private:
  using Base  = MultiBsplineBase<T>;
  using Alloc = OffloadAllocator<T>;
  ///use allocator
  BsplineAllocator<T, Alloc> myAllocator;

  typename Base::SplineType* createImpl(const Ugrid grid[3],
                                        const typename Base::BoundaryCondition bc[3],
                                        int num_splines) override
  {
    static_assert(std::is_same<T, typename Alloc::value_type>::value,
                  "MultiBsplineOffload and Alloc data types must agree!");
    if (getAlignedSize<T, Alloc::alignment>(num_splines) != num_splines)
      throw std::runtime_error("When creating the data space of MultiBsplineOffload, num_splines must be padded!\n");
    return myAllocator.allocateMultiBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], num_splines);
  }

public:
  MultiBsplineOffload() = default;

  ~MultiBsplineOffload() override
  {
    if (Base::spline_m != nullptr)
      myAllocator.destroy(Base::spline_m);
  }

  void finalize() override
  {
    if (auto* spline_m = Base::spline_m; spline_m != nullptr)
    {
      auto* coefs = spline_m->coefs;
      // attach pointers on the device to achieve deep copy
      PRAGMA_OFFLOAD("omp target map(always, to: spline_m[:1], coefs[:spline_m->coefs_size])")
      {
        spline_m->coefs = coefs;
      }
    }
  }
};

} // namespace qmcplusplus

#endif

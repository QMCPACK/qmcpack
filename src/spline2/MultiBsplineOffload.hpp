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
                                        int num_splines) override;

public:
  MultiBsplineOffload();
  ~MultiBsplineOffload() override;
  void finalize() override;
};

extern template class MultiBsplineOffload<float>;
extern template class MultiBsplineOffload<double>;
} // namespace qmcplusplus

#endif

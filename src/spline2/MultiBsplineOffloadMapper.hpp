//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEOFFLOADMAPPER_HPP

#include "MultiBsplineBase.hpp"
#include <vector>

namespace qmcplusplus
{
/** A mapper class to map host spline coeficients to devices and handle multi-walker evaluation.
 * @tparam T the precision of splines
 */
template<typename T>
class MultiBsplineOffloadMapper
{
  using HostBspline = MultiBsplineBase<T>;

  /// reference to a host spline object.
  const HostBspline& host_bsplines_;
  /// array of host coefficient pointers for all the blocks.
  std::vector<const T*> block_coefs_;

public:
  MultiBsplineOffloadMapper(const HostBspline& host_bsplines);

  ~MultiBsplineOffloadMapper();

  /// map host coefficients to devices
  virtual void mapToDevice();

  /// update device coeficients
  void updateToDevice();

  /** evaluate spline values
   * @param num_pos, number of electron positions
   * @param pos_arr, array of electron positions [num_pos, 3]
   * @param spline_v, result pointer
   * @param walker_stride, result distance between two positions
   */
  void mw_evaluate_v(int num_pos, T* pos_arr, T* spline_v, size_t walker_stride);
  /** evaluate spline value, gradients and hessian.
   * @param num_pos, number of electron positions
   * @param pos_arr, array of electron positions [num_pos, 3]
   * @param spline_vgh, result pointer
   * @param walker_stride, result distance between two positions
   * @param filed_stride, result distance of value, gradients and hessian fields for a given electron position.
   */
  void mw_evaluate_vgh(int num_pos, T* pos_arr, T* spline_vgh, size_t walker_stride, size_t field_stride);
};

extern template class MultiBsplineOffloadMapper<float>;
extern template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus

#endif

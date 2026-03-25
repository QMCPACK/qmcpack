//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBsplineBase.hpp
 *
 * define classes MultiBsplineBase
 * The evaluation functions are defined in MultiBsplineBaseEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP

#include <array>
#include <cstddef>
#include <vector>
#include "spline2/bspline_traits.hpp"
#include "spline2/MultiBsplineEval.hpp"

namespace qmcplusplus
{

/** container class to hold a 3D multi spline pointer
 * @tparam T the precision of splines
 *
 * This class contains a pointer to a C object, copy and assign of this class is forbidden.
 */
template<typename T>
class MultiBsplineBase
{
protected:
  ///define the einsplie object type
  using SplineType = typename bspline_traits<T, 3>::SplineType;
  ///define the real type
  using real_type = typename bspline_traits<T, 3>::real_type;
  ///define the boundary condition type
  using BoundaryCondition = typename bspline_traits<T, 3>::BCType;
  ///actual vector of einspline multi-bspline objects
  std::vector<SplineType*> spline_blocks;

  const std::vector<size_t> offsets_;

  MultiBsplineBase(const std::vector<size_t>& offsets) : offsets_(offsets) {}

  /** create BoundaryCondition
   * @tparam BCT boundary type
   */
  template<typename BCT>
  static auto createBoundaryCondition(const BCT& bc)
  {
    std::array<BoundaryCondition, 3> xyzBC;
    xyzBC[0].lCode = bc[0].lCode;
    xyzBC[1].lCode = bc[1].lCode;
    xyzBC[2].lCode = bc[2].lCode;
    xyzBC[0].rCode = bc[0].rCode;
    xyzBC[1].rCode = bc[1].rCode;
    xyzBC[2].rCode = bc[2].rCode;
    xyzBC[0].lVal  = static_cast<T>(bc[0].lVal);
    xyzBC[1].lVal  = static_cast<T>(bc[1].lVal);
    xyzBC[2].lVal  = static_cast<T>(bc[2].lVal);
    xyzBC[0].rVal  = static_cast<T>(bc[0].rVal);
    xyzBC[1].rVal  = static_cast<T>(bc[1].rVal);
    xyzBC[2].rVal  = static_cast<T>(bc[2].rVal);
    return xyzBC;
  }

  void setMetaData(SplineType& spline,
                   Ugrid x_grid,
                   Ugrid y_grid,
                   Ugrid z_grid,
                   const BoundaryCondition bc[3],
                   size_t num_splines,
                   size_t num_splines_padded)
  {
    auto& xBC          = bc[0];
    auto& yBC          = bc[1];
    auto& zBC          = bc[2];
    spline.spcode      = bspline_traits<T, 3>::spcode;
    spline.tcode       = bspline_traits<T, 3>::tcode;
    spline.xBC         = xBC;
    spline.yBC         = yBC;
    spline.zBC         = zBC;
    spline.num_splines = num_splines;

    // Setup internal variables
    int Mx = x_grid.num;
    int My = y_grid.num;
    int Mz = z_grid.num;
    int Nx, Ny, Nz;

    if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
      Nx = Mx + 3;
    else
      Nx = Mx + 2;
    x_grid.delta     = (x_grid.end - x_grid.start) / (double)(Nx - 3);
    x_grid.delta_inv = 1.0 / x_grid.delta;
    spline.x_grid    = x_grid;

    if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
      Ny = My + 3;
    else
      Ny = My + 2;
    y_grid.delta     = (y_grid.end - y_grid.start) / (double)(Ny - 3);
    y_grid.delta_inv = 1.0 / y_grid.delta;
    spline.y_grid    = y_grid;

    if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
      Nz = Mz + 3;
    else
      Nz = Mz + 2;
    z_grid.delta     = (z_grid.end - z_grid.start) / (double)(Nz - 3);
    z_grid.delta_inv = 1.0 / z_grid.delta;
    spline.z_grid    = z_grid;

    spline.x_stride = (size_t)Ny * (size_t)Nz * num_splines_padded;
    spline.y_stride = Nz * num_splines_padded;
    spline.z_stride = num_splines_padded;

    spline.coefs_size = (size_t)Nx * spline.x_stride;
  }

public:
  MultiBsplineBase()                                      = default;
  MultiBsplineBase(const MultiBsplineBase& in)            = delete;
  MultiBsplineBase& operator=(const MultiBsplineBase& in) = delete;

  virtual ~MultiBsplineBase() = default;

  SplineType* getSplinePtr()
  {
    if (spline_blocks.size() != 1)
      throw std::runtime_error("Bug! Cannot access splime_m. the number of spline_blocks is not 1.");
    return spline_blocks[0];
  }

  void flush_zero() const
  {
    for (auto spline_m : spline_blocks)
      std::fill(spline_m->coefs, spline_m->coefs + spline_m->coefs_size, T(0));
  }

  size_t num_splines() const
  {
    size_t num_splines = 0;
    for (auto spline_m : spline_blocks)
      num_splines += spline_m->num_splines;
    return num_splines;
  }

  size_t num_splines_padded() const
  {
    size_t num_splines_padded = 0;
    for (auto spline_m : spline_blocks)
      num_splines_padded += spline_m->z_stride;
    return num_splines_padded;
  }

  size_t sizeInByte() const
  {
    size_t num_T = 0;
    for (auto spline_m : spline_blocks)
      num_T += spline_m->coefs_size * sizeof(T);
    return num_T * sizeof(T);
  }

  /** copy a single spline to multi spline
   * @param single UBspline_3d_d
   * @param int index of single in multi
   */
  void set_spline(const typename bspline_traits<double, 3>::SingleSplineType& single, int i)
  {
    size_t iblock = 0;
    while (iblock < spline_blocks.size() && i >= offsets_[iblock + 1])
      iblock++;
    if (iblock == spline_blocks.size())
      throw std::runtime_error("Bug detected in MultiBsplineBase::set_spline i goes out of bound!");

    auto& multi(*spline_blocks[iblock]);

    if (single.x_grid.num != multi.x_grid.num || single.y_grid.num != multi.y_grid.num ||
        single.z_grid.num != multi.z_grid.num)
      throw std::runtime_error("Cannot copy a single spline to MultiSpline with a different grid!\n");

    intptr_t x_stride_in  = single.x_stride;
    intptr_t y_stride_in  = single.y_stride;
    intptr_t x_stride_out = multi.x_stride;
    intptr_t y_stride_out = multi.y_stride;
    intptr_t z_stride_out = multi.z_stride;
    const intptr_t istart = static_cast<intptr_t>(i - offsets_[iblock]);
    const intptr_t n0 = multi.x_grid.num + 3, n1 = multi.y_grid.num + 3, n2 = multi.z_grid.num + 3;
    for (intptr_t ix = 0; ix < n0; ++ix)
      for (intptr_t iy = 0; iy < n1; ++iy)
      {
        auto* __restrict__ out      = multi.coefs + ix * x_stride_out + iy * y_stride_out + istart;
        const auto* __restrict__ in = single.coefs + ix * x_stride_in + iy * y_stride_in;
        for (intptr_t iz = 0; iz < n2; ++iz)
          out[iz * z_stride_out] = in[iz];
      }
  }

  virtual void finalize() {};

  template<typename PT, typename VT>
  inline void evaluate_v(const PT& r, VT& psi)
  {
    for (size_t ib = 0; ib < spline_blocks.size(); ib++)
    {
      const auto* spline_m(spline_blocks[ib]);
      if (spline_m->num_splines == 0)
        continue;
      spline2::evaluate_v_impl(spline_m, r[0], r[1], r[2], psi.data() + offsets_[ib], 0, spline_m->num_splines);
    }
  }

  template<typename PT, typename VT, typename GT>
  inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, GT& lap)
  {
    for (size_t ib = 0; ib < spline_blocks.size(); ib++)
    {
      const auto* spline_m(spline_blocks[ib]);
      if (spline_m->num_splines == 0)
        continue;
      spline2::evaluate_vgl_impl(spline_m, r[0], r[1], r[2], psi.data() + offsets_[ib], grad.data() + offsets_[ib],
                                 lap.data() + offsets_[ib], psi.size(), 0, spline_m->num_splines);
    }
  }

  template<typename PT, typename VT, typename GT, typename HT>
  inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
  {
    for (size_t ib = 0; ib < spline_blocks.size(); ib++)
    {
      const auto* spline_m(spline_blocks[ib]);
      if (spline_m->num_splines == 0)
        continue;
      spline2::evaluate_vgh_impl(spline_m, r[0], r[1], r[2], psi.data() + offsets_[ib], grad.data() + offsets_[ib],
                                 hess.data() + offsets_[ib], psi.size(), 0, spline_m->num_splines);
    }
  }

  template<typename PT, typename VT, typename GT, typename HT, typename GHT>
  inline void evaluate_vghgh(const PT& r, VT& psi, GT& grad, HT& hess, GHT& ghess)
  {
    for (size_t ib = 0; ib < spline_blocks.size(); ib++)
    {
      const auto* spline_m(spline_blocks[ib]);
      if (spline_m->num_splines == 0)
        continue;
      spline2::evaluate_vghgh_impl(spline_m, r[0], r[1], r[2], psi.data() + offsets_[ib], grad.data() + offsets_[ib],
                                   hess.data() + offsets_[ib], ghess.data() + offsets_[ib], psi.size(), 0,
                                   spline_m->num_splines);
    }
  }
};

} // namespace qmcplusplus

#endif

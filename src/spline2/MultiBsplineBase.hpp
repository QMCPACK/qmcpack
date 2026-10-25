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
 * @brief Defines the MultiBsplineBase class for 3D multi-B-spline operations.
 *
 * Provides a C++ object-oriented interface around the einspline C library,
 * managing metadata, blocks of splines, and dispatching precision-mixed evaluations.
 * The inner evaluation functions are defined in MultiBsplineEval.hpp and
 * MultiBsplineValue.hpp, MultiBsplineVGLH.hpp, etc.
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEBASE_HPP

#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsSoA/VectorSoaContainer.h>
#include <CPU/SIMD/aligned_allocator.hpp>
#include "spline2/bspline_traits.hpp"

namespace qmcplusplus
{

/** @ingroup spline2
 * @brief Base container class for 3D multi-B-spline evaluation.
 *
 * MultiBsplineBase manages a collection (blocks) of 3D multi-B-spline objects
 * implemented in C (from the einspline library). It handles memory layout metadata,
 * boundary conditions, and dispatches evaluation calls (value, gradient, laplacian,
 * hessian) to optimized SIMD routines.
 *
 * @tparam T The storage precision of the internal spline coefficients (float or double).
 *           Note that spline evaluations can be performed in a different precision (VT).
 *
 * Copying and assignment are explicitly deleted to prevent unsafe aliasing of the
 * underlying C-style pointers.
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
  ///index offsets of spline_blocks.
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

  /** Configure internal einspline metadata including grid and stride information
   * @param spline reference to the internal spline object to configure
   * @param x_grid grid in x direction
   * @param y_grid grid in y direction
   * @param z_grid grid in z direction
   * @param bc array of boundary conditions
   * @param num_splines number of valid splines to store
   * @param num_splines_padded number of splines padded for SIMD alignment
   */
  void setMetaData(SplineType& spline,
                   Ugrid x_grid,
                   Ugrid y_grid,
                   Ugrid z_grid,
                   const BoundaryCondition bc[3],
                   size_t num_splines,
                   size_t num_splines_padded);

public:
  MultiBsplineBase(const MultiBsplineBase& in)            = delete;
  MultiBsplineBase& operator=(const MultiBsplineBase& in) = delete;

  virtual ~MultiBsplineBase() = default;

  /** Return the number of spline blocks currently stored
   * @return number of spline blocks
   */
  size_t getNumBlocks() const { return spline_blocks.size(); }
  /** Get the offsets array mapping global spline index to local block index
   * @return constant reference to offsets vector
   */
  const auto& getBlockOffsets() const { return offsets_; }

  /** Return pointer to the primary spline object (assumes a single block)
   * @return pointer to SplineType
   */
  SplineType* getSplinePtr();

  /** Access a specific spline block by index
   * @param iblock index of the block
   * @return reference to the requested spline block
   */
  SplineType& getBlock(size_t iblock) { return *spline_blocks[iblock]; }
  /** Const access to a specific spline block by index
   * @param iblock index of the block
   * @return const reference to the requested spline block
   */
  const SplineType& getBlock(size_t iblock) const { return *spline_blocks[iblock]; }

  /** Zero out the spline coefficients for a specific block
   * @param iblock block index to zero out
   */
  void flush_zero(size_t iblock = 0) const;

  /** Get the total number of splines across all blocks
   * @return total number of splines
   */
  size_t num_splines() const;

  /** Get the total padded number of splines (includes SIMD padding)
   * @return total padded number of splines
   */
  size_t num_splines_padded() const;

  /** Get the memory size occupied by the spline coefficients in bytes
   * @return size in bytes
   */
  size_t sizeInByte() const;

  /** Copy a single spline into the multi-spline data structure
   * @param single source UBspline_3d_d object
   * @param i destination index within the multi-spline structure
   */
  void setOneSpline(const UBspline_3d_d& single, int i);


  /** Evaluate spline values into single-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   */
  void evaluate_v(const TinyVector<float, 3>& r, Vector<float, aligned_allocator<float>>& psi);

  /** Evaluate spline values into double-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   */
  void evaluate_v(const TinyVector<double, 3>& r, Vector<double, aligned_allocator<double>>& psi);

  /** Evaluate spline values, gradients, and laplacians into single-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param lap output container for spline laplacians
   */
  void evaluate_vgl(const TinyVector<float, 3>& r,
                    Vector<float, aligned_allocator<float>>& psi,
                    VectorSoaContainer<float, 3>& grad,
                    VectorSoaContainer<float, 3>& lap);

  /** Evaluate spline values, gradients, and laplacians into double-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param lap output container for spline laplacians
   */
  void evaluate_vgl(const TinyVector<double, 3>& r,
                    Vector<double, aligned_allocator<double>>& psi,
                    VectorSoaContainer<double, 3>& grad,
                    VectorSoaContainer<double, 3>& lap);

  /** Evaluate spline values, gradients, and hessians into single-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   */
  void evaluate_vgh(const TinyVector<float, 3>& r,
                    Vector<float, aligned_allocator<float>>& psi,
                    VectorSoaContainer<float, 3>& grad,
                    VectorSoaContainer<float, 6>& hess);

  /** Evaluate spline values, gradients, and hessians into double-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   */
  void evaluate_vgh(const TinyVector<double, 3>& r,
                    Vector<double, aligned_allocator<double>>& psi,
                    VectorSoaContainer<double, 3>& grad,
                    VectorSoaContainer<double, 6>& hess);


  /** Evaluate spline values, gradients, hessians, and gradient-hessians into single-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   * @param ghess output container for spline gradient-hessians
   */
  void evaluate_vghgh(const TinyVector<float, 3>& r,
                      Vector<float, aligned_allocator<float>>& psi,
                      VectorSoaContainer<float, 3>& grad,
                      VectorSoaContainer<float, 6>& hess,
                      VectorSoaContainer<float, 10>& ghess);

  /** Evaluate spline values, gradients, hessians, and gradient-hessians into double-precision arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   * @param ghess output container for spline gradient-hessians
   */
  void evaluate_vghgh(const TinyVector<double, 3>& r,
                      Vector<double, aligned_allocator<double>>& psi,
                      VectorSoaContainer<double, 3>& grad,
                      VectorSoaContainer<double, 6>& hess,
                      VectorSoaContainer<double, 10>& ghess);

private:
  /** Evaluate the spline values at a given 3D position (implementation)
   * @tparam VT the precision of the output evaluation arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   */
  template<typename VT>
  void evaluate_v_impl(const TinyVector<VT, 3>& r, Vector<VT, aligned_allocator<VT>>& psi);

  /** Evaluate the spline values, gradients, and laplacians (implementation)
   * @tparam VT the precision of the output evaluation arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param lap output container for spline laplacians
   */
  template<typename VT>
  void evaluate_vgl_impl(const TinyVector<VT, 3>& r,
                         Vector<VT, aligned_allocator<VT>>& psi,
                         VectorSoaContainer<VT, 3>& grad,
                         VectorSoaContainer<VT, 3>& lap);

  /** Evaluate the spline values, gradients, and hessians (implementation)
   * @tparam VT the precision of the output evaluation arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   */
  template<typename VT>
  void evaluate_vgh_impl(const TinyVector<VT, 3>& r,
                         Vector<VT, aligned_allocator<VT>>& psi,
                         VectorSoaContainer<VT, 3>& grad,
                         VectorSoaContainer<VT, 6>& hess);

  /** Evaluate the spline values, gradients, hessians, and gradient-hessians (implementation)
   * @tparam VT the precision of the output evaluation arrays
   * @param r 3D position vector
   * @param psi output vector for spline values
   * @param grad output container for spline gradients
   * @param hess output container for spline hessians
   * @param ghess output container for spline gradient-hessians
   */
  template<typename VT>
  void evaluate_vghgh_impl(const TinyVector<VT, 3>& r,
                           Vector<VT, aligned_allocator<VT>>& psi,
                           VectorSoaContainer<VT, 3>& grad,
                           VectorSoaContainer<VT, 6>& hess,
                           VectorSoaContainer<VT, 10>& ghess);
};


extern template class MultiBsplineBase<float>;
extern template class MultiBsplineBase<double>;

} // namespace qmcplusplus

#endif

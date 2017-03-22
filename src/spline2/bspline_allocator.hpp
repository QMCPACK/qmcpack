//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file bspline_allocator.hpp
 * @brief Allocator and management classes
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H
#define QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H

#include <simd/allocator.hpp>
#include <simd/simd.hpp>
#include <spline2/bspline_traits.hpp>
#include "spline2/einspline_allocator.h"

namespace qmcplusplus { namespace einspline {

  class Allocator
  {
    ///Setting the allocation policy: default is using aligned allocator
    int Policy;
    public:
    ///constructor
    Allocator();
#if (__cplusplus >= 201103L)
    ///disable copy constructor
    Allocator(const Allocator&)=delete;
    ///disable assignement
    Allocator& operator=(const Allocator&)=delete;
#endif
    ///destructor
    ~Allocator();

    template<typename SplineType>
    void destroy(SplineType* spline)
    {
      einspline_free(spline->coefs);
      free(spline);
    }

    ///allocate a single multi-bspline
    multi_UBspline_3d_s* allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, int num_splines);

    ///allocate a double multi-bspline
    multi_UBspline_3d_d* allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, int num_splines);

    ///allocate a single bspline
    UBspline_3d_s* allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, float* data=nullptr);

    ///allocate a UBspline_3d_d
    UBspline_3d_d* allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, double* data=nullptr);

    /** allocate a multi_UBspline_3d_(s,d)
     * @tparam T datatype
     * @tparam ValT 3D container for start and end
     * @tparam IntT 3D container for ng
     */
    template<typename T, typename ValT, typename IntT>
      typename bspline_traits<T,3>::SplineType*
      createMultiBspline(T dummy, ValT& start , ValT& end, IntT& ng, bc_code bc, int num_splines);

    /** allocate a UBspline_3d_(s,d)
     * @tparam T datatype
     * @tparam ValT 3D container for start and end
     * @tparam IntT 3D container for ng
     */
    template< typename ValT, typename IntT, typename T>
      typename bspline_traits<T,3>::SingleSplineType*
      createUBspline(ValT& start , ValT& end, IntT& ng , bc_code bc,T* data=nullptr);

    /** set the data to a spline, interpolation is done
     * @param indata starting address of the input data
     * @param spline target MultiBsplineType
     * @param i the band index to copy to
     */
    void set(float* indata, multi_UBspline_3d_s* spline, int i);
    void set(double* indata, multi_UBspline_3d_d* spline, int i);
    /** set the data in double to multi_UBspline_3d_s */
    void set(double* indata, multi_UBspline_3d_s* spline, int i);

    /** copy a UBSpline_3d_X to multi_UBspline_3d_X at i-th band
     * @param single  UBspline_3d_X
     * @param multi target multi_UBspline_3d_X
     * @param i the band index to copy to
     * @param offset starting offset for AoSoA
     * @param N shape of AoSoA
     */
    template<typename UBT, typename MBT>
      void copy(UBT* single, MBT* multi, int i,  const int* offset, const int* N);

    /** copy double to single: only for testing */
    void copy(multi_UBspline_3d_d* in, multi_UBspline_3d_s* out);
  };

  template<typename T, typename ValT, typename IntT>
    typename bspline_traits<T,3>::SplineType*
    Allocator::createMultiBspline(T dummy, ValT& start , ValT& end, IntT& ng , bc_code bc, int num_splines)
    { 
      Ugrid x_grid, y_grid, z_grid;
      typename bspline_traits<T,3>::BCType xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return allocateMultiBspline(x_grid,y_grid,z_grid, xBC, yBC, zBC, num_splines);
    }

  template< typename ValT, typename IntT, typename T>
    typename bspline_traits<T,3>::SingleSplineType*
    Allocator::createUBspline(ValT& start , ValT& end, IntT& ng , bc_code bc, T* data)
    { 
      Ugrid x_grid, y_grid, z_grid;
      typename bspline_traits<T,3>::BCType xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return allocateUBspline(x_grid,y_grid,z_grid, xBC, yBC, zBC,data);
    }

  template<typename UBT, typename MBT>
    void Allocator::copy(UBT* single, MBT* multi, int i,  const int* offset, const int* N)
    {
      typedef typename bspline_type<MBT>::value_type out_type;
      typedef typename bspline_type<UBT>::value_type in_type;
      intptr_t x_stride_in=single->x_stride;
      intptr_t y_stride_in=single->y_stride;
      intptr_t x_stride_out=multi->x_stride;
      intptr_t y_stride_out=multi->y_stride;
      intptr_t z_stride_out=multi->z_stride;
      intptr_t offset0=static_cast<intptr_t>(offset[0]);
      intptr_t offset1=static_cast<intptr_t>(offset[1]);
      intptr_t offset2=static_cast<intptr_t>(offset[2]);
      const intptr_t istart=static_cast<intptr_t>(i);
      const intptr_t n0=N[0],n1=N[1],n2=N[2];
      for(intptr_t ix=0; ix<n0; ++ix)
        for(intptr_t iy=0; iy<n1; ++iy)
        {
          out_type* restrict out=multi->coefs+ix*x_stride_out+iy*y_stride_out+istart;
          const in_type* restrict in =single->coefs+(ix+offset0)*x_stride_in+(iy+offset1)*y_stride_in+offset2;
          for(intptr_t iz=0; iz<n2; ++iz)
          {
            out[iz*z_stride_out]=static_cast<out_type>(in[iz]);
          }
        }
    }

  void set(multi_UBspline_1d_d* spline, int i, UBspline_1d_d* spline_in,
       const int offset, const int N);

  void set(multi_UBspline_1d_s* spline, int i, UBspline_1d_d* spline_in,
       const int offset, const int N);

  /** create spline for double */
  multi_UBspline_1d_d* create(multi_UBspline_1d_d* s, Ugrid& grid, BCtype_d& bc, int num_splines);

  /** create spline for float */
  multi_UBspline_1d_s* create(multi_UBspline_1d_s* s, Ugrid& grid, BCtype_s& bc, int num_splines);

} }
#endif

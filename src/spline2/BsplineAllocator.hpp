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
 * @brief BsplineAllocator and management classes
 */
#ifndef QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H
#define QMCPLUSPLUS_EINSPLINE_BSPLINE_ALLOCATOR_H

#include <simd/allocator.hpp>
#include <simd/simd.hpp>
#include <spline2/bspline_traits.hpp>
#include "spline2/einspline_allocator.h"
#include "simd/allocator.hpp"

namespace qmcplusplus { namespace einspline {

  class BsplineAllocator
  {
    ///Setting the allocation policy: default is using aligned allocator
    int Policy;
    public:
    ///constructor
    BsplineAllocator();
    ///disable copy constructor
    BsplineAllocator(const BsplineAllocator&)=delete;
    ///disable assignement
    BsplineAllocator& operator=(const BsplineAllocator&)=delete;
    ///destructor
    ~BsplineAllocator();

    template<typename SplineType>
    void destroy(SplineType* spline)
    {
      einspline_free(spline->coefs);
      free(spline);
    }

    ///allocate a multi-bspline structure
    template<typename PRECISION, size_t ALIGN = QMC_CLINE>
    typename bspline_traits<PRECISION,3>::SplineType*
    allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                         typename bspline_traits<PRECISION,3>::BCType xBC,
                         typename bspline_traits<PRECISION,3>::BCType yBC,
                         typename bspline_traits<PRECISION,3>::BCType zBC,
                         int num_splines);

//    ///allocate a single-bspline structure
//    template<typename PRECISION>
//    bspline_traits<PRECISION,3>::SplineType
//    UBspline_3d_s* allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
//        bspline_traits<PRECISION,3>::BCType xBC,
//        bspline_traits<PRECISION,3>::BCType yBC,
//        bspline_traits<PRECISION,3>::BCType zBC,
//        PRECISION* data=nullptr);

    ///allocate a single bspline
    UBspline_3d_s* allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_s xBC, BCtype_s yBC, BCtype_s zBC, float* data=nullptr);

    ///allocate a UBspline_3d_d
    UBspline_3d_d* allocateUBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid, 
        BCtype_d xBC, BCtype_d yBC, BCtype_d zBC, double* data=nullptr);

    /** copy a UBSpline_3d_X to multi_UBspline_3d_X at i-th band
     * @param single  UBspline_3d_X
     * @param multi target multi_UBspline_3d_X
     * @param i the band index to copy to
     * @param offset starting offset for AoSoA
     * @param N shape of AoSoA
     */
    template<typename UBT, typename MBT>
      void copy(UBT* single, MBT* multi, int i,  const int* offset, const int* N);
  };

  template<typename PRECISION, size_t ALIGN>
    typename bspline_traits<PRECISION,3>::SplineType*
    BsplineAllocator::allocateMultiBspline(Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                         typename bspline_traits<PRECISION,3>::BCType xBC,
                         typename bspline_traits<PRECISION,3>::BCType yBC,
                         typename bspline_traits<PRECISION,3>::BCType zBC,
                         int num_splines)
    {
      using real_type = typename bspline_traits<PRECISION,3>::real_type;
      using SplineType = typename bspline_traits<PRECISION,3>::SplineType;
      // Create new spline
      SplineType* restrict spline = (SplineType*) malloc (sizeof(SplineType));
      spline->spcode = bspline_traits<PRECISION,3>::spcode;
      spline->tcode  = bspline_traits<PRECISION,3>::tcode;
      spline->xBC = xBC;
      spline->yBC = yBC;
      spline->zBC = zBC;
      spline->num_splines = num_splines;

      // Setup internal variables
      int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
      int Nx, Ny, Nz;

      if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
        Nx = Mx+3;
      else
        Nx = Mx+2;
      x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
      x_grid.delta_inv = 1.0/x_grid.delta;
      spline->x_grid   = x_grid;

      if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
        Ny = My+3;
      else
        Ny = My+2;
      y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
      y_grid.delta_inv = 1.0/y_grid.delta;
      spline->y_grid   = y_grid;

      if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
        Nz = Mz+3;
      else
        Nz = Mz+2;
      z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
      z_grid.delta_inv = 1.0/z_grid.delta;
      spline->z_grid   = z_grid;

      const int N = getAlignedSize<real_type, ALIGN>(num_splines);

      spline->x_stride = (size_t)Ny*(size_t)Nz*(size_t)N;
      spline->y_stride = Nz*N;
      spline->z_stride = N;

      spline->coefs_size=(size_t)Nx*spline->x_stride;
      spline->coefs=(real_type*)einspline_alloc(sizeof(real_type)*spline->coefs_size,ALIGN);

      return spline;
    }

  template<typename UBT, typename MBT>
    void BsplineAllocator::copy(UBT* single, MBT* multi, int i,  const int* offset, const int* N)
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

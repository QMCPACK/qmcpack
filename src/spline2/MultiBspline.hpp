//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBspline.hpp
 *
 * define classes MultiBspline and MultiBspline1D
 * The evaluation functions are defined in MultiBsplineEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#include "config.h"
#include <iostream>
#include <spline2/BsplineAllocator.hpp>
#include <spline2/MultiBsplineData.hpp>
#include <stdlib.h>

namespace qmcplusplus
{

  template<typename T>
    struct MultiBspline
    {

      ///define the einsplie object type
      using spliner_type=typename bspline_traits<T,3>::SplineType;
      ///define the real type
      using real_type=typename bspline_traits<T,3>::real_type;
      ///actual einspline multi-bspline object
      spliner_type* spline_m;
      ///use allocator
      einspline::BsplineAllocator myAllocator;

      MultiBspline():spline_m(nullptr) {}
      MultiBspline(const MultiBspline& in)=delete;
      MultiBspline& operator=(const MultiBspline& in)=delete;

      ~MultiBspline()
      {
        if(spline_m!=nullptr)
          myAllocator.destroy(spline_m);
      }

      /** create the einspline as used in the builder
       */
      template<typename GT, typename BCT>
      void create(GT& grid, BCT& bc, int num_splines)
      {
        if(getAlignedSize<T>(num_splines)!=num_splines)
          throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!\n");
        if(spline_m==nullptr)
        {
          typename bspline_traits<T,3>::BCType xBC, yBC, zBC;
          xBC.lCode=bc[0].lCode; yBC.lCode=bc[1].lCode; zBC.lCode=bc[2].lCode;
          xBC.rCode=bc[0].rCode; yBC.rCode=bc[1].rCode; zBC.rCode=bc[2].rCode;
          xBC.lVal=static_cast<T>(bc[0].lVal); yBC.lVal=static_cast<T>(bc[1].lVal); zBC.lVal=static_cast<T>(bc[2].lVal);
          xBC.rVal=static_cast<T>(bc[0].rVal); yBC.rVal=static_cast<T>(bc[1].rVal); zBC.rVal=static_cast<T>(bc[2].rVal);
          spline_m=myAllocator.allocateMultiBspline<T>(grid[0],grid[1],grid[2],xBC,yBC,zBC,num_splines);
        }
      }

      void flush_zero() const
      {
        if(spline_m!=nullptr) std::fill(spline_m->coefs, spline_m->coefs+spline_m->coefs_size, T(0));
      }

      int num_splines() const
      {
        return (spline_m==nullptr)?0:spline_m->num_splines;
      }

      size_t sizeInByte() const
      {
        return (spline_m==nullptr)?0:spline_m->coefs_size*sizeof(T);
      }

      /** copy a single spline to the big table
       * @param aSpline UBspline_3d_(d,s)
       * @param int index of aSpline
       * @param offset_ starting index for the case of multiple domains
       * @param base_ number of bases
       */
      template<typename SingleSpline>
      void copy_spline(SingleSpline* aSpline,int i)
      {
        if( aSpline->x_grid.num != spline_m->x_grid.num ||
            aSpline->y_grid.num != spline_m->y_grid.num ||
            aSpline->z_grid.num != spline_m->z_grid.num)
          throw std::runtime_error("Cannot copy a single spline to MultiSpline with a different grid!\n");

        const int BaseOffset[3] = {0, 0, 0};
        const int BaseN[3] = {spline_m->x_grid.num+3, spline_m->y_grid.num+3, spline_m->z_grid.num+3};
        myAllocator.copy(aSpline,spline_m,i,BaseOffset,BaseN);
      }

    };

}/** qmcplusplus namespace */

#endif


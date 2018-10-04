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
 * Master header file to define MultiBspline and MultiBspline1D
 *
 * The evaluation functions in MultiBspline and MultiBspline1D are memory
 * bandwith (BW) bound. Optimizations towards maximizing BW usage is always
 * needed. For this reason, with SIMD, memory alignment is required in order
 * to saturate the BW. The number of splines must be a multiple of aligned
 * size. The result vectors must be in Structure-of-Array datayout with
 * their starting address correctly aligned.
 *
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#include "config.h"
#include <iostream>
#include <spline2/bspline_allocator.hpp>
#include <spline2/MultiBsplineData.hpp>
#include <stdlib.h>

namespace qmcplusplus
{
  /** compute Trace(H*G)
   *
   * gg is symmetrized as
   *  gg[0]=GG(0,0)
   *  gg[1]=GG(0,1)+GG(1,0)
   *  gg[2]=GG(0,2)+GG(2,0)
   *  gg[3]=GG(1,1)
   *  gg[4]=GG(1,2)+GG(2,1)
   *  gg[5]=GG(2,2)
   */
  template<typename T>
    inline T SymTrace(T h00, T h01, T h02, T h11, T h12, T h22, const T* restrict gg)
    {
      return h00*gg[0]+h01*gg[1]+h02*gg[2]+h11*gg[3]+h12*gg[4]+h22*gg[5];
    }

  template<typename T>
    T v_m_v(T h00, T h01, T h02, T h11, T h12, T h22, T g1x, T g1y, T g1z, T g2x, T g2y, T g2z)
    {
      return g1x*g2x*h00+g1x*g2y*h01+g1x*g2z*h02
            +g1y*g2x*h01+g1y*g2y*h11+g1y*g2z*h12
            +g1z*g2x*h02+g1z*g2y*h12+g1z*g2z*h22;
    }

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
      einspline::Allocator myAllocator;

      MultiBspline():spline_m(nullptr) {}
      MultiBspline(const MultiBspline& in)=delete;
      MultiBspline& operator=(const MultiBspline& in)=delete;

      ~MultiBspline()
      {
        if(spline_m!=nullptr)
          myAllocator.destroy(spline_m);
      }

      template<typename RV, typename IV>
      void create(RV& start, RV& end, IV& ng, bc_code bc, int num_splines)
      {
        if(getAlignedSize<T>(num_splines)!=num_splines)
          throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!");
        if(spline_m==nullptr)
          spline_m=myAllocator.createMultiBspline(T(0),start,end,ng,bc,num_splines);
      }

      /** create the einspline as used in the builder
       */
      template<typename GT, typename BCT>
      void create(GT& grid, BCT& bc, int num_splines)
      {
        if(getAlignedSize<T>(num_splines)!=num_splines)
          throw std::runtime_error("When creating the data space of MultiBspline, num_splines must be padded!");
        if(spline_m==nullptr)
        {
          typename bspline_traits<T,3>::BCType xBC, yBC, zBC;
          xBC.lCode=bc[0].lCode; yBC.lCode=bc[1].lCode; zBC.lCode=bc[2].lCode;
          xBC.rCode=bc[0].rCode; yBC.rCode=bc[1].rCode; zBC.rCode=bc[2].rCode;
          xBC.lVal=static_cast<T>(bc[0].lVal); yBC.lVal=static_cast<T>(bc[1].lVal); zBC.lVal=static_cast<T>(bc[2].lVal);
          xBC.rVal=static_cast<T>(bc[0].rVal); yBC.rVal=static_cast<T>(bc[1].rVal); zBC.rVal=static_cast<T>(bc[2].rVal);
          spline_m=myAllocator.allocateMultiBspline(grid[0],grid[1],grid[2],xBC,yBC,zBC,num_splines);
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

      template<typename CT>
      inline void set(int i, CT& data)
      {
        myAllocator.set(data.data(),spline_m,i);
      }

      /** copy a single spline to the big table
       * @param aSpline UBspline_3d_(d,s)
       * @param int index of aSpline
       * @param offset_ starting index for the case of multiple domains
       * @param base_ number of bases
       */
      template<typename SingleSpline>
      void copy_spline(SingleSpline* aSpline,int i, const int* offset_, const int* base_)
      {
        myAllocator.copy(aSpline,spline_m,i,offset_,base_);
      }

    };

  template<typename T>
    struct MultiBspline1D
    {

      ///define the einsplie object type
      using spliner_type=typename bspline_traits<T,1>::SplineType;
      ///define the real type
      using real_type=typename bspline_traits<T,1>::real_type;
      ///actual einspline multi-bspline object
      spliner_type spline_m;
      ///use allocator
      //einspline::Allocator myAllocator;

      MultiBspline1D()
      {
        spline_m.coefs=nullptr;
        spline_m.num_splines=0;
        spline_m.coefs_size=0;
      }

      /** create the einspline as used in the builder
       */
      template<typename GT, typename BCT>
      void create(GT& grid, BCT& bc, int num_splines)
      {
        if(getAlignedSize<T>(num_splines)!=num_splines)
          throw std::runtime_error("When creating the data space of MultiBspline1D, num_splines must be padded!");
        spliner_type* temp_spline;
        temp_spline=einspline::create(temp_spline, grid, bc, num_splines);
        spline_m=*temp_spline;
        free(temp_spline);
      }

      void flush_zero() const
      {
        if(spline_m.coefs!=nullptr) std::fill(spline_m.coefs, spline_m.coefs+spline_m.coefs_size, T(0));
      }

      int num_splines() const
      {
        return spline_m.num_splines;
      }

      size_t sizeInByte() const
      {
        return (spline_m.coefs==nullptr)?0:spline_m.coefs_size*sizeof(T);
      }

      /** copy a single spline to the big table
       * @param aSpline UBspline_3d_(d,s)
       * @param int index of aSpline
       * @param offset_ starting index for the case of multiple domains
       * @param base_ number of bases
       */
      template<typename SingleSpline>
      void copy_spline(SingleSpline* aSpline,int i, const int offset_, const int base_)
      {
        einspline::set(&spline_m,i,aSpline,offset_,base_);
      }

      template<typename PT, typename VT>
        inline void evaluate(const PT& r, VT& psi) const
        {
          evaluate_v_impl(r,psi.data());
        }

      template<typename PT, typename VT, typename GT, typename LT>
        inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, LT& lap) const
        {
          evaluate_vgl_impl(r,psi.data(),grad.data(),lap.data());
        }

      //template<typename PT, typename VT, typename GT, typename HT>
      //  inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
      //  {
      //    evaluate_vgh_impl(r,psi.data(),grad.data(),hess.data());
      //  }

      /// compute values only.
      void evaluate_v_impl(T r, T* restrict vals) const;
      /// compute VGL.
      void evaluate_vgl_impl(T r, T* restrict vals, T* restrict grads, T* restrict lapl) const;
    };

}/** qmcplusplus namespace */

///include evaluate_vgl/vgh_impl
#ifdef BGQPX
#include <spline2/MultiBsplineVGLH_BGQ.hpp>
#else
#include <spline2/MultiBsplineVGLH.hpp>
#endif

///include evaluate_v/vgl/vgh_impl for 1D case
#include <spline2/MultiBspline1D.hpp>

#endif


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
  /** Coordinate transform for a 3rd rank symmetric tensor representing coordinate derivatives
   *  (hence t3_contract, for contraction with vectors).
   *
   * hijk are the symmetry inequivalent tensor elements, i,j,k range from 0 to 2 for x to z.
   * (gix,giy,giz) are vectors, labelled 1, 2 and 3.  g1 is contracted with the first tensor index,
   * g2 with the second, g3 with the third. 
   *
   * This would be easier with a for loop, but I'm sticking with the convention in this section.
   */
  template<typename T>
    T t3_contract(T h000, T h001, T h002, T h011, T h012, T h022, T h111, T h112, T h122, T h222,
                  T g1x, T g1y, T g1z, T g2x, T g2y, T g2z, T g3x, T g3y, T g3z) 
    { 
      return h000*(g1x*g2x*g3x)
            +h001*(g1x*g2x*g3y+g1x*g2y*g3x+g1y*g2x*g3x)
            +h002*(g1x*g2x*g3z+g1x*g2z*g3x+g1z*g2x*g3x)
            +h011*(g1x*g2y*g3y+g1y*g2x*g3y+g1y*g2y*g3x)
            +h012*(g1x*g2y*g3z+g1x*g2z*g3y+g1y*g2x*g3z+g1y*g2z*g3x+g1z*g2x*g3y+g1z*g2y*g3x)
            +h022*(g1x*g2z*g3z+g1z*g2x*g3z+g1z*g2z*g3x)
            +h111*(g1y*g2y*g3y)
            +h112*(g1y*g2y*g3z+g1y*g2z*g3y+g1z*g2y*g3y)
            +h122*(g1y*g2z*g3z+g1z*g2y*g3z+g1z*g2z*g3y)
            +h222*(g1z*g2z*g3z);  
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

///include evaluate_v/vgl/vgh_impl for 1D case
#include <spline2/MultiBspline1D.hpp>

#endif


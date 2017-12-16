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
 * Master header file to define MultiBspline
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#define QMCPLUSPLUS_MULTIEINSPLINE_COMMON_HPP
#include "config.h"
#include <iostream>
#include <spline2/bspline_allocator.hpp>
#include <stdlib.h>

namespace qmcplusplus
{
  template<typename T> struct SplineBound {};

  template<> struct SplineBound <double>
  {
    static inline void get(double x, double& dx, int& ind, int ng)
    {
      double ipart;
      dx=modf(x,&ipart);
      ind = std::min(std::max(int(0),static_cast<int>(ipart)),ng);
    }
  };

  template<> struct SplineBound <float>
  {
    static inline void get(float x, float& dx, int& ind, int ng)
    {
      float ipart;
      dx=modff(x,&ipart);
      ind = std::min(std::max(int(0),static_cast<int>(ipart)),ng);
    }
  };

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
    struct MultiBsplineData
    {
      static const T   A44[16];
      static const T  dA44[16];
      static const T d2A44[16];
      static const T d3A44[16];

      inline void compute_prefactors(T a[4], T tx) const
      {
        a[0] = ( ( A44[0]  * tx + A44[1] ) * tx + A44[2] ) * tx + A44[3];
        a[1] = ( ( A44[4]  * tx + A44[5] ) * tx + A44[6] ) * tx + A44[7];
        a[2] = ( ( A44[8]  * tx + A44[9] ) * tx + A44[10] ) * tx + A44[11];
        a[3] = ( ( A44[12] * tx + A44[13] ) * tx + A44[14] ) * tx + A44[15];
      }

      inline void compute_prefactors(T a[4], T da[4], T d2a[4], T tx) const
      {
        a[0] = ( ( A44[0]  * tx + A44[1] ) * tx + A44[2] ) * tx + A44[3];
        a[1] = ( ( A44[4]  * tx + A44[5] ) * tx + A44[6] ) * tx + A44[7];
        a[2] = ( ( A44[8]  * tx + A44[9] ) * tx + A44[10] ) * tx + A44[11];
        a[3] = ( ( A44[12] * tx + A44[13] ) * tx + A44[14] ) * tx + A44[15];
        da[0] = ( ( dA44[0]  * tx + dA44[1] ) * tx + dA44[2] ) * tx + dA44[3];
        da[1] = ( ( dA44[4]  * tx + dA44[5] ) * tx + dA44[6] ) * tx + dA44[7];
        da[2] = ( ( dA44[8]  * tx + dA44[9] ) * tx + dA44[10] ) * tx + dA44[11];
        da[3] = ( ( dA44[12] * tx + dA44[13] ) * tx + dA44[14] ) * tx + dA44[15];
        d2a[0] = ( ( d2A44[0]  * tx + d2A44[1] ) * tx + d2A44[2] ) * tx + d2A44[3];
        d2a[1] = ( ( d2A44[4]  * tx + d2A44[5] ) * tx + d2A44[6] ) * tx + d2A44[7];
        d2a[2] = ( ( d2A44[8]  * tx + d2A44[9] ) * tx + d2A44[10] ) * tx + d2A44[11];
        d2a[3] = ( ( d2A44[12] * tx + d2A44[13] ) * tx + d2A44[14] ) * tx + d2A44[15];
      }
    };

  template<typename T>
    struct MultiBspline: public MultiBsplineData<T>
    {

#if (__cplusplus< 201103L)
      ///define the einsplie object type
      typedef  typename bspline_traits<T,3>::SplineType spliner_type;
      ///define the real type
      typedef typename bspline_traits<T,3>::real_type real_type;
#else
      ///define the einsplie object type
      using spliner_type=typename bspline_traits<T,3>::SplineType;
      ///define the real type
      using real_type=typename bspline_traits<T,3>::real_type;
#endif
      ///set to true if create is invoked
      bool own_spline;
      ///actual einspline multi-bspline object
      spliner_type* spline_m;
      ///offset
      std::vector<int> offset;
      ///use allocator
      einspline::Allocator myAllocator;

      MultiBspline():own_spline(false),spline_m(nullptr) {}
#if (__cplusplus >= 201103L)
      MultiBspline(const MultiBspline& in)=delete;
      MultiBspline& operator=(const MultiBspline& in)=delete;
#endif

      //MultiBspline(const MultiBspline& in):own_spline(false),spline_m(in.splime_m),offset(in.offset)
      //{ }
      //MultiBspline& operator=(const MultiBspline& in)
      //{
      //}

      ~MultiBspline()
      {
        if(own_spline)
        {
          myAllocator.destroy(spline_m);
        }
      }

      template<typename RV, typename IV>
      void create(RV& start, RV& end, IV& ng, bc_code bc, int num_splines, int nteams=1)
      {
        if(spline_m==nullptr)
        {
          spline_m=myAllocator.createMultiBspline(T(0),start,end,ng,bc,num_splines);
          own_spline=true;
        }
        //should be refined to ensure alignment with minimal waste
        int nsb=num_splines/nteams;
        offset.resize(nteams+1);
        for(int i=0; i<nteams; ++i) offset[i]=i*nsb;
        offset[nteams]=num_splines;
      }

      /** create the einspline as used in the builder
       */
      template<typename GT, typename BCT>
      void create(GT& grid, BCT& bc, int num_splines, int nteams=1)
      {
        if(spline_m==nullptr)
        {
          typename bspline_traits<T,3>::BCType xBC, yBC, zBC;
          xBC.lCode=bc[0].lCode; yBC.lCode=bc[1].lCode; zBC.lCode=bc[2].lCode;
          xBC.rCode=bc[0].rCode; yBC.rCode=bc[1].rCode; zBC.rCode=bc[2].rCode;
          xBC.lVal=static_cast<T>(bc[0].lVal); yBC.lVal=static_cast<T>(bc[1].lVal); zBC.lVal=static_cast<T>(bc[2].lVal);
          xBC.rVal=static_cast<T>(bc[0].rVal); yBC.rVal=static_cast<T>(bc[1].rVal); zBC.rVal=static_cast<T>(bc[2].rVal);
          spline_m=myAllocator.allocateMultiBspline(grid[0],grid[1],grid[2],xBC,yBC,zBC,num_splines);
          own_spline=true;
        }
        //should be refined to ensure alignment with minimal waste
        int nsb=num_splines/nteams;
        offset.resize(nteams+1);
        for(int i=0; i<nteams; ++i) offset[i]=i*nsb;
        offset[nteams]=num_splines;
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

      /*
      void print(std::ostream& os)
      {
        std::copy(A44,A44+16,std::ostream_iterator<T>(std::cout," "));
        os << std::endl;
      }
      */

      template<typename PT, typename VT>
        void evaluate(const PT& r, VT& psi)
        {
          evaluate_v_impl(r[0],r[1],r[2],psi.data(),0,psi.size());
        }

      template<typename PT, typename VT>
        void evaluate(const PT& r, VT& psi, int ip)
        {
          const int first=offset[ip];
          evaluate_v_impl(r[0],r[1],r[2],psi.data()+first,first,offset[ip+1]);
        }

      template<typename PT, typename VT, typename GT>
        inline void evaluate(const PT& r, VT& psi, GT& grad)
        {
          //einspline::evaluate(spliner,r,psi,grad);
        }

      template<typename PT, typename VT, typename GT, typename LT>
        inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, LT& lap)
        {
          evaluate_vgl_impl(r[0],r[1],r[2],psi.data(),grad.data(),lap.data(),0,psi.size());
        }

      template<typename PT, typename VT, typename GT, typename LT>
        inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, LT& lap, int ip)
        {
          const int first=offset[ip];
          evaluate_vgl_impl(r[0],r[1],r[2],psi.data()+first,grad.data()+first,lap.data()+first,first,offset[ip+1]);
        }

      template<typename PT, typename VT, typename GT, typename HT>
        inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
        {
          evaluate_vgh_impl(r[0],r[1],r[2],psi.data(),grad.data(),hess.data(),0,psi.size());
        }

      template<typename PT, typename VT, typename GT, typename HT>
        inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess, int ip)
        {
          const int first=offset[ip];
          evaluate_vgh_impl(r[0],r[1],r[2],psi.data()+first,grad.data()+first,hess.data()+first,first,offset[ip+1]);
        }

      /** compute values vals[first,last)
       *
       * The base address for vals, grads and lapl are set by the callers, e.g., evaluate_vgh(r,psi,grad,hess,ip).
       */
      void evaluate_v_impl(T x, T y, T z, T* restrict vals, int first, int last) const;

      void evaluate_vgl_impl(T x, T y, T z, T* restrict vals, T* restrict grads, T* restrict lapl, int first, int last, size_t out_offset=0) const;

      void evaluate_vgh_impl(T x, T y, T z, T* restrict vals, T* restrict grads, T* restrict hess, int first, int last, size_t out_offset=0) const;
    };

  template<typename T>
    struct MultiBspline1D: public MultiBsplineData<T>
    {

#if (__cplusplus< 201103L)
      ///define the einsplie object type
      typedef  typename bspline_traits<T,1>::SplineType spliner_type;
      ///define the real type
      typedef typename bspline_traits<T,1>::real_type real_type;
#else
      ///define the einsplie object type
      using spliner_type=typename bspline_traits<T,1>::SplineType;
      ///define the real type
      using real_type=typename bspline_traits<T,1>::real_type;
#endif
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
      void create(GT& grid, BCT& bc, int num_splines, const bool ishandle=false)
      {
        spliner_type* temp_spline;
        temp_spline=einspline::create(temp_spline, grid, bc, num_splines);
        if(ishandle)
        {
          free(temp_spline->coefs);
          temp_spline->coefs=nullptr;
        }
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
        inline void evaluate(const PT& r, VT& psi)
        {
          //einspline::evaluate(&spline_m,r,psi);
          evaluate_v_impl(r,psi.data());
        }

      template<typename PT, typename VT, typename GT, typename LT>
        inline void evaluate_vgl(const PT& r, VT& psi, GT& grad, LT& lap)
        {
          //einspline::evaluate(&spline_m,r,psi,grad,lap);
          evaluate_vgl_impl(r,psi.data(),grad.data(),lap.data());
        }

      template<typename PT, typename VT, typename GT, typename HT>
        inline void evaluate_vgh(const PT& r, VT& psi, GT& grad, HT& hess)
        {
          //einspline::evaluate(&spline_m,r,psi,grad,hess);
          evaluate_vgl_impl(r,psi.data(),grad.data(),hess.data());
        }

      /// compute values only.
      void evaluate_v_impl(T r, T* restrict vals) const;
      /// compute VGL.
      void evaluate_vgl_impl(T r, T* restrict vals, T* restrict grads, T* restrict lapl) const;
    };

}/** qmcplusplus namespace */

///include evaluate_v_impl
#include <spline2/MultiBsplineValue.hpp>

/** choose vgl/vgh, default MultiBsplineStd.hpp based on Ye's BGQ version
 * Only used by tests
 */
#ifdef USE_EINSPLINE_UNROLLED
#include <spline2/MultiBsplineStd2.hpp>
#elif USE_EINSPLINE_STD4
#include <spline2/MultiBsplineStd4.hpp>
#elif USE_EINSPLINE_BASE
#include <spline2/MultiBsplineBase.hpp>
#elif USE_EINSPLINE_BLOCKED
#include <spline2/MultiBsplineStd5.hpp>
#else
#include <spline2/MultiBsplineStd.hpp>
#endif

#include <spline2/MultiBspline1D.hpp>

#endif


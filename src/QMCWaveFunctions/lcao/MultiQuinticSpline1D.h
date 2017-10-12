//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_MULTI_FUNCTOR_QUINTIC_SPLINE_SET_H
#define QMCPLUSPLUS_MULTI_FUNCTOR_QUINTIC_SPLINE_SET_H

#include "Numerics/OneDimGridBase.h"
#include "Numerics/NRSplineFunctions.h"
#include <simd/allocator.hpp>

namespace qmcplusplus
{

  template<typename T>
    struct LogGridLight
    {
      T lower_bound;
      T upper_bound;
      T OneOverLogDelta;
      double LogDelta;

      inline void set(T ri, T rf, int n)
      {
        lower_bound=ri;
        upper_bound=rf;
        //num_points=n;
        double ratio = rf/ri;
        double log_ratio = std::log(ratio);
        double dlog_ratio = log_ratio/static_cast<double>(n-1);
        LogDelta = dlog_ratio;
        OneOverLogDelta = 1.0/dlog_ratio;
      }

      inline int locate(T r) const
      {
        return static_cast<int>(std::log(r/lower_bound)*OneOverLogDelta);
      }

      inline T operator()(int i)
      {
        return static_cast<T>(lower_bound*std::exp(i*LogDelta));
      }

      //CHECK MIXED PRECISION SENSITIVITY
      inline T getCLForQuintic(T r, int& loc) const
      {
        loc=static_cast<int>(std::log(r/lower_bound)*OneOverLogDelta);
        //return r-static_cast<T>(lower_bound*std::exp(loc*LogDelta));
        return r-lower_bound*std::exp(loc*LogDelta);
      }
    };

  /** multivalue implementation for OneDimQuintic
   */
template<typename T>
struct MultiQuinticSpline1D
{
  typedef T value_type;
  typedef OneDimGridBase<T> grid_type;
  typedef Matrix<value_type, aligned_allocator<value_type> > coeff_type;

  ///number of splines
  size_t num_splines;
  ///order of spline
  size_t spline_order;
  ///default is false
  bool own_spline;
  ///maximum radius
  T r_max;

  ///will be the real grid
  LogGridLight<T> myGrid;

  /** Need to use shared_ptr
   *
   * Coeffs[6*spline_points][num_splines+padding]
   */
  coeff_type* Coeffs;
  aligned_vector<value_type> first_deriv;

  MultiQuinticSpline1D():own_spline(false),Coeffs(nullptr){}

  MultiQuinticSpline1D(const MultiQuinticSpline1D& in)=default;

  MultiQuinticSpline1D<T>* makeClone() const
  {
    MultiQuinticSpline1D<T>* myclone= new MultiQuinticSpline1D<T>(*this);
    myclone->own_spline=false;
    return myclone;
  }

  inline T rmax() const
  {
    return myGrid.upper_bound;
  }

  inline void evaluate(T r, T* restrict u) 
  {
    if(r<myGrid.lower_bound)
    {
      const value_type dr=r-myGrid.lower_bound;
      const value_type* restrict a=(*Coeffs)[0];
      for(size_t i=0; i<num_splines; ++i)
        u[i]=a[i]+first_deriv[i]*dr;
    }
#if 0
    else if(r>=r_max)
    {
      CONSTEXPR T czero(0);
      for(size_t i=0; i<num_splines; ++i)
        u[i]=czero;
    }
#endif
    else
    {
      int loc;
      const auto cL=myGrid.getCLForQuintic(r,loc);
      const size_t offset=loc*6;
      const value_type* restrict a=(*Coeffs)[offset+0];
      const value_type* restrict b=(*Coeffs)[offset+1];
      const value_type* restrict c=(*Coeffs)[offset+2];
      const value_type* restrict d=(*Coeffs)[offset+3];
      const value_type* restrict e=(*Coeffs)[offset+4];
      const value_type* restrict f=(*Coeffs)[offset+5];
      for(size_t i=0; i<num_splines; ++i)
        u[i]=a[i]+cL*(b[i]+cL*(c[i]+cL*(d[i]+cL*(e[i]+cL*f[i]))));
    }
  }

  inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u) 
  {
    constexpr value_type czero(0);

    if(r<myGrid.lower_bound)
    {
      const value_type dr=r-myGrid.lower_bound;
      const value_type* restrict a=(*Coeffs)[0];
      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]=a[i]+first_deriv[i]*dr;
        du[i]=first_deriv[i];
        d2u[i]=czero;
      }
    }
#if 0
    //should never come to this
    else if(r>=r_max)
    {
      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]=czero;
        du[i]=czero;
        du[i]=czero;
      }
    }
#endif
    else
    {
      int loc;
      const auto cL=myGrid.getCLForQuintic(r,loc);
      const size_t offset=loc*6;

      constexpr value_type ctwo(2);
      constexpr value_type cthree(3);
      constexpr value_type cfour(4);
      constexpr value_type cfive(5);
      constexpr value_type csix(6);
      constexpr value_type c12(12);
      constexpr value_type c20(20);

      const value_type* restrict a=(*Coeffs)[offset+0];
      const value_type* restrict b=(*Coeffs)[offset+1];
      const value_type* restrict c=(*Coeffs)[offset+2];
      const value_type* restrict d=(*Coeffs)[offset+3];
      const value_type* restrict e=(*Coeffs)[offset+4];
      const value_type* restrict f=(*Coeffs)[offset+5];

      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]  = a[i]+cL*(b[i]+cL*(c[i]+cL*(d[i]+cL*(e[i]+cL*f[i]))));
        du[i] = b[i]+cL*(ctwo*c[i]+cL*(cthree*d[i]+cL*(cfour*e[i]+cL*f[i]*cfive)));
        d2u[i]= ctwo*c[i]+cL*(csix*d[i]+cL*(c12*e[i]+cL*f[i]*c20));
      }
    }
  }

  /** compute upto 3rd derivatives */
  inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u, T* restrict d3u) 
  {
    constexpr value_type czero(0);

    if(r<myGrid.lower_bound)
    {
      const value_type dr=r-myGrid.lower_bound;
      const value_type* restrict a=(*Coeffs)[0];
      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]=a[i]+first_deriv[i]*dr;
        du[i]=first_deriv[i];
        d2u[i]=czero;
        d3u[i]=czero;
      }
    }
    else
    {
      int loc;
      const auto cL=myGrid.getCLForQuintic(r,loc);
      const size_t offset=loc*6;

      constexpr value_type ctwo(2);
      constexpr value_type cthree(3);
      constexpr value_type cfour(4);
      constexpr value_type cfive(5);
      constexpr value_type csix(6);
      constexpr value_type c12(12);
      constexpr value_type c20(20);
      constexpr value_type c24(24);
      constexpr value_type c60(60);

      const value_type* restrict a=(*Coeffs)[offset+0];
      const value_type* restrict b=(*Coeffs)[offset+1];
      const value_type* restrict c=(*Coeffs)[offset+2];
      const value_type* restrict d=(*Coeffs)[offset+3];
      const value_type* restrict e=(*Coeffs)[offset+4];
      const value_type* restrict f=(*Coeffs)[offset+5];

      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]  = a[i]+cL*(b[i]+cL*(c[i]+cL*(d[i]+cL*(e[i]+cL*f[i]))));
        du[i] = b[i]+cL*(ctwo*c[i]+cL*(cthree*d[i]+cL*(cfour*e[i]+cL*f[i]*cfive)));
        d2u[i]= ctwo*c[i]+cL*(csix*d[i]+cL*(c12*e[i]+cL*f[i]*c20));
        d3u[i] = csix*d+cL*(c24*e+cL*f*c60);
      }
    }
  }

  /** initialize grid and container 
   * @param ri minimum  grid point
   * @param rf maximum grid point
   * @param npts number of grid points
   * @param n number of splines
   * @param oreder 5=quintic and 3=cubic
   */
  template<typename GT>
  void initialize(GT* agrid, int norbs,int order=5)
  {
    myGrid.set(agrid->rmin(),agrid->rmax(),agrid->size());
    r_max=myGrid.upper_bound;
    if(Coeffs==nullptr && !own_spline)
    {
      spline_order=order;
      num_splines=norbs;
      Coeffs=new coeff_type((order+1)*agrid->size(),getAlignedSize<T>(norbs));
      first_deriv.resize(num_splines);
      own_spline=true;
    }
  }

  template<typename T1>
  void add_spline(int ispline, OneDimQuinticSpline<T1>& in)
  {
    first_deriv[ispline]=in.first_deriv;
    //if(spline_order==QUINTIC)
    {
      const T1* restrict A=in.m_Y.data();
      const T1* restrict B=in.B.data();
      const T1* restrict C=in.m_Y2.data();
      const T1* restrict D=in.D.data();
      const T1* restrict E=in.E.data();
      const T1* restrict F=in.F.data();
      value_type* restrict out=Coeffs->data();
      const size_t ncols=Coeffs->cols();
      const size_t num_points=in.size();
      for(size_t i=0; i<num_points; ++i)
      {
        out[(i*6+0)*ncols+ispline]=static_cast<T>(A[i]);
        out[(i*6+1)*ncols+ispline]=static_cast<T>(B[i]);
        out[(i*6+2)*ncols+ispline]=static_cast<T>(C[i]);
        out[(i*6+3)*ncols+ispline]=static_cast<T>(D[i]);
        out[(i*6+4)*ncols+ispline]=static_cast<T>(E[i]);
        out[(i*6+5)*ncols+ispline]=static_cast<T>(F[i]);
      }
    }
  }
};


}
#endif

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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_GRID_FUNCTOR_QUINTIC_SPLINE_SET_H
#define QMCPLUSPLUS_GRID_FUNCTOR_QUINTIC_SPLINE_SET_H

#include "Numerics/OneDimGridBase.h"
#include "Numerics/NRSplineFunctions.h"
#include <simd/allocator.hpp>

namespace qmcplusplus
{

  /** multivalue implementation for OneDimQuintic
   */
template<typename T>
struct MultiQuinticSpline1D
{
  typedef T value_type;
  typedef OneDimGridBase<T> grid_type;

  ///number of splines
  size_t num_splines;
  ///order of spline
  size_t spline_order;
  ///default is false
  bool own_spline;

  /** shared_ptr
   * Coeffs[6*spline_points][num_splines]
   */
  grid_type* m_grid;
  Matrix<value_type> *Coeffs;
  aligned_vector<value_type> first_deriv;

  MultiQuinticSpline1D():own_spline(false),m_grid(nullptr),Coeffs(nullptr){}

  MultiQuinticSpline1D(const MultiQuinticSpline1D& in)=default;

  MultiQuinticSpline1D<T>* makeClone() const
  {
    MultiQuinticSpline1D<T>* myclone= new MultiQuinticSpline1D<T>(*this);
    myclone->m_grid=m_grid->makeClone(); //wasting X, to remove this
    myclone->own_spline=false;
    return myclone;
  }

  inline T rmax() const
  {
    m_grid->rmax();
  }

  inline void evaluate(T r, T* restrict u) 
  {
    if(r<m_grid->rmin())
    {
      const value_type dr=r-m_grid->rmin();
      const value_type* restrict a=(*Coeffs)[0];
      for(size_t i=0; i<num_splines; ++i)
        u[i]=a[i]+first_deriv[i]*dr;
    }
    else
    {
      //m_grid->locate(r);
      m_grid->updateForQuintic(r,false);
      const size_t offset=m_grid->Loc*6; 
      const value_type* restrict a=(*Coeffs)[offset+0];
      const value_type* restrict b=(*Coeffs)[offset+1];
      const value_type* restrict c=(*Coeffs)[offset+2];
      const value_type* restrict d=(*Coeffs)[offset+3];
      const value_type* restrict e=(*Coeffs)[offset+4];
      const value_type* restrict f=(*Coeffs)[offset+5];
      const auto cL=m_grid->cL;
      for(size_t i=0; i<num_splines; ++i)
        u[i]=a[i]+cL*(b[i]+cL*(c[i]+cL*(d[i]+cL*(e[i]+cL*f[i]))));
    }
  }

  inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u) 
  {
    constexpr value_type czero(0);
    if(r<m_grid->rmin())
    {
      const T dr=r-m_grid->rmin();
      const value_type* restrict a=(*Coeffs)[0];
      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]=a[i]+first_deriv[i]*dr;
        du[i]=first_deriv[i];
        d2u[i]=czero;
      }
    }
    else
    {
      m_grid->updateForQuintic(r,true);
      constexpr value_type ctwo(2);
      constexpr value_type cthree(3);
      constexpr value_type cfour(4);
      constexpr value_type cfive(5);
      constexpr value_type csix(6);
      constexpr value_type c12(5);
      constexpr value_type c20(6);

      const size_t offset=m_grid->Loc*6; 
      const value_type* restrict a=(*Coeffs)[offset+0];
      const value_type* restrict b=(*Coeffs)[offset+1];
      const value_type* restrict c=(*Coeffs)[offset+2];
      const value_type* restrict d=(*Coeffs)[offset+3];
      const value_type* restrict e=(*Coeffs)[offset+4];
      const value_type* restrict f=(*Coeffs)[offset+5];
      const auto cR=m_grid->cR;
      const auto cL=m_grid->cL;
      for(size_t i=0; i<num_splines; ++i)
      {
        u[i]  = a[i]+cL*(b[i]+cL*(c[i]+cL*(d[i]+cL*(e[i]+cL*f[i]))));
        du[i] = b[i]+cL*(ctwo*c[i]+cL*(cthree*d[i]+cL*(cfour*e[i]+cL*f[i]*cfive)));
        d2u[i]= ctwo*c[i]+cL*(csix*d[i]+cL*(c12*e[i]+cL*f[i]*c20));
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
  void initialize(grid_type* agrid, int norbs,int order=5)
  {
    m_grid=agrid;
    //this is very shaky
    if(Coeffs==nullptr && !own_spline)
    {
      spline_order=order;
      num_splines=norbs;
      Coeffs=new Matrix<value_type>((order+1)*m_grid->size(),getAlignedSize<T>(norbs));
      first_deriv.resize(num_splines);
      own_spline=true;
    }
  }

  template<typename T1>
  void add_spline(int ispline, OneDimQuinticSpline<T1>& in)
  {
    const T1* restrict A=in.m_Y.data();
    const T1* restrict B=in.B.data();
    const T1* restrict C=in.m_Y2.data();
    const T1* restrict D=in.D.data();
    const T1* restrict E=in.E.data();
    const T1* restrict F=in.F.data();

    first_deriv[ispline]=in.first_deriv;

    value_type* restrict out=Coeffs->data();
    const size_t ncols=Coeffs->cols();
    const size_t num_points=in.size();
    for(size_t i=0; i<num_points; ++i)
    {
      out[(i*6+0)*ncols+ispline]=A[i];
      out[(i*6+1)*ncols+ispline]=B[i];
      out[(i*6+2)*ncols+ispline]=C[i];
      out[(i*6+3)*ncols+ispline]=D[i];
      out[(i*6+4)*ncols+ispline]=E[i];
      out[(i*6+5)*ncols+ispline]=F[i];
    }

  }
#if 0
  template<typename InType, typename VV>
  inline void add_spline(int ispline, VV& yin, InType yp1, InType ypn)
  {
    //if(ispline>num_splines) DO SOMETHING
    
    first_deriv[ispline]=yp1;
    const typename VV::value_type czero(0);
    int num_points=m_grid->size();

    aligned_vector<InType> B(num_points,czero),C(num_points,czero),D(num_points,czero),E(num_points,czero),F(num_points,czero);
    QuinticSplineSolve(num_points-1,m_grid->data(), yin.data(), B.data(), C.data(),D.data(),E.data(),F.data());

    value_type* restrict out=Coeffs->data();
    const size_t ncols=Coeffs->cols();
    for(size_t i=0; i<num_points; ++i)
    {
      //std::cout << m_grid->operator[](i) << " " << yin[i] << std::endl;
      out[(i*6+0)*ncols+ispline]=yin[i];
      out[(i*6+1)*ncols+ispline]=B[i];
      out[(i*6+2)*ncols+ispline]=C[i];
      out[(i*6+3)*ncols+ispline]=D[i];
      out[(i*6+4)*ncols+ispline]=E[i];
      out[(i*6+5)*ncols+ispline]=F[i];
    }

    first_deriv[ispline]=yp1;
  }
#endif

};


}
#endif

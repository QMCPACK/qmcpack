//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef GUARD_TRICUBICSPLINE_H
#define GUARD_TRICUBICSPLINE_H

#include "Numerics/Spline3D/Config.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/SetSplinePoint.h"
#include <fstream>
#include <blitz/array.h>

/// This is tri-cubic Splines with natural boundary conditions, i.e.
/// the second derivative vanishes at the boundary, not suitable for
/// functions like sines and cosines.

/// Each point of F contains:
/// 0) F(x,y,z)
/// 1) dF/dx
/// 2) dF/dy
/// 3) dF/dz
/// 4) d2F/dxdy
/// 5) d2F/dxdz
/// 6) d2F/dydz
/// 7) d3F/dxdydz

class TriCubicSpline
{

  /// functions which depend on the point where the interpolated value
  /// is required. t = (x - xi)/h
  inline double p1(double t)
  {
    return ((t-1.0)*(t-1.0)*(1.0+2.0*t));
  }
  inline double p2(double t)
  {
    return (t*t*(3.0-2.0*t));
  }
  inline double q1(double t)
  {
    return (t*(t-1.0)*(t-1.0));
  }
  inline double q2(double t)
  {
    return (t*t*(t-1.0));
  }
  inline double dp1(double t)
  {
    return (6.0*t*(t-1.0));
  }
  inline double dq1(double t)
  {
    return ((t-1.0)*(3.0*t-1.0));
  }
  inline double dp2(double t)
  {
    return (-dp1(t));
  }
  inline double dq2 (double t)
  {
    return ((3.0*t - 2.0)*t);
  }
  inline double d2p1(double t)
  {
    return (12.0*t-6.0);
  }
  inline double d2q1 (double t)
  {
    return (6.0*t - 4.0);
  }
  inline double d2p2 (double t)
  {
    return (-d2p1(t));
  }
  inline double d2q2 (double t)
  {
    return (6.0*t - 2.0);
  }

  // dim:     Dimension to calculate derivative w.r.t
  // source:  Function to differentiate
  // dest:    where to put result
  void UpdateX (int source, int dest);
  void UpdateY (int source, int dest);
  void UpdateZ (int source, int dest);

  void UpdateD2X (int source, int dest);
  void UpdateD2Y (int source, int dest);
  void UpdateD2Z (int source, int dest);

  void D2FDR2();
  double delsqf(const gridvec_t&, int);


  /// whether the first derivatives have been calculated using the
  /// m-relations.
  bool UpToDate;

  int h,k,l;
  double u,v,w;

  inline double dx(int dir, int i)
  {
    return m_grid->m_axis[dir].h(i);
  }

  /// function and derivatives at each point
  blitz::Array<blitz::TinyVector<double,8>,3> F;
  blitz::Array<blitz::TinyVector<double,8>,3> D2F;


public:

  typedef double value_type;

  /// points on the Grid3D
  int n_x, n_y, n_z;

  /// pointer to the Schroedinger Grid
  Grid3D* m_grid;

  /// pointer to the function to assign point r properties
  SetSplinePoint* m_set;

  /// constructor
  TriCubicSpline(SetSplinePoint* aset,
                 Grid3D* agrid)
  {
    m_grid = agrid;
    m_set = aset;
    n_x = agrid->n_x;
    n_y = agrid->n_y;
    n_z = agrid->n_z;
    F.resize(n_x,n_y,n_z);
    UpToDate = false;
  }

  TriCubicSpline(Grid3D* agrid)
  {
    m_grid = agrid;
    if(!m_set)
      m_set = new SetSplinePoint;
    n_x = agrid->n_x;
    n_y = agrid->n_y;
    n_z = agrid->n_z;
    F.resize(n_x,n_y,n_z);
    UpToDate = false;
  }


  inline void reset() { }

  inline double operator()(int ix, int iy, int iz) const
  {
    return (F(ix,iy,iz)[0]);
  }

  inline double& operator()(int ix, int iy, int iz)
  {
    UpToDate = false;
    return (F(ix,iy,iz)[0]);
  }

  inline void set_point(const posvec_t& r)
  {
    m_set->set_point(r,m_grid);
    return;
  }

  /// update the derivatives using the m-relations
  void Update(bool);

  inline double evaluate(const posvec_t& r)
  {
    if( !UpToDate )
      Update(false);              /// m-relations
    int ix = m_set->ix;
    int iy = m_set->iy;
    int iz = m_set->iz;
    double h = m_set->h;
    double k = m_set->k;
    double l = m_set->l;
    double u = m_set->u;
    double v = m_set->v;
    double w = m_set->w;
    double a0 = p1(u);
    double a1 = p2(u);
    double a2 = h*q1(u);
    double a3 = h*q2(u);
    register double b0 = p1(v);
    register double b1 = p2(v);
    register double b2 = k*q1(v);
    register double b3 = k*q2(v);
    register double c0 = p1(w);
    register double c1 = p2(w);
    register double c2 = l*q1(w);
    register double c3 = l*q2(w);
    double& Y000 = F(ix,iy,iz)[0];      //   F
    double& Y002 = F(ix,iy,iz)[3];      //  dF/dz
    double& Y001 = F(ix,iy,iz+1)[0];    //   F
    double& Y003 = F(ix,iy,iz+1)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz)[0];    //   F
    double& Y012 = F(ix,iy+1,iz)[3];    //  dF/dz
    double& Y011 = F(ix,iy+1,iz+1)[0];  //   F
    double& Y013 = F(ix,iy+1,iz+1)[3];  //  dF/dz
    double& Y021 = F(ix,iy,iz+1)[2];    //  dF/dy
    double& Y020 = F(ix,iy,iz)[2];      //  dF/dy
    double& Y022 = F(ix,iy,iz)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1)[6];  // d2F/dydz
    double& Y100 = F(ix+1,iy,iz)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1)[0];    //   F
    double& Y102 = F(ix+1,iy,iz)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz
    double& Y200 = F(ix,iy,iz)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz
    double& Y300 = F(ix+1,iy,iz)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
    double val =
      a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    return val;
  }

  inline double evaluate(const posvec_t& r,
                         posvec_t& gradf)
  {
    if( !UpToDate )
      Update(true);              /// m-relations
    int ix = m_set->ix;
    int iy = m_set->iy;
    int iz = m_set->iz;
    double h = m_set->h;
    double k = m_set->k;
    double l = m_set->l;
    double u = m_set->u;
    double v = m_set->v;
    double w = m_set->w;
    double hinv = m_set->hinv;
    double kinv = m_set->kinv;
    double linv = m_set->linv;
    double a0 = p1(u);
    double a1 = p2(u);
    double a2 = h*q1(u);
    double a3 = h*q2(u);
    double b0 = p1(v);
    double b1 = p2(v);
    double b2 = k*q1(v);
    double b3 = k*q2(v);
    double c0 = p1(w);
    double c1 = p2(w);
    double c2 = l*q1(w);
    double c3 = l*q2(w);
    double da0 = hinv*dp1(u);
    double da1 = hinv*dp2(u);
    double da2 = dq1(u);
    double da3 = dq2(u);
    double db0 = kinv*dp1(v);
    double db1 = kinv*dp2(v);
    double db2 = dq1(v);
    double db3 = dq2(v);
    double dc0 = linv*dp1(w);
    double dc1 = linv*dp2(w);
    double dc2 = dq1(w);
    double dc3 = dq2(w);
    double& Y000 = F(ix,iy,iz)[0];      //   F
    double& Y001 = F(ix,iy,iz+1)[0];    //   F
    double& Y002 = F(ix,iy,iz)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1)[0];  //   F
    double& Y012 = F(ix,iy+1,iz)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1)[6];  // d2F/dydz
    double& Y100 = F(ix+1,iy,iz)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1)[0];    //   F
    double& Y102 = F(ix+1,iy,iz)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz
    double& Y200 = F(ix,iy,iz)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz
    double& Y300 = F(ix+1,iy,iz)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
    double term0 =
      b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
      b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
      b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
      b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3);
    double term1 =
      b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
      b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
      b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
      b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3);
    double term2 =
      b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
      b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
      b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
      b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3);
    double term3 =
      b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
      b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
      b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
      b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3);
    double val = a0*term0 + a1*term1 + a2*term2 + a3*term3;
    gradf[0] = da0*term0 + da1*term1 + da2*term2 + da3*term3;
    gradf[1] =
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    gradf[2] =
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    return val;
  }

  inline double evaluate_all(const posvec_t& r,
                             posvec_t& gradf,
                             double& lapf)
  {
    if( !UpToDate )
      Update(false);              /// m-relations
    /// if the point is out of bounds, return small number
    if( m_set->ifout )
    {
      m_set->ifout = false;
      /*
      std::cout << r[0] << '\t' << m_grid->m_axis[0].m_start << '\t'
      << m_grid->m_axis[0].m_end << std::endl;
      std::cout << r[1] << '\t' << m_grid->m_axis[1].m_start << '\t'
      << m_grid->m_axis[1].m_end << std::endl;
      std::cout << r[2] << '\t' << m_grid->m_axis[2].m_start << '\t'
      << m_grid->m_axis[2].m_end << std::endl;
      */
      gradf[0] = 1e-40;
      gradf[1] = 1e-40;
      gradf[2] = 1e-40;
      lapf = 1e-40;
      return 1e-20;
    }
    //    std::cout << "pointr: " << r << std::endl;
    int ix = m_set->ix;
    int iy = m_set->iy;
    int iz = m_set->iz;
    double h = m_set->h;
    double k = m_set->k;
    double l = m_set->l;
    double u = m_set->u;
    double v = m_set->v;
    double w = m_set->w;
    double hinv = m_set->hinv;
    double kinv = m_set->kinv;
    double linv = m_set->linv;
    double a0 = p1(u);
    double a1 = p2(u);
    double a2 = h*q1(u);
    double a3 = h*q2(u);
    double b0 = p1(v);
    double b1 = p2(v);
    double b2 = k*q1(v);
    double b3 = k*q2(v);
    double c0 = p1(w);
    double c1 = p2(w);
    double c2 = l*q1(w);
    double c3 = l*q2(w);
    double da0 = hinv*dp1(u);
    double da1 = hinv*dp2(u);
    double da2 = dq1(u);
    double da3 = dq2(u);
    double db0 = kinv*dp1(v);
    double db1 = kinv*dp2(v);
    double db2 = dq1(v);
    double db3 = dq2(v);
    double dc0 = linv*dp1(w);
    double dc1 = linv*dp2(w);
    double dc2 = dq1(w);
    double dc3 = dq2(w);
    double d2a0 = hinv*hinv*d2p1(u);
    double d2a1 = hinv*hinv*d2p2(u);
    double d2a2 = hinv*d2q1(u);
    double d2a3 = hinv*d2q2(u);
    double d2b0 = kinv*kinv*d2p1(v);
    double d2b1 = kinv*kinv*d2p2(v);
    double d2b2 = kinv*d2q1(v);
    double d2b3 = kinv*d2q2(v);
    double d2c0 = linv*linv*d2p1(w);
    double d2c1 = linv*linv*d2p2(w);
    double d2c2 = linv*d2q1(w);
    double d2c3 = linv*d2q2(w);
    double& Y000 = F(ix,iy,iz)[0];      //   F
    double& Y001 = F(ix,iy,iz+1)[0];    //   F
    double& Y002 = F(ix,iy,iz)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1)[0];  //   F
    double& Y012 = F(ix,iy+1,iz)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1)[6];  // d2F/dydz
    double& Y100 = F(ix+1,iy,iz)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1)[0];    //   F
    double& Y102 = F(ix+1,iy,iz)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz
    double& Y200 = F(ix,iy,iz)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz
    double& Y300 = F(ix+1,iy,iz)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
    double term0 =
      b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
      b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
      b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
      b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3);
    double term1 =
      b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
      b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
      b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
      b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3);
    double term2 =
      b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
      b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
      b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
      b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3);
    double term3 =
      b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
      b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
      b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
      b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3);
    double val = a0*term0 + a1*term1 + a2*term2 + a3*term3;
    gradf[0] = da0*term0 + da1*term1 + da2*term2 + da3*term3;
    gradf[1] =
      a0*
      (db0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       db1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       db2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       db3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (db0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       db1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       db2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       db3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (db0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       db1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       db2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       db3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (db0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       db1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       db2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       db3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    gradf[2] =
      a0*
      (b0*(Y000*dc0+Y001*dc1+Y002*dc2+Y003*dc3) +
       b1*(Y010*dc0+Y011*dc1+Y012*dc2+Y013*dc3) +
       b2*(Y020*dc0+Y021*dc1+Y022*dc2+Y023*dc3) +
       b3*(Y030*dc0+Y031*dc1+Y032*dc2+Y033*dc3))+
      a1 *
      (b0*(Y100*dc0+Y101*dc1+Y102*dc2+Y103*dc3) +
       b1*(Y110*dc0+Y111*dc1+Y112*dc2+Y113*dc3) +
       b2*(Y120*dc0+Y121*dc1+Y122*dc2+Y123*dc3) +
       b3*(Y130*dc0+Y131*dc1+Y132*dc2+Y133*dc3))+
      a2 *
      (b0*(Y200*dc0+Y201*dc1+Y202*dc2+Y203*dc3) +
       b1*(Y210*dc0+Y211*dc1+Y212*dc2+Y213*dc3) +
       b2*(Y220*dc0+Y221*dc1+Y222*dc2+Y223*dc3) +
       b3*(Y230*dc0+Y231*dc1+Y232*dc2+Y233*dc3))+
      a3 *
      (b0*(Y300*dc0+Y301*dc1+Y302*dc2+Y303*dc3) +
       b1*(Y310*dc0+Y311*dc1+Y312*dc2+Y313*dc3) +
       b2*(Y320*dc0+Y321*dc1+Y322*dc2+Y323*dc3) +
       b3*(Y330*dc0+Y331*dc1+Y332*dc2+Y333*dc3));
    lapf =
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3))+
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3)) +
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
    return val;
  }

  inline double evaluate(const posvec_t& r,
                         posvec_t& gradf,
                         double& lapf)
  {
    double val = evaluate_all(r,gradf,lapf);
    //double val = evaluate(r,gradf);
    //    lapf = evaluateD2F(r);
    return val;
  }

  inline double evaluateD2F(const posvec_t& r)
  {
    if( !UpToDate )
      Update(false);              /// m-relations
    int ix = m_set->ix;
    int iy = m_set->iy;
    int iz = m_set->iz;
    double h = m_set->h;
    double k = m_set->k;
    double l = m_set->l;
    double u = m_set->u;
    double v = m_set->v;
    double w = m_set->w;
    double a0 = p1(u);
    double a1 = p2(u);
    double a2 = h*q1(u);
    double a3 = h*q2(u);
    register double b0 = p1(v);
    register double b1 = p2(v);
    register double b2 = k*q1(v);
    register double b3 = k*q2(v);
    register double c0 = p1(w);
    register double c1 = p2(w);
    register double c2 = l*q1(w);
    register double c3 = l*q2(w);
    double& Y000 = D2F(ix,iy,iz)[0];      //   D2F
    double& Y001 = D2F(ix,iy,iz+1)[0];    //   D2F
    double& Y002 = D2F(ix,iy,iz)[3];      //  dD2F/dz
    double& Y003 = D2F(ix,iy,iz+1)[3];    //  dD2F/dz
    double& Y010 = D2F(ix,iy+1,iz)[0];    //   D2F
    double& Y011 = D2F(ix,iy+1,iz+1)[0];  //   D2F
    double& Y012 = D2F(ix,iy+1,iz)[3];    //  dD2F/dz
    double& Y013 = D2F(ix,iy+1,iz+1)[3];  //  dD2F/dz
    double& Y020 = D2F(ix,iy,iz)[2];      //  dD2F/dy
    double& Y021 = D2F(ix,iy,iz+1)[2];    //  dD2F/dy
    double& Y022 = D2F(ix,iy,iz)[6];      // d2D2F/dydz
    double& Y023 = D2F(ix,iy,iz+1)[6];    // d2D2F/dydz
    double& Y030 = D2F(ix,iy+1,iz)[2];    //  dD2F/dy
    double& Y031 = D2F(ix,iy+1,iz+1)[2];  //  dD2F/dy
    double& Y032 = D2F(ix,iy+1,iz)[6];    // d2D2F/dydz
    double& Y033 = D2F(ix,iy+1,iz+1)[6];  // d2D2F/dydz
    double& Y100 = D2F(ix+1,iy,iz)[0];      //   D2F
    double& Y101 = D2F(ix+1,iy,iz+1)[0];    //   D2F
    double& Y102 = D2F(ix+1,iy,iz)[3];      //  dD2F/dz
    double& Y103 = D2F(ix+1,iy,iz+1)[3];    //  dD2F/dz
    double& Y110 = D2F(ix+1,iy+1,iz)[0];    //   D2F
    double& Y111 = D2F(ix+1,iy+1,iz+1)[0];  //   D2F
    double& Y112 = D2F(ix+1,iy+1,iz)[3];    //  dD2F/dz
    double& Y113 = D2F(ix+1,iy+1,iz+1)[3];  //  dD2F/dz
    double& Y120 = D2F(ix+1,iy,iz)[2];      //  dD2F/dy
    double& Y121 = D2F(ix+1,iy,iz+1)[2];    //  dD2F/dy
    double& Y122 = D2F(ix+1,iy,iz)[6];      // d2D2F/dydz
    double& Y123 = D2F(ix+1,iy,iz+1)[6];    // d2D2F/dydz
    double& Y130 = D2F(ix+1,iy+1,iz)[2];    //  dD2F/dy
    double& Y131 = D2F(ix+1,iy+1,iz+1)[2];  //  dD2F/dy
    double& Y132 = D2F(ix+1,iy+1,iz)[6];    // d2D2F/dydz
    double& Y133 = D2F(ix+1,iy+1,iz+1)[6];  // d2D2F/dydz
    double& Y200 = D2F(ix,iy,iz)[1];      //  dD2F/dx
    double& Y201 = D2F(ix,iy,iz+1)[1];    //  dD2F/dx
    double& Y202 = D2F(ix,iy,iz)[5];      // d2D2F/dxdz
    double& Y203 = D2F(ix,iy,iz+1)[5];    // d2D2F/dxdz
    double& Y210 = D2F(ix,iy+1,iz)[1];    //  dD2F/dx
    double& Y211 = D2F(ix,iy+1,iz+1)[1];  //  dD2F/dx
    double& Y212 = D2F(ix,iy+1,iz)[5];    // d2D2F/dxdz
    double& Y213 = D2F(ix,iy+1,iz+1)[5];  // d2D2F/dxdz
    double& Y220 = D2F(ix,iy,iz)[4];      // d2D2F/dxdy
    double& Y221 = D2F(ix,iy,iz+1)[4];    // d2D2F/dxdy
    double& Y222 = D2F(ix,iy,iz)[7];      // d3D2F/dxdydz
    double& Y223 = D2F(ix,iy,iz+1)[7];    // d3D2F/dxdydz
    double& Y230 = D2F(ix,iy+1,iz)[4];    // d2D2F/dxdy
    double& Y231 = D2F(ix,iy+1,iz+1)[4];  // d2D2F/dxdy
    double& Y232 = D2F(ix,iy+1,iz)[7];    // d3D2F/dxdydz
    double& Y233 = D2F(ix,iy+1,iz+1)[7];  // d3D2F/dxdydz
    double& Y300 = D2F(ix+1,iy,iz)[1];      //  dD2F/dx
    double& Y301 = D2F(ix+1,iy,iz+1)[1];    //  dD2F/dx
    double& Y302 = D2F(ix+1,iy,iz)[5];      // d2D2F/dxdz
    double& Y303 = D2F(ix+1,iy,iz+1)[5];    // d2D2F/dxdz
    double& Y310 = D2F(ix+1,iy+1,iz)[1];    //  dD2F/dx
    double& Y311 = D2F(ix+1,iy+1,iz+1)[1];  //  dD2F/dx
    double& Y312 = D2F(ix+1,iy+1,iz)[5];    // d2D2F/dxdz
    double& Y313 = D2F(ix+1,iy+1,iz+1)[5];  // d2D2F/dxdz
    double& Y320 = D2F(ix+1,iy,iz)[4];      // d2D2F/dxdy
    double& Y321 = D2F(ix+1,iy,iz+1)[4];    // d2D2F/dxdy
    double& Y322 = D2F(ix+1,iy,iz)[7];      // d3D2F/dxdydz
    double& Y323 = D2F(ix+1,iy,iz+1)[7];    // d3D2F/dxdydz
    double& Y330 = D2F(ix+1,iy+1,iz)[4];    // d2D2F/dxdy
    double& Y331 = D2F(ix+1,iy+1,iz+1)[4];  // d2D2F/dxdy
    double& Y332 = D2F(ix+1,iy+1,iz)[7];    // d3D2F/dxdydz
    double& Y333 = D2F(ix+1,iy+1,iz+1)[7];  // d3D2F/dxdydz
    double val =
      a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3));
    return val;
  }


  void read_data(const char* data_file, const double ufac)
  {
    std::ifstream infile(data_file,std::ios_base::in);
    for(int iz = 0; iz < n_z; iz++)
    {
      for(int iy = 0; iy < n_y; iy++)
      {
        for(int ix = 0; ix < n_x; ix++)
        {
          infile >> F(ix,iy,iz)[0];
          F(ix,iy,iz)[0] *= ufac;
        }
      }
    }
    return;
  }

  inline double laplacian(const posvec_t& r)
  {
    if( !UpToDate )
      Update(false);              /// m-relations
    m_set->set_point(r,m_grid);
    int ix = m_set->ix;
    int iy = m_set->iy;
    int iz = m_set->iz;
    double h = m_set->h;
    double k = m_set->k;
    double l = m_set->l;
    double u = m_set->u;
    double v = m_set->v;
    double w = m_set->w;
    double hinv = m_set->hinv;
    double kinv = m_set->kinv;
    double linv = m_set->linv;
    double a0 = p1(u);
    double a1 = p2(u);
    double a2 = h*q1(u);
    double a3 = h*q2(u);
    double b0 = p1(v);
    double b1 = p2(v);
    double b2 = k*q1(v);
    double b3 = k*q2(v);
    double c0 = p1(w);
    double c1 = p2(w);
    double c2 = l*q1(w);
    double c3 = l*q2(w);
    double d2a0 = hinv*hinv*d2p1(u);
    double d2a1 = hinv*hinv*d2p2(u);
    double d2a2 = hinv*d2q1(u);
    double d2a3 = hinv*d2q2(u);
    double d2b0 = kinv*kinv*d2p1(v);
    double d2b1 = kinv*kinv*d2p2(v);
    double d2b2 = kinv*d2q1(v);
    double d2b3 = kinv*d2q2(v);
    double d2c0 = linv*linv*d2p1(w);
    double d2c1 = linv*linv*d2p2(w);
    double d2c2 = linv*d2q1(w);
    double d2c3 = linv*d2q2(w);
    double& Y000 = F(ix,iy,iz)[0];      //   F
    double& Y001 = F(ix,iy,iz+1)[0];    //   F
    double& Y002 = F(ix,iy,iz)[3];      //  dF/dz
    double& Y003 = F(ix,iy,iz+1)[3];    //  dF/dz
    double& Y010 = F(ix,iy+1,iz)[0];    //   F
    double& Y011 = F(ix,iy+1,iz+1)[0];  //   F
    double& Y012 = F(ix,iy+1,iz)[3];    //  dF/dz
    double& Y013 = F(ix,iy+1,iz+1)[3];  //  dF/dz
    double& Y020 = F(ix,iy,iz)[2];      //  dF/dy
    double& Y021 = F(ix,iy,iz+1)[2];    //  dF/dy
    double& Y022 = F(ix,iy,iz)[6];      // d2F/dydz
    double& Y023 = F(ix,iy,iz+1)[6];    // d2F/dydz
    double& Y030 = F(ix,iy+1,iz)[2];    //  dF/dy
    double& Y031 = F(ix,iy+1,iz+1)[2];  //  dF/dy
    double& Y032 = F(ix,iy+1,iz)[6];    // d2F/dydz
    double& Y033 = F(ix,iy+1,iz+1)[6];  // d2F/dydz
    double& Y100 = F(ix+1,iy,iz)[0];      //   F
    double& Y101 = F(ix+1,iy,iz+1)[0];    //   F
    double& Y102 = F(ix+1,iy,iz)[3];      //  dF/dz
    double& Y103 = F(ix+1,iy,iz+1)[3];    //  dF/dz
    double& Y110 = F(ix+1,iy+1,iz)[0];    //   F
    double& Y111 = F(ix+1,iy+1,iz+1)[0];  //   F
    double& Y112 = F(ix+1,iy+1,iz)[3];    //  dF/dz
    double& Y113 = F(ix+1,iy+1,iz+1)[3];  //  dF/dz
    double& Y120 = F(ix+1,iy,iz)[2];      //  dF/dy
    double& Y121 = F(ix+1,iy,iz+1)[2];    //  dF/dy
    double& Y122 = F(ix+1,iy,iz)[6];      // d2F/dydz
    double& Y123 = F(ix+1,iy,iz+1)[6];    // d2F/dydz
    double& Y130 = F(ix+1,iy+1,iz)[2];    //  dF/dy
    double& Y131 = F(ix+1,iy+1,iz+1)[2];  //  dF/dy
    double& Y132 = F(ix+1,iy+1,iz)[6];    // d2F/dydz
    double& Y133 = F(ix+1,iy+1,iz+1)[6];  // d2F/dydz
    double& Y200 = F(ix,iy,iz)[1];      //  dF/dx
    double& Y201 = F(ix,iy,iz+1)[1];    //  dF/dx
    double& Y202 = F(ix,iy,iz)[5];      // d2F/dxdz
    double& Y203 = F(ix,iy,iz+1)[5];    // d2F/dxdz
    double& Y210 = F(ix,iy+1,iz)[1];    //  dF/dx
    double& Y211 = F(ix,iy+1,iz+1)[1];  //  dF/dx
    double& Y212 = F(ix,iy+1,iz)[5];    // d2F/dxdz
    double& Y213 = F(ix,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y220 = F(ix,iy,iz)[4];      // d2F/dxdy
    double& Y221 = F(ix,iy,iz+1)[4];    // d2F/dxdy
    double& Y222 = F(ix,iy,iz)[7];      // d3F/dxdydz
    double& Y223 = F(ix,iy,iz+1)[7];    // d3F/dxdydz
    double& Y230 = F(ix,iy+1,iz)[4];    // d2F/dxdy
    double& Y231 = F(ix,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y232 = F(ix,iy+1,iz)[7];    // d3F/dxdydz
    double& Y233 = F(ix,iy+1,iz+1)[7];  // d3F/dxdydz
    double& Y300 = F(ix+1,iy,iz)[1];      //  dF/dx
    double& Y301 = F(ix+1,iy,iz+1)[1];    //  dF/dx
    double& Y302 = F(ix+1,iy,iz)[5];      // d2F/dxdz
    double& Y303 = F(ix+1,iy,iz+1)[5];    // d2F/dxdz
    double& Y310 = F(ix+1,iy+1,iz)[1];    //  dF/dx
    double& Y311 = F(ix+1,iy+1,iz+1)[1];  //  dF/dx
    double& Y312 = F(ix+1,iy+1,iz)[5];    // d2F/dxdz
    double& Y313 = F(ix+1,iy+1,iz+1)[5];  // d2F/dxdz
    double& Y320 = F(ix+1,iy,iz)[4];      // d2F/dxdy
    double& Y321 = F(ix+1,iy,iz+1)[4];    // d2F/dxdy
    double& Y322 = F(ix+1,iy,iz)[7];      // d3F/dxdydz
    double& Y323 = F(ix+1,iy,iz+1)[7];    // d3F/dxdydz
    double& Y330 = F(ix+1,iy+1,iz)[4];    // d2F/dxdy
    double& Y331 = F(ix+1,iy+1,iz+1)[4];  // d2F/dxdy
    double& Y332 = F(ix+1,iy+1,iz)[7];    // d3F/dxdydz
    double& Y333 = F(ix+1,iy+1,iz+1)[7];  // d3F/dxdydz
    double lapf =
      d2a0*
      (b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      d2a1 *
      (b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      d2a2 *
      (b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      d2a3 *
      (b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3))+
      a0*
      (d2b0*(Y000*c0+Y001*c1+Y002*c2+Y003*c3) +
       d2b1*(Y010*c0+Y011*c1+Y012*c2+Y013*c3) +
       d2b2*(Y020*c0+Y021*c1+Y022*c2+Y023*c3) +
       d2b3*(Y030*c0+Y031*c1+Y032*c2+Y033*c3))+
      a1 *
      (d2b0*(Y100*c0+Y101*c1+Y102*c2+Y103*c3) +
       d2b1*(Y110*c0+Y111*c1+Y112*c2+Y113*c3) +
       d2b2*(Y120*c0+Y121*c1+Y122*c2+Y123*c3) +
       d2b3*(Y130*c0+Y131*c1+Y132*c2+Y133*c3))+
      a2 *
      (d2b0*(Y200*c0+Y201*c1+Y202*c2+Y203*c3) +
       d2b1*(Y210*c0+Y211*c1+Y212*c2+Y213*c3) +
       d2b2*(Y220*c0+Y221*c1+Y222*c2+Y223*c3) +
       d2b3*(Y230*c0+Y231*c1+Y232*c2+Y233*c3))+
      a3 *
      (d2b0*(Y300*c0+Y301*c1+Y302*c2+Y303*c3) +
       d2b1*(Y310*c0+Y311*c1+Y312*c2+Y313*c3) +
       d2b2*(Y320*c0+Y321*c1+Y322*c2+Y323*c3) +
       d2b3*(Y330*c0+Y331*c1+Y332*c2+Y333*c3)) +
      a0*
      (b0*(Y000*d2c0+Y001*d2c1+Y002*d2c2+Y003*d2c3) +
       b1*(Y010*d2c0+Y011*d2c1+Y012*d2c2+Y013*d2c3) +
       b2*(Y020*d2c0+Y021*d2c1+Y022*d2c2+Y023*d2c3) +
       b3*(Y030*d2c0+Y031*d2c1+Y032*d2c2+Y033*d2c3))+
      a1 *
      (b0*(Y100*d2c0+Y101*d2c1+Y102*d2c2+Y103*d2c3) +
       b1*(Y110*d2c0+Y111*d2c1+Y112*d2c2+Y113*d2c3) +
       b2*(Y120*d2c0+Y121*d2c1+Y122*d2c2+Y123*d2c3) +
       b3*(Y130*d2c0+Y131*d2c1+Y132*d2c2+Y133*d2c3))+
      a2 *
      (b0*(Y200*d2c0+Y201*d2c1+Y202*d2c2+Y203*d2c3) +
       b1*(Y210*d2c0+Y211*d2c1+Y212*d2c2+Y213*d2c3) +
       b2*(Y220*d2c0+Y221*d2c1+Y222*d2c2+Y223*d2c3) +
       b3*(Y230*d2c0+Y231*d2c1+Y232*d2c2+Y233*d2c3))+
      a3 *
      (b0*(Y300*d2c0+Y301*d2c1+Y302*d2c2+Y303*d2c3) +
       b1*(Y310*d2c0+Y311*d2c1+Y312*d2c2+Y313*d2c3) +
       b2*(Y320*d2c0+Y321*d2c1+Y322*d2c2+Y323*d2c3) +
       b3*(Y330*d2c0+Y331*d2c1+Y332*d2c2+Y333*d2c3));
    return lapf;
  }


};
#endif

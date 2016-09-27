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
    
    



#include "Numerics/Spline3D/TriCubicSpline.h"
#include <iostream>


void TriCubicSpline::Update(bool iflap)
{
  UpdateX(0, 1);  // Do dF/dx
  UpdateY(0, 2);  // Do dF/dy
  UpdateZ(0, 3);  // Do dF/dz
  UpdateY(1, 4);  // Do d2F/dxdy
  UpdateZ(1, 5);  // Do d2F/dxdz
  UpdateZ(2, 6);  // Do d2F/dydz
  UpdateZ(4, 7);  // Do d3F/dxdydz
  if(iflap)
  {
    D2F.resize(n_x,n_y,n_z);
    /// calculate D2F at all points.
    D2FDR2();
    UpdateD2X(0, 1);    // Do d (D2F) /dx
    UpdateD2Y(0, 2);    // Do d (D2F) /dy
    UpdateD2Z(0, 3);    // Do d (D2F) /dy
    UpdateD2Y(1, 4);    // Do d2 (D2F) /dxdy
    UpdateD2Z(1, 5);    // Do d2 (D2F) /dxdz
    UpdateD2Z(2, 6);    // Do d2 (D2F) /dydz
    UpdateD2Z(4, 7);    // Do d3 (D2F) /dxdydz
  }
  UpToDate = true;
}

void TriCubicSpline::UpdateX(int source, int dest)
{
  blitz::Array<double,1> mu(n_x);
  /// Loop over all y and z
  for(int iy = 0; iy < n_y; iy++)
  {
    for(int iz = 0; iz < n_z; iz++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      F(0,iy,iz)[dest] = 1.5 *
                         ( F(1,iy,iz)[source] - F(0,iy,iz)[source] ) / dx(0,0);
      F(n_x-1,iy,iz)[dest] = 1.5 *
                             (F(n_x-1,iy,iz)[source] - F(n_x-2,iy,iz)[source]) / dx(0,n_x-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_x - 1; j++)
      {
        double lambda = 0.5*dx(0,j)/(dx(0,j)+dx(0,j-1));
        mu(j) = 0.5 - lambda;
        F(j,iy,iz)[dest] =
          3.0*(lambda*(F(j,iy,iz)[source]-F(j-1,iy,iz)[source])/dx(0,j-1)+
               mu(j) *(F(j+1,iy,iz)[source]-F(j,iy,iz)[source])/dx(0,j) );
        double cj = 1.0 - lambda * mu(j-1);
        F(j,iy,iz)[dest] -= lambda * F(j-1,iy,iz)[dest];
        mu(j) /= cj;
        F(j,iy,iz)[dest] /= cj;
      }
      /// last element
      int j = n_x - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(j,iy,iz)[dest] =
        3.0*(lambda*(F(j,iy,iz)[source]-F(j-1,iy,iz)[source])/dx(0,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      F(j,iy,iz)[dest] -= lambda * F(j-1,iy,iz)[dest];
      mu(j) /= cj;
      F(j,iy,iz)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_x-2; j>=0; j--)
        F(j,iy,iz)[dest] -= mu(j) * F(j+1,iy,iz)[dest];
    }
  }
}

void TriCubicSpline::UpdateY(int source, int dest)
{
  blitz::Array<double,1> mu(n_y);
  /// Loop over all x and z
  for(int ix = 0; ix < n_x; ix++)
  {
    for(int iz = 0; iz < n_z; iz++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      F(ix,0,iz)[dest] = 1.5 *
                         ( F(ix,1,iz)[source] - F(ix,0,iz)[source] ) / dx(1,0);
      F(ix,n_y-1,iz)[dest] = 1.5 *
                             (F(ix,n_y-1,iz)[source] - F(ix,n_y-2,iz)[source]) / dx(1,n_y-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_y - 1; j++)
      {
        double lambda = 0.5*dx(1,j)/(dx(1,j)+dx(1,j-1));
        mu(j) = 0.5 - lambda;
        F(ix,j,iz)[dest] =
          3.0*(lambda*(F(ix,j,iz)[source]-F(ix,j-1,iz)[source])/dx(1,j-1)+
               mu(j) *(F(ix,j+1,iz)[source]-F(ix,j,iz)[source])/dx(1,j) );
        double cj = 1.0 - lambda * mu(j-1);
        F(ix,j,iz)[dest] -= lambda * F(ix,j-1,iz)[dest];
        mu(j) /= cj;
        F(ix,j,iz)[dest] /= cj;
      }
      /// last element
      int j = n_y - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,j,iz)[dest] =
        3.0*(lambda*(F(ix,j,iz)[source]-F(ix,j-1,iz)[source])/dx(1,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,j,iz)[dest] -= lambda * F(ix,j-1,iz)[dest];
      mu(j) /= cj;
      F(ix,j,iz)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_y-2; j>=0; j--)
        F(ix,j,iz)[dest] -= mu(j) * F(ix,j+1,iz)[dest];
    }
  }
}

void TriCubicSpline::UpdateZ(int source, int dest)
{
  blitz::Array<double,1> mu(n_z);
  /// Loop over all x and y
  for(int ix = 0; ix < n_x; ix++)
  {
    for(int iy = 0; iy < n_y; iy++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      F(ix,iy,0)[dest] = 1.5 *
                         ( F(ix,iy,1)[source] - F(ix,iy,0)[source] ) / dx(2,0);
      F(ix,iy,n_z-1)[dest] = 1.5 *
                             (F(ix,iy,n_z-1)[source] - F(ix,iy,n_z-2)[source]) / dx(2,n_z-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_z - 1; j++)
      {
        double lambda = 0.5*dx(2,j)/(dx(2,j)+dx(2,j-1));
        mu(j) = 0.5 - lambda;
        F(ix,iy,j)[dest] =
          3.0*(lambda*(F(ix,iy,j)[source]-F(ix,iy,j-1)[source])/dx(2,j-1)+
               mu(j) *(F(ix,iy,j+1)[source]-F(ix,iy,j)[source])/dx(2,j) );
        double cj = 1.0 - lambda * mu(j-1);
        F(ix,iy,j)[dest] -= lambda * F(ix,iy,j-1)[dest];
        mu(j) /= cj;
        F(ix,iy,j)[dest] /= cj;
      }
      /// last element
      int j = n_z - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      F(ix,iy,j)[dest] =
        3.0*(lambda*(F(ix,iy,j)[source]-F(ix,iy,j-1)[source])/dx(2,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      F(ix,iy,j)[dest] -= lambda * F(ix,iy,j-1)[dest];
      mu(j) /= cj;
      F(ix,iy,j)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_z-2; j>=0; j--)
        F(ix,iy,j)[dest] -= mu(j) * F(ix,iy,j+1)[dest];
    }
  }
  return;
}

void TriCubicSpline::UpdateD2X(int source, int dest)
{
  blitz::Array<double,1> mu(n_x);
  /// Loop over all y and z
  for(int iy = 0; iy < n_y; iy++)
  {
    for(int iz = 0; iz < n_z; iz++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      D2F(0,iy,iz)[dest] = 1.5 *
                           ( D2F(1,iy,iz)[source] - D2F(0,iy,iz)[source] ) / dx(0,0);
      D2F(n_x-1,iy,iz)[dest] = 1.5 *
                               (D2F(n_x-1,iy,iz)[source] - D2F(n_x-2,iy,iz)[source]) / dx(0,n_x-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_x - 1; j++)
      {
        double lambda = 0.5*dx(0,j)/(dx(0,j)+dx(0,j-1));
        mu(j) = 0.5 - lambda;
        D2F(j,iy,iz)[dest] =
          3.0*(lambda*(D2F(j,iy,iz)[source]-D2F(j-1,iy,iz)[source])/dx(0,j-1)+
               mu(j) *(D2F(j+1,iy,iz)[source]-D2F(j,iy,iz)[source])/dx(0,j) );
        double cj = 1.0 - lambda * mu(j-1);
        D2F(j,iy,iz)[dest] -= lambda * D2F(j-1,iy,iz)[dest];
        mu(j) /= cj;
        D2F(j,iy,iz)[dest] /= cj;
      }
      /// last element
      int j = n_x - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      D2F(j,iy,iz)[dest] =
        3.0*(lambda*(D2F(j,iy,iz)[source]-D2F(j-1,iy,iz)[source])/dx(0,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      D2F(j,iy,iz)[dest] -= lambda * D2F(j-1,iy,iz)[dest];
      mu(j) /= cj;
      D2F(j,iy,iz)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_x-2; j>=0; j--)
        D2F(j,iy,iz)[dest] -= mu(j) * D2F(j+1,iy,iz)[dest];
    }
  }
}

void TriCubicSpline::UpdateD2Y(int source, int dest)
{
  blitz::Array<double,1> mu(n_y);
  /// Loop over all x and z
  for(int ix = 0; ix < n_x; ix++)
  {
    for(int iz = 0; iz < n_z; iz++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      D2F(ix,0,iz)[dest] = 1.5 *
                           ( D2F(ix,1,iz)[source] - D2F(ix,0,iz)[source] ) / dx(1,0);
      D2F(ix,n_y-1,iz)[dest] = 1.5 *
                               (D2F(ix,n_y-1,iz)[source] - D2F(ix,n_y-2,iz)[source]) / dx(1,n_y-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_y - 1; j++)
      {
        double lambda = 0.5*dx(1,j)/(dx(1,j)+dx(1,j-1));
        mu(j) = 0.5 - lambda;
        D2F(ix,j,iz)[dest] =
          3.0*(lambda*(D2F(ix,j,iz)[source]-D2F(ix,j-1,iz)[source])/dx(1,j-1)+
               mu(j) *(D2F(ix,j+1,iz)[source]-D2F(ix,j,iz)[source])/dx(1,j) );
        double cj = 1.0 - lambda * mu(j-1);
        D2F(ix,j,iz)[dest] -= lambda * D2F(ix,j-1,iz)[dest];
        mu(j) /= cj;
        D2F(ix,j,iz)[dest] /= cj;
      }
      /// last element
      int j = n_y - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      D2F(ix,j,iz)[dest] =
        3.0*(lambda*(D2F(ix,j,iz)[source]-D2F(ix,j-1,iz)[source])/dx(1,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      D2F(ix,j,iz)[dest] -= lambda * D2F(ix,j-1,iz)[dest];
      mu(j) /= cj;
      D2F(ix,j,iz)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_y-2; j>=0; j--)
        D2F(ix,j,iz)[dest] -= mu(j) * D2F(ix,j+1,iz)[dest];
    }
  }
}

void TriCubicSpline::UpdateD2Z(int source, int dest)
{
  blitz::Array<double,1> mu(n_z);
  /// Loop over all x and y
  for(int ix = 0; ix < n_x; ix++)
  {
    for(int iy = 0; iy < n_y; iy++)
    {
      /// set up tridiagonal set of equations : initialise RHS
      D2F(ix,iy,0)[dest] = 1.5 *
                           ( D2F(ix,iy,1)[source] - D2F(ix,iy,0)[source] ) / dx(2,0);
      D2F(ix,iy,n_z-1)[dest] = 1.5 *
                               (D2F(ix,iy,n_z-1)[source] - D2F(ix,iy,n_z-2)[source]) / dx(2,n_z-2);
      mu(0) = 0.5;
      /// solve tri-diagonal set of equations. First eliminate lower elements
      for(int j = 1; j < n_z - 1; j++)
      {
        double lambda = 0.5*dx(2,j)/(dx(2,j)+dx(2,j-1));
        mu(j) = 0.5 - lambda;
        D2F(ix,iy,j)[dest] =
          3.0*(lambda*(D2F(ix,iy,j)[source]-D2F(ix,iy,j-1)[source])/dx(2,j-1)+
               mu(j) *(D2F(ix,iy,j+1)[source]-D2F(ix,iy,j)[source])/dx(2,j) );
        double cj = 1.0 - lambda * mu(j-1);
        D2F(ix,iy,j)[dest] -= lambda * D2F(ix,iy,j-1)[dest];
        mu(j) /= cj;
        D2F(ix,iy,j)[dest] /= cj;
      }
      /// last element
      int j = n_z - 1;
      double lambda = 0.5;
      mu(j) = 0.5 - lambda;
      D2F(ix,iy,j)[dest] =
        3.0*(lambda*(D2F(ix,iy,j)[source]-D2F(ix,iy,j-1)[source])/dx(2,j-1) );
      double cj = 1.0 - lambda * mu(j-1);
      D2F(ix,iy,j)[dest] -= lambda * D2F(ix,iy,j-1)[dest];
      mu(j) /= cj;
      D2F(ix,iy,j)[dest] /= cj;
      ///  now last d/dx is correct, so proceed upward, back substituting
      for (j=n_z-2; j>=0; j--)
        D2F(ix,iy,j)[dest] -= mu(j) * D2F(ix,iy,j+1)[dest];
    }
  }
}

void TriCubicSpline::D2FDR2()
{
  D2F = 0.0;
  for(int ix = 1; ix < n_x - 1; ix++)
  {
    for(int iy = 1; iy < n_y - 1; iy++)
    {
      for(int iz = 1; iz < n_z - 1; iz++)
      {
        posvec_t r = m_grid->ptr(ix,iy,iz);
        D2F(ix,iy,iz)[0] = laplacian(r);
      }
    }
  }
  return;
}

double TriCubicSpline::delsqf(const gridvec_t& P0, int dir)
{
  double numer, denom;
  gridvec_t Pi,Pf;
  for(int idim = 0; idim < 3; idim++)
  {
    Pi[idim] = P0[idim];
    Pf[idim] = P0[idim];
  }
  Pi[dir] -= 1;
  Pf[dir] += 1;
  double f0 = F(P0[0],P0[1],P0[2])[0];
  double fi = F(Pi[0],Pi[1],Pi[2])[0];
  double ff = F(Pf[0],Pf[1],Pf[2])[0];
  double hi = dx(dir,Pi[dir]);
  double h0 = dx(dir,P0[dir]);
  numer = hi*(ff-f0) + h0*(fi-f0);
  denom = hi*h0*(hi+h0);
  return 2.0*numer/(denom+1e-30);
}


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
    
    



#include "Numerics/Spline3D/Grid1D.h"
#include <iostream>



void Grid1D::init(int nsections,
                  int npts,
                  double xi,
                  const std::vector<int>& nrho,
                  const std::vector<double>& dh)
{
  m_size = npts;
  m_coord.resize(m_size);
  m_sections = nsections;
  m_start = xi;
  /// assign and initialise arrays
  m_h.resize(m_size);
  m_M.resize(m_size+1);
  m_M[0] = 0;
  for(int isec = 0; isec < m_sections; isec++)
  {
    m_h[isec] = dh[isec];
    m_M[isec+1] = m_M[isec] + nrho[isec];
  }
  /// compute the coordinates
  int ipt = 0;
  m_coord[0] = m_start;
  for(int isec = 0; isec < m_size; isec++)
  {
    for(int ix = m_M[isec]; ix < m_M[isec+1]; ix++)
    {
      ipt++;
      m_coord[ipt] = m_coord[ipt-1] + m_h[isec];
    }
  }
  m_end = m_coord[m_size-1];
  return;
}

void Grid1D::init(int xi,
                  int xf,
                  const Grid1D& agrid)
{
  m_size = xf - xi + 1;
  m_coord.resize(m_size);
  int si, sf;   /// the indices of the initial and final sections
  bool flag1 = true, flag2 = true;
  for(int isec = 0; isec < agrid.m_sections; isec++)
  {
    if(flag1 && xi <= agrid.m_M[isec+1])
    {
      si = isec;
      flag1 = false;
    }
    if(flag2 && xf <= agrid.m_M[isec+1])
    {
      sf = isec;
      flag2 = false;
    }
  }
  m_sections = sf - si + 1;
  /// assign arrays
  m_M.resize(m_size+1);
  m_h.resize(m_size);
  m_M[0] = 0;
  for(int isec = 1; isec < m_sections; isec++)
    m_M[isec] = agrid.m_M[isec+si] - xi;
  m_M[m_sections] = xf - xi;
  /// the grid spacings
  for(int isec = 0; isec < m_sections; isec++)
    m_h[isec] = agrid.m_h[si+isec];
  /// compute the coordinates
  int ipt = 0;
  m_coord[0] = agrid.m_coord[xi];
  for(int isec = 0; isec < m_size; isec++)
  {
    for(int ix = m_M[isec]; ix < m_M[isec+1]; ix++)
    {
      ipt++;
      m_coord[ipt] = m_coord[ipt-1] + m_h[isec];
    }
  }
  m_start = m_coord[0];
  m_end = m_coord[m_size-1];
  return;
}

int Grid1D::xl(double x)
{
  /// the start and end points
  int a = 0;
  int b = m_size-1;
  while( b - a > 1 )
  {
    int m = a + ( b - a )/2 ; // divide by 2^{1}
    if( x >= m_coord[m] )
    {
      a = m;
    }
    else
    {
      b = m;
    }
  }
  return a;
}

int Grid1D::xn(double x)
{
  int i = xl(x);       /// get the lowest Grid point
  if( x - m_coord[i] > 0.5*h(i) )
    i += 1;
  return i;
}

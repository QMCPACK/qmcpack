//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef GUARD_UGRID1D_H
#define GUARD_UGRID1D_H


#include <vector>

class uGrid1D
{

public:

  /// total number of Grid points : 0 to m_size - 1
  int m_size;

  /// coordinate of first and last grid points
  double m_start;
  double m_end;

  /// Grid width
  double m_h;

  /// coordinates
  std::vector<double> m_coord;

  /// Constructor
  uGrid1D() {}

  /// initialise with  intervals
  void init(int npts, double dx, double x0)
  {
    m_size = npts;
    m_h = dx;
    m_start = x0;
    m_coord.resize(m_size);
    m_coord[0] = m_start;
    for(int i = 1; i < m_size; i++)
      m_coord[i] = m_coord[i-1] + m_h;
    m_end = m_coord[m_size-1];
  }

  /// initialise with start and end.
  void init(double xi, double xf, int npts)
  {
    m_size = npts;
    m_start = xi;
    m_end = xf;
    m_h = ( m_end - m_start )/ ( m_size - 1 );
    m_coord.resize(m_size);
    m_coord[0] = m_start;
    for(int i = 1; i < m_size; i++)
      m_coord[i] = m_coord[i-1] + m_h;
  }

  /// initialise from part of a supplied Grid
  void init(int ix, int fx, const uGrid1D& agrid)
  {
    m_size = fx - ix + 1;
    m_h = agrid.m_h;
    m_start = agrid.m_coord[ix];
    m_end = agrid.m_coord[fx];
    m_coord.resize(m_size);
    m_coord[0] = m_start;
    for(int i = 1; i < m_size; i++)
      m_coord[i] = m_coord[i-1] + m_h;
  }

  /// the lowest Grid-point for coordinate x
  int xl(double x)
  {
    return int((x - m_start)/m_h);
  }

  /// the nearest Grid-point for coordinate x
  int xn(double x)
  {
    int ipt = xl(x);
    if( x - m_coord[ipt] > 0.5 * m_h )
    {
      return ipt + 1;
    }
    return ipt;
  }

  /// the Grid width for coordinate x
  inline double h(double x)
  {
    return m_h;
  }

  /// the Grid width for some Grid point i
  inline double h(int i)
  {
    return m_h;
  }

  /// the size of the Grid :
  inline int Size()
  {
    return m_size;
  }

};
#endif

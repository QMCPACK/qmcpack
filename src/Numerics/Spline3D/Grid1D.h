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
    
    



#ifndef GUARD_GRID1D_H
#define GUARD_GRID1D_H


#include <vector>

class Grid1D
{

public:

  /// number of different Grid density sections (regions)
  int m_sections;

  /// total number of Grid points : 0 to m_size - 1
  int m_size;

  /// coordinate of first and last grid points
  double m_start;
  double m_end;

  /// coordinates
  std::vector<double> m_coord;

  /// the interval points
  std::vector<int> m_M;

  /// Grid widths
  std::vector<double> m_h;

  /// Constructor
  Grid1D() {}

  /// initialise with the entire data from input
  void init(int, int, double, const std::vector<int>&,
            const std::vector<double>&);

  /// initialise from part of a supplied Grid
  void init(int, int, const Grid1D&);


  /// the lowest Grid-point for coordinate x
  int xl(double);

  /// the nearest Grid-point for coordinate x
  int xn(double);

  /// the Grid width for coordinate x
  inline double h(double x)
  {
    return h(xl(x));
  }

  /// the Grid width for some Grid point i
  inline double h(int i)
  {
    return m_coord[i+1] - m_coord[i];
  }

  /// the size of the Grid :
  inline int Size()
  {
    return m_size;
  }

};
#endif

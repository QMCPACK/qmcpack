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


/** @file OneDimNumGridFunctor.h
 * @brief Definition of OneDimNumGridFunctor
 */
#ifndef QMCPLUSPLUS_ONEDIMNUMGRIDFUNCTOR_H
#define QMCPLUSPLUS_ONEDIMNUMGRIDFUNCTOR_H
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
/** adaptor class to handle a temporary OneDimCubicSpline on a numerical grid.
 */
template<class T>
struct OneDimNumGridFunctor
{
  using GridType = NumericalGrid<T>;
  using FuncType = OneDimCubicSpline<T>;

  FuncType myFunc;

  inline OneDimNumGridFunctor() { myFunc.m_grid = std::make_unique<GridType>(); }

  inline T r(int i) const { return (*myFunc.m_grid)[i]; }
  inline T operator()(int i) const { return myFunc(i); }

  inline T rmax() const { return myFunc.m_grid->rmax(); }

  inline T rmin() const { return myFunc.m_grid->rmin(); }

  inline T splint(T r) { return myFunc.splint(r); }

  void put(int npoints, std::istream& fin)
  {
    auto& myGrid = dynamic_cast<GridType&>(*myFunc.m_grid);
    myGrid.resize(npoints);
    myFunc.resize(npoints);
    for (int j = 0; j < npoints; j++)
    {
      fin >> myGrid(j) >> myFunc(j);
    }
    myGrid.set(myGrid(0), myGrid(npoints - 1), npoints);
    T yprime_i = (myFunc(1) - myFunc(0)) / (myGrid(1) - myGrid(0));
    myFunc.spline(0, yprime_i, npoints - 1, 0.0);
  }

  bool put(xmlNodePtr cur) { return true; }
};
} // namespace qmcplusplus
#endif

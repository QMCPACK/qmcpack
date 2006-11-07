//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_GRID_FUNCTOR_LINEAR_SPLINE_H
#define QMCPLUSPLUS_GRID_FUNCTOR_LINEAR_SPLINE_H

#include "Numerics/OneDimGridFunctor.h"

namespace qmcplusplus {
/** Perform One-Dimensional Cubic Spline Interpolation with fixed first-derivatives at the ends
 *
 * Using m-relationship and the first-order derivaties.
 * Each funtor checks the bounds r_min and r_max.
 */
template <class Td, 
	  class Tg = Td, 
	  class CTd= Vector<Td>,
	  class CTg= Vector<Tg> >
class OneDimLinearSpline: public OneDimGridFunctor<Td,Tg,CTd,CTg> {

public:

  typedef OneDimGridFunctor<Td,Tg,CTd,CTg> base_type;
  typedef typename base_type::value_type  value_type;
  typedef typename base_type::point_type  point_type;
  typedef typename base_type::data_type data_type;
  typedef typename base_type::grid_type grid_type;

  using base_type::GridManager;
  using base_type::m_grid;
  using base_type::Y;
  using base_type::dY;
  using base_type::d2Y;
  using base_type::m_Y;
  using base_type::FirstAddress;

  data_type m_Y1;
  int First,Last;
  point_type r_min, r_max;

  OneDimLinearSpline(grid_type* gt = 0): base_type(gt){ }

  template<class VV>
  OneDimLinearSpline(grid_type* gt, const VV& nv, bool pbc=true): base_type(gt)
  {
    m_Y.resize(nv.size());
    std::copy(nv.begin(), nv.end(), m_Y.data());
  }

  /** evaluate the value at r
   * @param r value on a grid
   * @return value obtained by cubic-spline
   */
  inline value_type splint(point_type r) {
    if(r>=r_max) return 0.0;
    int k = m_grid->getIndex(r);
    //m_grid->locate(r);
    //int k=m_grid->currentIndex();
    point_type dr=r-(*m_grid)[k];
    return m_Y[k]+m_Y1[k]*dr;
  }

  /** evaluate the value at r
   * @param r value on a grid
   * @param du first derivative (assigned)
   * @param d2u second derivative (assigned)
   * @return value obtained by cubic-spline
   */
  inline value_type 
  splint(point_type r, value_type& du, value_type& d2u) {
    cerr << "  OneDimLinearSpline cannot be used for derivates." << endl;
    return 0.0;
  }

  /** evaluate the spline coefficients 
   * @param imin index of the first valid grid
   * @param yp1 first derivative at imin grid point
   * @param imax index of the last valid grid
   * @param ypn first derivative at imax grid point
   *
   * Use m-relation to evalaute the spline coefficients on [imin,imax] grid points.
   */
  inline 
  void spline(int imin, value_type yp1, int imax, value_type ypn) {
    int npts(imax-imin+1);
    First=imin; Last=imax;
    m_Y1.resize(npts);
    r_min=m_grid->r(imin);
    r_max=m_grid->r(imax);
    for(int i=0; i<imax-1; i++)
    {
      m_Y1[i]= (m_Y[i+1]-m_Y[i])/((*m_grid)[i+1]-(*m_grid)[i]);
    }
  }

  inline void spline() {
    spline(0,0.0,this->size()-1,0.0);
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

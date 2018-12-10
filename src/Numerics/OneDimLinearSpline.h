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
    
    



#ifndef QMCPLUSPLUS_GRID_FUNCTOR_LINEAR_SPLINE_H
#define QMCPLUSPLUS_GRID_FUNCTOR_LINEAR_SPLINE_H

#include "Numerics/OneDimGridFunctor.h"
//#define USE_MEMORYSAVEMODE

namespace qmcplusplus
{
/** Perform One-Dimensional linear spline Interpolation.
 *
 * Only valid with linear grid.
 * @todo Have to prevent OneDimLinearSpline<T> being used with other than
 * linear grid!!
 */
template <class Td,
         class Tg = Td,
         class CTd= Vector<Td>,
         class CTg= Vector<Tg> >
class OneDimLinearSpline: public OneDimGridFunctor<Td,Tg,CTd,CTg>
{

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
  //using base_type::FirstAddress;

  int First;
  int Last;
  value_type ConstValue;
  point_type r_min;
  point_type r_max;
  point_type delta;
  point_type delta_inv;
  data_type m_Y1;

  OneDimLinearSpline(grid_type* gt = 0): base_type(gt),r_min(0),r_max(0)
  {
    if(gt)
    {
      r_min=gt->rmin();
      r_max=gt->rmax();
    }
  }

  OneDimLinearSpline(point_type ri, point_type rf): base_type(0), r_min(ri), r_max(rf)
  {
  }

  template<class VV>
  OneDimLinearSpline(grid_type* gt, const VV& nv, bool pbc=false): base_type(gt)
  {
    if(gt)
    {
      r_min=gt->rmin();
      r_max=gt->rmax();
    }
    assign(nv.begin(), nv.end());
  }

  OneDimLinearSpline<Td,Tg,CTd,CTg>* makeClone() const
  {
    return new OneDimLinearSpline<Td,Tg,CTd,CTg>(*this);
  }

  template<class IT>
  void assign(IT d_first, IT d_last)
  {
    m_Y.resize(d_last-d_first);
    copy(d_first,d_last,m_Y.data());
    delta=(r_max-r_min)/static_cast<point_type>(m_Y.size()-1);
    //r_min=m_grid->rmin();
    //r_max=m_grid->rmax();
    //delta=m_grid->dh();
    delta_inv=1.0/delta;
  }


  inline point_type rmax() const
  {
    return r_max;
  }

  /** evaluate the value at r
   * @param r value on a grid
   * @return value obtained by cubic-spline
   *
   * Performance may be tunned: define USE_MEMORYSAVEMODE
   * to evaluate the coefficients instead of using aux. arrays
   */
  inline value_type splint(point_type r)
  {
    if(r>=r_max)
      return ConstValue;
    int k = static_cast<int>((r-r_min)*delta_inv);
#if defined(USE_MEMORYSAVEMODE)
    return m_Y[k]+(m_Y[k+1]-m_Y[k])*(r*delta_inv-k);
#else
    return m_Y[k]+m_Y1[k]*(r-(*m_grid)[k]);
#endif
  }

  //template<class IT1, class IT2>
  //void assign(IT1 g_first, IT1 g_last, IT2 d_first, IT2 d_last)
  //{
  //  if(m_grid ==0)
  //  {
  //    NumericalGrid<Td> *agrid=new NumericalGrid<Td>;
  //    agrid->assign(g_first,g_last);
  //    m_grid=agrid;
  //  }
  //  assign(d_first,d_last);
  //}


//  /** evaluate the value at r using a binary search on a grid
//   * @param r distance
//   */
//  inline value_type splintNG(point_type r) const
//  {
//    if(r>=r_max) return ConstValue;
//    int k=m_grid->getIndex(r);
//    //int k = static_cast<int>((r-r_min)*delta_inv);
//#if defined(USE_MEMORYSAVEMODE)
//    return m_Y[k]+(m_Y[k+1]-m_Y[k])*(r*delta_inv-k);
//#else
//    return m_Y[k]+m_Y1[k]*(r-(*m_grid)[k]);
//#endif
//  }

  /** evaluate the index and the linear coefficient
   * @param r distance
   * @param k return index
   * @param rfrac (r-floor(r/delta))/delta
   */
  inline void locate(point_type r, int &k, point_type& rfrac)
  {
    k=static_cast<int>((r-r_min)*delta_inv);
    rfrac=r*delta_inv-k;
  }

  /** evaluate the value at r=(k + rfrac)*delta
   */
  inline value_type f(int k, point_type rfrac)
  {
    return m_Y[k]*(1.0-rfrac)+m_Y[k+1]*rfrac;
  }

  /** evaluate the value at r
   * @param r value on a grid
   * @param du first derivative (assigned)
   * @param d2u second derivative (assigned)
   * @return value obtained by cubic-spline
   */
  inline value_type
  splint(point_type r, value_type& du, value_type& d2u)
  {
    std::cerr << "  OneDimLinearSpline cannot be used for derivates." << std::endl;
    return 0.0;
  }

  /** evaluate the spline coefficients
   * @param imin index of the first valid grid
   * @param yp1 first derivative at imin grid point
   * @param imax index of the last valid grid
   * @param ypn first derivative at imax grid point
   *
   */
  inline
  void spline(int imin, value_type yp1, int imax, value_type ypn)
  {
    int npts(imax-imin+1);
    First=imin;
    Last=imax;
    m_Y1.resize(npts);
    //r_min=m_grid->r(imin);
    //r_max=m_grid->r(imax);
    for(int i=imin; i<imax-1; i++)
    {
      m_Y1[i]= (m_Y[i+1]-m_Y[i])/((*m_grid)[i+1]-(*m_grid)[i]);
      //m_Y1[i]= (m_Y[i+1]-m_Y[i])*delta_inv;
    }
    m_Y1[imax]=0.0;
    ConstValue=m_Y[imax];
  }

  inline void spline()
  {
    spline(0,0.0,this->size()-1,0.0);
  }
};

}
#endif

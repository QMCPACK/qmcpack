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


#ifndef QMCPLUSPLUS_GRID_FUNCTOR_QUINTIC_SPLINE_H
#define QMCPLUSPLUS_GRID_FUNCTOR_QUINTIC_SPLINE_H

#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/SplineSolvers.h"

namespace qmcplusplus
{
/*
 * Perform One-Dimensional Quintic Spline Interpolation.
 */

template<class Td, class Tg = Td, class CTd = Vector<Td>, class CTg = Vector<Tg>>
class OneDimQuinticSpline : public OneDimGridFunctor<Td, Tg, CTd, CTg>
{
public:
  using base_type  = OneDimGridFunctor<Td, Tg, CTd, CTg>;
  using value_type = typename base_type::value_type;
  using point_type = typename base_type::point_type;
  using data_type  = typename base_type::data_type;
  using grid_type  = typename base_type::grid_type;

  using base_type::d2Y;
  using base_type::dY;
  using base_type::m_grid;
  using base_type::m_Y;
  using base_type::Y;

  data_type m_Y2;
  data_type B;
  data_type D;
  data_type E;
  data_type F;

  using base_type::NumNodes;

  point_type r_min;
  point_type r_max;
  value_type first_deriv;
  value_type last_deriv;
  value_type ConstValue;

  OneDimQuinticSpline(std::unique_ptr<grid_type> gt = std::unique_ptr<grid_type>()) : base_type(std::move(gt)) {}

  template<class VV>
  OneDimQuinticSpline(std::unique_ptr<grid_type> gt, const VV& nv)
      : base_type(std::move(gt)), first_deriv(0.0), last_deriv(0.0)
  {
    int n = nv.size();
    m_Y.resize(nv.size());
    m_Y2.resize(nv.size());
    std::copy(nv.begin(), nv.end(), m_Y.data());
    B.resize(n);
    D.resize(n);
    E.resize(n);
    F.resize(n);
  }

  void set(Vector<Td>& data)
  {
    int n = data.size();
    m_Y.resize(n);
    m_Y = data;
    m_Y2.resize(n);
    B.resize(n);
    D.resize(n);
    E.resize(n);
    F.resize(n);
  }

  OneDimQuinticSpline<Td, Tg, CTd, CTg>* makeClone() const { return new OneDimQuinticSpline<Td, Tg, CTd, CTg>(*this); }

  OneDimQuinticSpline<Td, Tg, CTd, CTg>(const OneDimQuinticSpline<Td, Tg, CTd, CTg>& a)
      : OneDimGridFunctor<Td, Tg, CTd, CTg>(a)
  {
    m_Y2.resize(a.m_Y2.size());
    m_Y2        = a.m_Y2;
    ConstValue  = a.ConstValue;
    r_min       = a.r_min;
    r_max       = a.r_max;
    first_deriv = a.first_deriv;
    last_deriv  = a.last_deriv;
    B.resize(a.B.size());
    B = a.B;
    D.resize(a.D.size());
    D = a.D;
    E.resize(a.E.size());
    E = a.E;
    F.resize(a.F.size());
    F = a.F;
  }

private:
  inline Td quinticInterpolate(Td cL, Td a, Td b, Td c, Td d, Td e, Td f) const
  {
    return a + cL * (b + cL * (c + cL * (d + cL * (e + cL * f))));
  }

  inline Td quinticInterpolateSecondDeriv(Td cL, Td a, Td b, Td c, Td d, Td e, Td f, Td& du, Td& d2u) const
  {
    du  = b + cL * (2.0 * c + cL * (3.0 * d + cL * (4.0 * e + cL * f * 5.0)));
    d2u = 2.0 * c + cL * (6.0 * d + cL * (12.0 * e + cL * f * 20.0));
    return a + cL * (b + cL * (c + cL * (d + cL * (e + cL * f))));
  }

  inline Td quinticInterpolateThirdDeriv(Td cL, Td a, Td b, Td c, Td d, Td e, Td f, Td& du, Td& d2u, Td& d3u) const
  {
    du  = b + cL * (2.0 * c + cL * (3.0 * d + cL * (4.0 * e + cL * f * 5.0)));
    d2u = 2.0 * c + cL * (6.0 * d + cL * (12.0 * e + cL * f * 20.0));
    d3u = 6.0 * d + cL * (24.0 * e + cL * f * 60.0);
    return a + cL * (b + cL * (c + cL * (d + cL * (e + cL * f))));
  }

public:
  inline value_type splint(point_type r) const override
  {
    if (r < r_min)
    {
      return m_Y[0] + first_deriv * (r - r_min);
    }
    else if (r >= r_max)
    {
      return ConstValue;
    }
    Td cL;
    int Loc = m_grid->getIndexAndDistanceFromGridPoint(r, cL);
    return quinticInterpolate(cL, m_Y[Loc], B[Loc], m_Y2[Loc], D[Loc], E[Loc], F[Loc]);
  }


  inline value_type splint(point_type r, value_type& du, value_type& d2u) const override
  {
    if (r < r_min)
    {
      return m_Y[0] + first_deriv * (r - r_min);
    }
    else if (r >= r_max)
    {
      return ConstValue;
    }
    Td cL;
    int Loc = m_grid->getIndexAndDistanceFromGridPoint(r, cL);
    return quinticInterpolateSecondDeriv(cL, m_Y[Loc], B[Loc], m_Y2[Loc], D[Loc], E[Loc], F[Loc], du, d2u);
  }

  inline value_type splint(point_type r, value_type& du, value_type& d2u, value_type& d3u) const
  {
    if (r < r_min)
    {
      return m_Y[0] + first_deriv * (r - r_min);
    }
    else if (r >= r_max)
    {
      return ConstValue;
    }
    Td cL;
    int Loc = m_grid->getIndexAndDistanceFromGridPoint(r, cL);
    return quinticInterpolateThirdDeriv(cL, m_Y[Loc], B[Loc], m_Y2[Loc], D[Loc], E[Loc], F[Loc], du, d2u, d3u);
  }

  inline void spline(int imin, value_type yp1, int imax, value_type ypn) override
  {
    first_deriv = yp1;
    last_deriv  = ypn;
    r_min       = m_grid->r(imin);
    r_max       = m_grid->r(imax);
    int npts(this->size());
    m_Y2.resize(npts);
    B.resize(npts);
    D.resize(npts);
    E.resize(npts);
    F.resize(npts);
    m_Y2 = 0.0;
    B    = 0.0;
    D    = 0.0;
    E    = 0.0;
    F    = 0.0;
    QuinticSplineSolve(npts - imin, m_grid->data() + imin, m_Y.data() + imin, B.data() + imin, m_Y2.data() + imin,
                       D.data() + imin, E.data() + imin, F.data() + imin);
    ConstValue = m_Y[imax];
  }

  inline void spline() override { spline(0, 0.0, m_grid->size() - 1, 0.0); }
};


} // namespace qmcplusplus
#endif

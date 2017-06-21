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
    
    


#ifndef OHMMS_NUCLEAR_LINEAR_TRANSFORMFUNCTION_H
#define OHMMS_NUCLEAR_LINEAR_TRANSFORMFUNCTION_H

#include <numeric>
#include "Numerics/OneDimGridFunctor.h"

/**@ingroup NumerovTransform
 *\brief Transform function for the radial Schroedinger equation on a
 *linear grid with a Nuclear potential \f$V(r) = -\frac{Z}{r}.\f$
 *
 *Implements a transformation for the Numerov algorithm:
 \f[
 \frac{d^2u_{nl}}{dr^2}+k^2(r) u_{nl} = 0,
 \f]
 *where
 \f[
 k^2(r)=2[\varepsilon-\frac{L(L+1)}{2r^2}-V(r)].
 \f]
*/
template<class SourceType>
struct NuclearLinearTransform
{

  typedef SourceType Source_t;
  typedef typename SourceType::value_type value_type;

  ///source function on a linear grid
  Source_t& V;

  ///number of nodes for given quantum numbers
  int NumNodes;

  ///reference energy \f$\varepsilon\f$
  value_type E;

  ///Nuclear charge
  value_type Z;

  ///parameter to determine the cusp condition
  value_type CuspParam;

  ///principal quantum number
  value_type N;

  ///angular momentum quantum number
  value_type L;

  ///LL \f$=l(l+1)/(2m_{eff})\f$
  value_type LL;

  /** Constructor for \f$u_{nl}\f$
   * \param v the source function on a grid
   * \param n the principal quantum number
   * \param l the anuglar quantum number
   * \param z the charge
   * \param meff the effective mass
   * \note For the Nuclear potential the number of nodes for
   the radial function with quantum numbers
   \f$n\f$ and \f$l\f$ is \f$n-l-1\f$.
  */
  NuclearLinearTransform(Source_t& v, int n, int l, value_type z,
                         value_type meff=1.0):
    V(v), E(0.0), N(n), CuspParam(z), Z(z)
  {
    NumNodes = v.getNumOfNodes();
    L = static_cast<value_type>(l);
    LL = 0.5*L*(L+1)/meff;
  }

  ///returns the starting valid index
  inline int first() const
  {
    return 0;
  }

  /**
   *\param i index of the first value of the linear grid
   *\param z0 \f$u_{nl}(r_0)\f$
   *\param z1 \f$u_{nl}(r_0 + dr)\f$
   *\return the derivate of the target function before transformation
   *\brief  Assign the first and second values using the cusp condition.
   *
   For \f$l=0\f$ in the limit \f$r\rightarrow 0\f$
   \f[
   u(r) \approx r(1-Zr),
   \f]
   while for positive l the condition is
   \f[
   u(r) \approx r^{l+1}.
   \f]
  */
  value_type setCusp(int i, value_type& z0, value_type& z1)
  {
    value_type r0 = V.r(i);
    value_type dr = V.dr(i);
    value_type deriv;
    if(L<std::numeric_limits<value_type>::epsilon())
    {
      z0 = pow(r0,L+1);
      deriv = static_cast<value_type>(L+1)*pow(r0,L);
    }
    else
    {
      z0 = r0*(1-Z*r0);
      deriv = 1.0-2*Z*r0;
    }
    z1 = z0 + deriv*dr;
    return deriv;
  }

  ///returns the number of nodes
  inline int nodes() const
  {
    return NumNodes;
  }

  ///return \f$k^2(r_i)=2[\varepsilon-L(L+1)/2r_i^2-V(r_i)]\f$
  inline value_type k2(int i)
  {
    value_type rinv = 1.0/V.r(i);
    return 2.0*(E-LL*rinv*rinv-V(i));
  }

  /*!\fn value_type convert(value_type y, value_type r)
   *\param y the value of \f$u_{nl}(r)\f$
   *\param r the grid value
   *\return The same value: no transformation is needed.
   */
  inline value_type convert(value_type y, value_type r)
  {
    return y;
  }

  ///reset the reference energy and turning point
  inline void reset(value_type e)
  {
    E = e;
  }
};
#endif

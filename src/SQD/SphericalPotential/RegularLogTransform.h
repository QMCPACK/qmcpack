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
    
    


#ifndef OHMMS_REGULAR_LOG_TRANSFORMFUNCTION_H
#define OHMMS_REGULAR_LOG_TRANSFORMFUNCTION_H

#include "Numerics/OneDimGridFunctor.h"

/**@ingroup NumerovTransform
 *\brief Transform function for the radial Schroedinger equation on a log grid
 *with a regular potential.
 *
 *For regular potentials, the asymptotic approximation for \f$u_{l}\f$
 *near the origin follows\f[\lim_{r\rightarrow 0} u(r) \approx r^{(l+1)}.\f]
 *An example for Regular potentials is a Harmonic potential,
 \f[
 V(r) = \frac{1}{2}\omega^2r^2.
 \f]
 *
 *Since the Numerov algorithm can only be integrated on a uniform grid,
 *it is necessary to make the following change of variables
 *\f[
 Z_{nl}(r) = \frac{u_{nl}(r)}{\sqrt{r}} \mbox{ and } y = log(r).
 \f]
 *The Transform function for the Numerov algorithm becomes
 \f[
 \frac{d^2Z_{nl}(y)}{dy^2} + k'^2(r)Z_{nl}(y) = 0,
 \f]
 where
 \f[
 k'^2(r) = 2r^2[\varepsilon -\frac{(L+\frac{1}{2})^2}{2r^2}-V(r)]
 \f]
*/
template<class SourceType>
struct RegularLogTransform
{

  typedef SourceType Source_t;
  typedef typename SourceType::value_type value_type;

  Source_t& V;

  ///number of nodes for given quantum numbers
  int NumNodes;

  ///reference energy \f$\varepsilon\f$
  value_type E;

  ///parameter to determine the cusp condition
  value_type CuspParam;

  ///principal quantum number
  value_type N;

  ///angular momentum quantum number
  value_type L;

///\f[LL=(l+\frac{1}{2})^2/(2m_{eff})\f]
  value_type LL;

  /** Constructor for \f$Z_{nl}\f$
   * \param v the source function on a grid
   * \param n the principal quantum number
   * \param l the anuglar quantum number
   * \param z the charge
   * \param meff the effective mass
   * \note For the Harmonic potential \f$V(r) =  \frac{1}{2}\omega^2r^2\f$
   * the number of nodes for the radial function with quantum numbers
   * \f$n\f$ and \f$l\f$ is \f$n\f$.
   */
  RegularLogTransform(Source_t& v, int n, int l,
                      value_type cuspparam,
                      value_type meff=1.0):
    V(v), E(0.0), N(n), CuspParam(cuspparam)
  {
    NumNodes = v.getNumOfNodes();
    L = static_cast<value_type>(l);
    LL = 0.5*(L+0.5)*(L+0.5)/meff;
  }

  ///returns the starting valid index
  inline int first() const
  {
    return 1;
  }

  /**
    *\param i index of the first value of the log grid
    *\param z0 \f$Z_{nl}(r_0)\f$
    *\param z1 \f$Z_{nl}(r_0 + dr)\f$
    *\return the derivate of the target function before transformation
    *\brief  Assign the first and second values using the cusp condition.
    *
    In the limit \f$r\rightarrow 0 \f$ \f[ Z_{nl}(r)
    \approx r^{(l+\frac{1}{2})}. \f]
    */
  value_type setCusp(int i, value_type& z0, value_type& z1)
  {
    value_type r0 = V.r(i);
    value_type dr = V.dr(i);
    z0 = pow(r0,L+0.5);
    z1 = z0 + (L+0.5)*pow(r0,L-0.5)*dr;
    return static_cast<value_type>(L+1)*pow(r0,L);
  }

  ///returns the number of nodes
  inline int nodes() const
  {
    return NumNodes;
  }

  /*return \f$ k'^2(r_i) = 2r_i^2[\varepsilon
    -(L+\frac{1}{2})^2/2r_i^2-V(r_i)] \f$
  */
  inline value_type k2(int i)
  {
    value_type rsq = pow(V.r(i),2.0);
    return 2.0*(rsq*(E-V(i))-LL); //2.0*rsq*(E-LL/rsq-V(i));
  }

  /*!\fn value_type convert(value_type y, value_type r)
   *\param y the value \f$Z_{nl}(r)\f$
   *\param r the grid value
   *\return  \f$u_{nl}(r) = \sqrt{r}Z_{nl}(r)\f$
   */
  inline value_type convert(value_type y, value_type r)
  {
    return y*sqrt(r);
  }

  ///reset the reference energy and turning point
  inline void reset(value_type e)
  {
    E = e;
  }
};

#endif

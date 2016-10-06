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
    
    



#ifndef OHMMS_UOVERRN_H
#define OHMMS_UOVERRN_H
#include <math.h>

/**class UOverRN
 * \brief The purpose of this class is to serve as
 a function wrapper for the template parameter Rin.
 *
 * Rin is a radial function (grid or analytic) which
 must include the member function
 <ul>
 <li> R = evaluate(r,rinv,dRdr,d2Rdr2)
 </ul>
 UOverRN returns
 \f[ U(r) = \frac{R(r)}{r^N} \f] and the derivatives
 \f[
 U' = \frac{R'}{r^N}-N\frac{R}{r^{N+1}}
 \f]
 \f[
 U'' = \frac{R''}{r^N} - \frac{2NR'}{r^{N+1}}
 + \frac{N(N+1)R}{r^{N+2}}
 \f]
 *
*/

template<class Rin>
struct UOverRN
{

  ///the inverse exponent
  int X;
  typedef typename Rin::value_type value_type;
  typedef typename Rin::point_type point_type;

  ///storage
  value_type Y, dY, d2Y;
  ///the radial function
  Rin& U;

  UOverRN(int x, Rin& u): X(x), U(u) { }

  ~UOverRN() { }

  inline void setgrid(value_type r)
  {
    U.setgrid(r);
  }

  inline value_type f(value_type r) const
  {
    value_type rinv = 1.0/r;
    evaluate(r,rinv);
    return Y;
  }


  inline value_type df(value_type r) const
  {
    value_type rinv = 1.0/r;
    evaluate(r,rinv);
    return dY;
  }

  inline void evaluate(value_type r, value_type rinv)
  {
    Y = evaluate(r,rinv,dY,d2Y);
  }

  inline value_type evaluate(value_type r, value_type rinv,
                             value_type & drnl,
                             value_type& d2rnl)
  {
    if(X)
    {
      value_type rinv_x = pow(rinv,X);
      value_type unl, dunl, d2unl;
      unl = U.evaluate(r,rinv,dunl,d2unl);
      drnl = rinv_x*(dunl-static_cast<value_type>(X)*unl*rinv);
      d2rnl = rinv_x*(d2unl-2.0*static_cast<value_type>(X)*dunl*rinv+
                      static_cast<value_type>(X)*static_cast<value_type>(X+1)*unl*rinv*rinv);
      return unl*rinv_x;
    }
    else
    {
      return U.evaluate(r,rinv,drnl,d2rnl);
    }
  }
};

#endif

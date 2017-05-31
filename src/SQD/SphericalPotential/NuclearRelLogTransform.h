//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_NUCLEAR_RELATIVISTIC_LOG_TRANSFORMFUNCTION_H
#define OHMMS_NUCLEAR_RELATIVISTIC_LOG_TRANSFORMFUNCTION_H

#include <numeric>
#include "Numerics/OneDimGridFunctor.h"

//#define PRINT_DEBUG
/**@ingroup NumerovTransform
 *\brief Transform function for the radial Schroedinger equation on a
 *log grid with the scalar Nuclear relativistic potential.
 *
 *The Nuclear scalar relativistic radial equation with
 *\f$V(r) = -Z/r\f$
 \f[
 -\frac{1}{2M}\frac{d^2u_{nl}}{dr^2}
 +\left[ V+\frac{l(l+1)}{2Mr^2}
 +\frac{\alpha^2}{4M^2r}\frac{dV}{dr}\right] u_{nl}
 -\frac{\alpha^2}{4M^2}\frac{dV}{dr}\frac{du_{nl}}{dr}
 = \varepsilon' u_{nl},
 \f]
 where
 \f[
 \varepsilon'=\varepsilon-\frac{m}{\alpha^2}
 \f]
 and
 \f[
 M(r)=m+\frac{\alpha^2}{2}\left(\varepsilon'-V(r)\right).
 \f]
 *
 To convert this equation to a form that is compatable with the Numerov
 algoritm the following change of variables is necessary:
 \f[ Z_{nl}(r)=u_{nl}(r)/\sqrt{M}, \f] this leads to
 \f[
 \frac{d^2Z_{nl}}{dr^2}=\left[ \frac{l(l+1)}{r^2}
 + 2M(V-\varepsilon')+\frac{\alpha^2}{2Mr}\frac{dV}{dr}
 +\frac{3\alpha^4}{16M^2} \left( \frac{dV}{dr} \right)^2
 + \frac{\alpha^2}{4M}\frac{d^2V}{dr^2}\right]Z_{nl}.
 \f]
 The Numerov algorithm can only be integrated on a uniform grid.  To
 solve the equation on the logarithmic grid it is necessary to make
 another change of variables
 \f[ \Psi_{nl}(r) = \frac{Z_{nl}(r)}{\sqrt{r}} \mbox{ and } y = \log(r) \f]
 with the final result
 \f[
 \frac{d^2\Psi_{nl}}{dy^2} = \left[ \frac{(l+\frac{1}{2})^2}{r^2}
 + 2M(V-\varepsilon')+\frac{\alpha^2}{2Mr}\frac{dV}{dr}
 + \frac{3\alpha^4}{16M^2} \left( \frac{dV}{dr} \right)^2
 + \frac{\alpha^2}{4M}\frac{d^2V}{dr^2}\right]r^2 \Psi_{nl}.
 \f]
 *The Transform function for the Numerov algorithm becomes:
 \f[
 \frac{d^2\Psi_{nl}(y)}{dy^2} + k'^2(r)\Psi_{nl}(y) = 0,
 \f] where
 \f[
 k'^2(r) =  r^2 \left[ 2M(V-\varepsilon')
 + \frac{(l+\frac{1}{2})^2}{r^2} + \frac{\alpha^2}{2Mr}\frac{dV}{dr}
 + \frac{3\alpha^4}{16M^2} \left( \frac{dV}{dr} \right)^2
 + \frac{\alpha^2}{4M}\frac{d^2V}{dr^2} \right]
 \f]
 *
 See Koelling, D.D. and Harmon, B.N. J. Phys. C: Solid State Phys.,
 \textbf{10}, 3107, (1977).
*/
template<class SourceType>
struct NuclearRelLogTransform
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

  /**
   *Fine-Structure constant \f$\alpha = 7.297352568e-3 \approx 1/137\f$
   *For reference see http://physics.nist.gov/cgi-bin/cuu/Value?alph
   */
  value_type alpha;

  ///\f$\alpha^2\f$
  value_type alpha2;

  /**
   * \param v the source function on a grid
   * \param n the principal quantum number
   * \param l the anuglar momentum quantum number
   * \param z the charge
   * \param meff the effective mass
   * \note For the Nuclear Scalar Relativistic potential the
   number of nodes for the radial function with quantum numbers
   \f$n\f$ and \f$l\f$ is \f$n-l-1\f$, the same as in the
   nonrelativistic case.
   * \brief Constructor for \f$\Psi_{nl}\f$
   *
   Calculate the boundary conditions for the One-Dimensional Cubic
   *Spline for the Nuclear Potential \f$V(r) = -Z/r.\f$
   */
  NuclearRelLogTransform(Source_t& v, int n, int l, value_type z,
                         value_type meff=1.0):
    V(v), E(0.0), N(n), CuspParam(z), alpha(7.297352568e-3)
  {
    NumNodes = v.getNumOfNodes();
    L = static_cast<value_type>(l);
    LL = 0.5*(L+0.5)*(L+0.5)/meff;
    alpha2 = alpha*alpha;
    int imin = 0;
    int imax = V.size()-1;
    //calculate first derivative at the boundaries
    value_type yp1 = CuspParam/(V.r(imin)*V.r(imin));
    value_type ypn = CuspParam/(V.r(imax)*V.r(imax));
    V.spline(imin,yp1,imax,ypn);
    //cubic-spline test
#ifdef PRINT_DEBUG
    std::string fname("check_spline");
    fname.append(".dat");
    std::ofstream dfile(fname.c_str());
    for(int ig=imin+1; ig<V.size(); ig++)
    {
      value_type _r = V.r(ig),y,dy,d2y;
      //   V.setgrid(_r);
      y = V.evaluate(_r,1.0/_r,dy,d2y);
      dfile << std::setw(15) <<  setprecision(12) << _r << std::setw(20) << V(ig)
            << std::setw(20) << dV[ig] << std::setw(20) << dy << std::setw(20) << d2V[ig]
            << std::setw(20) << d2y << std::endl;
    }
#endif
  }

  ///returns the starting valid index
  inline int first() const
  {
    return 1;
  }

  /**
   *\param i index of the first value of the log grid
   *\param z0 \f$\Psi_{nl}(r_0)\f$
   *\param z1 \f$\Psi_{nl}(r_0 + dr)\f$
   *\return the derivate of the target function before transformation
   *\brief Assign the first and second values using the cusp condition.
   *
   For \f$l=0\f$ in the limit \f$r\rightarrow 0\f$
   \f[
   \Psi(r) \approx \sqrt{r}(1-Zr),
   \f]
   while for positive \f$l\f$ the condition is
   \f[
   \Psi(r) \approx r^{l+\frac{1}{2}}.
   \f]
  */
  value_type setCusp(int i, value_type& z0, value_type& z1)
  {
    value_type r0 = V.r(i);
    value_type dr = V.dr(i);
    value_type deriv;
    if(L)
    {
      z0 = pow(r0,L+0.5);
      deriv = (L+0.5)*pow(r0,L-0.5);
    }
    else
    {
      value_type r0h = sqrt(r0);
      z0 = r0h*(1-CuspParam*r0);
      deriv = 0.5*(1/r0h - 3.0*CuspParam*r0h);
    }
    z1 = z0 + deriv*dr;
    return deriv;
  }

  ///returns the number of nodes
  inline int nodes() const
  {
    return NumNodes;
  }


  /**
   *@brief return \f$ k'^2(r_i) \f$
   *
   \f[ = r_i^2 \left[ 2M(r_i)(V-\varepsilon')
   + \frac{(l+\frac{1}{2})^2}{r_i^2}
   + \left. \frac{\alpha^2}{2M(r_i)r_i}\frac{dV}{dr} \right|_{r=r_i}
   + \frac{3\alpha^4}{16M(r_i)^2}
   \left. \left( \frac{dV}{dr} \right)^2 \right|_{r=r_i}
   + \frac{\alpha^2}{4M(r_i)}\left. \frac{d^2V}{dr^2} \right|_{r=r_i}
   \right] \f]
  */
  inline value_type k2(int i)
  {
    value_type r0 = V.r(i);
    value_type rinv = 1.0/r0;
    value_type rsq = pow(r0,2.0);
    //    V.setgrid(r0);
    V.evaluate(r0,rinv);
    value_type M = 1.0 + 0.5*alpha2*(E-V.Y);
    value_type Minv = 1.0/M;
    return 2.0*(rsq*(M*(E-V.Y) - 0.5*Minv*alpha2*
                     (0.5*rinv*V.dY + 0.1875*alpha2*Minv*V.dY*V.dY
                      + 0.25*V.d2Y)) - LL);
    //return 2.0*(rsq*(E-V.Y)-LL);
  }


  /*!\fn value_type convert(value_type y, value_type r)
   *\param y the value \f$\Psi_{nl}(r)\f$
   *\param r the grid value
   *\return  \f$u_{nl}(r) = \sqrt{rM}\Psi_{nl}\f$
   */
  inline value_type convert(value_type y, value_type r)
  {
    //    int i = V.setgrid(r);
    value_type rinv = 1.0/r;
    V.evaluate(r,rinv);
    return y*sqrt(r)*sqrt(1.0 + 0.5*alpha2*(E-V.Y));
  }

  ///reset the reference energy and turning point
  inline void reset(value_type e)
  {
    E = e;
  }
};
/**@}*/

#endif

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
    
    



#ifndef QMCPLUSPLUS_SphericalTensor_h_
#define QMCPLUSPLUS_SphericalTensor_h_

//#include <valarray>
//#include "Point.h"

/**Class SphericalTensor
 * @author John Shumway
 * @author Jeongnim Kim
 * \brief evaluates the Real Spherical Harmonics
 \f[
 r^l \Re (Y_l^m(\theta,\phi)).
 \f]
 *
 A list of the real spherical harmonics for the \textit{s},
 \textit{p} and \textit{d} states
 \f[ \Re (Y_0^0) =  \sqrt{\frac{1}{4\pi}} = s \f]
 \f[ r\Re (Y_1^{-1}) = -\sqrt{\frac{3}{8\pi}}y = p_y \f]
 \f[ r\Re (Y_1^0) = \sqrt{\frac{3}{4\pi}}z = p_z \f]
 \f[ r\Re (Y_1^1) = -\sqrt{\frac{3}{8\pi}}x = p_x \f]
 \f[ r^2\Re (Y_2^{-2}) = \sqrt{\frac{15}{8\pi}}xy = d_{xy} \f]
 \f[ r^2\Re (Y_2^{-1}) = -\sqrt{\frac{15}{8\pi}}yz = d_{yz} \f]
 \f[ r^2\Re (Y_2^0) = \sqrt{\frac{5}{16\pi}}(3z^2-r^2) = d_{z^2} \f]
 \f[ r^2\Re (Y_2^1) = -\sqrt{\frac{15}{8\pi}}xz = d_{xz} \f]
 \f[ r^2\Re (Y_2^2) = \sqrt{\frac{15}{32\pi}}(x^2-y^2) = d_{x^2-y^2} \f]
 *
 The template parameter T is the value_type, e.g. double, and the
 template parameter Point_t is a vector type which must have the
 operator[] defined.
 */

template<class T, class Point_t>
class SphericalTensor
{
public :

  typedef T value_type;

  ///constructor
  explicit SphericalTensor(const int lmax);

  ///makes a table of \f$ r^{l} \Re (Y_l^m) \f$ and their gradients up to lmax.
  void evaluate(const Point_t& p);

  ///returns the index \f$ l(l+1)+m \f$
  inline int index(int l, int m) const
  {
    return (l*(l+1))+m;
  }

  ///returns the value of \f$ r^{l} \Re (Y_l^m) \f$ given l,m
  inline value_type getYlm(int l, int m) const
  {
    return Ylm[index(l,m)];
  }

  ///returns the gradient of \f$ r^{l} \Re (Y_l^m) \f$ given l,m
  inline Point_t getGradYlm(int l, int m) const
  {
    return gradYlm[index(l,m)];
  }

  ///returns the value of \f$ r^{l} \Re (Y_l^m) \f$ given index lm
  inline value_type getYlm(int lm) const
  {
    return Ylm[lm];
  }

  ///returns the gradient of \f$ r^{l} \Re (Y_l^m) \f$ given index lm
  inline Point_t getGradYlm(int lm) const
  {
    return gradYlm[lm];
  }

  inline int size() const
  {
    return Ylm.size();
  }

  inline int lmax() const
  {
    return Lmax;
  }

private :

  int Lmax;
  std::vector<value_type> Ylm;
  std::vector<Point_t> gradYlm;
};

template<class T, class Point_t>
SphericalTensor<T, Point_t>::SphericalTensor(const int lmax) : Lmax(lmax)
{
  Ylm.resize((lmax+1)*(lmax+1));
  gradYlm.resize((lmax+1)*(lmax+1));
}


template<class T, class Point_t>
void SphericalTensor<T,Point_t>::evaluate(const Point_t& p)
{
  value_type x=p[0], y=p[1], z=p[2];
  const value_type pi = 4.0*atan(1.0);
  const value_type pi4 = 4.0*pi;
  /*  Calculate r, cos(theta), sin(theta), cos(phi), sin(phi) from input
      coordinates. Check here the coordinate singularity at cos(theta) = +-1.
      This also takes care of r=0 case. */
  value_type cphi,sphi,ctheta;
  value_type r2xy=x*x+y*y;
  value_type r=sqrt(r2xy+z*z);
//@todo use tolerance
  if (r2xy<1e-16)
  {
    cphi = 0.0;
    sphi = 1.0;
    ctheta = (z<0)?-1.0:1.0;
  }
  else
  {
    ctheta = z/r;
    value_type rxyi = 1.0/sqrt(r2xy);
    cphi = x*rxyi;
    sphi = y*rxyi;
  }
  value_type stheta = sqrt(1.0-ctheta*ctheta);
  /* Now to calculate the associated legendre functions P_lm from the
     recursion relation from l=0 to lmax. Conventions of J.D. Jackson,
     Classical Electrodynamics are used. */
  Ylm[0] = 1.0;
// calculate P_ll and P_l,l-1
  value_type fac = 1.0;
  int j = -1;
  for (int l=1; l<=Lmax; l++)
  {
    j += 2;
    fac *= -j*stheta;
    int ll=index(l,l);
    int l1=index(l,l-1);
    int l2=index(l-1,l-1);
    Ylm[ll] = fac;
    Ylm[l1] = j*ctheta*Ylm[l2];
  }
// Use recurence to get other plm's //
  for (int m=0; m<Lmax-1; m++)
  {
    int j = 2*m+1;
    for (int l=m+2; l<=Lmax; l++)
    {
      j += 2;
      int lm=index(l,m);
      int l1=index(l-1,m);
      int l2=index(l-2,m);
      Ylm[lm] = (ctheta*j*Ylm[l1]-(l+m-1)*Ylm[l2])/(l-m);
    }
  }
// Now to calculate r^l Y_lm. //
  value_type sphim,cphim,fac2,temp;
  Ylm[0] = 1.0/sqrt(pi4);
  value_type rpow = 1.0;
  for (int l=1; l<=Lmax; l++)
  {
    rpow *= r;
    fac = rpow*sqrt((2*l+1)/pi4);
    int l0=index(l,0);
    Ylm[l0] *= fac;
    cphim = 1.0;
    sphim = 0.0;
    for (int m=1; m<=l; m++)
    {
      fac2 = (l+m)*(l+1-m);
      fac = fac/sqrt(fac2);
      temp = cphim*cphi-sphim*sphi;
      sphim = sphim*cphi+cphim*sphi;
      cphim = temp;
      int lm = index(l,m);
      temp = fac*Ylm[lm];
      Ylm[lm] = temp*cphim;
      lm = index(l,-m);
      Ylm[lm] = temp*sphim;
    }
  }
// Calculating Gradient now//
  for (int l=0; l<Lmax+1; l++)
  {
    for (int m=-l; m<l+1; m++)
    {
      int lm = index(l-1,0);
      value_type fac = ((value_type) (2*l+1))/(2*l-1);
      value_type gx,gy,gz,dpr,dpi,dmr,dmi;
      int ma = std::abs(m);
      value_type cp = sqrt(fac*(l-ma-1)*(l-ma));
      value_type cm = sqrt(fac*(l+ma-1)*(l+ma));
      value_type c0 = sqrt(fac*(l-ma)*(l+ma));
      gz = (l > ma) ? c0*Ylm[lm+m]:0.0;
      if (l > ma+1)
      {
        dpr = cp*Ylm[lm+ma+1];
        dpi = cp*Ylm[lm-ma-1];
      }
      else
      {
        dpr = 0.0;
        dpi = 0.0;
      }
      if (l > 1)
      {
        switch (ma)
        {
        case 0:
          dmr = -cm*Ylm[lm+1];
          dmi = cm*Ylm[lm-1];
          break;
        case 1:
          dmr = cm*Ylm[lm];
          dmi = 0.0;
          break;
        default:
          dmr = cm*Ylm[lm+ma-1];
          dmi = cm*Ylm[lm-ma+1];
        }
      }
      else
      {
        dmr = (l==1) ? cm*Ylm[lm]:0.0;
        dmi = 0.0;
      }
      if (m < 0)
      {
        gx = 0.5*(dpi-dmi);
        gy = -0.5*(dpr+dmr);
      }
      else
      {
        gx = 0.5*(dpr-dmr);
        gy = 0.5*(dpi+dmi);
      }
      lm = index(l,m);
      gradYlm[lm] = Point_t(gx,gy,gz);
    }
  }
}

#endif

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    

#ifndef QMCPLUSPLUS_SOA_SPHERICAL_CARTESIAN_TENSOR_H
#define QMCPLUSPLUS_SOA_SPHERICAL_CARTESIAN_TENSOR_H

#include "OhmmsSoA/VectorSoaContainer.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{

/** SoaSphericalTensor that evaluates the Real Spherical Harmonics
 *
 * The template parameters
 * - T, the value_type, e.g. double
 * - Point_t, a vector type to provide xyz coordinate.
 * Point_t must have the operator[] defined, e.g., TinyVector\<double,3\>.
 *
 * Real Spherical Harmonics Ylm\f$=r^l S_l^m(x,y,z) \f$ is stored
 * in an array ordered as [0,-1 0 1,-2 -1 0 1 2, -Lmax,-Lmax+1,..., Lmax-1,Lmax]
 * where Lmax is the maximum angular momentum of a center.
 * All the data members, e.g, Ylm and pre-calculated factors,
 * can be accessed by index(l,m) which returns the
 * locator of the combination for l and m.
 */
template<typename T>
struct SoaSphericalTensor
{
  ///maximum angular momentum for the center
  int Lmax;
  /// Normalization factors
  aligned_vector<T> NormFactor;
  ///pre-evaluated factor \f$1/\sqrt{(l+m)\times(l+1-m)}\f$
  aligned_vector<T> FactorLM;
  ///pre-evaluated factor \f$\sqrt{(2l+1)/(4\pi)}\f$
  aligned_vector<T> FactorL;
  ///pre-evaluated factor \f$(2l+1)/(2l-1)\f$
  aligned_vector<T> Factor2L;
  ///composite 
  VectorSoaContainer<T,5> cYlm;

  explicit SoaSphericalTensor(const int l_max, bool addsign=false);

#if (__cplusplus >= 201103L)
  SoaSphericalTensor(const SoaSphericalTensor& rhs)=default;
#endif

  ///compute Ylm
  void evaluate_bare(T x, T y, T z, T* Ylm) const;

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* Ylm) const
  {
    evaluate_bare(x,y,z,Ylm);
    for (int i=0, nl=cYlm.size(); i<nl; i++)
      Ylm[i]*= NormFactor[i];
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z)
  {
    T* restrict Ylm=cYlm.data(0);
    evaluate_bare(x,y,z,Ylm);
    for (int i=0, nl=cYlm.size(); i<nl; i++)
      Ylm[i]*= NormFactor[i];
  }

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGL(T x, T y, T z);

  ///returns the index/locator for (\f$l,m\f$) combo, \f$ l(l+1)+m \f$
  inline int index(int l, int m) const
  {
    return (l*(l+1))+m;
  }

  /** return the starting address of the component
   *
   * component=0(V), 1(dx), 2(dy), 3(dz), 4(Lap)
   */
  inline const T* operator[](size_t component) const
  {
    return cYlm.data(component);
  }

  inline size_t size() const
  {
    return cYlm.size();
  }

  inline int lmax() const
  {
    return Lmax;
  }

};

/** constructor
 * @param l_max maximum angular momentum
 * @param addsign flag to determine what convention to use
 *
 * Evaluate all the constants and prefactors.
 * The spherical harmonics is defined as
 * \f[ Y_l^m (\theta,\phi) = \sqrt{\frac{(2l+1)(l-m)!}{4\pi(l+m)!}} P_l^m(\cos\theta)e^{im\phi}\f]
 * Note that the data member Ylm is a misnomer and should not be confused with "spherical harmonics"
 * \f$Y_l^m\f$.
 - When addsign == true, e.g., Gaussian packages
 \f{eqnarray*}
 S_l^m &=& (-1)^m \sqrt{2}\Re(Y_l^{|m|}), \;\;\;m > 0 \\
 &=& Y_l^0, \;\;\;m = 0 \\
 &=& (-1)^m \sqrt{2}\Im(Y_l^{|m|}),\;\;\;m < 0
 \f}
 - When addsign == false, e.g., SIESTA package,
 \f{eqnarray*}
 S_l^m &=& \sqrt{2}\Re(Y_l^{|m|}), \;\;\;m > 0 \\
 &=& Y_l^0, \;\;\;m = 0 \\
 &=&\sqrt{2}\Im(Y_l^{|m|}),\;\;\;m < 0
 \f}
 */
template<typename T>
inline SoaSphericalTensor<T>::SoaSphericalTensor(const int l_max, bool addsign) : Lmax(l_max)
{
  CONSTEXPR T czero(0);
  CONSTEXPR T cone(1);
  const T pi = 4.0*std::atan(1.0);
  const int ntot = (Lmax+1)*(Lmax+1);
  cYlm.resize(ntot);
  cYlm=czero; 
  NormFactor.resize(ntot,cone);
  const T sqrt2 = std::sqrt(2.0);
  if(addsign)
  {
    for (int l=0; l<=Lmax; l++)
    {
      NormFactor[index(l,0)]=cone;
      for (int m=1; m<=l; m++)
      {
        NormFactor[index(l,m)]=std::pow(-cone,m)*sqrt2;
        NormFactor[index(l,-m)]=std::pow(-cone,-m)*sqrt2;
      }
    }
  }
  else
  {
    for (int l=0; l<=Lmax; l++)
    {
      for (int m=1; m<=l; m++)
      {
        NormFactor[index(l,m)]=sqrt2;
        NormFactor[index(l,-m)]=sqrt2;
      }
    }
  }
  FactorL.resize(Lmax+1);
  const T omega = 1.0/std::sqrt(16.0*std::atan(1.0));
  for(int l=1; l<=Lmax; l++)
    FactorL[l] = std::sqrt(static_cast<T>(2*l+1))*omega;
  Factor2L.resize(Lmax+1);
  for(int l=1; l<=Lmax; l++)
    Factor2L[l] = static_cast<T>(2*l+1)/static_cast<T>(2*l-1);
  FactorLM.resize(ntot);
  for(int l=1; l<=Lmax; l++)
    for(int m=1; m<=l; m++)
    {
      T fac2 = 1.0/std::sqrt(static_cast<T>((l+m)*(l+1-m)));
      FactorLM[index(l,m)]=fac2;
      FactorLM[index(l,-m)]=fac2;
    }
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluate_bare(T x, T y, T z, T* restrict Ylm) const
{
  CONSTEXPR T czero(0);
  CONSTEXPR T cone(1);
  const T pi = 4.0*std::atan(1.0);
  const T omega = 1.0/std::sqrt(4.0*pi);
  const T sqrt2 = std::sqrt(2.0);
  CONSTEXPR T eps2 = std::numeric_limits<T>::epsilon()*std::numeric_limits<T>::epsilon();

  /*  Calculate r, cos(theta), sin(theta), cos(phi), sin(phi) from input
      coordinates. Check here the coordinate singularity at cos(theta) = +-1.
      This also takes care of r=0 case. */
  T cphi,sphi,ctheta;
  T r2xy=x*x+y*y;
  T r=std::sqrt(r2xy+z*z);
  if (r2xy<eps2)
  {
    cphi = czero;
    sphi = cone;
    ctheta = (z<czero)?-cone:cone;
  }
  else
  {
    ctheta = z/r;
    //protect ctheta, when ctheta is slightly >1 or <-1
    if(ctheta>cone) ctheta=cone;
    if(ctheta<-cone) ctheta=-cone;
    T rxyi = cone/std::sqrt(r2xy);
    cphi = x*rxyi;
    sphi = y*rxyi;
  }
  T stheta = std::sqrt(cone-ctheta*ctheta);
  /* Now to calculate the associated legendre functions P_lm from the
     recursion relation from l=0 to Lmax. Conventions of J.D. Jackson,
     Classical Electrodynamics are used. */
  Ylm[0] = cone;
  // calculate P_ll and P_l,l-1
  T fac = cone;
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
  T sphim,cphim,fac2,temp;
  Ylm[0] = omega; //1.0/sqrt(pi4);
  T rpow = 1.0;
  for (int l=1; l<=Lmax; l++)
  {
    rpow *= r;
    //fac = rpow*sqrt(static_cast<T>(2*l+1))*omega;//rpow*sqrt((2*l+1)/pi4);
    //FactorL[l] = sqrt(2*l+1)/sqrt(4*pi)
    fac = rpow*FactorL[l];
    int l0=index(l,0);
    Ylm[l0] *= fac;
    cphim = cone;
    sphim = czero;
    for (int m=1; m<=l; m++)
    {
      temp = cphim*cphi-sphim*sphi;
      sphim = sphim*cphi+cphim*sphi;
      cphim = temp;
      int lm = index(l,m);
      fac *= FactorLM[lm];
      temp = fac*Ylm[lm];
      Ylm[lm] = temp*cphim;
      lm = index(l,-m);
      Ylm[lm] = temp*sphim;
    }
  }
  //for (int i=0; i<Ylm.size(); i++)
  //  Ylm[i]*= NormFactor[i];
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGL(T x, T y, T z)
{
  T* restrict Ylm=cYlm.data(0);
  evaluate_bare(x,y,z,Ylm);

  CONSTEXPR T czero(0);
  CONSTEXPR T ahalf(0.5);
  T* restrict gYlmX=cYlm.data(1);
  T* restrict gYlmY=cYlm.data(2);
  T* restrict gYlmZ=cYlm.data(3);

  // Calculating Gradient now//
  for (int l=1; l<=Lmax; l++)
  {
    //T fac = ((T) (2*l+1))/(2*l-1);
    T fac = Factor2L[l];
    for (int m=-l; m<=l; m++)
    {
      int lm = index(l-1,0);
      T gx,gy,gz,dpr,dpi,dmr,dmi;
      const int ma = std::abs(m);
      const T cp = std::sqrt(fac*(l-ma-1)*(l-ma));
      const T cm = std::sqrt(fac*(l+ma-1)*(l+ma));
      const T c0 = std::sqrt(fac*(l-ma)*(l+ma));
      gz = (l > ma) ? c0*Ylm[lm+m]:czero;
      if (l > ma+1)
      {
        dpr = cp*Ylm[lm+ma+1];
        dpi = cp*Ylm[lm-ma-1];
      }
      else
      {
        dpr = czero;
        dpi = czero;
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
          dmi = czero;
          break;
        default:
          dmr = cm*Ylm[lm+ma-1];
          dmi = cm*Ylm[lm-ma+1];
        }
      }
      else
      {
        dmr = cm*Ylm[lm];
        dmi = czero;
        //dmr = (l==1) ? cm*Ylm[lm]:0.0;
        //dmi = 0.0;
      }
      if (m < 0)
      {
        gx = ahalf*(dpi-dmi);
        gy = -ahalf*(dpr+dmr);
      }
      else
      {
        gx = ahalf*(dpr-dmr);
        gy = ahalf*(dpi+dmi);
      }
      lm = index(l,m);
      if(ma)
      {
        gYlmX[lm]=NormFactor[lm]*gx;
        gYlmY[lm]=NormFactor[lm]*gy;
        gYlmZ[lm]=NormFactor[lm]*gz;
      }
      else
      {
        gYlmX[lm]=gx;
        gYlmY[lm]=gy;
        gYlmZ[lm]=gz;
      }
    }
  }
  for (int i=0; i<cYlm.size(); i++)
    Ylm[i]*= NormFactor[i];
//for (int i=0; i<Ylm.size(); i++) gradYlm[i]*= NormFactor[i];
}

}
#endif

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

#include <stdexcept>
#include <limits>
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
  ///whether to multiply by (-1)^m (see comment below)
  bool Addsign;
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
  VectorSoaContainer<T, 5> cYlm;

  explicit SoaSphericalTensor(const int l_max, bool addsign = false);

  SoaSphericalTensor(const SoaSphericalTensor& rhs) = default;

  ///compute Ylm
  void evaluate_bare(T x, T y, T z, T* Ylm) const;

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* Ylm) const { evaluate_bare(x, y, z, Ylm); }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z)
  {
    T* restrict Ylm = cYlm.data(0);
    evaluate_bare(x, y, z, Ylm);
  }

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGL(T x, T y, T z);

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGH(T x, T y, T z);

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGHGH(T x, T y, T z);

  ///returns the index/locator for (\f$l,m\f$) combo, \f$ l(l+1)+m \f$
  inline int index(int l, int m) const { return (l * (l + 1)) + m; }

  /** return the starting address of the component
   *
   * component=0(V), 1(dx), 2(dy), 3(dz), 4(Lap)
   */
  inline const T* operator[](size_t component) const { return cYlm.data(component); }

  inline size_t size() const { return cYlm.size(); }

  inline int lmax() const { return Lmax; }
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
inline SoaSphericalTensor<T>::SoaSphericalTensor(const int l_max, bool addsign) : Lmax(l_max), Addsign(addsign)
{
  constexpr T czero(0);
  constexpr T cone(1);
  const int ntot = (Lmax + 1) * (Lmax + 1);
  cYlm.resize(ntot);
  cYlm = czero;
  NormFactor.resize(ntot, cone);
  const T sqrt2 = std::sqrt(2.0);
  if (addsign)
  {
    for (int l = 0; l <= Lmax; l++)
    {
      NormFactor[index(l, 0)] = cone;
      for (int m = 1; m <= l; m++)
      {
        NormFactor[index(l, m)]  = std::pow(-cone, m) * sqrt2;
        NormFactor[index(l, -m)] = std::pow(-cone, -m) * sqrt2;
      }
    }
  }
  else
  {
    for (int l = 0; l <= Lmax; l++)
    {
      for (int m = 1; m <= l; m++)
      {
        NormFactor[index(l, m)]  = sqrt2;
        NormFactor[index(l, -m)] = sqrt2;
      }
    }
  }
  FactorL.resize(Lmax + 1);
  const T omega = 1.0 / std::sqrt(16.0 * std::atan(1.0));
  for (int l = 1; l <= Lmax; l++)
    FactorL[l] = std::sqrt(static_cast<T>(2 * l + 1)) * omega;
  Factor2L.resize(Lmax + 1);
  for (int l = 1; l <= Lmax; l++)
    Factor2L[l] = static_cast<T>(2 * l + 1) / static_cast<T>(2 * l - 1);
  FactorLM.resize(ntot);
  for (int l = 1; l <= Lmax; l++)
    for (int m = 1; m <= l; m++)
    {
      T fac2                 = 1.0 / std::sqrt(static_cast<T>((l + m) * (l + 1 - m)));
      FactorLM[index(l, m)]  = fac2;
      FactorLM[index(l, -m)] = fac2;
    }
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluate_bare(T x, T y, T z, T* restrict Ylm) const
{
  SHEval(x, y, z, Ylm, Lmax, Addsign);
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGL(T x, T y, T z)
{
  const T sqrt2inv = std::sqrt(0.5);
  T* restrict Ylm  = cYlm.data(0);
  // remove sign
  if (Addsign)
  {
    evaluate_bare(-x, -y, z, Ylm);
  }
  else
  {
    evaluate_bare(x, y, z, Ylm);
  }
  // remove norm factor for m!=0
  // might be better to just leave it in, add it for m==0, and then remove at the end
  for (int l = 1; l <= Lmax; l++)
  {
    int ll = index(l, 0);
    for (int m = 1; m <= l; m++)
    {
      Ylm[ll + m] *= sqrt2inv;
      Ylm[ll - m] *= sqrt2inv;
    }
  }

  constexpr T czero(0);
  constexpr T ahalf(0.5);
  T* restrict gYlmX = cYlm.data(1);
  T* restrict gYlmY = cYlm.data(2);
  T* restrict gYlmZ = cYlm.data(3);

  // Calculating Gradient now//
  for (int l = 1; l <= Lmax; l++)
  {
    //T fac = ((T) (2*l+1))/(2*l-1);
    T fac = Factor2L[l];
    for (int m = -l; m <= l; m++)
    {
      int lm = index(l - 1, 0);
      T gx, gy, gz, dpr, dpi, dmr, dmi;
      const int ma = std::abs(m);
      const T cp   = std::sqrt(fac * (l - ma - 1) * (l - ma));
      const T cm   = std::sqrt(fac * (l + ma - 1) * (l + ma));
      const T c0   = std::sqrt(fac * (l - ma) * (l + ma));
      gz           = (l > ma) ? c0 * Ylm[lm + m] : czero;
      if (l > ma + 1)
      {
        dpr = cp * Ylm[lm + ma + 1];
        dpi = cp * Ylm[lm - ma - 1];
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
          dmr = -cm * Ylm[lm + 1];
          dmi = cm * Ylm[lm - 1];
          break;
        case 1:
          dmr = cm * Ylm[lm];
          dmi = czero;
          break;
        default:
          dmr = cm * Ylm[lm + ma - 1];
          dmi = cm * Ylm[lm - ma + 1];
        }
      }
      else
      {
        dmr = cm * Ylm[lm];
        dmi = czero;
        //dmr = (l==1) ? cm*Ylm[lm]:0.0;
        //dmi = 0.0;
      }
      if (m < 0)
      {
        gx = ahalf * (dpi - dmi);
        gy = -ahalf * (dpr + dmr);
      }
      else
      {
        gx = ahalf * (dpr - dmr);
        gy = ahalf * (dpi + dmi);
      }
      lm = index(l, m);
      if (ma)
      {
        gYlmX[lm] = NormFactor[lm] * gx;
        gYlmY[lm] = NormFactor[lm] * gy;
        gYlmZ[lm] = NormFactor[lm] * gz;
      }
      else
      {
        gYlmX[lm] = gx;
        gYlmY[lm] = gy;
        gYlmZ[lm] = gz;
      }
    }
  }
  for (int i = 0; i < cYlm.size(); i++)
    Ylm[i] *= NormFactor[i];
  //for (int i=0; i<Ylm.size(); i++) gradYlm[i]*= NormFactor[i];
}
// Lmax == 0
template<typename T>
void SHEval0(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;

  pSH[0] = 0.28209479177387814346198748050032;
  fC0    = fX;
  fS0    = fY;
}

// Lmax == 1
template<typename T>
void SHEval1(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;

  pSH[0] = 0.28209479177387814346198748050032;
  pSH[2] = 0.48860251190291992159616188406979 * fZ;
  fC0    = fX;
  fS0    = fY;

  fTmpB  = -0.48860251190291992156905682975765;
  pSH[3] = fTmpB * fC0;
  pSH[1] = fTmpB * fS0;
}

// Lmax == 2
template<typename T>
void SHEval2(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2  = fZ * fZ;
  T fR2  = fX * fX + fY * fY + fZ * fZ;
  pSH[0] = 0.28209479177387814346198748050032;
  pSH[2] = 0.48860251190291992159616188406979 * fZ;
  pSH[6] = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  fC0    = fX;
  fS0    = fY;

  fTmpA  = -0.48860251190291992156905682975765;
  pSH[3] = fTmpA * fC0;
  pSH[1] = fTmpA * fS0;
  fTmpB  = -1.0925484305920790705553974353492 * fZ;
  pSH[7] = fTmpB * fC0;
  pSH[5] = fTmpB * fS0;
  fC1    = fX * fC0 - fY * fS0;
  fS1    = fX * fS0 + fY * fC0;

  fTmpC  = 0.54627421529603953527769871767461;
  pSH[8] = fTmpC * fC1;
  pSH[4] = fTmpC * fS1;
}

// Lmax == 3
template<typename T>
void SHEval3(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpC   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpC * fC0;
  pSH[9]  = fTmpC * fS0;
}

// Lmax == 4
template<typename T>
void SHEval4(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpC   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpC * fC1;
  pSH[16] = fTmpC * fS1;
}

// Lmax == 5
template<typename T>
void SHEval5(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30] = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fTmpB   = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31] = fTmpB * fC0;
  pSH[29] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fTmpA   = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32] = fTmpA * fC1;
  pSH[28] = fTmpA * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fTmpC   = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33] = fTmpC * fC0;
  pSH[27] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpA * fC1;
  pSH[16] = fTmpA * fS1;
  fTmpB   = 2.075662314881041278823506357476 * fZ;
  pSH[34] = fTmpB * fC1;
  pSH[26] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpC   = -0.65638205684017010278471018769331;
  pSH[35] = fTmpC * fC0;
  pSH[25] = fTmpC * fS0;
}

// Lmax == 6
template<typename T>
void SHEval6(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30] = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42] = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fTmpB   = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31] = fTmpB * fC0;
  pSH[29] = fTmpB * fS0;
  fTmpC   = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43] = fTmpC * fC0;
  pSH[41] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fTmpA   = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32] = fTmpA * fC1;
  pSH[28] = fTmpA * fS1;
  fTmpB   = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44] = fTmpB * fC1;
  pSH[40] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fTmpC   = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33] = fTmpC * fC0;
  pSH[27] = fTmpC * fS0;
  fTmpA   = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45] = fTmpA * fC0;
  pSH[39] = fTmpA * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpA * fC1;
  pSH[16] = fTmpA * fS1;
  fTmpB   = 2.075662314881041278823506357476 * fZ;
  pSH[34] = fTmpB * fC1;
  pSH[26] = fTmpB * fS1;
  fTmpC   = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46] = fTmpC * fC1;
  pSH[38] = fTmpC * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.65638205684017010278471018769331;
  pSH[35] = fTmpA * fC0;
  pSH[25] = fTmpA * fS0;
  fTmpB   = -2.3666191622317520320151890134142 * fZ;
  pSH[47] = fTmpB * fC0;
  pSH[37] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpC   = 0.68318410519191432195978963548555;
  pSH[48] = fTmpC * fC1;
  pSH[36] = fTmpC * fS1;
}

// Lmax == 7
template<typename T>
void SHEval7(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30] = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42] = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56] = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fTmpB   = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31] = fTmpB * fC0;
  pSH[29] = fTmpB * fS0;
  fTmpC   = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43] = fTmpC * fC0;
  pSH[41] = fTmpC * fS0;
  fTmpA   = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57] = fTmpA * fC0;
  pSH[55] = fTmpA * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fTmpA   = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32] = fTmpA * fC1;
  pSH[28] = fTmpA * fS1;
  fTmpB   = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44] = fTmpB * fC1;
  pSH[40] = fTmpB * fS1;
  fTmpC   = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58] = fTmpC * fC1;
  pSH[54] = fTmpC * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fTmpC   = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33] = fTmpC * fC0;
  pSH[27] = fTmpC * fS0;
  fTmpA   = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45] = fTmpA * fC0;
  pSH[39] = fTmpA * fS0;
  fTmpB   = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59] = fTmpB * fC0;
  pSH[53] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpA * fC1;
  pSH[16] = fTmpA * fS1;
  fTmpB   = 2.075662314881041278823506357476 * fZ;
  pSH[34] = fTmpB * fC1;
  pSH[26] = fTmpB * fS1;
  fTmpC   = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46] = fTmpC * fC1;
  pSH[38] = fTmpC * fS1;
  fTmpA   = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60] = fTmpA * fC1;
  pSH[52] = fTmpA * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.65638205684017010278471018769331;
  pSH[35] = fTmpA * fC0;
  pSH[25] = fTmpA * fS0;
  fTmpB   = -2.3666191622317520320151890134142 * fZ;
  pSH[47] = fTmpB * fC0;
  pSH[37] = fTmpB * fS0;
  fTmpC   = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61] = fTmpC * fC0;
  pSH[51] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.68318410519191432195978963548555;
  pSH[48] = fTmpA * fC1;
  pSH[36] = fTmpA * fS1;
  fTmpB   = 2.6459606618019002217956359146456 * fZ;
  pSH[62] = fTmpB * fC1;
  pSH[50] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpC   = -0.70716273252459617822346729654193;
  pSH[63] = fTmpC * fC0;
  pSH[49] = fTmpC * fS0;
}

// Lmax == 8
template<typename T>
void SHEval8(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30] = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42] = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56] = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72] = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fTmpB   = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31] = fTmpB * fC0;
  pSH[29] = fTmpB * fS0;
  fTmpC   = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43] = fTmpC * fC0;
  pSH[41] = fTmpC * fS0;
  fTmpA   = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57] = fTmpA * fC0;
  pSH[55] = fTmpA * fS0;
  fTmpB   = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73] = fTmpB * fC0;
  pSH[71] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fTmpA   = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32] = fTmpA * fC1;
  pSH[28] = fTmpA * fS1;
  fTmpB   = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44] = fTmpB * fC1;
  pSH[40] = fTmpB * fS1;
  fTmpC   = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58] = fTmpC * fC1;
  pSH[54] = fTmpC * fS1;
  fTmpA   = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74] = fTmpA * fC1;
  pSH[70] = fTmpA * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fTmpC   = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33] = fTmpC * fC0;
  pSH[27] = fTmpC * fS0;
  fTmpA   = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45] = fTmpA * fC0;
  pSH[39] = fTmpA * fS0;
  fTmpB   = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59] = fTmpB * fC0;
  pSH[53] = fTmpB * fS0;
  fTmpC   = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75] = fTmpC * fC0;
  pSH[69] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpA * fC1;
  pSH[16] = fTmpA * fS1;
  fTmpB   = 2.075662314881041278823506357476 * fZ;
  pSH[34] = fTmpB * fC1;
  pSH[26] = fTmpB * fS1;
  fTmpC   = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46] = fTmpC * fC1;
  pSH[38] = fTmpC * fS1;
  fTmpA   = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60] = fTmpA * fC1;
  pSH[52] = fTmpA * fS1;
  fTmpB   = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76] = fTmpB * fC1;
  pSH[68] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.65638205684017010278471018769331;
  pSH[35] = fTmpA * fC0;
  pSH[25] = fTmpA * fS0;
  fTmpB   = -2.3666191622317520320151890134142 * fZ;
  pSH[47] = fTmpB * fC0;
  pSH[37] = fTmpB * fS0;
  fTmpC   = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61] = fTmpC * fC0;
  pSH[51] = fTmpC * fS0;
  fTmpA   = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77] = fTmpA * fC0;
  pSH[67] = fTmpA * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.68318410519191432195978963548555;
  pSH[48] = fTmpA * fC1;
  pSH[36] = fTmpA * fS1;
  fTmpB   = 2.6459606618019002217956359146456 * fZ;
  pSH[62] = fTmpB * fC1;
  pSH[50] = fTmpB * fS1;
  fTmpC   = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78] = fTmpC * fC1;
  pSH[66] = fTmpC * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.70716273252459617822346729654193;
  pSH[63] = fTmpA * fC0;
  pSH[49] = fTmpA * fS0;
  fTmpB   = -2.9157066406993194756895604324853 * fZ;
  pSH[79] = fTmpB * fC0;
  pSH[65] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpC   = 0.72892666017482986886817999949706;
  pSH[80] = fTmpC * fC1;
  pSH[64] = fTmpC * fS1;
}

// Lmax == 9
template<typename T>
void SHEval9(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2   = fZ * fZ;
  T fR2   = fX * fX + fY * fY + fZ * fZ;
  pSH[0]  = 0.28209479177387814346198748050032;
  pSH[2]  = 0.48860251190291992159616188406979 * fZ;
  pSH[6]  = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12] = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20] = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30] = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42] = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56] = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72] = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90] = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  fC0     = fX;
  fS0     = fY;

  fTmpA   = -0.48860251190291992156905682975765;
  pSH[3]  = fTmpA * fC0;
  pSH[1]  = fTmpA * fS0;
  fTmpB   = -1.0925484305920790705553974353492 * fZ;
  pSH[7]  = fTmpB * fC0;
  pSH[5]  = fTmpB * fS0;
  fTmpC   = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13] = fTmpC * fC0;
  pSH[11] = fTmpC * fS0;
  fTmpA   = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21] = fTmpA * fC0;
  pSH[19] = fTmpA * fS0;
  fTmpB   = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31] = fTmpB * fC0;
  pSH[29] = fTmpB * fS0;
  fTmpC   = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43] = fTmpC * fC0;
  pSH[41] = fTmpC * fS0;
  fTmpA   = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57] = fTmpA * fC0;
  pSH[55] = fTmpA * fS0;
  fTmpB   = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73] = fTmpB * fC0;
  pSH[71] = fTmpB * fS0;
  fTmpC   = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91] = fTmpC * fC0;
  pSH[89] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.54627421529603953527769871767461;
  pSH[8]  = fTmpA * fC1;
  pSH[4]  = fTmpA * fS1;
  fTmpB   = 1.4453057213202770277206063442854 * fZ;
  pSH[14] = fTmpB * fC1;
  pSH[10] = fTmpB * fS1;
  fTmpC   = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22] = fTmpC * fC1;
  pSH[18] = fTmpC * fS1;
  fTmpA   = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32] = fTmpA * fC1;
  pSH[28] = fTmpA * fS1;
  fTmpB   = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44] = fTmpB * fC1;
  pSH[40] = fTmpB * fS1;
  fTmpC   = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58] = fTmpC * fC1;
  pSH[54] = fTmpC * fS1;
  fTmpA   = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74] = fTmpA * fC1;
  pSH[70] = fTmpA * fS1;
  fTmpB   = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92] = fTmpB * fC1;
  pSH[88] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.59004358992664351034589803601804;
  pSH[15] = fTmpA * fC0;
  pSH[9]  = fTmpA * fS0;
  fTmpB   = -1.7701307697799305311461143253027 * fZ;
  pSH[23] = fTmpB * fC0;
  pSH[17] = fTmpB * fS0;
  fTmpC   = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33] = fTmpC * fC0;
  pSH[27] = fTmpC * fS0;
  fTmpA   = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45] = fTmpA * fC0;
  pSH[39] = fTmpA * fS0;
  fTmpB   = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59] = fTmpB * fC0;
  pSH[53] = fTmpB * fS0;
  fTmpC   = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75] = fTmpC * fC0;
  pSH[69] = fTmpC * fS0;
  fTmpA   = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93] = fTmpA * fC0;
  pSH[87] = fTmpA * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.62583573544917613459435609679637;
  pSH[24] = fTmpA * fC1;
  pSH[16] = fTmpA * fS1;
  fTmpB   = 2.075662314881041278823506357476 * fZ;
  pSH[34] = fTmpB * fC1;
  pSH[26] = fTmpB * fS1;
  fTmpC   = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46] = fTmpC * fC1;
  pSH[38] = fTmpC * fS1;
  fTmpA   = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60] = fTmpA * fC1;
  pSH[52] = fTmpA * fS1;
  fTmpB   = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76] = fTmpB * fC1;
  pSH[68] = fTmpB * fS1;
  fTmpC   = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94] = fTmpC * fC1;
  pSH[86] = fTmpC * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.65638205684017010278471018769331;
  pSH[35] = fTmpA * fC0;
  pSH[25] = fTmpA * fS0;
  fTmpB   = -2.3666191622317520320151890134142 * fZ;
  pSH[47] = fTmpB * fC0;
  pSH[37] = fTmpB * fS0;
  fTmpC   = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61] = fTmpC * fC0;
  pSH[51] = fTmpC * fS0;
  fTmpA   = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77] = fTmpA * fC0;
  pSH[67] = fTmpA * fS0;
  fTmpB   = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95] = fTmpB * fC0;
  pSH[85] = fTmpB * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.68318410519191432195978963548555;
  pSH[48] = fTmpA * fC1;
  pSH[36] = fTmpA * fS1;
  fTmpB   = 2.6459606618019002217956359146456 * fZ;
  pSH[62] = fTmpB * fC1;
  pSH[50] = fTmpB * fS1;
  fTmpC   = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78] = fTmpC * fC1;
  pSH[66] = fTmpC * fS1;
  fTmpA   = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96] = fTmpA * fC1;
  pSH[84] = fTmpA * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpA   = -0.70716273252459617822346729654193;
  pSH[63] = fTmpA * fC0;
  pSH[49] = fTmpA * fS0;
  fTmpB   = -2.9157066406993194756895604324853 * fZ;
  pSH[79] = fTmpB * fC0;
  pSH[65] = fTmpB * fS0;
  fTmpC   = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97] = fTmpC * fC0;
  pSH[83] = fTmpC * fS0;
  fC1     = fX * fC0 - fY * fS0;
  fS1     = fX * fS0 + fY * fC0;

  fTmpA   = 0.72892666017482986886817999949706;
  pSH[80] = fTmpA * fC1;
  pSH[64] = fTmpA * fS1;
  fTmpB   = 3.1773176489546974771236570456168 * fZ;
  pSH[98] = fTmpB * fC1;
  pSH[82] = fTmpB * fS1;
  fC0     = fX * fC1 - fY * fS1;
  fS0     = fX * fS1 + fY * fC1;

  fTmpC   = -0.74890095185318829654197436695995;
  pSH[99] = fTmpC * fC0;
  pSH[81] = fTmpC * fS0;
}

// Lmax == 10
template<typename T>
void SHEval10(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpC    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpC * fC1;
  pSH[100] = fTmpC * fS1;
}

// Lmax == 11
template<typename T>
void SHEval11(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpC    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpC * fC0;
  pSH[121] = fTmpC * fS0;
}

// Lmax == 12
template<typename T>
void SHEval12(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpC    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpC * fC1;
  pSH[144] = fTmpC * fS1;
}

// Lmax == 13
template<typename T>
void SHEval13(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpC    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpC * fC0;
  pSH[169] = fTmpC * fS0;
}

// Lmax == 14
template<typename T>
void SHEval14(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpC    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpC * fC1;
  pSH[196] = fTmpC * fS1;
}

// Lmax == 15
template<typename T>
void SHEval15(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  pSH[240] = 1.9988885800753266486651932298813 * fZ * pSH[210] + -1.0000823011400101476917751108786 * pSH[182] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fTmpC    = 2.003345416333103836325699176335 * fZ * fTmpB + -0.99975195336341716702754575663015 * fTmpA * fR2;
  pSH[241] = fTmpC * fC0;
  pSH[239] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fTmpB    = 2.0168969490698876096713976213692 * fZ * fTmpA + -0.99874229606879180774856724633892 * fTmpC * fR2;
  pSH[242] = fTmpB * fC1;
  pSH[238] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fTmpA    = 2.0401071141087266541148947940343 * fZ * fTmpC + -0.99699479851094886096287209231726 * fTmpB * fR2;
  pSH[243] = fTmpA * fC0;
  pSH[237] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fTmpC    = 2.0739902137422357166710723541669 * fZ * fTmpB + -0.99440219512922986143561854266437 * fTmpA * fR2;
  pSH[244] = fTmpC * fC1;
  pSH[236] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fTmpB    = 2.120141504711419020139800961644 * fZ * fTmpA + -0.99079092984678996129579292562184 * fTmpC * fR2;
  pSH[245] = fTmpB * fC0;
  pSH[235] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fTmpA    = 2.180966243804281491170185547368 * fZ * fTmpC + -0.98588907503509976264916003785288 * fTmpB * fR2;
  pSH[246] = fTmpA * fC1;
  pSH[234] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fTmpC    = 2.2600784378986817490121002949266 * fZ * fTmpB + -0.97926740294193723591464201261303 * fTmpA * fR2;
  pSH[247] = fTmpC * fC0;
  pSH[233] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fTmpB    = 2.3630173363047971158562576832196 * fZ * fTmpA + -0.9702261872276653011690043804105 * fTmpC * fR2;
  pSH[248] = fTmpB * fC1;
  pSH[232] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fTmpA    = 2.4986107250941583107772814287273 * fZ * fTmpC + -0.95756141751469769717681687626332 * fTmpB * fR2;
  pSH[249] = fTmpA * fC0;
  pSH[231] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fTmpC    = 2.6817904466978772499811956020466 * fZ * fTmpB + -0.93903023262181382105036678287213 * fTmpA * fR2;
  pSH[250] = fTmpC * fC1;
  pSH[230] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fTmpB    = 2.9401072717216916591054243212966 * fZ * fTmpA + -0.90994035683194807477211507595882 * fTmpC * fR2;
  pSH[251] = fTmpB * fC0;
  pSH[229] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fTmpA    = fZ * (53.794072123058085929669935865149 * fZ2 + -5.564904012730146820250864969637 * fR2);
  pSH[252] = fTmpA * fC1;
  pSH[228] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fTmpC    = -17.608243388844820554936521084244 * fZ2 + 0.60718080651189036396706694143077 * fR2;
  pSH[253] = fTmpC * fC0;
  pSH[227] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpA * fC1;
  pSH[196] = fTmpA * fS1;
  fTmpB    = 4.6241512566300120266882256458985 * fZ;
  pSH[254] = fTmpB * fC1;
  pSH[226] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpC    = -0.84425065085737263610244154876661;
  pSH[255] = fTmpC * fC0;
  pSH[225] = fTmpC * fS0;
}

// Lmax == 16
template<typename T>
void SHEval16(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  pSH[240] = 1.9988885800753266486651932298813 * fZ * pSH[210] + -1.0000823011400101476917751108786 * pSH[182] * fR2;
  pSH[272] = 1.9990231989649344736563116309291 * fZ * pSH[240] + -1.0000673468701305762361061790777 * pSH[210] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fTmpC    = 2.003345416333103836325699176335 * fZ * fTmpB + -0.99975195336341716702754575663015 * fTmpA * fR2;
  pSH[241] = fTmpC * fC0;
  pSH[239] = fTmpC * fS0;
  fTmpA    = 2.0029390170153341292971771459008 * fZ * fTmpC + -0.99979713966725040625964371354684 * fTmpB * fR2;
  pSH[273] = fTmpA * fC0;
  pSH[271] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fTmpB    = 2.0168969490698876096713976213692 * fZ * fTmpA + -0.99874229606879180774856724633892 * fTmpC * fR2;
  pSH[242] = fTmpB * fC1;
  pSH[238] = fTmpB * fS1;
  fTmpC    = 2.014825999813336120320209077228 * fZ * fTmpB + -0.99897320026315349013367253516726 * fTmpA * fR2;
  pSH[274] = fTmpC * fC1;
  pSH[270] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fTmpA    = 2.0401071141087266541148947940343 * fZ * fTmpC + -0.99699479851094886096287209231726 * fTmpB * fR2;
  pSH[243] = fTmpA * fC0;
  pSH[237] = fTmpA * fS0;
  fTmpB    = 2.0351168037383750113438612983074 * fZ * fTmpA + -0.99755389786357772280477387849551 * fTmpC * fR2;
  pSH[275] = fTmpB * fC0;
  pSH[269] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fTmpC    = 2.0739902137422357166710723541669 * fZ * fTmpB + -0.99440219512922986143561854266437 * fTmpA * fR2;
  pSH[244] = fTmpC * fC1;
  pSH[236] = fTmpC * fS1;
  fTmpA    = 2.064582282206257818736941378468 * fZ * fTmpC + -0.99546384960081245812197822675493 * fTmpB * fR2;
  pSH[276] = fTmpA * fC1;
  pSH[268] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fTmpB    = 2.120141504711419020139800961644 * fZ * fTmpA + -0.99079092984678996129579292562184 * fTmpC * fR2;
  pSH[245] = fTmpB * fC0;
  pSH[235] = fTmpB * fS0;
  fTmpC    = 2.1044171232366050514495797729708 * fZ * fTmpB + -0.99258333397093026685296598965458 * fTmpA * fR2;
  pSH[277] = fTmpC * fC0;
  pSH[267] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fTmpA    = 2.180966243804281491170185547368 * fZ * fTmpC + -0.98588907503509976264916003785288 * fTmpB * fR2;
  pSH[246] = fTmpA * fC1;
  pSH[234] = fTmpA * fS1;
  fTmpB    = 2.1563858652847824675605897803976 * fZ * fTmpA + -0.98872959240459255315664616192706 * fTmpC * fR2;
  pSH[278] = fTmpB * fC1;
  pSH[266] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fTmpC    = 2.2600784378986817490121002949266 * fZ * fTmpB + -0.97926740294193723591464201261303 * fTmpA * fR2;
  pSH[247] = fTmpC * fC0;
  pSH[233] = fTmpC * fS0;
  fTmpA    = 2.2230674720995866292406334396858 * fZ * fTmpC + -0.98362403482177094963569141672366 * fTmpB * fR2;
  pSH[279] = fTmpA * fC0;
  pSH[265] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fTmpB    = 2.3630173363047971158562576832196 * fZ * fTmpA + -0.9702261872276653011690043804105 * fTmpC * fR2;
  pSH[248] = fTmpB * fC1;
  pSH[232] = fTmpB * fS1;
  fTmpC    = 2.3082731640774234846327089831775 * fZ * fTmpB + -0.9768329366922967098513588823927 * fTmpA * fR2;
  pSH[280] = fTmpC * fC1;
  pSH[264] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fTmpA    = 2.4986107250941583107772814287273 * fZ * fTmpC + -0.95756141751469769717681687626332 * fTmpB * fR2;
  pSH[249] = fTmpA * fC0;
  pSH[231] = fTmpA * fS0;
  fTmpB    = 2.4177911997760033447311955878689 * fZ * fTmpA + -0.96765421499777267555627777162464 * fTmpC * fR2;
  pSH[281] = fTmpB * fC0;
  pSH[263] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fTmpC    = 2.6817904466978772499811956020466 * fZ * fTmpB + -0.93903023262181382105036678287213 * fTmpA * fR2;
  pSH[250] = fTmpC * fC1;
  pSH[230] = fTmpC * fS1;
  fTmpA    = 2.5607991541103546086991654684439 * fZ * fTmpC + -0.95488413617980451758614560131555 * fTmpB * fR2;
  pSH[282] = fTmpA * fC1;
  pSH[262] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fTmpB    = 2.9401072717216916591054243212966 * fZ * fTmpA + -0.90994035683194807477211507595882 * fTmpC * fR2;
  pSH[251] = fTmpB * fC0;
  pSH[229] = fTmpB * fS0;
  fTmpC    = 2.7527763762750104249103083597916 * fZ * fTmpB + -0.93628433314374189411904980673285 * fTmpA * fR2;
  pSH[283] = fTmpC * fC0;
  pSH[261] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fTmpA    = fZ * (53.794072123058085929669935865149 * fZ2 + -5.564904012730146820250864969637 * fR2);
  pSH[252] = fTmpA * fC1;
  pSH[228] = fTmpA * fS1;
  fTmpB    = 3.0222389997200041803718933985934 * fZ * fTmpA + -0.90717582656041188333062227910908 * fTmpC * fR2;
  pSH[284] = fTmpB * fC1;
  pSH[260] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fTmpC    = -17.608243388844820554936521084244 * fZ2 + 0.60718080651189036396706694143077 * fR2;
  pSH[253] = fTmpC * fC0;
  pSH[227] = fTmpC * fS0;
  fTmpA    = fZ * (-60.380154942952015812568378194669 * fZ2 + 5.8432408009308402399399617888065 * fR2);
  pSH[285] = fTmpA * fC0;
  pSH[259] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpA * fC1;
  pSH[196] = fTmpA * fS1;
  fTmpB    = 4.6241512566300120266882256458985 * fZ;
  pSH[254] = fTmpB * fC1;
  pSH[226] = fTmpB * fS1;
  fTmpC    = 19.093881509360250419912730102112 * fZ2 + -0.61593166159226614257433257693108 * fR2;
  pSH[286] = fTmpC * fC1;
  pSH[258] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.84425065085737263610244154876661;
  pSH[255] = fTmpA * fC0;
  pSH[225] = fTmpA * fS0;
  fTmpB    = -4.849850753230681765781895364853 * fZ;
  pSH[287] = fTmpB * fC0;
  pSH[257] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpC    = 0.85734058883802509643716829867977;
  pSH[288] = fTmpC * fC1;
  pSH[256] = fTmpC * fS1;
}

// Lmax == 17
template<typename T>
void SHEval17(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  pSH[240] = 1.9988885800753266486651932298813 * fZ * pSH[210] + -1.0000823011400101476917751108786 * pSH[182] * fR2;
  pSH[272] = 1.9990231989649344736563116309291 * fZ * pSH[240] + -1.0000673468701305762361061790777 * pSH[210] * fR2;
  pSH[306] = 1.9991347609372268759927310233238 * fZ * pSH[272] + -1.0000558082429209263942634922095 * pSH[240] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fTmpC    = 2.003345416333103836325699176335 * fZ * fTmpB + -0.99975195336341716702754575663015 * fTmpA * fR2;
  pSH[241] = fTmpC * fC0;
  pSH[239] = fTmpC * fS0;
  fTmpA    = 2.0029390170153341292971771459008 * fZ * fTmpC + -0.99979713966725040625964371354684 * fTmpB * fR2;
  pSH[273] = fTmpA * fC0;
  pSH[271] = fTmpA * fS0;
  fTmpB    = 2.0026024734496526301299329508865 * fZ * fTmpA + -0.99983197513113354919136316345529 * fTmpC * fR2;
  pSH[307] = fTmpB * fC0;
  pSH[305] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fTmpB    = 2.0168969490698876096713976213692 * fZ * fTmpA + -0.99874229606879180774856724633892 * fTmpC * fR2;
  pSH[242] = fTmpB * fC1;
  pSH[238] = fTmpB * fS1;
  fTmpC    = 2.014825999813336120320209077228 * fZ * fTmpB + -0.99897320026315349013367253516726 * fTmpA * fR2;
  pSH[274] = fTmpC * fC1;
  pSH[270] = fTmpC * fS1;
  fTmpA    = 2.0131148946216081338771858311176 * fZ * fTmpC + -0.9991507429465936452315545646119 * fTmpB * fR2;
  pSH[308] = fTmpA * fC1;
  pSH[304] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fTmpA    = 2.0401071141087266541148947940343 * fZ * fTmpC + -0.99699479851094886096287209231726 * fTmpB * fR2;
  pSH[243] = fTmpA * fC0;
  pSH[237] = fTmpA * fS0;
  fTmpB    = 2.0351168037383750113438612983074 * fZ * fTmpA + -0.99755389786357772280477387849551 * fTmpC * fR2;
  pSH[275] = fTmpB * fC0;
  pSH[269] = fTmpB * fS0;
  fTmpC    = 2.0310096011589900901091187979119 * fZ * fTmpB + -0.99798183447169211033548480438427 * fTmpA * fR2;
  pSH[309] = fTmpC * fC0;
  pSH[303] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fTmpC    = 2.0739902137422357166710723541669 * fZ * fTmpB + -0.99440219512922986143561854266437 * fTmpA * fR2;
  pSH[244] = fTmpC * fC1;
  pSH[236] = fTmpC * fS1;
  fTmpA    = 2.064582282206257818736941378468 * fZ * fTmpC + -0.99546384960081245812197822675493 * fTmpB * fR2;
  pSH[276] = fTmpA * fC1;
  pSH[268] = fTmpA * fS1;
  fTmpB    = 2.0568833780186057910190078334978 * fZ * fTmpA + -0.99627096277343579159169878467495 * fTmpC * fR2;
  pSH[310] = fTmpB * fC1;
  pSH[302] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fTmpB    = 2.120141504711419020139800961644 * fZ * fTmpA + -0.99079092984678996129579292562184 * fTmpC * fR2;
  pSH[245] = fTmpB * fC0;
  pSH[235] = fTmpB * fS0;
  fTmpC    = 2.1044171232366050514495797729708 * fZ * fTmpB + -0.99258333397093026685296598965458 * fTmpA * fR2;
  pSH[277] = fTmpC * fC0;
  pSH[267] = fTmpC * fS0;
  fTmpA    = 2.0916500663351888698003600008946 * fZ * fTmpC + -0.99393320993236341956249269014023 * fTmpB * fR2;
  pSH[311] = fTmpA * fC0;
  pSH[301] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fTmpA    = 2.180966243804281491170185547368 * fZ * fTmpC + -0.98588907503509976264916003785288 * fTmpB * fR2;
  pSH[246] = fTmpA * fC1;
  pSH[234] = fTmpA * fS1;
  fTmpB    = 2.1563858652847824675605897803976 * fZ * fTmpA + -0.98872959240459255315664616192706 * fTmpC * fR2;
  pSH[278] = fTmpB * fC1;
  pSH[266] = fTmpB * fS1;
  fTmpC    = 2.1366369348357590876238965016398 * fZ * fTmpB + -0.99084165280112553362114671817729 * fTmpA * fR2;
  pSH[312] = fTmpC * fC1;
  pSH[300] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fTmpC    = 2.2600784378986817490121002949266 * fZ * fTmpB + -0.97926740294193723591464201261303 * fTmpA * fR2;
  pSH[247] = fTmpC * fC0;
  pSH[233] = fTmpC * fS0;
  fTmpA    = 2.2230674720995866292406334396858 * fZ * fTmpC + -0.98362403482177094963569141672366 * fTmpB * fR2;
  pSH[279] = fTmpA * fC0;
  pSH[265] = fTmpA * fS0;
  fTmpB    = 2.1937410968480305151814130359966 * fZ * fTmpA + -0.98680814882156560015943197461397 * fTmpC * fR2;
  pSH[313] = fTmpB * fC0;
  pSH[299] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fTmpB    = 2.3630173363047971158562576832196 * fZ * fTmpA + -0.9702261872276653011690043804105 * fTmpC * fR2;
  pSH[248] = fTmpB * fC1;
  pSH[232] = fTmpB * fS1;
  fTmpC    = 2.3082731640774234846327089831775 * fZ * fTmpB + -0.9768329366922967098513588823927 * fTmpA * fR2;
  pSH[280] = fTmpC * fC1;
  pSH[264] = fTmpC * fS1;
  fTmpA    = 2.2656860623955237929363221160983 * fZ * fTmpC + -0.98155023315928857430062368094603 * fTmpB * fR2;
  pSH[314] = fTmpA * fC1;
  pSH[298] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fTmpA    = 2.4986107250941583107772814287273 * fZ * fTmpC + -0.95756141751469769717681687626332 * fTmpB * fR2;
  pSH[249] = fTmpA * fC0;
  pSH[231] = fTmpA * fS0;
  fTmpB    = 2.4177911997760033447311955878689 * fZ * fTmpA + -0.96765421499777267555627777162464 * fTmpC * fR2;
  pSH[281] = fTmpB * fC0;
  pSH[263] = fTmpB * fS0;
  fTmpC    = 2.3564559438666820498364806724112 * fZ * fTmpB + -0.97463169858712211832086833029898 * fTmpA * fR2;
  pSH[315] = fTmpC * fC0;
  pSH[297] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fTmpC    = 2.6817904466978772499811956020466 * fZ * fTmpB + -0.93903023262181382105036678287213 * fTmpA * fR2;
  pSH[250] = fTmpC * fC1;
  pSH[230] = fTmpC * fS1;
  fTmpA    = 2.5607991541103546086991654684439 * fZ * fTmpC + -0.95488413617980451758614560131555 * fTmpB * fR2;
  pSH[282] = fTmpA * fC1;
  pSH[262] = fTmpA * fS1;
  fTmpB    = 2.4720661623652209828994052998041 * fZ * fTmpA + -0.96534949193391143146009830688925 * fTmpC * fR2;
  pSH[316] = fTmpB * fC1;
  pSH[296] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fTmpB    = 2.9401072717216916591054243212966 * fZ * fTmpA + -0.90994035683194807477211507595882 * fTmpC * fR2;
  pSH[251] = fTmpB * fC0;
  pSH[229] = fTmpB * fS0;
  fTmpC    = 2.7527763762750104249103083597916 * fZ * fTmpB + -0.93628433314374189411904980673285 * fTmpA * fR2;
  pSH[283] = fTmpC * fC0;
  pSH[261] = fTmpC * fS0;
  fTmpA    = 2.6220221204253788677401154627589 * fZ * fTmpC + -0.95250095250142875237531203680419 * fTmpB * fR2;
  pSH[317] = fTmpA * fC0;
  pSH[295] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fTmpA    = fZ * (53.794072123058085929669935865149 * fZ2 + -5.564904012730146820250864969637 * fR2);
  pSH[252] = fTmpA * fC1;
  pSH[228] = fTmpA * fS1;
  fTmpB    = 3.0222389997200041803718933985934 * fZ * fTmpA + -0.90717582656041188333062227910908 * fTmpC * fR2;
  pSH[284] = fTmpB * fC1;
  pSH[260] = fTmpB * fS1;
  fTmpC    = 2.8223247937435036367635754483985 * fZ * fTmpB + -0.93385228435109810061348634135925 * fTmpA * fR2;
  pSH[318] = fTmpC * fC1;
  pSH[294] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fTmpC    = -17.608243388844820554936521084244 * fZ2 + 0.60718080651189036396706694143077 * fR2;
  pSH[253] = fTmpC * fC0;
  pSH[227] = fTmpC * fS0;
  fTmpA    = fZ * (-60.380154942952015812568378194669 * fZ2 + 5.8432408009308402399399617888065 * fR2);
  pSH[285] = fTmpA * fC0;
  pSH[259] = fTmpA * fS0;
  fTmpB    = 3.1024184114977141491897166813985 * fZ * fTmpA + -0.90473663963430495837036646178397 * fTmpC * fR2;
  pSH[319] = fTmpB * fC0;
  pSH[293] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpA * fC1;
  pSH[196] = fTmpA * fS1;
  fTmpB    = 4.6241512566300120266882256458985 * fZ;
  pSH[254] = fTmpB * fC1;
  pSH[226] = fTmpB * fS1;
  fTmpC    = 19.093881509360250419912730102112 * fZ2 + -0.61593166159226614257433257693108 * fR2;
  pSH[286] = fTmpC * fC1;
  pSH[258] = fTmpC * fS1;
  fTmpA    = fZ * (67.288948373844056319303952307109 * fZ2 + -6.1171771248949142110035159802806 * fR2);
  pSH[320] = fTmpA * fC1;
  pSH[292] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.84425065085737263610244154876661;
  pSH[255] = fTmpA * fC0;
  pSH[225] = fTmpA * fS0;
  fTmpB    = -4.849850753230681765781895364853 * fZ;
  pSH[287] = fTmpB * fC0;
  pSH[257] = fTmpB * fS0;
  fTmpC    = -20.602948605549728459257474710853 * fZ2 + 0.62433177592574934725022650638948 * fR2;
  pSH[321] = fTmpC * fC0;
  pSH[291] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.85734058883802509643716829867977;
  pSH[288] = fTmpA * fC1;
  pSH[256] = fTmpA * fS1;
  fTmpB    = 5.0720953248553606105934743464303 * fZ;
  pSH[322] = fTmpB * fC1;
  pSH[290] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpC    = -0.86985717192062808221161840371849;
  pSH[323] = fTmpC * fC0;
  pSH[289] = fTmpC * fS0;
}

// Lmax == 18
template<typename T>
void SHEval18(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  pSH[240] = 1.9988885800753266486651932298813 * fZ * pSH[210] + -1.0000823011400101476917751108786 * pSH[182] * fR2;
  pSH[272] = 1.9990231989649344736563116309291 * fZ * pSH[240] + -1.0000673468701305762361061790777 * pSH[210] * fR2;
  pSH[306] = 1.9991347609372268759927310233238 * fZ * pSH[272] + -1.0000558082429209263942634922095 * pSH[240] * fR2;
  pSH[342] = 1.9992282461607312885921994283223 * fZ * pSH[306] + -1.0000467628422711158057978320102 * pSH[272] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fTmpC    = 2.003345416333103836325699176335 * fZ * fTmpB + -0.99975195336341716702754575663015 * fTmpA * fR2;
  pSH[241] = fTmpC * fC0;
  pSH[239] = fTmpC * fS0;
  fTmpA    = 2.0029390170153341292971771459008 * fZ * fTmpC + -0.99979713966725040625964371354684 * fTmpB * fR2;
  pSH[273] = fTmpA * fC0;
  pSH[271] = fTmpA * fS0;
  fTmpB    = 2.0026024734496526301299329508865 * fZ * fTmpA + -0.99983197513113354919136316345529 * fTmpC * fR2;
  pSH[307] = fTmpB * fC0;
  pSH[305] = fTmpB * fS0;
  fTmpC    = 2.0023206350873464509643184783272 * fZ * fTmpB + -0.99985926394976398475303303037265 * fTmpA * fR2;
  pSH[343] = fTmpC * fC0;
  pSH[341] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fTmpB    = 2.0168969490698876096713976213692 * fZ * fTmpA + -0.99874229606879180774856724633892 * fTmpC * fR2;
  pSH[242] = fTmpB * fC1;
  pSH[238] = fTmpB * fS1;
  fTmpC    = 2.014825999813336120320209077228 * fZ * fTmpB + -0.99897320026315349013367253516726 * fTmpA * fR2;
  pSH[274] = fTmpC * fC1;
  pSH[270] = fTmpC * fS1;
  fTmpA    = 2.0131148946216081338771858311176 * fZ * fTmpC + -0.9991507429465936452315545646119 * fTmpB * fR2;
  pSH[308] = fTmpA * fC1;
  pSH[304] = fTmpA * fS1;
  fTmpB    = 2.0116846174288851485396217855239 * fZ * fTmpA + -0.99928952033659667242206786630376 * fTmpC * fR2;
  pSH[344] = fTmpB * fC1;
  pSH[340] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fTmpA    = 2.0401071141087266541148947940343 * fZ * fTmpC + -0.99699479851094886096287209231726 * fTmpB * fR2;
  pSH[243] = fTmpA * fC0;
  pSH[237] = fTmpA * fS0;
  fTmpB    = 2.0351168037383750113438612983074 * fZ * fTmpA + -0.99755389786357772280477387849551 * fTmpC * fR2;
  pSH[275] = fTmpB * fC0;
  pSH[269] = fTmpB * fS0;
  fTmpC    = 2.0310096011589900901091187979119 * fZ * fTmpB + -0.99798183447169211033548480438427 * fTmpA * fR2;
  pSH[309] = fTmpC * fC0;
  pSH[303] = fTmpC * fS0;
  fTmpA    = 2.0275875100994065630428259128237 * fZ * fTmpC + -0.9983150788368352763742577526962 * fTmpB * fR2;
  pSH[345] = fTmpA * fC0;
  pSH[339] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fTmpC    = 2.0739902137422357166710723541669 * fZ * fTmpB + -0.99440219512922986143561854266437 * fTmpA * fR2;
  pSH[244] = fTmpC * fC1;
  pSH[236] = fTmpC * fS1;
  fTmpA    = 2.064582282206257818736941378468 * fZ * fTmpC + -0.99546384960081245812197822675493 * fTmpB * fR2;
  pSH[276] = fTmpA * fC1;
  pSH[268] = fTmpA * fS1;
  fTmpB    = 2.0568833780186057910190078334978 * fZ * fTmpA + -0.99627096277343579159169878467495 * fTmpC * fR2;
  pSH[310] = fTmpB * fC1;
  pSH[302] = fTmpB * fS1;
  fTmpC    = 2.0504988306618110688143985509413 * fZ * fTmpB + -0.99689600906642312800210598000561 * fTmpA * fR2;
  pSH[346] = fTmpC * fC1;
  pSH[338] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fTmpB    = 2.120141504711419020139800961644 * fZ * fTmpA + -0.99079092984678996129579292562184 * fTmpC * fR2;
  pSH[245] = fTmpB * fC0;
  pSH[235] = fTmpB * fS0;
  fTmpC    = 2.1044171232366050514495797729708 * fZ * fTmpB + -0.99258333397093026685296598965458 * fTmpA * fR2;
  pSH[277] = fTmpC * fC0;
  pSH[267] = fTmpC * fS0;
  fTmpA    = 2.0916500663351888698003600008946 * fZ * fTmpC + -0.99393320993236341956249269014023 * fTmpB * fR2;
  pSH[311] = fTmpA * fC0;
  pSH[301] = fTmpA * fS0;
  fTmpB    = 2.081130384894172334750081510002 * fZ * fTmpA + -0.99497063031224519067405656636005 * fTmpC * fR2;
  pSH[347] = fTmpB * fC0;
  pSH[337] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fTmpA    = 2.180966243804281491170185547368 * fZ * fTmpC + -0.98588907503509976264916003785288 * fTmpB * fR2;
  pSH[246] = fTmpA * fC1;
  pSH[234] = fTmpA * fS1;
  fTmpB    = 2.1563858652847824675605897803976 * fZ * fTmpA + -0.98872959240459255315664616192706 * fTmpC * fR2;
  pSH[278] = fTmpB * fC1;
  pSH[266] = fTmpB * fS1;
  fTmpC    = 2.1366369348357590876238965016398 * fZ * fTmpB + -0.99084165280112553362114671817729 * fTmpA * fR2;
  pSH[312] = fTmpC * fC1;
  pSH[300] = fTmpC * fS1;
  fTmpA    = 2.1205017749999120852885670096555 * fZ * fTmpC + -0.99244833805276922518608107015581 * fTmpB * fR2;
  pSH[348] = fTmpA * fC1;
  pSH[336] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fTmpC    = 2.2600784378986817490121002949266 * fZ * fTmpB + -0.97926740294193723591464201261303 * fTmpA * fR2;
  pSH[247] = fTmpC * fC0;
  pSH[233] = fTmpC * fS0;
  fTmpA    = 2.2230674720995866292406334396858 * fZ * fTmpC + -0.98362403482177094963569141672366 * fTmpB * fR2;
  pSH[279] = fTmpA * fC0;
  pSH[265] = fTmpA * fS0;
  fTmpB    = 2.1937410968480305151814130359966 * fZ * fTmpA + -0.98680814882156560015943197461397 * fTmpC * fR2;
  pSH[313] = fTmpB * fC0;
  pSH[299] = fTmpB * fS0;
  fTmpC    = 2.1700439878239586252457044013298 * fZ * fTmpB + -0.98919785518075951594044528669691 * fTmpA * fR2;
  pSH[349] = fTmpC * fC0;
  pSH[335] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fTmpB    = 2.3630173363047971158562576832196 * fZ * fTmpA + -0.9702261872276653011690043804105 * fTmpC * fR2;
  pSH[248] = fTmpB * fC1;
  pSH[232] = fTmpB * fS1;
  fTmpC    = 2.3082731640774234846327089831775 * fZ * fTmpB + -0.9768329366922967098513588823927 * fTmpA * fR2;
  pSH[280] = fTmpC * fC1;
  pSH[264] = fTmpC * fS1;
  fTmpA    = 2.2656860623955237929363221160983 * fZ * fTmpC + -0.98155023315928857430062368094603 * fTmpB * fR2;
  pSH[314] = fTmpA * fC1;
  pSH[298] = fTmpA * fS1;
  fTmpB    = 2.2317637040621551359261681701796 * fZ * fTmpA + -0.98502777640009739969521526670171 * fTmpC * fR2;
  pSH[350] = fTmpB * fC1;
  pSH[334] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fTmpA    = 2.4986107250941583107772814287273 * fZ * fTmpC + -0.95756141751469769717681687626332 * fTmpB * fR2;
  pSH[249] = fTmpA * fC0;
  pSH[231] = fTmpA * fS0;
  fTmpB    = 2.4177911997760033447311955878689 * fZ * fTmpA + -0.96765421499777267555627777162464 * fTmpC * fR2;
  pSH[281] = fTmpB * fC0;
  pSH[263] = fTmpB * fS0;
  fTmpC    = 2.3564559438666820498364806724112 * fZ * fTmpB + -0.97463169858712211832086833029898 * fTmpA * fR2;
  pSH[315] = fTmpC * fC0;
  pSH[297] = fTmpC * fS0;
  fTmpA    = 2.3085099321848032225944907791515 * fZ * fTmpC + -0.97965333839290678338363252408705 * fTmpB * fR2;
  pSH[351] = fTmpA * fC0;
  pSH[333] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fTmpC    = 2.6817904466978772499811956020466 * fZ * fTmpB + -0.93903023262181382105036678287213 * fTmpA * fR2;
  pSH[250] = fTmpC * fC1;
  pSH[230] = fTmpC * fS1;
  fTmpA    = 2.5607991541103546086991654684439 * fZ * fTmpC + -0.95488413617980451758614560131555 * fTmpB * fR2;
  pSH[282] = fTmpA * fC1;
  pSH[262] = fTmpA * fS1;
  fTmpB    = 2.4720661623652209828994052998041 * fZ * fTmpA + -0.96534949193391143146009830688925 * fTmpC * fR2;
  pSH[316] = fTmpB * fC1;
  pSH[296] = fTmpB * fS1;
  fTmpC    = 2.4044230077089180938696572065183 * fZ * fTmpB + -0.97263699666048446728717005727027 * fTmpA * fR2;
  pSH[352] = fTmpC * fC1;
  pSH[332] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fTmpB    = 2.9401072717216916591054243212966 * fZ * fTmpA + -0.90994035683194807477211507595882 * fTmpC * fR2;
  pSH[251] = fTmpB * fC0;
  pSH[229] = fTmpB * fS0;
  fTmpC    = 2.7527763762750104249103083597916 * fZ * fTmpB + -0.93628433314374189411904980673285 * fTmpA * fR2;
  pSH[283] = fTmpC * fC0;
  pSH[261] = fTmpC * fS0;
  fTmpA    = 2.6220221204253788677401154627589 * fZ * fTmpC + -0.95250095250142875237531203680419 * fTmpB * fR2;
  pSH[317] = fTmpA * fC0;
  pSH[295] = fTmpA * fS0;
  fTmpB    = 2.5257296658248257995461188984976 * fZ * fTmpA + -0.9632754987646972096742604974029 * fTmpC * fR2;
  pSH[353] = fTmpB * fC0;
  pSH[331] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fTmpA    = fZ * (53.794072123058085929669935865149 * fZ2 + -5.564904012730146820250864969637 * fR2);
  pSH[252] = fTmpA * fC1;
  pSH[228] = fTmpA * fS1;
  fTmpB    = 3.0222389997200041803718933985934 * fZ * fTmpA + -0.90717582656041188333062227910908 * fTmpC * fR2;
  pSH[284] = fTmpB * fC1;
  pSH[260] = fTmpB * fS1;
  fTmpC    = 2.8223247937435036367635754483985 * fZ * fTmpB + -0.93385228435109810061348634135925 * fTmpA * fR2;
  pSH[318] = fTmpC * fC1;
  pSH[294] = fTmpC * fS1;
  fTmpA    = 2.6822461565718468646472155691995 * fZ * fTmpC + -0.95036764107299718119109543934542 * fTmpB * fR2;
  pSH[354] = fTmpA * fC1;
  pSH[330] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fTmpC    = -17.608243388844820554936521084244 * fZ2 + 0.60718080651189036396706694143077 * fR2;
  pSH[253] = fTmpC * fC0;
  pSH[227] = fTmpC * fS0;
  fTmpA    = fZ * (-60.380154942952015812568378194669 * fZ2 + 5.8432408009308402399399617888065 * fR2);
  pSH[285] = fTmpA * fC0;
  pSH[259] = fTmpA * fS0;
  fTmpB    = 3.1024184114977141491897166813985 * fZ * fTmpA + -0.90473663963430495837036646178397 * fTmpC * fR2;
  pSH[319] = fTmpB * fC0;
  pSH[293] = fTmpB * fS0;
  fTmpC    = 2.890473786367456292181049581913 * fZ * fTmpB + -0.93168406158731500391750879330743 * fTmpA * fR2;
  pSH[355] = fTmpC * fC0;
  pSH[329] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpA * fC1;
  pSH[196] = fTmpA * fS1;
  fTmpB    = 4.6241512566300120266882256458985 * fZ;
  pSH[254] = fTmpB * fC1;
  pSH[226] = fTmpB * fS1;
  fTmpC    = 19.093881509360250419912730102112 * fZ2 + -0.61593166159226614257433257693108 * fR2;
  pSH[286] = fTmpC * fC1;
  pSH[258] = fTmpC * fS1;
  fTmpA    = fZ * (67.288948373844056319303952307109 * fZ2 + -6.1171771248949142110035159802806 * fR2);
  pSH[320] = fTmpA * fC1;
  pSH[292] = fTmpA * fS1;
  fTmpB    = 3.1807526624998681277160100799861 * fZ * fTmpA + -0.90256893466271141007219169782871 * fTmpC * fR2;
  pSH[356] = fTmpB * fC1;
  pSH[328] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.84425065085737263610244154876661;
  pSH[255] = fTmpA * fC0;
  pSH[225] = fTmpA * fS0;
  fTmpB    = -4.849850753230681765781895364853 * fZ;
  pSH[287] = fTmpB * fC0;
  pSH[257] = fTmpB * fS0;
  fTmpC    = -20.602948605549728459257474710853 * fZ2 + 0.62433177592574934725022650638948 * fR2;
  pSH[321] = fTmpC * fC0;
  pSH[291] = fTmpC * fS0;
  fTmpA    = fZ * (-74.51550792153200047235328540296 * fZ2 + 6.3870435361313143262512737052816 * fR2);
  pSH[357] = fTmpA * fC0;
  pSH[327] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.85734058883802509643716829867977;
  pSH[288] = fTmpA * fC1;
  pSH[256] = fTmpA * fS1;
  fTmpB    = 5.0720953248553606105934743464303 * fZ;
  pSH[322] = fTmpB * fC1;
  pSH[290] = fTmpB * fS1;
  fTmpC    = 22.134404124649146408318478584931 * fZ2 + -0.63241154641854704026400144090125 * fR2;
  pSH[358] = fTmpC * fC1;
  pSH[326] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.86985717192062808221161840371849;
  pSH[323] = fTmpA * fC0;
  pSH[289] = fTmpA * fS0;
  fTmpB    = -5.2911346120699731671033205770982 * fZ;
  pSH[359] = fTmpB * fC0;
  pSH[325] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpC    = 0.88185576867832886122002683526588;
  pSH[360] = fTmpC * fC1;
  pSH[324] = fTmpC * fS1;
}

// Lmax == 19
template<typename T>
void SHEval19(const T fX, const T fY, const T fZ, T* pSH)
{
  T fC0, fC1, fS0, fS1, fTmpA, fTmpB, fTmpC;
  T fZ2    = fZ * fZ;
  T fR2    = fX * fX + fY * fY + fZ * fZ;
  pSH[0]   = 0.28209479177387814346198748050032;
  pSH[2]   = 0.48860251190291992159616188406979 * fZ;
  pSH[6]   = 0.94617469575756001814222789780828 * fZ2 + -0.31539156525252000603837428116538 * fR2;
  pSH[12]  = fZ * (1.8658816629505769570714773797349 * fZ2 + -1.1195289977703461743296226016398 * fR2);
  pSH[20]  = 1.9843134832984429428027334241236 * fZ * pSH[12] + -1.0062305898749053633539629615257 * pSH[6] * fR2;
  pSH[30]  = 1.9899748742132399095816630563149 * fZ * pSH[20] + -1.0028530728448139497073834935925 * pSH[12] * fR2;
  pSH[42]  = 1.9930434571835663370075950040494 * fZ * pSH[30] + -1.0015420209622192480305655215567 * pSH[20] * fR2;
  pSH[56]  = 1.9948914348241344530240221066819 * fZ * pSH[42] + -1.0009272139219581055366928290518 * pSH[30] * fR2;
  pSH[72]  = 1.9960899278339139998987572521827 * fZ * pSH[56] + -1.000600781069514794846542216078 * pSH[42] * fR2;
  pSH[90]  = 1.9969111950679364953821492978392 * fZ * pSH[72] + -1.0004114379931337589987871972141 * pSH[56] * fR2;
  pSH[110] = 1.9974984355438178915695748849579 * fZ * pSH[90] + -1.0002940744071803443092719132501 * pSH[72] * fR2;
  pSH[132] = 1.9979328159850827789905530762482 * fZ * pSH[110] + -1.0002174622185106380497371381111 * pSH[90] * fR2;
  pSH[156] = 1.9982631347136331421672841845982 * fZ * pSH[132] + -1.0001653302482984141553307155803 * pSH[110] * fR2;
  pSH[182] = 1.9985201625794737999687947227478 * fZ * pSH[156] + -1.0001286256356210523695698944024 * pSH[132] * fR2;
  pSH[210] = 1.9987240828047460813876590179916 * fZ * pSH[182] + -1.000102035610693605705533160144 * pSH[156] * fR2;
  pSH[240] = 1.9988885800753266486651932298813 * fZ * pSH[210] + -1.0000823011400101476917751108786 * pSH[182] * fR2;
  pSH[272] = 1.9990231989649344736563116309291 * fZ * pSH[240] + -1.0000673468701305762361061790777 * pSH[210] * fR2;
  pSH[306] = 1.9991347609372268759927310233238 * fZ * pSH[272] + -1.0000558082429209263942634922095 * pSH[240] * fR2;
  pSH[342] = 1.9992282461607312885921994283223 * fZ * pSH[306] + -1.0000467628422711158057978320102 * pSH[272] * fR2;
  pSH[380] = 1.9993073592865872622516276724269 * fZ * pSH[342] + -1.0000395718327849260524675667483 * pSH[306] * fR2;
  fC0      = fX;
  fS0      = fY;

  fTmpA    = -0.48860251190291992156905682975765;
  pSH[3]   = fTmpA * fC0;
  pSH[1]   = fTmpA * fS0;
  fTmpB    = -1.0925484305920790705553974353492 * fZ;
  pSH[7]   = fTmpB * fC0;
  pSH[5]   = fTmpB * fS0;
  fTmpC    = -2.2852289973223286806327386733173 * fZ2 + 0.45704579946446573612112672380103 * fR2;
  pSH[13]  = fTmpC * fC0;
  pSH[11]  = fTmpC * fS0;
  fTmpA    = fZ * (-4.6833258049010241751662630971254 * fZ2 + 2.0071396306718675037975702091231 * fR2);
  pSH[21]  = fTmpA * fC0;
  pSH[19]  = fTmpA * fS0;
  fTmpB    = 2.0310096011589900901091187979119 * fZ * fTmpA + -0.99103120896511485331606405857485 * fTmpC * fR2;
  pSH[31]  = fTmpB * fC0;
  pSH[29]  = fTmpB * fS0;
  fTmpC    = 2.021314989237027776210198215523 * fZ * fTmpB + -0.99522670305623857703314349976154 * fTmpA * fR2;
  pSH[43]  = fTmpC * fC0;
  pSH[41]  = fTmpC * fS0;
  fTmpA    = 2.0155644370746374129804712183045 * fZ * fTmpC + -0.99715504402183205234040316855548 * fTmpB * fR2;
  pSH[57]  = fTmpA * fC0;
  pSH[55]  = fTmpA * fS0;
  fTmpB    = 2.0118695404073912316820355039582 * fZ * fTmpA + -0.99816681789017427597768619684793 * fTmpC * fR2;
  pSH[73]  = fTmpB * fC0;
  pSH[71]  = fTmpB * fS0;
  fTmpC    = 2.009353129741011949457515917139 * fZ * fTmpB + -0.99874921777190894578478744247896 * fTmpA * fR2;
  pSH[91]  = fTmpC * fC0;
  pSH[89]  = fTmpC * fS0;
  fTmpA    = 2.0075614636426527856265245031153 * fZ * fTmpC + -0.99910833687128449445659025829336 * fTmpB * fR2;
  pSH[111] = fTmpA * fC0;
  pSH[109] = fTmpA * fS0;
  fTmpB    = 2.0062402647738879432354891507728 * fZ * fTmpA + -0.99934188870792151407600883983307 * fTmpC * fR2;
  pSH[133] = fTmpB * fC0;
  pSH[131] = fTmpB * fS0;
  fTmpC    = 2.0052378963551982949536922617995 * fZ * fTmpB + -0.99950037468777319166679182216306 * fTmpA * fR2;
  pSH[157] = fTmpC * fC0;
  pSH[155] = fTmpC * fS0;
  fTmpA    = 2.0044593143431828851340481545407 * fZ * fTmpC + -0.99961172586383361700476321565212 * fTmpB * fR2;
  pSH[183] = fTmpA * fC0;
  pSH[181] = fTmpA * fS0;
  fTmpB    = 2.0038424627162224563904635576961 * fZ * fTmpA + -0.9996922603404586649177357426943 * fTmpC * fR2;
  pSH[211] = fTmpB * fC0;
  pSH[209] = fTmpB * fS0;
  fTmpC    = 2.003345416333103836325699176335 * fZ * fTmpB + -0.99975195336341716702754575663015 * fTmpA * fR2;
  pSH[241] = fTmpC * fC0;
  pSH[239] = fTmpC * fS0;
  fTmpA    = 2.0029390170153341292971771459008 * fZ * fTmpC + -0.99979713966725040625964371354684 * fTmpB * fR2;
  pSH[273] = fTmpA * fC0;
  pSH[271] = fTmpA * fS0;
  fTmpB    = 2.0026024734496526301299329508865 * fZ * fTmpA + -0.99983197513113354919136316345529 * fTmpC * fR2;
  pSH[307] = fTmpB * fC0;
  pSH[305] = fTmpB * fS0;
  fTmpC    = 2.0023206350873464509643184783272 * fZ * fTmpB + -0.99985926394976398475303303037265 * fTmpA * fR2;
  pSH[343] = fTmpC * fC0;
  pSH[341] = fTmpC * fS0;
  fTmpA    = 2.0020822493926999834177454795636 * fZ * fTmpC + -0.99988094529394086349338363617356 * fTmpB * fR2;
  pSH[381] = fTmpA * fC0;
  pSH[379] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.54627421529603953527769871767461;
  pSH[8]   = fTmpA * fC1;
  pSH[4]   = fTmpA * fS1;
  fTmpB    = 1.4453057213202770277206063442854 * fZ;
  pSH[14]  = fTmpB * fC1;
  pSH[10]  = fTmpB * fS1;
  fTmpC    = 3.3116114351514600632267470992076 * fZ2 + -0.473087347878780009044008894592 * fR2;
  pSH[22]  = fTmpC * fC1;
  pSH[18]  = fTmpC * fS1;
  fTmpA    = fZ * (7.1903051774599856323273716629529 * fZ2 + -2.3967683924866618773701770761519 * fR2);
  pSH[32]  = fTmpA * fC1;
  pSH[28]  = fTmpA * fS1;
  fTmpB    = 2.1139418156609703623206863998263 * fZ * fTmpA + -0.97361012046232688413316402886721 * fTmpC * fR2;
  pSH[44]  = fTmpB * fC1;
  pSH[40]  = fTmpB * fS1;
  fTmpC    = 2.0816659994661327352447055538676 * fZ * fTmpB + -0.98473192783466186186287424919605 * fTmpA * fR2;
  pSH[58]  = fTmpC * fC1;
  pSH[54]  = fTmpC * fS1;
  fTmpA    = 2.0615528128088302748256027685514 * fZ * fTmpC + -0.99033793766028713591656473802516 * fTmpB * fR2;
  pSH[74]  = fTmpA * fC1;
  pSH[70]  = fTmpA * fS1;
  fTmpB    = 2.0481223583578191105563498508602 * fZ * fTmpA + -0.99348527267040401402898100458039 * fTmpC * fR2;
  pSH[92]  = fTmpB * fC1;
  pSH[88]  = fTmpB * fS1;
  fTmpC    = 2.0386883037875113094776480249237 * fZ * fTmpB + -0.99539380324041186210446210957947 * fTmpA * fR2;
  pSH[112] = fTmpC * fC1;
  pSH[108] = fTmpC * fS1;
  fTmpA    = 2.0317984959648750084224011480671 * fZ * fTmpC + -0.99662047022596034228682226885354 * fTmpB * fR2;
  pSH[134] = fTmpA * fC1;
  pSH[130] = fTmpA * fS1;
  fTmpB    = 2.0266087084444439302393509150235 * fZ * fTmpA + -0.99744571741206722651201105334096 * fTmpC * fR2;
  pSH[158] = fTmpB * fC1;
  pSH[154] = fTmpB * fS1;
  fTmpC    = 2.0225995873897262588344408973384 * fZ * fTmpB + -0.99802175869569072528620160000834 * fTmpA * fR2;
  pSH[184] = fTmpC * fC1;
  pSH[180] = fTmpC * fS1;
  fTmpA    = 2.0194368026754390117814136340613 * fZ * fTmpC + -0.99843627738579290858836681743504 * fTmpB * fR2;
  pSH[212] = fTmpA * fC1;
  pSH[208] = fTmpA * fS1;
  fTmpB    = 2.0168969490698876096713976213692 * fZ * fTmpA + -0.99874229606879180774856724633892 * fTmpC * fR2;
  pSH[242] = fTmpB * fC1;
  pSH[238] = fTmpB * fS1;
  fTmpC    = 2.014825999813336120320209077228 * fZ * fTmpB + -0.99897320026315349013367253516726 * fTmpA * fR2;
  pSH[274] = fTmpC * fC1;
  pSH[270] = fTmpC * fS1;
  fTmpA    = 2.0131148946216081338771858311176 * fZ * fTmpC + -0.9991507429465936452315545646119 * fTmpB * fR2;
  pSH[308] = fTmpA * fC1;
  pSH[304] = fTmpA * fS1;
  fTmpB    = 2.0116846174288851485396217855239 * fZ * fTmpA + -0.99928952033659667242206786630376 * fTmpC * fR2;
  pSH[344] = fTmpB * fC1;
  pSH[340] = fTmpB * fS1;
  fTmpC    = 2.0104767610501468005738956446038 * fZ * fTmpB + -0.99939957965166423677029136629635 * fTmpA * fR2;
  pSH[382] = fTmpC * fC1;
  pSH[378] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.59004358992664351034589803601804;
  pSH[15]  = fTmpA * fC0;
  pSH[9]   = fTmpA * fS0;
  fTmpB    = -1.7701307697799305311461143253027 * fZ;
  pSH[23]  = fTmpB * fC0;
  pSH[17]  = fTmpB * fS0;
  fTmpC    = -4.4031446949172534889235808286401 * fZ2 + 0.48923829943525038764311728411993 * fR2;
  pSH[33]  = fTmpC * fC0;
  pSH[27]  = fTmpC * fS0;
  fTmpA    = fZ * (-10.13325785466415849059296228063 * fZ2 + 2.7636157785447704976315719260782 * fR2);
  pSH[45]  = fTmpA * fC0;
  pSH[39]  = fTmpA * fS0;
  fTmpB    = 2.2079402165819617136124225487137 * fZ * fTmpA + -0.95940322360024695439372974248293 * fTmpC * fR2;
  pSH[59]  = fTmpB * fC0;
  pSH[53]  = fTmpB * fS0;
  fTmpC    = 2.1532216876958202239969453195556 * fZ * fTmpB + -0.97521738656001772939082086755214 * fTmpA * fR2;
  pSH[75]  = fTmpC * fC0;
  pSH[69]  = fTmpC * fS0;
  fTmpA    = 2.1180441711898057371380593716381 * fZ * fTmpC + -0.98366284497920962828255986298842 * fTmpB * fR2;
  pSH[93]  = fTmpA * fC0;
  pSH[87]  = fTmpA * fS0;
  fTmpB    = 2.093947321356338375558311937219 * fZ * fTmpA + -0.98862306548596150408952237809146 * fTmpC * fR2;
  pSH[113] = fTmpB * fC0;
  pSH[107] = fTmpB * fS0;
  fTmpC    = 2.0766559657295187131878511088701 * fZ * fTmpB + -0.99174222032690902703735980061595 * fTmpA * fR2;
  pSH[135] = fTmpC * fC0;
  pSH[129] = fTmpC * fS0;
  fTmpA    = 2.0637972912229677745644257358393 * fZ * fTmpC + -0.99380798999990653172769555778743 * fTmpB * fR2;
  pSH[159] = fTmpA * fC0;
  pSH[153] = fTmpA * fS0;
  fTmpB    = 2.0539595906443729254452906785033 * fZ * fTmpA + -0.99523320404555565076671827529076 * fTmpC * fR2;
  pSH[185] = fTmpB * fC0;
  pSH[179] = fTmpB * fS0;
  fTmpC    = 2.0462565272714634686482965131304 * fZ * fTmpB + -0.99624965193668058786695407302858 * fTmpA * fR2;
  pSH[213] = fTmpC * fC0;
  pSH[207] = fTmpC * fS0;
  fTmpA    = 2.0401071141087266541148947940343 * fZ * fTmpC + -0.99699479851094886096287209231726 * fTmpB * fR2;
  pSH[243] = fTmpA * fC0;
  pSH[237] = fTmpA * fS0;
  fTmpB    = 2.0351168037383750113438612983074 * fZ * fTmpA + -0.99755389786357772280477387849551 * fTmpC * fR2;
  pSH[275] = fTmpB * fC0;
  pSH[269] = fTmpB * fS0;
  fTmpC    = 2.0310096011589900901091187979119 * fZ * fTmpB + -0.99798183447169211033548480438427 * fTmpA * fR2;
  pSH[309] = fTmpC * fC0;
  pSH[303] = fTmpC * fS0;
  fTmpA    = 2.0275875100994065630428259128237 * fZ * fTmpC + -0.9983150788368352763742577526962 * fTmpB * fR2;
  pSH[345] = fTmpA * fC0;
  pSH[339] = fTmpA * fS0;
  fTmpB    = 2.02470536577098500605371989014 * fZ * fTmpA + -0.99857853517341885088258457425781 * fTmpC * fR2;
  pSH[383] = fTmpB * fC0;
  pSH[377] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.62583573544917613459435609679637;
  pSH[24]  = fTmpA * fC1;
  pSH[16]  = fTmpA * fS1;
  fTmpB    = 2.075662314881041278823506357476 * fZ;
  pSH[34]  = fTmpB * fC1;
  pSH[26]  = fTmpB * fS1;
  fTmpC    = 5.5502139080159657515307902730939 * fZ2 + -0.50456490072872415923992822639477 * fR2;
  pSH[46]  = fTmpC * fC1;
  pSH[38]  = fTmpC * fS1;
  fTmpA    = fZ * (13.491805046726768313111732844334 * fZ2 + -3.1134934723215619183436797534625 * fR2);
  pSH[60]  = fTmpA * fC1;
  pSH[52]  = fTmpA * fS1;
  fTmpB    = 2.3048861143232218272932504410377 * fZ * fTmpA + -0.94817638735546538522732870624132 * fTmpC * fR2;
  pSH[76]  = fTmpB * fC1;
  pSH[68]  = fTmpB * fS1;
  fTmpC    = 2.2291771507062351975400615877732 * fZ * fTmpB + -0.96715283972318221404618210357285 * fTmpA * fR2;
  pSH[94]  = fTmpC * fC1;
  pSH[86]  = fTmpC * fS1;
  fTmpA    = 2.179449471770336776293985892039 * fZ * fTmpC + -0.97769236109380361197549944018981 * fTmpB * fR2;
  pSH[114] = fTmpA * fC1;
  pSH[106] = fTmpA * fS1;
  fTmpB    = 2.1447610589527216608599080593933 * fZ * fTmpA + -0.9840838646332836542503230692347 * fTmpC * fR2;
  pSH[136] = fTmpB * fC1;
  pSH[128] = fTmpB * fS1;
  fTmpC    = 2.119478119726646293515329166901 * fZ * fTmpB + -0.98821176880261854123983777942186 * fTmpA * fR2;
  pSH[160] = fTmpC * fC1;
  pSH[152] = fTmpC * fS1;
  fTmpA    = 2.1004201260420147052629391559719 * fZ * fTmpC + -0.99100816681840077734160290856558 * fTmpB * fR2;
  pSH[186] = fTmpA * fC1;
  pSH[178] = fTmpA * fS1;
  fTmpB    = 2.0856653614614210205148447929702 * fZ * fTmpA + -0.99297532698451274258593171606613 * fTmpC * fR2;
  pSH[214] = fTmpB * fC1;
  pSH[206] = fTmpB * fS1;
  fTmpC    = 2.0739902137422357166710723541669 * fZ * fTmpB + -0.99440219512922986143561854266437 * fTmpA * fR2;
  pSH[244] = fTmpC * fC1;
  pSH[236] = fTmpC * fS1;
  fTmpA    = 2.064582282206257818736941378468 * fZ * fTmpC + -0.99546384960081245812197822675493 * fTmpB * fR2;
  pSH[276] = fTmpA * fC1;
  pSH[268] = fTmpA * fS1;
  fTmpB    = 2.0568833780186057910190078334978 * fZ * fTmpA + -0.99627096277343579159169878467495 * fTmpC * fR2;
  pSH[310] = fTmpB * fC1;
  pSH[302] = fTmpB * fS1;
  fTmpC    = 2.0504988306618110688143985509413 * fZ * fTmpB + -0.99689600906642312800210598000561 * fTmpA * fR2;
  pSH[346] = fTmpC * fC1;
  pSH[338] = fTmpC * fS1;
  fTmpA    = 2.045142707894041782679811491974 * fZ * fTmpC + -0.99738789279580297795612525701969 * fTmpB * fR2;
  pSH[384] = fTmpA * fC1;
  pSH[376] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.65638205684017010278471018769331;
  pSH[35]  = fTmpA * fC0;
  pSH[25]  = fTmpA * fS0;
  fTmpB    = -2.3666191622317520320151890134142 * fZ;
  pSH[47]  = fTmpB * fC0;
  pSH[37]  = fTmpB * fS0;
  fTmpC    = -6.7459025233633841565558664221669 * fZ2 + 0.518915578720260319705876589369 * fR2;
  pSH[61]  = fTmpC * fC0;
  pSH[51]  = fTmpC * fS0;
  fTmpA    = fZ * (-17.249553110490540087582078676576 * fZ2 + 3.4499106220981080176465199960134 * fR2);
  pSH[77]  = fTmpA * fC0;
  pSH[67]  = fTmpA * fS0;
  fTmpB    = 2.4016363469220611496363071424298 * fZ * fTmpA + -0.93922460420437088487706847605985 * fTmpC * fR2;
  pSH[95]  = fTmpB * fC0;
  pSH[85]  = fTmpB * fS0;
  fTmpC    = 2.306512518934159177829562592521 * fZ * fTmpB + -0.9603920767980494892571141640758 * fTmpA * fR2;
  pSH[115] = fTmpC * fC0;
  pSH[105] = fTmpC * fS0;
  fTmpA    = 2.2430448056157950944330264908544 * fZ * fTmpC + -0.9724832565193738675448330288642 * fTmpB * fR2;
  pSH[137] = fTmpA * fC0;
  pSH[127] = fTmpA * fS0;
  fTmpB    = 2.1981657747106435412363933945556 * fZ * fTmpA + -0.97999191510005049935479876088706 * fTmpC * fR2;
  pSH[161] = fTmpB * fC0;
  pSH[151] = fTmpB * fS0;
  fTmpC    = 2.1650635094610966170213667281175 * fZ * fTmpB + -0.98494096049061433725130276783943 * fTmpA * fR2;
  pSH[187] = fTmpC * fC0;
  pSH[177] = fTmpC * fS0;
  fTmpA    = 2.1398475105532759976776496779749 * fZ * fTmpC + -0.98835322899414756493696732064791 * fTmpB * fR2;
  pSH[215] = fTmpA * fC0;
  pSH[205] = fTmpA * fS0;
  fTmpB    = 2.120141504711419020139800961644 * fZ * fTmpA + -0.99079092984678996129579292562184 * fTmpC * fR2;
  pSH[245] = fTmpB * fC0;
  pSH[235] = fTmpB * fS0;
  fTmpC    = 2.1044171232366050514495797729708 * fZ * fTmpB + -0.99258333397093026685296598965458 * fTmpA * fR2;
  pSH[277] = fTmpC * fC0;
  pSH[267] = fTmpC * fS0;
  fTmpA    = 2.0916500663351888698003600008946 * fZ * fTmpC + -0.99393320993236341956249269014023 * fTmpB * fR2;
  pSH[311] = fTmpA * fC0;
  pSH[301] = fTmpA * fS0;
  fTmpB    = 2.081130384894172334750081510002 * fZ * fTmpA + -0.99497063031224519067405656636005 * fTmpC * fR2;
  pSH[347] = fTmpB * fC0;
  pSH[337] = fTmpB * fS0;
  fTmpC    = 2.0723520109148583403634730215614 * fZ * fTmpB + -0.99578192022804934251112296550446 * fTmpA * fR2;
  pSH[385] = fTmpC * fC0;
  pSH[375] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.68318410519191432195978963548555;
  pSH[48]  = fTmpA * fC1;
  pSH[36]  = fTmpA * fS1;
  fTmpB    = 2.6459606618019002217956359146456 * fZ;
  pSH[62]  = fTmpB * fC1;
  pSH[50]  = fTmpB * fS1;
  fTmpC    = 7.9849914908931386142015851348219 * fZ2 + -0.53233276605954257430178971910451 * fR2;
  pSH[78]  = fTmpC * fC1;
  pSH[66]  = fTmpC * fS1;
  fTmpA    = fZ * (21.392890190908636255037733597817 * fZ2 + -3.7752159160427005153270324511183 * fR2);
  pSH[96]  = fTmpA * fC1;
  pSH[84]  = fTmpA * fS1;
  fTmpB    = 2.4968730444297723642180231173882 * fZ * fTmpA + -0.93196897827695329100811116873615 * fTmpC * fR2;
  pSH[116] = fTmpB * fC1;
  pSH[104] = fTmpB * fS1;
  fTmpC    = 2.3837686425440851888980786643657 * fZ * fTmpB + -0.95470158078801415946516503718833 * fTmpA * fR2;
  pSH[138] = fTmpC * fC1;
  pSH[126] = fTmpC * fS1;
  fTmpA    = 2.307395517477243014484514227469 * fZ * fTmpC + -0.96796118394051333064806788564205 * fTmpB * fR2;
  pSH[162] = fTmpA * fC1;
  pSH[150] = fTmpA * fS1;
  fTmpB    = 2.2528177844479149152714242410056 * fZ * fTmpA + -0.97634660697919710712960536524996 * fTmpC * fR2;
  pSH[188] = fTmpB * fC1;
  pSH[176] = fTmpB * fS1;
  fTmpC    = 2.2121821805628938751361184378297 * fZ * fTmpB + -0.98196232106939826308469529414502 * fTmpA * fR2;
  pSH[216] = fTmpC * fC1;
  pSH[204] = fTmpC * fS1;
  fTmpA    = 2.180966243804281491170185547368 * fZ * fTmpC + -0.98588907503509976264916003785288 * fTmpB * fR2;
  pSH[246] = fTmpA * fC1;
  pSH[234] = fTmpA * fS1;
  fTmpB    = 2.1563858652847824675605897803976 * fZ * fTmpA + -0.98872959240459255315664616192706 * fTmpC * fR2;
  pSH[278] = fTmpB * fC1;
  pSH[266] = fTmpB * fS1;
  fTmpC    = 2.1366369348357590876238965016398 * fZ * fTmpB + -0.99084165280112553362114671817729 * fTmpA * fR2;
  pSH[312] = fTmpC * fC1;
  pSH[300] = fTmpC * fS1;
  fTmpA    = 2.1205017749999120852885670096555 * fZ * fTmpC + -0.99244833805276922518608107015581 * fTmpB * fR2;
  pSH[348] = fTmpA * fC1;
  pSH[336] = fTmpA * fS1;
  fTmpB    = 2.1071307505705477694686600376173 * fZ * fTmpA + -0.99369440545299007368088353708835 * fTmpC * fR2;
  pSH[386] = fTmpB * fC1;
  pSH[374] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.70716273252459617822346729654193;
  pSH[63]  = fTmpA * fC0;
  pSH[49]  = fTmpA * fS0;
  fTmpB    = -2.9157066406993194756895604324853 * fZ;
  pSH[79]  = fTmpB * fC0;
  pSH[65]  = fTmpB * fS0;
  fTmpC    = -9.263393182848904240309084734406 * fZ2 + 0.54490548134405319058090610973011 * fR2;
  pSH[97]  = fTmpC * fC0;
  pSH[83]  = fTmpC * fS0;
  fTmpA    = fZ * (-25.91024131336630202497584019028 * fZ2 + 4.0910907336894161089556332111528 * fR2);
  pSH[117] = fTmpA * fC0;
  pSH[103] = fTmpA * fS0;
  fTmpB    = 2.5900450446533421888090087392698 * fZ * fTmpA + -0.92598927658525138716496052926352 * fTmpC * fR2;
  pSH[139] = fTmpB * fC0;
  pSH[125] = fTmpB * fS0;
  fTmpC    = 2.4602096615832091356258076730867 * fZ * fTmpB + -0.94987138029195529654495969151817 * fTmpA * fR2;
  pSH[163] = fTmpC * fC0;
  pSH[149] = fTmpC * fS0;
  fTmpA    = 2.3717082451262844991924511051096 * fZ * fTmpC + -0.96402688037572713822742284661693 * fTmpB * fR2;
  pSH[189] = fTmpA * fC0;
  pSH[175] = fTmpA * fS0;
  fTmpB    = 2.3079277744862160132166550852162 * fZ * fTmpA + -0.97310779233865149528189680827595 * fTmpC * fR2;
  pSH[217] = fTmpB * fC0;
  pSH[203] = fTmpB * fS0;
  fTmpC    = 2.2600784378986817490121002949266 * fZ * fTmpB + -0.97926740294193723591464201261303 * fTmpA * fR2;
  pSH[247] = fTmpC * fC0;
  pSH[233] = fTmpC * fS0;
  fTmpA    = 2.2230674720995866292406334396858 * fZ * fTmpC + -0.98362403482177094963569141672366 * fTmpB * fR2;
  pSH[279] = fTmpA * fC0;
  pSH[265] = fTmpA * fS0;
  fTmpB    = 2.1937410968480305151814130359966 * fZ * fTmpA + -0.98680814882156560015943197461397 * fTmpC * fR2;
  pSH[313] = fTmpB * fC0;
  pSH[299] = fTmpB * fS0;
  fTmpC    = 2.1700439878239586252457044013298 * fZ * fTmpB + -0.98919785518075951594044528669691 * fTmpA * fR2;
  pSH[349] = fTmpC * fC0;
  pSH[335] = fTmpC * fS0;
  fTmpA    = 2.1505813167606566930090822298283 * fZ * fTmpC + -0.99103120896511485326185394995058 * fTmpB * fR2;
  pSH[387] = fTmpA * fC0;
  pSH[373] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.72892666017482986886817999949706;
  pSH[80]  = fTmpA * fC1;
  pSH[64]  = fTmpA * fS1;
  fTmpB    = 3.1773176489546974771236570456168 * fZ;
  pSH[98]  = fTmpB * fC1;
  pSH[82]  = fTmpB * fS1;
  fTmpC    = 10.577811721687949635117842461796 * fZ2 + -0.55672693272041840183974800715383 * fR2;
  pSH[118] = fTmpC * fC1;
  pSH[102] = fTmpC * fS1;
  fTmpA    = fZ * (30.791579703357485663442472123563 * fZ2 + -4.3987971004796408091251647132225 * fR2);
  pSH[140] = fTmpA * fC1;
  pSH[124] = fTmpA * fS1;
  fTmpB    = 2.6809513236909020762917255087388 * fZ * fTmpA + -0.92098549701625905379792982885512 * fTmpC * fR2;
  pSH[164] = fTmpB * fC1;
  pSH[148] = fTmpB * fS1;
  fTmpC    = 2.5354627641855497323982587820623 * fZ * fTmpB + -0.94573248748692077369614600312886 * fTmpA * fR2;
  pSH[190] = fTmpC * fC1;
  pSH[174] = fTmpC * fS1;
  fTmpA    = 2.4355324226579661401302645540312 * fZ * fTmpC + -0.96058694178469484717007229046537 * fTmpB * fR2;
  pSH[218] = fTmpA * fC1;
  pSH[202] = fTmpA * fS1;
  fTmpB    = 2.3630173363047971158562576832196 * fZ * fTmpA + -0.9702261872276653011690043804105 * fTmpC * fR2;
  pSH[248] = fTmpB * fC1;
  pSH[232] = fTmpB * fS1;
  fTmpC    = 2.3082731640774234846327089831775 * fZ * fTmpB + -0.9768329366922967098513588823927 * fTmpA * fR2;
  pSH[280] = fTmpC * fC1;
  pSH[264] = fTmpC * fS1;
  fTmpA    = 2.2656860623955237929363221160983 * fZ * fTmpC + -0.98155023315928857430062368094603 * fTmpB * fR2;
  pSH[314] = fTmpA * fC1;
  pSH[298] = fTmpA * fS1;
  fTmpB    = 2.2317637040621551359261681701796 * fZ * fTmpA + -0.98502777640009739969521526670171 * fTmpC * fR2;
  pSH[350] = fTmpB * fC1;
  pSH[334] = fTmpB * fS1;
  fTmpC    = 2.2042200113840402628164610865369 * fZ * fTmpB + -0.98765832931686222317859552566333 * fTmpA * fR2;
  pSH[388] = fTmpC * fC1;
  pSH[372] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.74890095185318829654197436695995;
  pSH[99]  = fTmpA * fC0;
  pSH[81]  = fTmpA * fS0;
  fTmpB    = -3.4318952998917144348034469203412 * fZ;
  pSH[119] = fTmpB * fC0;
  pSH[101] = fTmpB * fS0;
  fTmpC    = -11.925527539452185581299314964809 * fZ2 + 0.56788226378343740856618343526563 * fR2;
  pSH[141] = fTmpC * fC0;
  pSH[123] = fTmpC * fS0;
  fTmpA    = fZ * (-36.028090689310769891007257825777 * fZ2 + 4.6993161768666221601679910957472 * fR2);
  pSH[165] = fTmpA * fC0;
  pSH[147] = fTmpA * fS0;
  fTmpB    = 2.7695585470349864862758815231558 * fZ * fTmpA + -0.91674152287482094271639163074461 * fTmpC * fR2;
  pSH[191] = fTmpB * fC0;
  pSH[173] = fTmpB * fS0;
  fTmpC    = 2.6093477445855914596296865060054 * fZ * fTmpB + -0.94215294613615865838493479422766 * fTmpA * fR2;
  pSH[219] = fTmpC * fC0;
  pSH[201] = fTmpC * fS0;
  fTmpA    = 2.4986107250941583107772814287273 * fZ * fTmpC + -0.95756141751469769717681687626332 * fTmpB * fR2;
  pSH[249] = fTmpA * fC0;
  pSH[231] = fTmpA * fS0;
  fTmpB    = 2.4177911997760033447311955878689 * fZ * fTmpA + -0.96765421499777267555627777162464 * fTmpC * fR2;
  pSH[281] = fTmpB * fC0;
  pSH[263] = fTmpB * fS0;
  fTmpC    = 2.3564559438666820498364806724112 * fZ * fTmpB + -0.97463169858712211832086833029898 * fTmpA * fR2;
  pSH[315] = fTmpC * fC0;
  pSH[297] = fTmpC * fS0;
  fTmpA    = 2.3085099321848032225944907791515 * fZ * fTmpC + -0.97965333839290678338363252408705 * fTmpB * fR2;
  pSH[351] = fTmpA * fC0;
  pSH[333] = fTmpA * fS0;
  fTmpB    = 2.2701478869385202613399854509879 * fZ * fTmpA + -0.98338233476432278864099584270164 * fTmpC * fR2;
  pSH[389] = fTmpB * fC0;
  pSH[371] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.76739511822199001257619258020704;
  pSH[120] = fTmpA * fC1;
  pSH[100] = fTmpA * fS1;
  fTmpB    = 3.6802976988053108635141202897856 * fZ;
  pSH[142] = fTmpB * fC1;
  pSH[122] = fTmpB * fS1;
  fTmpC    = 13.304254200257634744956891648116 * fZ2 + -0.57844583479381020635432669729781 * fR2;
  pSH[166] = fTmpC * fC1;
  pSH[146] = fTmpC * fS1;
  fTmpA    = fZ * (41.611931535496447638611261510277 * fZ2 + -4.9934317842595737162690594512782 * fR2);
  pSH[192] = fTmpA * fC1;
  pSH[172] = fTmpA * fS1;
  fTmpB    = 2.855914914698965607160394131192 * fZ * fTmpA + -0.91309911838748371456196684103901 * fTmpC * fR2;
  pSH[220] = fTmpB * fC1;
  pSH[200] = fTmpB * fS1;
  fTmpC    = 2.6817904466978772499811956020466 * fZ * fTmpB + -0.93903023262181382105036678287213 * fTmpA * fR2;
  pSH[250] = fTmpC * fC1;
  pSH[230] = fTmpC * fS1;
  fTmpA    = 2.5607991541103546086991654684439 * fZ * fTmpC + -0.95488413617980451758614560131555 * fTmpB * fR2;
  pSH[282] = fTmpA * fC1;
  pSH[262] = fTmpA * fS1;
  fTmpB    = 2.4720661623652209828994052998041 * fZ * fTmpA + -0.96534949193391143146009830688925 * fTmpC * fR2;
  pSH[316] = fTmpB * fC1;
  pSH[296] = fTmpB * fS1;
  fTmpC    = 2.4044230077089180938696572065183 * fZ * fTmpB + -0.97263699666048446728717005727027 * fTmpA * fR2;
  pSH[352] = fTmpC * fC1;
  pSH[332] = fTmpC * fS1;
  fTmpA    = 2.3513263559497452387294508246995 * fZ * fTmpC + -0.97791709213023767974758138077362 * fTmpB * fR2;
  pSH[390] = fTmpA * fC1;
  pSH[370] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.78464210578719688376396826368087;
  pSH[143] = fTmpA * fC0;
  pSH[121] = fTmpA * fS0;
  fTmpB    = -3.9232105289359844187656312097801 * fZ;
  pSH[167] = fTmpB * fC0;
  pSH[145] = fTmpB * fS0;
  fTmpC    = -14.712039483509941569829015950432 * fZ2 + 0.58848157934039766278231861629244 * fR2;
  pSH[193] = fTmpC * fC0;
  pSH[171] = fTmpC * fS0;
  fTmpA    = fZ * (-47.536054360662613679777699360329 * fZ2 + 5.2817838178514015198307396392607 * fR2);
  pSH[221] = fTmpA * fC0;
  pSH[199] = fTmpA * fS0;
  fTmpB    = 2.9401072717216916591054243212966 * fZ * fTmpA + -0.90994035683194807477211507595882 * fTmpC * fR2;
  pSH[251] = fTmpB * fC0;
  pSH[229] = fTmpB * fS0;
  fTmpC    = 2.7527763762750104249103083597916 * fZ * fTmpB + -0.93628433314374189411904980673285 * fTmpA * fR2;
  pSH[283] = fTmpC * fC0;
  pSH[261] = fTmpC * fS0;
  fTmpA    = 2.6220221204253788677401154627589 * fZ * fTmpC + -0.95250095250142875237531203680419 * fTmpB * fR2;
  pSH[317] = fTmpA * fC0;
  pSH[295] = fTmpA * fS0;
  fTmpB    = 2.5257296658248257995461188984976 * fZ * fTmpA + -0.9632754987646972096742604974029 * fTmpC * fR2;
  pSH[353] = fTmpB * fC0;
  pSH[331] = fTmpB * fS0;
  fTmpC    = 2.4520399670478456538667139108512 * fZ * fTmpB + -0.97082439194737994584208373716194 * fTmpA * fR2;
  pSH[391] = fTmpC * fC0;
  pSH[369] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.80082199578397171406147353467375;
  pSH[168] = fTmpA * fC1;
  pSH[144] = fTmpA * fS1;
  fTmpB    = 4.1611931535496447637743899772289 * fZ;
  pSH[194] = fTmpB * fC1;
  pSH[170] = fTmpB * fS1;
  fTmpC    = 16.147194793928202583357944810416 * fZ2 + -0.59804425162697046606434872484392 * fR2;
  pSH[222] = fTmpC * fC1;
  pSH[198] = fTmpC * fS1;
  fTmpA    = fZ * (53.794072123058085929669935865149 * fZ2 + -5.564904012730146820250864969637 * fR2);
  pSH[252] = fTmpA * fC1;
  pSH[228] = fTmpA * fS1;
  fTmpB    = 3.0222389997200041803718933985934 * fZ * fTmpA + -0.90717582656041188333062227910908 * fTmpC * fR2;
  pSH[284] = fTmpB * fC1;
  pSH[260] = fTmpB * fS1;
  fTmpC    = 2.8223247937435036367635754483985 * fZ * fTmpB + -0.93385228435109810061348634135925 * fTmpA * fR2;
  pSH[318] = fTmpC * fC1;
  pSH[294] = fTmpC * fS1;
  fTmpA    = 2.6822461565718468646472155691995 * fZ * fTmpC + -0.95036764107299718119109543934542 * fTmpB * fR2;
  pSH[354] = fTmpA * fC1;
  pSH[330] = fTmpA * fS1;
  fTmpB    = 2.5787147157554005437538752198989 * fZ * fTmpA + -0.96140121570767059499011339407382 * fTmpC * fR2;
  pSH[392] = fTmpB * fC1;
  pSH[368] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.81607711883762830927784709400541;
  pSH[195] = fTmpA * fC0;
  pSH[169] = fTmpA * fS0;
  fTmpB    = -4.3947097802721183800074566949689 * fZ;
  pSH[223] = fTmpB * fC0;
  pSH[197] = fTmpB * fS0;
  fTmpC    = -17.608243388844820554936521084244 * fZ2 + 0.60718080651189036396706694143077 * fR2;
  pSH[253] = fTmpC * fC0;
  pSH[227] = fTmpC * fS0;
  fTmpA    = fZ * (-60.380154942952015812568378194669 * fZ2 + 5.8432408009308402399399617888065 * fR2);
  pSH[285] = fTmpA * fC0;
  pSH[259] = fTmpA * fS0;
  fTmpB    = 3.1024184114977141491897166813985 * fZ * fTmpA + -0.90473663963430495837036646178397 * fTmpC * fR2;
  pSH[319] = fTmpB * fC0;
  pSH[293] = fTmpB * fS0;
  fTmpC    = 2.890473786367456292181049581913 * fZ * fTmpB + -0.93168406158731500391750879330743 * fTmpA * fR2;
  pSH[355] = fTmpC * fC0;
  pSH[329] = fTmpC * fS0;
  fTmpA    = 2.7414640249326636017076358475819 * fZ * fTmpC + -0.94844798034925006024672966553624 * fTmpB * fR2;
  pSH[393] = fTmpA * fC0;
  pSH[367] = fTmpA * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.83052208306452400299628099911153;
  pSH[224] = fTmpA * fC1;
  pSH[196] = fTmpA * fS1;
  fTmpB    = 4.6241512566300120266882256458985 * fZ;
  pSH[254] = fTmpB * fC1;
  pSH[226] = fTmpB * fS1;
  fTmpC    = 19.093881509360250419912730102112 * fZ2 + -0.61593166159226614257433257693108 * fR2;
  pSH[286] = fTmpC * fC1;
  pSH[258] = fTmpC * fS1;
  fTmpA    = fZ * (67.288948373844056319303952307109 * fZ2 + -6.1171771248949142110035159802806 * fR2);
  pSH[320] = fTmpA * fC1;
  pSH[292] = fTmpA * fS1;
  fTmpB    = 3.1807526624998681277160100799861 * fZ * fTmpA + -0.90256893466271141007219169782871 * fTmpC * fR2;
  pSH[356] = fTmpB * fC1;
  pSH[328] = fTmpB * fS1;
  fTmpC    = 2.9572714696920445983426700697905 * fZ * fTmpB + -0.92973952503676233431045838884188 * fTmpA * fR2;
  pSH[394] = fTmpC * fC1;
  pSH[366] = fTmpC * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.84425065085737263610244154876661;
  pSH[255] = fTmpA * fC0;
  pSH[225] = fTmpA * fS0;
  fTmpB    = -4.849850753230681765781895364853 * fZ;
  pSH[287] = fTmpB * fC0;
  pSH[257] = fTmpB * fS0;
  fTmpC    = -20.602948605549728459257474710853 * fZ2 + 0.62433177592574934725022650638948 * fR2;
  pSH[321] = fTmpC * fC0;
  pSH[291] = fTmpC * fS0;
  fTmpA    = fZ * (-74.51550792153200047235328540296 * fZ2 + 6.3870435361313143262512737052816 * fR2);
  pSH[357] = fTmpA * fC0;
  pSH[327] = fTmpA * fS0;
  fTmpB    = 3.2573446421352252764733203882486 * fZ * fTmpA + -0.90063003157873466783230401166982 * fTmpC * fR2;
  pSH[395] = fTmpB * fC0;
  pSH[365] = fTmpB * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.85734058883802509643716829867977;
  pSH[288] = fTmpA * fC1;
  pSH[256] = fTmpA * fS1;
  fTmpB    = 5.0720953248553606105934743464303 * fZ;
  pSH[322] = fTmpB * fC1;
  pSH[290] = fTmpB * fS1;
  fTmpC    = 22.134404124649146408318478584931 * fZ2 + -0.63241154641854704026400144090125 * fR2;
  pSH[358] = fTmpC * fC1;
  pSH[326] = fTmpC * fS1;
  fTmpA    = fZ * (82.055245832745453651857481247589 * fZ2 + -6.6531280404928746204443190670474 * fR2);
  pSH[396] = fTmpA * fC1;
  pSH[364] = fTmpA * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpA    = -0.86985717192062808221161840371849;
  pSH[323] = fTmpA * fC0;
  pSH[289] = fTmpA * fS0;
  fTmpB    = -5.2911346120699731671033205770982 * fZ;
  pSH[359] = fTmpB * fC0;
  pSH[325] = fTmpB * fS0;
  fTmpC    = -23.68730913497825269972696382581 * fZ2 + 0.64019754418860142430522386369773 * fR2;
  pSH[397] = fTmpC * fC0;
  pSH[363] = fTmpC * fS0;
  fC1      = fX * fC0 - fY * fS0;
  fS1      = fX * fS0 + fY * fC0;

  fTmpA    = 0.88185576867832886122002683526588;
  pSH[360] = fTmpA * fC1;
  pSH[324] = fTmpA * fS1;
  fTmpB    = 5.5071875102722446008139678408355 * fZ;
  pSH[398] = fTmpB * fC1;
  pSH[362] = fTmpB * fS1;
  fC0      = fX * fC1 - fY * fS1;
  fS0      = fX * fS1 + fY * fC1;

  fTmpC    = -0.89338378434994943305767073349344;
  pSH[399] = fTmpC * fC0;
  pSH[361] = fTmpC * fS0;
}

// auxiliary function
template<typename T>
void SHEval(const T fX0, const T fY0, const T fZ, T* pSH, const int l_max, bool addsign = false)
{
  T fX, fY;
  if (addsign)
  {
    fX = -1.0L * fX0;
    fY = -1.0L * fY0;
  }
  else
  {
    fX = fX0;
    fY = fY0;
  }
  switch (l_max)
  {
  case 0:
    SHEval0(fX, fY, fZ, pSH);
    break;

  case 1:
    SHEval1(fX, fY, fZ, pSH);
    break;

  case 2:
    SHEval2(fX, fY, fZ, pSH);
    break;

  case 3:
    SHEval3(fX, fY, fZ, pSH);
    break;

  case 4:
    SHEval4(fX, fY, fZ, pSH);
    break;

  case 5:
    SHEval5(fX, fY, fZ, pSH);
    break;

  case 6:
    SHEval6(fX, fY, fZ, pSH);
    break;

  case 7:
    SHEval7(fX, fY, fZ, pSH);
    break;

  case 8:
    SHEval8(fX, fY, fZ, pSH);
    break;

  case 9:
    SHEval9(fX, fY, fZ, pSH);
    break;

  case 10:
    SHEval10(fX, fY, fZ, pSH);
    break;

  case 11:
    SHEval11(fX, fY, fZ, pSH);
    break;

  case 12:
    SHEval12(fX, fY, fZ, pSH);
    break;

  case 13:
    SHEval13(fX, fY, fZ, pSH);
    break;

  case 14:
    SHEval14(fX, fY, fZ, pSH);
    break;

  case 15:
    SHEval15(fX, fY, fZ, pSH);
    break;

  case 16:
    SHEval16(fX, fY, fZ, pSH);
    break;

  case 17:
    SHEval17(fX, fY, fZ, pSH);
    break;

  case 18:
    SHEval18(fX, fY, fZ, pSH);
    break;

  case 19:
    SHEval19(fX, fY, fZ, pSH);
    break;
  }
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGH(T x, T y, T z)
{
  throw std::runtime_error("SoaSphericalTensor<T>::evaluateVGH(x,y,z):  Not implemented\n");
}

template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGHGH(T x, T y, T z)
{
  throw std::runtime_error("SoaSphericalTensor<T>::evaluateVGHGH(x,y,z):  Not implemented\n");
}

extern template struct SoaSphericalTensor<float>;
extern template struct SoaSphericalTensor<double>;
} // namespace qmcplusplus
#endif

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
#include "OhmmsPETE/OhmmsArray.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

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
class SoaSphericalTensor
{
private:
  using OffloadVector  = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadArray2D = Array<T, 2, OffloadPinnedAllocator<T>>;
  using OffloadArray3D = Array<T, 3, OffloadPinnedAllocator<T>>;
  using OffloadArray4D = Array<T, 4, OffloadPinnedAllocator<T>>;
  ///maximum angular momentum for the center
  int Lmax;
  /// Normalization factors
  const std::shared_ptr<OffloadVector> norm_factor_ptr_;
  ///pre-evaluated factor \f$1/\sqrt{(l+m)\times(l+1-m)}\f$
  const std::shared_ptr<OffloadVector> factorLM_ptr_;
  ///pre-evaluated factor \f$\sqrt{(2l+1)/(4\pi)}\f$
  const std::shared_ptr<OffloadVector> factorL_ptr_;
  ///pre-evaluated factor \f$(2l+1)/(2l-1)\f$
  const std::shared_ptr<OffloadVector> factor2L_ptr_;
  /// norm_factor reference
  OffloadVector& norm_factor_;
  /// factorLM reference
  OffloadVector& factorLM_;
  /// factorL reference
  OffloadVector& factorL_;
  /// factor2L reference
  OffloadVector& factor2L_;
  ///composite
  VectorSoaContainer<T, 5> cYlm;

public:
  explicit SoaSphericalTensor(const int l_max, bool addsign = false);

  SoaSphericalTensor(const SoaSphericalTensor& rhs) = default;

  ///compute Ylm for single position
  static void evaluate_bare(T x, T y, T z, T* Ylm, int lmax, const T* factorL, const T* factorLM);
  ///compute Ylm_vgl for single position
  static void evaluateVGL_impl(const T x,
                               const T y,
                               const T z,
                               T* restrict Ylm_vgl,
                               int lmax,
                               const T* factorL,
                               const T* factorLM,
                               const T* factor2L,
                               const T* normfactor,
                               size_t offset);

  ///compute Ylm
  inline void evaluateV(T x, T y, T z, T* Ylm) const
  {
    evaluate_bare(x, y, z, Ylm, Lmax, factorL_.data(), factorLM_.data());
    for (int i = 0, nl = cYlm.size(); i < nl; i++)
      Ylm[i] *= norm_factor_[i];
  }

  /**
   * @brief evaluate V for multiple electrons and multiple pbc images
   * 
   * @param [in] xyz electron positions [Nelec, Npbc, 3(x,y,z)]
   * @param [out] Ylm Spherical tensor elements [Nelec, Npbc, Nlm]
  */
  inline void batched_evaluateV(const OffloadArray3D& xyz, OffloadArray3D& Ylm) const
  {
    const size_t nElec = xyz.size(0);
    const size_t Npbc  = xyz.size(1); // number of PBC images
    assert(xyz.size(2) == 3);

    assert(Ylm.size(0) == nElec);
    assert(Ylm.size(1) == Npbc);
    const size_t Nlm = Ylm.size(2);

    size_t nR = nElec * Npbc; // total number of positions to evaluate

    auto* xyz_ptr          = xyz.device_data();
    auto* Ylm_ptr          = Ylm.device_data();
    auto* factorLM__ptr    = factorLM_.device_data();
    auto* factorL__ptr     = factorL_.device_data();
    auto* norm_factor__ptr = norm_factor_.device_data();

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    is_device_ptr(factorLM__ptr, factorL__ptr, norm_factor__ptr, xyz_ptr, Ylm_ptr)")
    for (uint32_t ir = 0; ir < nR; ir++)
    {
      evaluate_bare(xyz_ptr[0 + 3 * ir], xyz_ptr[1 + 3 * ir], xyz_ptr[2 + 3 * ir], Ylm_ptr + (ir * Nlm), Lmax,
                    factorL__ptr, factorLM__ptr);
      for (int i = 0; i < Nlm; i++)
        Ylm_ptr[ir * Nlm + i] *= norm_factor__ptr[i];
    }
  }

  /**
   * @brief evaluate VGL for multiple electrons and multiple pbc images
   * 
   * when offload is enabled, xyz is assumed to be up to date on the device before entering the function
   * Ylm_vgl will be up to date on the device (but not host) when this function exits
   * 
   * @param [in] xyz electron positions [Nelec, Npbc, 3(x,y,z)]
   * @param [out] Ylm_vgl Spherical tensor elements [5(v, gx, gy, gz, lapl), Nelec, Npbc, Nlm]
  */
  inline void batched_evaluateVGL(const OffloadArray3D& xyz, OffloadArray4D& Ylm_vgl) const
  {
    const size_t nElec = xyz.size(0);
    const size_t Npbc  = xyz.size(1); // number of PBC images
    assert(xyz.size(2) == 3);

    assert(Ylm_vgl.size(0) == 5);
    assert(Ylm_vgl.size(1) == nElec);
    assert(Ylm_vgl.size(2) == Npbc);
    const size_t Nlm = Ylm_vgl.size(3);
    assert(norm_factor_.size() == Nlm);

    const size_t nR     = nElec * Npbc; // total number of positions to evaluate
    const size_t offset = Nlm * nR;     // stride for v/gx/gy/gz/l

    auto* xyz_ptr          = xyz.device_data();
    auto* Ylm_vgl_ptr      = Ylm_vgl.device_data();
    auto* factorLM__ptr    = factorLM_.device_data();
    auto* factorL__ptr     = factorL_.device_data();
    auto* factor2L__ptr    = factor2L_.device_data();
    auto* norm_factor__ptr = norm_factor_.device_data();

    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    is_device_ptr(factorLM__ptr, factorL__ptr, norm_factor__ptr, factor2L__ptr, xyz_ptr, Ylm_vgl_ptr)")
    for (uint32_t ir = 0; ir < nR; ir++)
      evaluateVGL_impl(xyz_ptr[0 + 3 * ir], xyz_ptr[1 + 3 * ir], xyz_ptr[2 + 3 * ir], Ylm_vgl_ptr + (ir * Nlm), Lmax,
                       factorL__ptr, factorLM__ptr, factor2L__ptr, norm_factor__ptr, offset);
  }

  ///compute Ylm
  inline void evaluateV(T x, T y, T z)
  {
    T* restrict Ylm = cYlm.data(0);
    evaluate_bare(x, y, z, Ylm, Lmax, factorL_.data(), factorLM_.data());
    for (int i = 0, nl = cYlm.size(); i < nl; i++)
      Ylm[i] *= norm_factor_[i];
  }

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGL(T x, T y, T z);

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGH(T x, T y, T z);

  ///makes a table of \f$ r^l S_l^m \f$ and their gradients up to Lmax.
  void evaluateVGHGH(T x, T y, T z);

  ///returns the index/locator for (\f$l,m\f$) combo, \f$ l(l+1)+m \f$
  static inline int index(int l, int m) { return (l * (l + 1)) + m; }

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
inline SoaSphericalTensor<T>::SoaSphericalTensor(const int l_max, bool addsign)
    : Lmax(l_max),
      norm_factor_ptr_(std::make_shared<OffloadVector>()),
      factorLM_ptr_(std::make_shared<OffloadVector>()),
      factorL_ptr_(std::make_shared<OffloadVector>()),
      factor2L_ptr_(std::make_shared<OffloadVector>()),
      norm_factor_(*norm_factor_ptr_),
      factorLM_(*factorLM_ptr_),
      factorL_(*factorL_ptr_),
      factor2L_(*factor2L_ptr_)
{
  constexpr T czero(0);
  constexpr T cone(1);
  const int ntot = (Lmax + 1) * (Lmax + 1);
  cYlm.resize(ntot);
  norm_factor_.resize(ntot, cone);
  const T sqrt2 = std::sqrt(2.0);
  if (addsign)
  {
    for (int l = 0; l <= Lmax; l++)
    {
      norm_factor_[index(l, 0)] = cone;
      for (int m = 1; m <= l; m++)
      {
        norm_factor_[index(l, m)]  = std::pow(-cone, m) * sqrt2;
        norm_factor_[index(l, -m)] = std::pow(-cone, -m) * sqrt2;
      }
    }
  }
  else
  {
    for (int l = 0; l <= Lmax; l++)
    {
      for (int m = 1; m <= l; m++)
      {
        norm_factor_[index(l, m)]  = sqrt2;
        norm_factor_[index(l, -m)] = sqrt2;
      }
    }
  }
  factorL_.resize(Lmax + 1);
  const T omega = 1.0 / std::sqrt(16.0 * std::atan(1.0));
  for (int l = 1; l <= Lmax; l++)
    factorL_[l] = std::sqrt(static_cast<T>(2 * l + 1)) * omega;
  factor2L_.resize(Lmax + 1);
  for (int l = 1; l <= Lmax; l++)
    factor2L_[l] = static_cast<T>(2 * l + 1) / static_cast<T>(2 * l - 1);
  factorLM_.resize(ntot);
  for (int l = 1; l <= Lmax; l++)
    for (int m = 1; m <= l; m++)
    {
      T fac2                  = 1.0 / std::sqrt(static_cast<T>((l + m) * (l + 1 - m)));
      factorLM_[index(l, m)]  = fac2;
      factorLM_[index(l, -m)] = fac2;
    }
  norm_factor_.updateTo();
  factorLM_.updateTo();
  factorL_.updateTo();
  factor2L_.updateTo();
}

PRAGMA_OFFLOAD("omp declare target")
template<typename T>
inline void SoaSphericalTensor<T>::evaluate_bare(
    T x, T y, T z,
    T* restrict Ylm,
    int lmax,
    const T* restrict factorL,
    const T* restrict factorLM)
{
  constexpr T czero(0);
  constexpr T cone(1);
  // Instead of std::atan:
  constexpr T myPi = T(3.14159265358979323846);
  const T omega    = cone / std::sqrt(T(4) * myPi);

  // Epsilon^2 for near-axis detection
  constexpr T eps2 = std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon();

  // 1) Compute r, ctheta, cphi, sphi
  const T r2xy = x*x + y*y;
  const T r    = std::sqrt(r2xy + z*z);

  T ctheta, cphi, sphi;
  if (r2xy < eps2)
  {
    // near z-axis
    cphi   = czero;
    sphi   = cone;
    ctheta = (z < czero) ? -cone : cone;
  }
  else
  {
    const T rInv  = cone / r; // safe since r > 0 if r2xy>=eps2 or z!=0
    ctheta = (z * rInv);
    // clamp ctheta to [-1,1]
    ctheta = (ctheta > cone)  ? cone  : ctheta;
    ctheta = (ctheta < -cone) ? -cone : ctheta;

    // reciprocal for x,y plane
    const T rxyi = cone / std::sqrt(r2xy);
    cphi = x * rxyi;
    sphi = y * rxyi;
  }

  // stheta = sqrt(1 - ctheta^2)
  const T stheta = std::sqrt(cone - ctheta*ctheta);

  // 2) Compute Associated Legendre polynomials P_l^m
  //    - Ylm[0] = 1
  Ylm[0] = cone;

  // (a) P_ll and P_l,l-1
  T fac = cone;
  int j  = -1;
  for (int l = 1; l <= lmax; l++)
  {
    j   += 2;         // j: 1, 3, 5, ...
    fac *= -T(j) * stheta;
    const int ll = index(l, l);
    const int l1 = index(l, l-1);
    const int l2 = index(l-1, l-1);

    Ylm[ll] = fac;
    Ylm[l1] = T(j) * ctheta * Ylm[l2];
  }

  // (b) Recur for other P_l^m
  // Jackson's convention: For m in [0..lmax-1], l in [m+2..lmax]
  for (int m = 0; m < lmax - 1; m++)
  {
    int j = 2*m + 1; // starts at 1 for m=0, 3 for m=1, etc.
    for (int l = m + 2; l <= lmax; l++)
    {
      j += 2;
      const int lm  = index(l,   m);
      const int l1  = index(l-1, m);
      const int l2  = index(l-2, m);

      Ylm[lm] = (ctheta * T(j) * Ylm[l1] - (l + m - 1)*Ylm[l2]) / T(l - m);
    }
  }

  // 3) Multiply by r^l and incorporate azimuthal part
  //    Y_0^0 => 1/sqrt(4*pi)
  Ylm[0] = omega;

  T rpow = cone;
  for (int l = 1; l <= lmax; l++)
  {
    rpow *= r;
    // factorL[l] is sqrt( (2*l+1)/(4*pi) ), for instance
    T scale_l = rpow * factorL[l];

    // multiply the m=0 term
    const int l0 = index(l,0);
    Ylm[l0] *= scale_l;

    // Build up cphi^m, sphi^m iteratively
    T cphim = cone; // cphi^0
    T sphim = czero;

    for (int m = 1; m <= l; m++)
    {
      // rotate (cphim, sphim) by (cphi, sphi)
      T tmpC  = cphim*cphi - sphim*sphi;
      T tmpS  = sphim*cphi + cphim*sphi;
      cphim   = tmpC;
      sphim   = tmpS;

      const int lmPos = index(l, m);
      // factorLM[lmPos] might be e.g. sqrt( (l-m+1)(l+m) ) or something.
      scale_l *= factorLM[lmPos];

      // multiply P_l^m by scale_l
      T tmpVal = Ylm[lmPos] * scale_l;

      // store real part => cphim, imaginary => sphi => Ylm(l, -m)
      Ylm[lmPos]        = tmpVal * cphim;
      Ylm[ index(l, -m) ] = tmpVal * sphim;
    }
  }
}
PRAGMA_OFFLOAD("omp end declare target")



PRAGMA_OFFLOAD("omp declare target")
template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGL_impl(const T x,
                                                    const T y,
                                                    const T z,
                                                    T* restrict Ylm_vgl,
                                                    int lmax,
                                                    const T* factorL,
                                                    const T* factorLM,
                                                    const T* factor2L,
                                                    const T* normfactor,
                                                    size_t offset)
{
  // 1) Evaluate the base spherical harmonics [0..Nlm-1] in Ylm_vgl
  evaluate_bare(x, y, z, Ylm_vgl, lmax, factorL, factorLM);

  // 2) Basic parameters
  const int Nlm = (lmax + 1) * (lmax + 1);

  constexpr T czero(0);
  constexpr T ahalf(0.5);

  // 3) Identify gradient / laplacian blocks
  T* restrict gYlmX = Ylm_vgl + offset * 1;
  T* restrict gYlmY = Ylm_vgl + offset * 2;
  T* restrict gYlmZ = Ylm_vgl + offset * 3;
  T* restrict lYlm  = Ylm_vgl + offset * 4;

  // For the l=0 component, just set gradient/laplacian to 0
  gYlmX[0] = czero;
  gYlmY[0] = czero;
  gYlmZ[0] = czero;
  lYlm[0]  = czero;

  // 4) Compute gradients for l=1..lmax
  //    we store `lmBasePrev = index(l-1,0)` once per l instead of re-calling it inside the loop.
  for (int l = 1; l <= lmax; l++)
  {
    // factor2L[l] often = ( (2*l + 1) / (2*l - 1) ) or similar
    T fac = factor2L[l];

    // The base offset for shell (l-1) with m=0
    const int lmBasePrev = index(l - 1, 0);

    for (int m = -l; m <= l; m++)
    {
      T gx, gy, gz, dpr, dpi, dmr, dmi;
      const int ma = std::abs(m);

      const T cp = std::sqrt(fac * (l - ma - 1) * (l - ma));
      const T cm = std::sqrt(fac * (l + ma - 1) * (l + ma));
      const T c0 = std::sqrt(fac * (l - ma) * (l + ma));

      // if l>ma => gz = c0 * Ylm[(l-1,m)] = c0 * Ylm[lmBasePrev + m]
      gz = (l > ma) ? (c0 * Ylm_vgl[lmBasePrev + m]) : czero;

      // dpr, dpi if (l>ma+1)
      if (l > ma + 1)
      {
        dpr = cp * Ylm_vgl[lmBasePrev + (ma + 1)];
        dpi = cp * Ylm_vgl[lmBasePrev - (ma + 1)];
      }
      else
      {
        dpr = czero;
        dpi = czero;
      }

      // dmr, dmi if (l>1)
      if (l > 1)
      {
        switch (ma)
        {
        case 0: // ma=0
          dmr = -cm * Ylm_vgl[lmBasePrev + 1];
          dmi =  cm * Ylm_vgl[lmBasePrev - 1];
          break;
        case 1: // ma=1
          dmr =  cm * Ylm_vgl[lmBasePrev];
          dmi =  czero;
          break;
        default: // ma >= 2
          dmr =  cm * Ylm_vgl[lmBasePrev + (ma - 1)];
          dmi =  cm * Ylm_vgl[lmBasePrev - (ma - 1)];
          break;
        }
      }
      else
      {
        dmr = cm * Ylm_vgl[lmBasePrev];
        dmi = czero;
      }

      // combine them into (gx, gy)
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

      // Store results in the shell l => index(l,m)
      const int lmCurr = index(l, m);

      if (ma > 0) // multiply gradient by normfactor if ma != 0
      {
        const T nf = normfactor[lmCurr];
        gYlmX[lmCurr] = nf * gx;
        gYlmY[lmCurr] = nf * gy;
        gYlmZ[lmCurr] = nf * gz;
      }
      else
      {
        // if ma=0 => skip the normfactor for gradient
        gYlmX[lmCurr] = gx;
        gYlmY[lmCurr] = gy;
        gYlmZ[lmCurr] = gz;
      }
    } // end for m
  } // end for l

  // 5) Multiply the base spherical harmonics Ylm by normfactor for all i,
  //    and set lYlm[i] = 0
  for (int i = 0; i < Nlm; i++)
  {
    Ylm_vgl[i] *= normfactor[i];
    lYlm[i] = czero;
  }
}
PRAGMA_OFFLOAD("omp end declare target")


template<typename T>
inline void SoaSphericalTensor<T>::evaluateVGL(T x, T y, T z)
{
  evaluateVGL_impl(x, y, z, cYlm.data(), Lmax, factorL_.data(), factorLM_.data(), factor2L_.data(), norm_factor_.data(),
                   cYlm.capacity());
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

extern template class SoaSphericalTensor<float>;
extern template class SoaSphericalTensor<double>;
} // namespace qmcplusplus
#endif

//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim           //
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////
/** helper functions for EinsplineSetBuilder
 */
#ifndef QMCPLUSPLUS_EINSPLINEBUILDER_HELPER_H
#define QMCPLUSPLUS_EINSPLINEBUILDER_HELPER_H
#include <complex>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "CPU/SIMD/simd.hpp"
#include "CPU/BLAS.hpp"
#include "CPU/math.hpp"

namespace qmcplusplus
{
/** unpack packed cG to fftbox
   * @param cG packed vector
   * @param gvecs g-coordinate for cG[i]
   * @param maxg  fft grid
   * @param fftbox unpacked data to be transformed
   */
template<typename T>
inline void unpack4fftw(const Vector<std::complex<T>>& cG,
                        const std::vector<TinyVector<int, 3>>& gvecs,
                        const TinyVector<int, 3>& maxg,
                        Array<std::complex<T>, 3>& fftbox)
{
  fftbox                   = std::complex<T>();
  const int upper_bound[3] = {(maxg[0] - 1) / 2, (maxg[1] - 1) / 2, (maxg[2] - 1) / 2};
  const int lower_bound[3] = {upper_bound[0] - maxg[0] + 1, upper_bound[1] - maxg[1] + 1, upper_bound[2] - maxg[2] + 1};
  //only coefficient indices between [lower_bound,upper_bound] are taken for FFT.
  //this is rather unsafe
  //#pragma omp parallel for
  for (int iG = 0; iG < cG.size(); iG++)
  {
    if (gvecs[iG][0] > upper_bound[0] || gvecs[iG][0] < lower_bound[0] || gvecs[iG][1] > upper_bound[1] ||
        gvecs[iG][1] < lower_bound[1] || gvecs[iG][2] > upper_bound[2] || gvecs[iG][2] < lower_bound[2])
    {
      //std::cout << "Warning: cG out of bound "
      //          << "x " << gvecs[iG][0]    << " y " << gvecs[iG][1]    << " z " << gvecs[iG][2] << std::endl
      //          << "xu " << upper_bound[0] << " yu " << upper_bound[1] << " zu " << upper_bound[2] << std::endl
      //          << "xd " << lower_bound[0] << " yd " << lower_bound[1] << " zd " << lower_bound[2] << std::endl;
      continue;
    }
    fftbox((gvecs[iG][0] + maxg[0]) % maxg[0], (gvecs[iG][1] + maxg[1]) % maxg[1], (gvecs[iG][2] + maxg[2]) % maxg[2]) =
        cG[iG];
  }
}

/** rotate the state after 3dfft
   *
   * First, add the eikr phase factor.
   * Then, rotate the phase of the orbitals so that neither
   * the real part nor the imaginary part are very near
   * zero.  This sometimes happens in crystals with high
   * symmetry at special k-points.
   */
template<typename T, typename T1, typename T2>
inline void fix_phase_rotate_c2r(Array<std::complex<T>, 3>& in,
                                 Array<T1, 3>& out,
                                 const TinyVector<T2, 3>& twist,
                                 T& phase_r,
                                 T& phase_i)
{
  const T two_pi = -2.0 * M_PI;
  const int nx   = in.size(0);
  const int ny   = in.size(1);
  const int nz   = in.size(2);
  T nx_i         = 1.0 / static_cast<T>(nx);
  T ny_i         = 1.0 / static_cast<T>(ny);
  T nz_i         = 1.0 / static_cast<T>(nz);

  T rNorm = 0.0, iNorm = 0.0, riNorm = 0.0;
#pragma omp parallel for reduction(+ : rNorm, iNorm, riNorm)
  for (int ix = 0; ix < nx; ix++)
  {
    T s, c;
    std::complex<T>* restrict in_ptr = in.data() + ix * ny * nz;
    T rux                            = static_cast<T>(ix) * nx_i * twist[0];
    for (int iy = 0; iy < ny; iy++)
    {
      T ruy = static_cast<T>(iy) * ny_i * twist[1];
      for (int iz = 0; iz < nz; iz++)
      {
        T ruz = static_cast<T>(iz) * nz_i * twist[2];
        qmcplusplus::sincos(-two_pi * (rux + ruy + ruz), &s, &c);
        std::complex<T> eikr(c, s);
        *in_ptr *= eikr;
        rNorm += in_ptr->real() * in_ptr->real();
        iNorm += in_ptr->imag() * in_ptr->imag();
        riNorm += in_ptr->real() * in_ptr->imag();
        ++in_ptr;
      }
    }
  }

  const T x   = (rNorm - iNorm) / riNorm;
  const T y   = 1.0 / std::sqrt(x * x + 4.0);
  const T phs = std::sqrt(0.5 - y);
  phase_r     = phs;
  phase_i     = (x < 0) ? std::sqrt(1.0 - phs * phs) : -std::sqrt(1.0 - phs * phs);

#pragma omp parallel for
  for (int ix = 0; ix < nx; ix++)
  {
    const std::complex<T>* restrict in_ptr = in.data() + ix * ny * nz;
    T1* restrict out_ptr                   = out.data() + ix * ny * nz;
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++)
      {
        *out_ptr = static_cast<T1>(phase_r * in_ptr->real() - phase_i * in_ptr->imag());
        ++in_ptr;
        ++out_ptr;
      }
  }
}

template<typename T, typename T1, typename T2>
inline void fix_phase_rotate_c2c(const Array<std::complex<T>, 3>& in,
                                 Array<std::complex<T1>, 3>& out,
                                 const TinyVector<T2, 3>& twist)
{
  const int nx = in.size(0);
  const int ny = in.size(1);
  const int nz = in.size(2);
  T phase_r, phase_i;

  compute_phase(in, twist, phase_r, phase_i);

#pragma omp parallel for
  for (int ix = 0; ix < nx; ++ix)
  {
    const std::complex<T>* restrict in_ptr = in.data() + ix * ny * nz;
    std::complex<T1>* restrict out_ptr     = out.data() + ix * ny * nz;
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz)
      {
        *out_ptr = std::complex<T1>(static_cast<T1>(phase_r * in_ptr->real() - phase_i * in_ptr->imag()),
                                    static_cast<T1>(phase_i * in_ptr->real() + phase_r * in_ptr->imag()));
        ++out_ptr;
        ++in_ptr;
      }
  }
}

template<typename T, typename T1, typename T2>
inline void fix_phase_rotate_c2c(const Array<std::complex<T>, 3>& in,
                                 Array<T1, 3>& out_r,
                                 Array<T1, 3>& out_i,
                                 const TinyVector<T2, 3>& twist,
                                 T& phase_r,
                                 T& phase_i)
{
  const int nx = in.size(0);
  const int ny = in.size(1);
  const int nz = in.size(2);

  compute_phase(in, twist, phase_r, phase_i);

#pragma omp parallel for
  for (size_t ix = 0; ix < nx; ++ix)
    for (size_t iy = 0; iy < ny; ++iy)
    {
      const size_t offset                    = ix * ny * nz + iy * nz;
      const std::complex<T>* restrict in_ptr = in.data() + offset;
      T1* restrict r_ptr                     = out_r.data() + offset;
      T1* restrict i_ptr                     = out_i.data() + offset;
      for (size_t iz = 0; iz < nz; ++iz)
      {
        r_ptr[iz] = static_cast<T1>(phase_r * in_ptr[iz].real() - phase_i * in_ptr[iz].imag());
        i_ptr[iz] = static_cast<T1>(phase_i * in_ptr[iz].real() + phase_r * in_ptr[iz].imag());
      }
    }
}

/** Split FFTs into real/imaginary components. 
  * @param in ffts
  * @param out_r real component
  * @param out_i imaginary components
  */
template<typename T, typename T1>
inline void split_real_components_c2c(const Array<std::complex<T>, 3>& in, Array<T1, 3>& out_r, Array<T1, 3>& out_i)
{
  const int nx = in.size(0);
  const int ny = in.size(1);
  const int nz = in.size(2);

#pragma omp parallel for
  for (size_t ix = 0; ix < nx; ++ix)
    for (size_t iy = 0; iy < ny; ++iy)
    {
      const size_t offset                    = ix * ny * nz + iy * nz;
      const std::complex<T>* restrict in_ptr = in.data() + offset;
      T1* restrict r_ptr                     = out_r.data() + offset;
      T1* restrict i_ptr                     = out_i.data() + offset;
      for (size_t iz = 0; iz < nz; ++iz)
      {
        r_ptr[iz] = static_cast<T1>(in_ptr[iz].real());
        i_ptr[iz] = static_cast<T1>(in_ptr[iz].imag());
      }
    }
}

/** Compute the norm of an orbital.
   * @param cG the plane wave coefficients
   * @return norm of the orbital
   */
template<typename T>
inline T compute_norm(const Vector<std::complex<T>>& cG)
{
  T total_norm2(0);
#pragma omp parallel for reduction(+ : total_norm2)
  for (size_t ig = 0; ig < cG.size(); ++ig)
    total_norm2 += cG[ig].real() * cG[ig].real() + cG[ig].imag() * cG[ig].imag();
  return std::sqrt(total_norm2);
}

/** Compute the phase factor for rotation. The algorithm aims at balanced real and imaginary parts.
   * @param in the real space orbital value on a 3D grid
   * @param twist k-point in reduced coordinates
   * @param phase_r output real part of the phase
   * @param phase_i output imaginary part of the phase
   */
template<typename T, typename T2>
inline void compute_phase(const Array<std::complex<T>, 3>& in, const TinyVector<T2, 3>& twist, T& phase_r, T& phase_i)
{
  const T two_pi  = -2.0 * M_PI;
  const size_t nx = in.size(0);
  const size_t ny = in.size(1);
  const size_t nz = in.size(2);

  const T nx_i = 1.0 / static_cast<T>(nx);
  const T ny_i = 1.0 / static_cast<T>(ny);
  const T nz_i = 1.0 / static_cast<T>(nz);

  T rNorm = 0.0, iNorm = 0.0, riNorm = 0.0;
#pragma omp parallel for reduction(+ : rNorm, iNorm, riNorm)
  for (size_t ix = 0; ix < nx; ++ix)
  {
    for (size_t iy = 0; iy < ny; ++iy)
    {
      const T rux = static_cast<T>(ix) * nx_i * twist[0];
      T s, c;
      T rsum = 0, isum = 0, risum = 0.0;
      const T ruy                            = static_cast<T>(iy) * ny_i * twist[1];
      const std::complex<T>* restrict in_ptr = in.data() + ix * ny * nz + iy * nz;
      for (size_t iz = 0; iz < nz; ++iz)
      {
        const T ruz = static_cast<T>(iz) * nz_i * twist[2];
        qmcplusplus::sincos(-two_pi * (rux + ruy + ruz), &s, &c);
        const T re = c * in_ptr[iz].real() - s * in_ptr[iz].imag();
        const T im = s * in_ptr[iz].real() + c * in_ptr[iz].imag();
        rsum += re * re;
        isum += im * im;
        risum += re * im;
      }
      rNorm += rsum;
      iNorm += isum;
      riNorm += risum;
    }
  }

  const T x   = (rNorm - iNorm) / riNorm;
  const T y   = 1.0 / std::sqrt(x * x + 4.0);
  const T phs = std::sqrt(0.5 - y);
  phase_r     = phs;
  phase_i     = (x < 0) ? std::sqrt(1.0 - phs * phs) : -std::sqrt(1.0 - phs * phs);
}

/** rotate the state after 3dfft
   *
   */
template<typename T>
inline void fix_phase_rotate(const Array<std::complex<T>, 3>& e2pi, Array<std::complex<T>, 3>& in, Array<T, 3>& out)
{
  const int nx = e2pi.size(0);
  const int ny = e2pi.size(1);
  const int nz = e2pi.size(2);
  T rNorm = 0.0, iNorm = 0.0;
  //#pragma omp parallel for reduction(+:rNorm,iNorm)
  for (int ix = 0; ix < nx; ix++)
  {
    T rpart = 0.0, ipart = 0.0;
    const std::complex<T>* restrict p_ptr = e2pi.data() + ix * ny * nz;
    std::complex<T>* restrict in_ptr      = in.data() + ix * ny * nz;
    for (int iyz = 0; iyz < ny * nz; ++iyz)
    {
      in_ptr[iyz] *= p_ptr[iyz];
      rpart += in_ptr[iyz].real() * in_ptr[iyz].real();
      ipart += in_ptr[iyz].imag() * in_ptr[iyz].imag();
    }
    rNorm += rpart;
    iNorm += ipart;
  }

  //#pragma omp parallel
  {
    T arg = std::atan2(iNorm, rNorm);
    T phase_i, phase_r;
    qmcplusplus::sincos(0.125 * M_PI - 0.5 * arg, &phase_i, &phase_r);
    //#pragma omp for
    for (int ix = 0; ix < nx; ix++)
    {
      const std::complex<T>* restrict in_ptr = in.data() + ix * ny * nz;
      T* restrict out_ptr                    = out.data() + ix * ny * nz;
      for (int iyz = 0; iyz < ny * nz; iyz++)
        out_ptr[iyz] = phase_r * in_ptr[iyz].real() - phase_i * in_ptr[iyz].imag();
    }
  }
}

} // namespace qmcplusplus

#endif

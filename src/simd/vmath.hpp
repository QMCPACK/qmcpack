//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file vmath.h
 *
 * Define vectorized math functions. 
 * HAVE_MKL_VML 
 * HAVE_MASSV 
 */
#ifndef QMCPLUSPLUS_VECTORIZED_STDMATH_HPP
#define QMCPLUSPLUS_VECTORIZED_STDMATH_HPP

#include <cmath>
#if defined(HAVE_MKL_VML)
#include <mkl_vml_functions.h>
#elif defined(HAVE_MASSV)
#include <mass.h>
#include <massv.h>
#endif

namespace qmcplusplus {

  namespace simd 
  {
    /**  mod on an array
     * out[i]=in[i]-floor(in[i])
     */
    template<typename T, typename SIZET>
      inline void remainder(const T* restrict in, T* restrict out, SIZET n)
      {
        for(SIZET i=0; i<n; ++i) out[i]=in[i]-std::floor(in[i]);
      }

    template<typename T, typename SIZET>
      inline void remainder(T* restrict inout, SIZET n)
      {
        for(SIZET i=0; i<n; ++i) inout[i]-=std::floor(inout[i]);
      }

    template<typename T, typename SIZET>
      inline void sqrt(T* restrict inout, SIZET n)
      {
        for(SIZET i=0; i<n; ++i) inout[i]=std::sqrt(inout[i]);
      }

    template<typename T>
      inline void sqrt(const T* restrict in, T* restrict out, int n)
      {
        for(int i=0; i<n; ++i) out[i]=std::sqrt(in[i]);
      }
    template<typename T>
      inline void inv(const T* restrict in, T* restrict out, int n)
      {
        for(int i=0; i<n; ++i) out[i]=1.0/in[i];
      }

#if defined(HAVE_MKL_VML)
    inline void sqrt(const double* restrict in, double* restrict out, int n)
    {
      vdSqrt(n,in,out);
    }

    inline void sqrt(const float* restrict in, float* restrict out, int n)
    {
      vsSqrt(n,in,out);
    }

    inline void inv(const double* restrict in, double* restrict out, int n)
    {
      vdInv(n,in,out);
    }

    inline void inv(const float* restrict in, float* restrict out, int n)
    {
      vsInv(n,in,out);
    }

    inline void inv_sqrt(const double* restrict in, double* restrict out, int n)
    {
      vdInvSqrt(n,in,out);
    }

    inline void inv_sqrt(const float* restrict in, float* restrict out, int n)
    {
      vsInvSqrt(n,in,out);
    }

#elif defined(HAVE_MASSV)
    inline void sqrt(double* restrict in, double* restrict out, int n)
    {
      vsqrt(out,in,&n);
    }
    inline void sqrt(float* restrict in, float* restrict out, int n)
    {
      vssqrt(out,in,&n);
    }
    inline void inv_sqrt(double* restrict in, double* restrict out, int n)
    {
      vrsqrt(out,in,&n);
    }

    inline void inv_sqrt(float* restrict in, float* restrict out, int n)
    {
      vsrsqrt(out,in,&n);
    }
#endif

    template<typename T>
      inline void add(int n, const T* restrict in, T* restrict out)
      {
        for(int i=0; i<n; ++i) out[i]+=in[i];
      }

    template<typename T>
      inline void get_phase(int n, const T* restrict kpts, const T* restrict xyz, T* restrict phi)
      {
        T x=xyz[0]; T y=xyz[1]; T z=xyz[2];
        for(int i=0; i<n; ++i)
          phi[i]=x*kpts[i*3]+y*kpts[i*3+1]+z*kpts[i*3+2];
      }

    template<typename AT, typename BT, typename CT>
      inline void get_phase(const AT& kpts, const BT& pos, CT& phase)
      {
        const char transa = 'T';
        const char transb = 'N';
        const double zone(1.0);
        const double zero(0.0);
        dgemm(transa, transb
            , phase.cols(), phase.rows(), 3
            , zone
            , &(kpts[0][0]), 3
            , &(pos[0][0]), 3
            , zero, phase.data(), phase.rows());
      }
  }
}
#endif

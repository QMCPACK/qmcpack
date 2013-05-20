//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_E2IPHI_H
#define QMCPLUSPLUS_E2IPHI_H

#include <config.h>
#include <vector>
#include <complex>
#include <config/stdlib/math.h>

// #if defined(HAVE_ACML)
// extern "C"
// {
// #include <acml_mv.h>
// }

#if defined(HAVE_AMDLIBM)
namespace std
{
#include <amdlibm.h>
}
using namespace std;

/* additional translations for AMD libm */

#undef vrda_sincos
#define vrda_sincos amd_vrda_sincos

#undef vrsa_sincosf
#define vrsa_sincosf amd_vrsa_sincosf

#undef vrsa_sincos
#define vrsa_sincos amd_vrsa_sincosf

inline void
eval_e2iphi(int n, double* restrict phi, double* restrict c, double *restrict s)
{
  vrda_sincos(n, phi, s, c);
}

inline void
eval_e2iphi(int n, float* restrict phi, float* restrict c, float *restrict s)
{
  vrsa_sincos(n, phi, s, c);
}

inline void
eval_e2iphi(int n, double* restrict phi, std::complex<double>* restrict z)
{
  double c[n], s[n];
  vrda_sincos(n, phi, s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (int n, float* restrict phi, std::complex<float>* restrict z)
{
  float c[n], s[n];
  vrsa_sincosf(n, phi, s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
}
#elif defined(HAVE_MASSV)
#include <massv.h>
inline void
eval_e2iphi(int n, double* restrict phi, double* restrict c, double *restrict s)
{
  vsincos(s,c,phi,&n);
}
inline void
eval_e2iphi(int n, float* restrict phi, float* restrict c, float *restrict s)
{
  vssincos(s,c,phi,&n);
}

inline void
eval_e2iphi(int n, double* restrict phi, std::complex<double>* restrict z)
{
  double s[n],c[n];
  vsincos(s,c,phi,&n);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (int n, float* restrict phi, std::complex<float>* restrict z)
{
  float s[n],c[n];
  vssincos(s,c,phi,&n);
  for (int i=0; i<n; i++)
    z[i] = std::complex<float>(c[i],s[i]);
}
#elif defined(HAVE_MKL_VML)
#include <mkl_vml_functions.h>
inline void
eval_e2iphi(int n, double* restrict phi, double* restrict c, double *restrict s)
{
  vdSinCos(n,phi,s,c);
}

inline void
eval_e2iphi(int n, float* restrict phi, float* restrict c, float *restrict s)
{
  vsSinCos(n,phi,s,c);
}
inline void
eval_e2iphi(int n, const double* restrict phi, std::complex<double>* restrict z)
{
  vzCIS (n,phi,(MKL_Complex16*)(z));
}

inline void
eval_e2iphi (int n, const float* restrict phi, std::complex<float>* restrict z)
{
  vcCIS (n,phi,(MKL_Complex8*)(z));
}
#else/* generic case */
template<typename T>
inline void
eval_e2iphi (int n, const T* restrict phi, T* restrict phase_r, T* restrict phase_i)
{
  for (int i=0; i<n; i++)
    sincos(*phi++,phase_i++,phase_r++);
}
template<typename T>
inline void
eval_e2iphi (int n, const T* restrict phi, std::complex<T>* restrict z)
{
  T s,c;
  for (int i=0; i<n; i++)
  {
    sincos(phi[i],&s,&c);
    z[i]=std::complex<T>(c,s);
  }
}
#endif


template<typename T>
inline void
eval_e2iphi(std::vector<T>& phi, std::vector<std::complex<T> >& z)
{
  eval_e2iphi(phi.size(),&phi[0],&z[0]);
}

//template<typename T>
//inline void
//eval_e2iphi(qmcplusplus::Vector<T>& phi, qmcplusplus::Vector<std::complex<T> >& z)
//{
//  eval_e2iphi(phi.size(),phi.data(),z.data());
//}

#endif

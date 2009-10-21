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

#include <vector>
#include <complex>
#include <config.h>
#include <config/stdlib/math.h>

#if defined(HAVE_ACML)
#include <acml.h>
inline void
eval_e2iphi (const std::vector<double> &phi, std::vector<std::complex<double> > &z)
{
  int n = phi.size();
  double c[n], s[n];
  vrda_sincos(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (const std::vector<float> &phi, std::vector<std::complex<float> > &z)
{
  int n = phi.size();
  float c[n], s[n];
  vrsa_sincosf(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<float>(c[i],s[i]);
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<double> &phi, APPNAMESPACE::Vector<std::complex<double> > &z)
{
  int n = phi.size();
  double c[n], s[n];
  vrda_sincos(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (const std::vector<float> &phi, std::vector<std::complex<float> > &z)
{
  int n = phi.size();
  float c[n], s[n];
  vrsa_sincosf(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<float>(c[i],s[i]);
}
#elif defined(HAVE_MASSV)
#include <massv.h>
inline void
eval_e2iphi (const std::vector<double> &phi, std::vector<std::complex<double> > &z)
{
  int n = phi.size();
  double s[n],c[n];
  vsincos(s,c,&(phi[0]),&n);
  for (int i=0; i<n; i++) z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (const std::vector<float> &phi, std::vector<std::complex<float> > &z)
{
  int n = phi.size();
  float s[n],c[n];
  vssincos(s,c,&(phi[0]),&n);
  for (int i=0; i<n; i++) z[i] = std::complex<float>(c[i],s[i]);
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<double> &phi, APPNAMESPACE::Vector<std::complex<double> > &z)
{
  int n = phi.size();
  double s[n],c[n];
  vsincos(s,c,&(phi[0]),&n);
  for (int i=0; i<n; i++) z[i] = std::complex<double>(c[i],s[i]);
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<float> &phi, APPNAMESPACE::Vector<std::complex<float> > &z)
{
  int n = phi.size();
  float s[n],c[n];
  vssincos(s,c,&(phi[0]),&n);
  for (int i=0; i<n; i++) z[i] = std::complex<float>(c[i],s[i]);
}

#elif defined(HAVE_MKL_VML)
#include <mkl_vml_functions.h>
inline void
eval_e2iphi (const std::vector<double> &phi, std::vector<std::complex<double> > &z)
{
  vzCIS (phi.size(), &(phi[0]), (MKL_Complex16*) &(z[0]));
}

inline void
eval_e2iphi (const std::vector<float> &phi, std::vector<std::complex<float> > &z)
{
  vcCIS (phi.size(), &(phi[0]), (MKL_Complex8*) &(z[0]));
}
inline void
eval_e2iphi (const APPNAMESPACE::Vector<double> &phi, APPNAMESPACE::Vector<std::complex<double> > &z)
{
  vzCIS (phi.size(), phi.data(),(MKL_Complex16*)z.data());
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<float> &phi, APPNAMESPACE::Vector<std::complex<float> > &z)
{
  vcCIS (phi.size(), phi.data(),(MKL_Complex8*)z.data());
}

#else
inline void
eval_e2iphi (const std::vector<double> &phi, std::vector<std::complex<double> > &z)
{
  int n = phi.size();
  for (int i=0; i<n; i++) sincos (phi[i], &(z[i].imag()), &(z[i].real()));
}

inline void
eval_e2iphi (const std::vector<float> &phi, std::vector<std::complex<float> > &z)
{
  int n = phi.size();
#if defined(__APPLE__) && !defined(__INTEL_COMPILER)
  for (int i=0; i<n; i++) sincos (phi[i], &(z[i].imag()), &(z[i].real()));
#else
  for (int i=0; i<n; i++) sincosf (phi[i], &(z[i].imag()), &(z[i].real()));
#endif
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<double> &phi, APPNAMESPACE::Vector<std::complex<double> > &z)
{
  int n = phi.size();
  for (int i=0; i<n; i++) sincos (phi[i], &(z[i].imag()), &(z[i].real()));
}

inline void
eval_e2iphi (const APPNAMESPACE::Vector<float> &phi, APPNAMESPACE::Vector<std::complex<float> > &z)
{
  int n = phi.size();
#if defined(__APPLE__) && !defined(__INTEL_COMPILER)
  for (int i=0; i<n; i++) sincos (phi[i], &(z[i].imag()), &(z[i].real()));
#else
  for (int i=0; i<n; i++) sincosf (phi[i], &(z[i].imag()), &(z[i].real()));
#endif
}
#endif/* generic case */
#endif

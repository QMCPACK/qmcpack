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

#ifdef HAVE_ACML
#include <acml.h>
#endif

#ifdef HAVE_MKL
#include <mkl.h>
#endif


inline void
eval_e2iphi (std::vector<double> phi, std::vector<std::complex<double> > z)
{
#ifdef HAVE_ACML
  int n = phi.size();
  double c[n], s[n];
  vrda_sincos(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<double>(c[i],s[i]);
#elifdef HAVE_MKL
  vzCIS (n, &(phi[0]), (double*) &(z[0]));
#else
  for (int i=0; i<n; i++)
    sincos (phi[i], &(z[i].imag()), &(z[i].real()));
#endif
}


inline void
eval_e2iphi (std::vector<float> phi, std::vector<std::complex<float> > z)
{
#ifdef HAVE_ACML
  int n = phi.size();
  float c[n], s[n];
  vrsa_sincosf(n, &(phi[0]), s, c);
  for (int i=0; i<n; i++)
    z[i] = std::complex<float>(c[i],s[i]);
#elifdef HAVE_MKL
  vcCIS (n, &(phi[0]), (float*) &(z[0]));
#else
  for (int i=0; i<n; i++)
    sincos (phi[i], &(z[i].imag()), &(z[i].real()));
#endif
}


#endif

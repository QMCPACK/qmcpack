//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Kenneth Esler 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef BSPLINE_ONEDIM_SOLVERS
#define BSPLINE_ONEDIM_SOLVERS
#include <vector>
#include <complex>
/** dummy declaration to be specialized
 *
 * Specializations for double and complex<double>
 */
template<typename T> 
struct SolvePeriodicInterp1D
{ };

template<typename T>
inline void FuncSolvePeriodicInterp1D(const std::vector<T>& data, std::vector<T>& p)
{
  const T ratio = 0.25;
  int N = data.size();
  std::vector<T> d(N), gamma(N), mu(N);
  //d = 1.5*data;
  for(int i=0; i<N; i++) d[i]=1.5*data[i];
  // First, eliminate leading coefficients
  gamma [0] = ratio;
  mu[0] = ratio;
  mu[N-1] = ratio;
  gamma[N-1] = 1.0;
  for (int row=1; row <(N-1); row++) {
    double diag = 1.0- mu[row-1]*ratio;
    double diagInv = 1.0/diag;
    gamma[row] = -ratio*gamma[row-1]*diagInv;
    mu[row] = diagInv*ratio;
    d[row]  = diagInv*(d[row]-ratio*d[row-1]);
    // Last row
    d[N-1] -= mu[N-1] * d[row-1];
    gamma[N-1] -= mu[N-1]*gamma[row-1];
    mu[N-1] = -mu[N-1]*mu[row-1];
  }
  // Last row:  gamma(N-1) hold diagonal element
  mu[N-1] += ratio;
  gamma[N-1] -= mu[N-1]*(mu[N-2]+gamma[N-2]);
  d[N-1] -= mu[N-1] * d[N-2];
  p[N-1] = d[N-1]/gamma[N-1];

  // Now go back upward, back substituting
  for (int row=N-2; row>=0; row--) 
    p[row] = d[row] - mu[row]*p[row+1] - gamma[row]*p[N-1];
}

template<> 
struct SolvePeriodicInterp1D<double>
{
  static inline void apply(const std::vector<double>& data, std::vector<double>& p)
  {
    FuncSolvePeriodicInterp1D(data,p);
  }

  template<typename CT>
  static inline void apply(const CT& data, CT& p, int N)
  {
    const double ratio = 0.25;
    CT d(N), gamma(N), mu(N);
    for(int i=0; i<N; i++) d[i]=1.5*data[i];
    // First, eliminate leading coefficients
    gamma [0] = ratio;
    mu[0] = ratio;
    mu[N-1] = ratio;
    gamma[N-1] = 1.0;
    for (int row=1; row <(N-1); row++) {
      double diag = 1.0- mu[row-1]*ratio;
      double diagInv = 1.0/diag;
      gamma[row] = -ratio*gamma[row-1]*diagInv;
      mu[row] = diagInv*ratio;
      d[row]  = diagInv*(d[row]-ratio*d[row-1]);
      // Last row
      d[N-1] -= mu[N-1] * d[row-1];
      gamma[N-1] -= mu[N-1]*gamma[row-1];
      mu[N-1] = -mu[N-1]*mu[row-1];
    }
    // Last row:  gamma(N-1) hold diagonal element
    mu[N-1] += ratio;
    gamma[N-1] -= mu[N-1]*(mu[N-2]+gamma[N-2]);
    d[N-1] -= mu[N-1] * d[N-2];
    p[N] = d[N-1]/gamma[N-1];

    // Now go back upward, back substituting
    for (int row=N-2; row>=0; row--)
      p[row+1] = d[row] - mu[row]*p[row+2] - gamma[row]*p[N];

  }
};

template<>
struct SolvePeriodicInterp1D<std::complex<double> >
{

  typedef std::complex<double> value_type;
  typedef double real_type;

  static inline void apply(const std::vector<value_type>& data, std::vector<value_type>& p)
  {
    int N = data.size();
    std::vector<real_type> dataReal(N), dataImag(N), pReal(N), pImag(N);
    for (int i=0; i<N; i++) {
      dataReal[i] = data[i].real();
      dataImag[i] = data[i].imag();
    }
    SolvePeriodicInterp1D<double>::apply(dataReal, pReal);
    SolvePeriodicInterp1D<double>::apply(dataImag, pImag);
    for (int i=0; i<N; i++) 
      p[i] = value_type(pReal[i], pImag[i]);
  }
};


template<typename T> 
struct SolveFirstDerivInterp1D
{ };

template<>
struct SolveFirstDerivInterp1D<double>
{
  template<class CT>
  static inline void apply(const CT& data, CT& p, int N, double* bcLower, double* bcUpper)
  {
    const double ratio=0.25;
    CT d(N+2), mu(N+2,ratio);
    for(int i=1; i<=N; i++) d[i]=1.5*data[i-1];
    double al=0.25*bcLower[0];
    double bl=0.25*bcLower[1];
    double cl=0.25*bcLower[2];
    double dl=1.5*bcLower[3];
    double ar=0.25*bcUpper[0];
    double br=0.25*bcUpper[1];
    double cr=0.25*bcUpper[2];
    double dr=1.5*bcUpper[3];

    // First, eliminate leading coefficients
    double alInv = 1.0/al;
    bl *= alInv;
    cl *= alInv;
    dl *= alInv;

    d[0] = dl;  
    mu[0] = bl;
    mu[1] = ratio - ratio*cl;

    for (int row=1; row <=N; row++) {
      double diag = 1.0- mu[row-1]*ratio;
      double diagInv = 1.0/diag;
      mu[row] *= diagInv;
      d[row]  = diagInv*(d[row]-ratio*d[row-1]);
    }

    br -= ar*mu[N-1];
    dr -= ar*d[N-1];
    cr -= br*mu[N];
    dr -= br*d[N];
    p[N+1] = dr/cr; 

    // Now go back upward, back substituting
    for (int row=N; row>=1; row--) 
      p[row] = d[row] - mu[row]*p[row+1];

    // And do 0th row
    p[0] = dl -bl*p[1] - cl*p[2];
  }
};

template<>
struct SolveFirstDerivInterp1D<std::complex<double> >
{

  typedef std::complex<double> value_type;
  typedef double real_type;

  template<class CT>
  static inline void apply(const CT& data, CT& p, int N, double* bcLower, double* bcUpper)
  {
    std::vector<real_type> dataReal(N), dataImag(N), pReal(N), pImag(N);
    for (int i=0; i<N; i++) {
      dataReal[i] = data[i].real();
      dataImag[i] = data[i].imag();
    }
    SolveFirstDerivInterp1D<double>::apply(dataReal, pReal, N, bcLower, bcUpper);
    SolveFirstDerivInterp1D<double>::apply(dataImag, pImag, N, bcLower, bcUpper);
    for (int i=0; i<N; i++) 
      p[i] = value_type(pReal[i], pImag[i]);
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

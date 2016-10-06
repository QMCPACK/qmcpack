//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef EXP_FIT_CLASS_H
#define EXP_FIT_CLASS_H

#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus
{
template<int M> class ComplexExpFitClass;

template<int M>
class ExpFitClass
{
private:
  TinyVector<double,M> Coefs, dCoefs, d2Coefs;
  double sign;
public:
  inline void Fit (std::vector<double> &r, std::vector<double> &u);
  inline void FitCusp(std::vector<double> &r, std::vector<double> &u, double cusp);
  inline void eval (double r, double &u);
  inline void eval (double r, double &u, double &du, double &d2u);
  friend class ComplexExpFitClass<M>;
};

template<int M> void
ExpFitClass<M>::FitCusp (std::vector<double> &r, std::vector<double> &u, double cusp)
{
  int N = r.size();
  sign = u[0] < 0.0 ? -1.0 : 1.0;
  if (r.size() != u.size())
    app_error() << "Different number of rows of basis functions than"
                << " of data points in LinFit.  Exitting.\n";
  std::vector<TinyVector<double,M-1> > F(N);
  std::vector<double> log_u(N);
  for (int i=0; i<N; i++)
  {
    log_u[i] = std::log (sign*u[i]) - cusp * r[i];
    double r2jp1 = r[i]*r[i];
    F[i][0] = 1.0;
    for (int j=1; j<M-1; j++)
    {
      F[i][j] = r2jp1;
      r2jp1 *= r[i];
    }
  }
  // Next, construct alpha matrix
  Matrix<double> alpha(M-1,M-1), alphaInv(M-1,M-1);
  alpha = 0.0;
  for (int j=0; j<M-1; j++)
    for (int k=0; k<M-1; k++)
    {
      alpha(k,j) = 0.0;
      for (int i=0; i<N; i++)
        alpha(k,j) += F[i][j] * F[i][k];
    }
  // Next, construct beta vector
  TinyVector<double,M-1> beta;
  beta = 0.0;
  for (int k=0; k<M-1; k++)
    for (int i=0; i<N; i++)
      beta[k] += log_u[i]*F[i][k];
  // Now, invert alpha
  for (int i=0; i<M-1; i++)
    for (int j=0; j<M-1; j++)
      alphaInv(i,j) = alpha(i,j);
  double det = invert_matrix(alphaInv);
  TinyVector<double,M-1> c;
  for (int i=0; i<M-1; i++)
  {
    c[i] = 0.0;
    for (int j=0; j<M-1; j++)
      c[i] += alphaInv(i,j) * beta[j];
  }
  Coefs[0] = c[0];
  Coefs[1] = cusp;
  for (int i=2; i<M; i++)
    Coefs[i] = c[i-1];
  dCoefs  = 0.0;
  d2Coefs = 0.0;
  for (int i=0; i<M-1; i++)
    dCoefs[i] = (double)(i+1) * Coefs[i+1];
  for (int i=0; i<M-2; i++)
    d2Coefs[i] = (double)(i+1) * dCoefs[i+1];
}

template<int M> void
ExpFitClass<M>::Fit (std::vector<double> &r, std::vector<double> &u)
{
  int N = r.size();
  sign = u[0] < 0.0 ? -1.0 : 1.0;
  if (r.size() != u.size())
    app_error() << "Different number of rows of basis functions than"
                << " of data points in LinFit.  Exitting.\n";
  std::vector<TinyVector<double,M> > F(N);
  std::vector<double> log_u(N);
  for (int i=0; i<N; i++)
  {
    log_u[i] = std::log (sign*u[i]);
    double r2j=1.0;
    for (int j=0; j<M; j++)
    {
      F[i][j] = r2j;
      r2j *= r[i];
    }
  }
  // Next, construct alpha matrix
  Matrix<double> alpha(M,M), alphaInv(M,M);
  alpha = 0.0;
  for (int j=0; j<M; j++)
    for (int k=0; k<M; k++)
    {
      alpha(k,j) = 0.0;
      for (int i=0; i<N; i++)
        alpha(k,j) += F[i][j] * F[i][k];
    }
  // Next, construct beta vector
  TinyVector<double,M> beta;
  beta = 0.0;
  for (int k=0; k<M; k++)
    for (int i=0; i<N; i++)
      beta[k] += log_u[i]*F[i][k];
  // Now, invert alpha
  for (int i=0; i<M; i++)
    for (int j=0; j<M; j++)
      alphaInv(i,j) = alpha(i,j);
  double det = invert_matrix(alphaInv);
  for (int i=0; i<M; i++)
  {
    Coefs[i] = 0.0;
    for (int j=0; j<M; j++)
      Coefs[i] += alphaInv(i,j) * beta[j];
  }
  dCoefs  = 0.0;
  d2Coefs = 0.0;
  for (int i=0; i<M-1; i++)
    dCoefs[i] = (double)(i+1) * Coefs[i+1];
  for (int i=0; i<M-2; i++)
    d2Coefs[i] = (double)(i+1) * dCoefs[i+1];
}



template<int M> void
ExpFitClass<M>::eval (double r, double &u)
{
  double r2j = 1.0;
  double poly = 0.0;
  for (int j=0; j<M; j++)
  {
    poly += Coefs[j] * r2j;
    r2j *= r;
  }
  u = sign*std::exp(poly);
}

template<int M> void
ExpFitClass<M>::eval (double r, double &u, double &du, double &d2u)
{
  double r2j = 1.0;
  double P=0.0, dP=0.0, d2P=0.0;
  for (int j=0; j<M; j++)
  {
    P   +=   Coefs[j] * r2j;
    dP  +=  dCoefs[j] * r2j;
    d2P += d2Coefs[j] * r2j;
    r2j *= r;
  }
  u   = sign*std::exp (P);
  du  = dP * u;
  d2u = (d2P + dP*dP)*u;
}

template<int M>
class ComplexExpFitClass
{
private:
  TinyVector<std::complex<double>,M> Coefs, dCoefs, d2Coefs;
  double realSign, imagSign;
public:
  inline void Fit (std::vector<double> &r, std::vector<std::complex<double> > &u);
  inline void FitCusp(std::vector<double> &r, std::vector<std::complex<double> > &u, double cusp);
  inline void eval (double r, std::complex<double> &u);
  inline void eval (double r, std::complex<double> &u, std::complex<double> &du, std::complex<double> &d2u);
};


template<int M> void
ComplexExpFitClass<M>::FitCusp (std::vector<double> &r, std::vector<std::complex<double> > &u, double cusp)
{
  int N = u.size();
  ExpFitClass<M> realFit, imagFit;
  std::vector<double> realVals(N), imagVals(N);
  for (int i=0; i<N; i++)
  {
    realVals[i] = real(u[i]);
    imagVals[i] = imag(u[i]);
  }
  realFit.FitCusp (r, realVals, cusp);
  imagFit.FitCusp (r, imagVals, cusp);
  realSign = realFit.sign;
  imagSign = imagFit.sign;
  for (int i=0; i<M; i++)
  {
    Coefs[i]   = std::complex<double>(realFit.Coefs[i]  ,   imagFit.Coefs[i]);
    dCoefs[i]  = std::complex<double>(realFit.dCoefs[i] ,  imagFit.dCoefs[i]);
    d2Coefs[i] = std::complex<double>(realFit.d2Coefs[i], imagFit.d2Coefs[i]);
  }
}

template<int M> void
ComplexExpFitClass<M>::Fit (std::vector<double> &r, std::vector<std::complex<double> > &u)
{
  int N = u.size();
  ExpFitClass<M> realFit, imagFit;
  std::vector<double> realVals(N), imagVals(N);
  for (int i=0; i<N; i++)
  {
    realVals[i] = real(u[i]);
    imagVals[i] = imag(u[i]);
  }
  realFit.Fit (r, realVals);
  imagFit.Fit (r, imagVals);
  realSign = realFit.sign;
  imagSign = imagFit.sign;
  for (int i=0; i<M; i++)
  {
    Coefs[i]   = std::complex<double>(realFit.Coefs[i]  ,   imagFit.Coefs[i]);
    dCoefs[i]  = std::complex<double>(realFit.dCoefs[i] ,  imagFit.dCoefs[i]);
    d2Coefs[i] = std::complex<double>(realFit.d2Coefs[i], imagFit.d2Coefs[i]);
  }
}


template<int M> void
ComplexExpFitClass<M>::eval (double r, std::complex<double> &u)
{
  double r2j = 1.0;
  std::complex<double> P = std::complex<double>();
  for (int j=0; j<M; j++)
  {
    P += Coefs[j] * r2j;
    r2j *= r;
  }
  u = std::complex<double>(realSign*std::exp(P.real()), imagSign*std::exp(P.imag()));
}

template<int M> void
ComplexExpFitClass<M>::eval (double r, std::complex<double> &u,
                             std::complex<double> &du, std::complex<double> &d2u)
{
  double r2j = 1.0;
  std::complex<double> P, dP, d2P;
  P = dP = d2P = std::complex<double>();
  for (int j=0; j<M; j++)
  {
    P   +=   Coefs[j] * r2j;
    dP  +=  dCoefs[j] * r2j;
    d2P += d2Coefs[j] * r2j;
    r2j *= r;
  }
  u= std::complex<double>(realSign * std::exp (P.real()),imagSign * std::exp (P.imag()));
  du= std::complex<double>(dP.real() * u.real(),dP.imag() * u.imag());
  d2u= std::complex<double>((d2P.real() + dP.real()*dP.real())*u.real(),(d2P.imag() + dP.imag()*dP.imag())*u.imag());
  //u.real()   = realSign * std::exp (P.real());
  //du.real()  = dP.real() * u.real();
  //d2u.real() = (d2P.real() + dP.real()*dP.real())*u.real();
  //u.imag()   = imagSign * std::exp (P.imag());
  //du.imag()  = dP.imag() * u.imag();
  //d2u.imag() = (d2P.imag() + dP.imag()*dP.imag())*u.imag();
}

}

#endif

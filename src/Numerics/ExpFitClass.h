//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Ken Esler                            //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: esler@uiuc.edu                                     //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef EXP_FIT_CLASS_H
#define EXP_FIT_CLASS_H

#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "Numerics/DeterminantOperators.h"


namespace qmcplusplus {
  template<int M> class ComplexExpFitClass; 

  template<int M>
  class ExpFitClass 
  {
  private:
    TinyVector<double,M> Coefs, dCoefs, d2Coefs;
  public:
    inline void Fit (vector<double> &r, vector<double> &u);
    inline void FitCusp(vector<double> &r, vector<double> &u, double cusp);
    inline void eval (double r, double &u);
    inline void eval (double r, double &u, double &du, double &d2u);
    friend class ComplexExpFitClass<M>;
  };

  template<int M> void 
  ExpFitClass<M>::FitCusp (vector<double> &r, vector<double> &u, double cusp)
  {
    int N = r.size();

    if (r.size() != u.size()) 
      app_error() << "Different number of rows of basis functions than"
		  << " of data points in LinFit.  Exitting.\n";
    vector<TinyVector<double,M-1> > F(N);
    vector<double> log_u(N);
    for (int i=0; i<N; i++) {
      log_u[i] = std::log (u[i]) - cusp * r[i];
      double r2jp1 = r[i]*r[i];
      F[i][0] = 1.0;
      for (int j=1; j<M-1; j++) {
	F[i][j] = r2jp1;
	r2jp1 *= r[i];
      }
    }
      
    // Next, construct alpha matrix
    Matrix<double> alpha(M-1,M-1), alphaInv(M-1,M-1);
    alpha = 0.0;
    for (int j=0; j<M-1; j++)
      for (int k=0; k<M-1; k++) {
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

    for (int i=0; i<M-1; i++) {
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
  ExpFitClass<M>::Fit (vector<double> &r, vector<double> &u)
  {
    int N = r.size();

    if (r.size() != u.size()) 
      app_error() << "Different number of rows of basis functions than"
		  << " of data points in LinFit.  Exitting.\n";
    vector<TinyVector<double,M> > F(N);
    vector<double> log_u(N);
    for (int i=0; i<N; i++) {
      log_u[i] = std::log (u[i]);
      double r2j 1.0;
      for (int j=0; j<M; j++) {
	F[i][j] = r2j;
	r2j *= r[i];
      }
    }
      
    // Next, construct alpha matrix
    Matrix<double> alpha(M,M), alphaInv(M,M);
    alpha = 0.0;
    for (int j=0; j<M; j++)
      for (int k=0; k<M; k++) {
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
	    
    for (int i=0; i<M; i++) {
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
    for (int j=0; j<M; j++) {
      poly += Coefs[j] * r2j;
      rj2 *= r;
    }
    return std::exp(poly);
  }

  template<int M> void
  ExpFitClass<M>::eval (double r, double &u, double &du, double &d2u) {
    double r2j = 1.0;
    double P=0.0, dP=0.0, d2P=0.0;
    for (int j=0; j<M; j++) {
      P   +=   Coefs[j] * r2j;
      dP  +=  dCoefs[j] * r2j;      
      d2P += d2Coefs[j] * r2j;
      r2j *= r;
    }
    u   = std::exp (P);
    du  = dP * u;
    d2u = (d2P + dP*dP)*u;
  }

  template<int M>
  class ComplexExpFitClass
  {
  private:
    TinyVector<complex<double>,M> Coefs, dCoefs, d2Coefs;
  public:
    inline void Fit (vector<double> &r, vector<complex<double> > &u);
    inline void FitCusp(vector<double> &r, vector<complex<double> > &u, double cusp);
    inline void eval (double r, complex<double> &u);
    inline void eval (double r, complex<double> &u, complex<double> &du, complex<double> &d2u);
  };


  template<int M> void
  ComplexExpFitClass<M>::FitCusp (vector<double> &r, vector<complex<double> > &u, double cusp)
  {
    int N = u.size();
    ExpFitClass<M> realFit, imagFit;
    vector<double> realVals(N), imagVals[N];
    for (int i=0; i<N; i++) {
      realVals[i] = real(u[i]);
      iamgVals[i] = imag(u[i]);
    }
    realFit.FitCusp (r, realVals, cusp);
    imagFit.FitCusp (r, imagVals, cusp);
    for (int i=0; i<M; i++) {
      Coefs[i]   = complex<double>(realFit.Coefs[i]  ,   imag.Coefs[i]);
      dCoefs[i]  = complex<double>(realFit.dCoefs[i] ,  imag.dCoefs[i]);
      d2Coefs[i] = complex<double>(realFit.d2Coefs[i], imag.d2Coefs[i]);
    }
  }

  template<int M> void
  ComplexExpFitClass<M>::Fit (vector<double> &r, vector<complex<double> > &u)
  {
    int N = u.size();
    ExpFitClass<M> realFit, imagFit;
    vector<double> realVals(N), imagVals[N];
    for (int i=0; i<N; i++) {
      realVals[i] = real(u[i]);
      iamgVals[i] = imag(u[i]);
    }
    realFit.Fit (r, realVals);
    imagFit.Fit (r, imagVals);
    for (int i=0; i<M; i++) {
      Coefs[i]   = complex<double>(realFit.Coefs[i]  ,   imag.Coefs[i]);
      dCoefs[i]  = complex<double>(realFit.dCoefs[i] ,  imag.dCoefs[i]);
      d2Coefs[i] = complex<double>(realFit.d2Coefs[i], imag.d2Coefs[i]);
    }
  }


}

#endif

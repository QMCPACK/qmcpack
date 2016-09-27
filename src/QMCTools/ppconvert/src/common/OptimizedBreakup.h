//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef OPTIMIZED_BREAKUP_H
#define OPTIMIZED_BREAKUP_H

#include "blitz/array.h"

using namespace blitz;

class BasisClass
{
public:
  //protected:
  double r_c;
  TinyVector<double,3> Box;
  double Omega;
public:
  /// Set the cutoff radius
  virtual void Set_rc(double rc) = 0;
  inline double Get_rc() { return r_c; }
  inline void SetBox (TinyVector<double,3> box);
  inline TinyVector<double,3> GetBox ();
  /// Returns the number of basis elements
  virtual int NumElements() = 0;
  /// Returns the basis element n evaluated in real space at r
  virtual double h(int n, double r) = 0;
  /// Returns the basis element n evaluated in k space at k
  virtual double c(int n, double k) = 0;
  virtual double dc_dk(int n, double k) = 0;
  double c_numerical (int n, double k);
  /// This returns the coefficent of the nth basis function
  //virtual double  Get_t(int n) const     = 0;
  /// This sets the coefficent of the nth basis function
  //virtual void    Set_t(int n, double t) = 0;
  /// This returns the linear combination of the basis functions with
  /// coefficients t_n
  //virtual double f (double r) = 0;
  BasisClass() : r_c (0.0)
  { /* do nothing */ }
};


class OptimizedBreakupClass
{
private:
  BasisClass &Basis;
  void Addk(double k, double degeneracy=1.0);
public:
  // First element is |k|, second is degeneracy of the point.
  Array<TinyVector<double,2>,1> kpoints;
  void SetkVecs(double kc, double kcont, double kMax);
  /// kc is the k-space cutoff for the Ewald sum.  
  /// kMax is largest k we use in determining the error in the breakup.  
  /// t is the set of coefficients of the breakup.
  /// inFit is a boolean array telling whether t_n should be optimized
  /// or left at its initial value.  Returns chi-squared for the breakup.
  double DoBreakup (const Array<double,1> &Vk, Array<double,1> &t, 
		    const Array<bool,1> &adjust);
  /// Same as above, but we assume that all t's are adjusted.
  double DoBreakup (const Array<double,1> &Vk, Array<double,1> &t);
  OptimizedBreakupClass (BasisClass &basis) : Basis(basis)
  { /* Do nothing */ }
};


inline void BasisClass::SetBox (TinyVector<double,3> box)
{
  Box = box;
  Omega = box[0]*box[1]*box[2];
}

inline TinyVector<double,3> BasisClass::GetBox ()
{
  return Box;
}



/// Locally Piecewise Quintic Hermite Interpolant
class LPQHI_BasisClass : public BasisClass
{
  /// public is HACK
  //private:
public:
  int NumKnots;
  double delta, deltaInv;
  TinyMatrix<double,3,6> S;
  /// The following are helpers to calculate the Fourier tranform of
  /// the basis functions
  inline std::complex<double> Eplus(int i, double k, int n);
  inline std::complex<double> Eminus(int i, double k, int n);
  inline double Dplus(int i, double k, int n);
  inline double Dminus(int i, double k, int n);
  inline std::complex<double> dEplus_dk(int i, double k, int n);
  inline std::complex<double> dEminus_dk(int i, double k, int n);
  inline double dDplus_dk(int i, double k, int n);
  inline double dDminus_dk(int i, double k, int n);
  Array<double,1> tvec;
public:
  inline double GetDelta() { return delta; }
  // n must be at least 2;
  void SetNumKnots(int n);
  void Set_rc(double rc);
  int NumElements();
  double h(int n, double r);
  double c(int n, double k);
  double dc_dk(int n, double k);
  LPQHI_BasisClass() : NumKnots(0), delta(0.0) 
  { 
    S = 
      1.0,   0.0,   0.0, -10.0,  15.0,  -6.0,
      0.0,   1.0,   0.0,  -6.0,   8.0,  -3.0,
      0.0,   0.0,   0.5,  -1.5,   1.5,  -0.5;
  }
};

inline std::complex<double> LPQHI_BasisClass::Eplus(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0);

  if (n == 0) {
    std::complex<double> e1(cos(k*delta)-1.0, sin(k*delta));
    std::complex<double> e2(cos(k*delta*i),   sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    std::complex<double> t1, t2;
    double sign = 1.0;
    t1 = std::complex<double>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    t2=-(double)n/delta*Eplus(i,k,n-1);;
    return (-(eye/k)*(t1+t2));
  }
}

inline std::complex<double> LPQHI_BasisClass::dEplus_dk(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0);

  if (n == 0) {
    std::complex<double>  e1(cos(k*delta)-1.0,        sin(k*delta));
    std::complex<double> de1(-delta*sin(k*delta),     delta*cos(k*delta));
    std::complex<double>  e2(cos(k*delta*i),          sin(k*delta*i));
    std::complex<double> de2(-delta*i*sin(k*delta*i), delta*i*cos(k*delta*i));
    return ((eye/(k*k))*e1*e2-(eye/k)*(de1*e2+e1*de2));
  }
  else {
    std::complex<double> t1, t2, dt1, dt2;
    double sign = 1.0;
    t1  = std::complex<double>(cos(k*delta*(i+1)),sin(k*delta*(i+1)));
    dt1 = std::complex<double>(-delta*(i+1)*sin(k*delta*(i+1)),
			  delta*(i+1)*cos(k*delta*(i+1)));
    t2  = -(double)n/delta*Eplus(i,k,n-1);
    dt2 = -(double)n/delta*dEplus_dk(i,k,n-1);
    return ((eye/(k*k))*(t1+t2) -(eye/k)*(dt1+dt2));
  }
}

inline std::complex<double> LPQHI_BasisClass::Eminus(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0);

  if (n == 0) {
    std::complex<double> e1(cos(k*delta)-1.0, -sin(k*delta));
    std::complex<double> e2(cos(k*delta*i),    sin(k*delta*i));
    return (-(eye/k)*e1*e2);
  }
  else {
    std::complex<double> t1, t2;
    double sign = (n & 1) ? -1.0 : 1.0;
    t1 = sign*
      std::complex<double> (cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    t2=-(double)n/delta*Eminus(i,k,n-1);
    return (-(eye/k)*(t1+t2));
  }
}

inline std::complex<double> LPQHI_BasisClass::dEminus_dk(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0);

  if (n == 0) {
    std::complex<double>  e1(cos(k*delta)-1.0,        -sin(k*delta));
    std::complex<double> de1(-delta*sin(k*delta),     -delta*cos(k*delta));
    std::complex<double>  e2(cos(k*delta*i),          sin(k*delta*i));
    std::complex<double> de2(-delta*i*sin(k*delta*i), delta*i*cos(k*delta*i));
    return ((eye/(k*k))*e1*e2-(eye/k)*(de1*e2+e1*de2));
  }
  else {
    std::complex<double> t1, t2, dt1, dt2;
    double sign = (n & 1) ? -1.0 : 1.0;
    t1  = sign * std::complex<double>(cos(k*delta*(i-1)),sin(k*delta*(i-1)));
    dt1 = sign * std::complex<double>(-delta*(i-1)*sin(k*delta*(i-1)),
			  delta*(i-1)*cos(k*delta*(i-1)));
    t2  = -(double)n/delta*Eminus(i,k,n-1);
    dt2 = -(double)n/delta*dEminus_dk(i,k,n-1);
    return ((eye/(k*k))*(t1+t2) -(eye/k)*(dt1+dt2));
  }
}



inline double LPQHI_BasisClass::Dplus(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0); 
  std::complex<double> Z1 = Eplus(i,k,n+1);
  std::complex<double> Z2 = Eplus(i,k,n);
  return 4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}

inline double LPQHI_BasisClass::dDplus_dk(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0); 
  std::complex<double> Z1 = Eplus(i,k,n+1);
  std::complex<double> Z2 = Eplus(i,k,n);
  std::complex<double> dZ1 = dEplus_dk(i,k,n+1);
  std::complex<double> dZ2 = dEplus_dk(i,k,n);
  return -4.0*M_PI/(k*k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag()) +
    4.0*M_PI/(k*Omega)*(delta*dZ1.imag() + i*delta*dZ2.imag());
}

inline double LPQHI_BasisClass::Dminus(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0); 
  std::complex<double> Z1 = Eminus(i,k,n+1);
  std::complex<double> Z2 = Eminus(i,k,n);
  return -4.0*M_PI/(k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag());
}

inline double LPQHI_BasisClass::dDminus_dk(int i, double k, int n)
{
  std::complex<double> eye(0.0, 1.0); 
  std::complex<double> Z1  = Eminus(i,k,n+1);
  std::complex<double> dZ1 = dEminus_dk(i,k,n+1);
  std::complex<double> Z2  = Eminus(i,k,n);
  std::complex<double> dZ2 = dEminus_dk(i,k,n);
  return 4.0*M_PI/(k*k*Omega)*(delta* Z1.imag() + i*delta*Z2.imag())
    -4.0*M_PI/(k*Omega)*(delta* dZ1.imag() + i*delta*dZ2.imag());
}


#endif

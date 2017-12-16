//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_SLATERTYPEORBITAL_H
#define QMCPLUSPLUS_SLATERTYPEORBITAL_H
#include "Numerics/OptimizableFunctorBase.h"
#include <cmath>

/** class to evaluate the normalization factors for the Slater-Type orbitals
 */
template<class T>
struct STONorm
{

  std::vector<T> Factorial;

  explicit STONorm(int nmax=1)
  {
    set(nmax);
  }

  inline void set(int nmax)
  {
    int n = 2*nmax+2;
    Factorial.resize(n+1);
    Factorial[0] = 1.0;
    for (int i=1; i<n+1; i++)
      Factorial[i] = Factorial[i-1]*static_cast<T>(i);
  }

  inline T operator()(int n, T screen)
  {
    return
      1.0/sqrt(Factorial[2*n+2]/pow(2.0*screen,2*n+3));
    //return
    //  1.0/sqrt(Factorial[2*n+2]*4.0*(4.0*atan(1.0))/pow(2.0*screen,2*n+3));
  }

};



/** Generic Slater-Type Orbital
 *
 * This class evalaute \f$ \frac{\chi_{n,\xi}(r)}{r^l} \f$
 * where a normalized STO is defined as
 * \f$\chi_{n,\xi}(r) =  C r^{n-1} \exp{-\xi r}\f$
 * C is a contraction factor.
 * The physical principal quantum number n has to be a positive integer,
 * where n-1 is the number of nodes.
 */
template<class T>
struct GenericSTO: public OptimizableFunctorBase
{

  typedef T value_type;

  int ID;
  ///Principal number
  int N;
  ///N-l-1
  int Power;
  real_type Z;
  real_type Norm;
  real_type Y, dY, d2Y, d3Y;

  GenericSTO(): N(-1), Power(0), Z(1.0), Norm(1.0) { }

  /** constructor with a known contraction factor
   */
  explicit GenericSTO(int power, real_type z, real_type norm=1.0):
    N(-1), Power(power), Z(z), Norm(norm) { }

  /** constructor with a set of quantum numbers
   * @param n principal quantum number
   * @param l angular quantum number
   * @param z exponent
   *
   * Power = n-l-1
   * Contraction factor is the normalization factor evaluated based on N and Z.
   */
  explicit GenericSTO(int n, int l, real_type z) :
    N(n), Power(n-l-1), Z(z)
  {
    reset();
  }


  OptimizableFunctorBase* makeClone() const
  {
    return new GenericSTO<T>(*this);
  }

  inline void reset()
  {
    if(N>0)
    {
      STONorm<T> anorm(N);
      Norm = anorm(N-1,Z);
    }
  }

  inline void setgrid(real_type r) { }

  inline real_type f(real_type r)
  {
    return exp(-Z*r)*Norm*pow(r,Power);
  }

  inline real_type df(real_type r)
  {
    real_type rnl = exp(-Z*r)*Norm;
    if(Power == 0)
    {
      return  -Z*rnl;
    }
    else
    {
      return rnl*pow(r,Power)*(Power/r-Z);
    }
  }

  /** return the value only
   * @param r distance
   * @param rinv inverse of r
   */
  inline real_type evaluate(real_type r, real_type rinv)
  {
    return Y = Norm*pow(r,Power)*exp(-Z*r);
  }

  inline void evaluateAll(real_type r, real_type rinv)
  {
    Y = evaluate(r,rinv,dY,d2Y);
  }

  inline real_type evaluate(real_type r, real_type rinv, real_type& drnl, real_type& d2rnl)
  {
    real_type rnl = Norm*exp(-Z*r);
    if(Power == 0)
    {
      drnl = -Z*rnl;
      d2rnl = rnl*Z*Z;
    }
    else
    {
      rnl *= pow(r,Power);
      real_type x = Power*rinv-Z;
      drnl = rnl*x;
      d2rnl = rnl*(x*x-Power*rinv*rinv);
    }
    return rnl;
  }

  inline real_type evaluate(real_type r, real_type rinv, real_type& drnl, real_type& d2rnl, real_type& d3rnl)
  {
    real_type rnl = Norm*exp(-Z*r);
    if(Power == 0)
    {
      drnl = -Z*rnl;
      d2rnl = rnl*Z*Z;
      d3rnl = -d2rnl*Z;
    }
    else
    {
      rnl *= pow(r,Power);
      real_type x = Power*rinv-Z;
      drnl = rnl*x;
      d2rnl = rnl*(x*x-Power*rinv*rinv);
      d3rnl = rnl*(x*x*x+(2.0*rinv-3.0*x)*Power*rinv*rinv);
    }
    return rnl;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void checkInVariables(opt_variables_type& active) { }
  void checkOutVariables(const opt_variables_type& active) { }
  void resetParameters(const opt_variables_type& active)
  {
    //DO NOTHING FOR NOW
  }

};

/**class for Slater-type orbitals,
 *
 *@f[
 *\Psi_{n,l,m}({\bf R}) = N r^{n-1} \exp{-Zr} Y_{lm}(\theta,\phi)
 *@f]
 */
template<class T>
struct RadialSTO
{
  typedef T real_type;
  int NminusOne;
  T Z;
  T Norm;
  T Y, dY, d2Y;
  RadialSTO(): NminusOne(0), Z(1.0), Norm(1.0) { }
  RadialSTO(int n, double z, double norm=1.0):
    NminusOne(n-1), Z(z),Norm(norm) { }

  inline void setgrid(T r) { }

  inline T f(T r) const
  {
    return pow(r,NminusOne)*exp(-Z*r)*Norm;
  }

  inline T df(T r) const
  {
    T rnl = pow(r,NminusOne)*exp(-Z*r)*Norm;
    return  (NminusOne/r-Z)*rnl;
  }

  inline T evaluate(T r)
  {
    return pow(r,NminusOne)*exp(-Z*r)*Norm;
  }

  inline void evaluateAll(T r, T rinv)
  {
    Y = evaluate(r,rinv,dY,d2Y);
  }

  inline T evaluate(T r, T rinv, T& drnl, T& d2rnl)
  {
    T rnl = pow(r,NminusOne)*exp(-Z*r)*Norm;
    T x = NminusOne*rinv-Z;
    drnl = rnl*x;
    d2rnl = rnl*(x*x-NminusOne*rinv*rinv);
    return rnl;
  }
};



#endif

//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file PadeFunctors.h
 * @brief Functors which implement Pade functions
 */
#ifndef QMCPLUSPLUS_PADEFUNCTORS_H
#define QMCPLUSPLUS_PADEFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include <cmath>

namespace qmcplusplus {

  /** Implements a Pade Function \f$u[r]=A*r/(1+B*r)\f$
   * 
   * Similar to PadeJastrow with a scale.
   */
  template<class T>
    struct PadeFunctor:public OptimizableFunctorBase<T> {

      typedef typename OptimizableFunctorBase<T>::real_type real_type;

      ///input A
      real_type A;
      ///input B
      real_type B0;
      ///input scaling, default=1.0
      real_type Scale;
      ///B=B0*Scale
      real_type B; 
      ///AB=A*B
      real_type AB;
      ///B2=2*B
      real_type B2;
      ///AoverB=A/B 
      real_type AoverB;

      ///constructor
      explicit PadeFunctor(real_type a, real_type b, real_type s=1.0): Scale(s) {
        reset(a,b);
      }

      inline void reset() {
        B = B0*Scale;
        AB = A*B; B2=2.0*B;
        AoverB=A/B;
      }

      void reset(real_type a, real_type b, real_type s=1.0) {
        A=a; B0=b; Scale=s;
        reset();
      }

      inline real_type evaluate(real_type r) {
        return A*r/(1.0+B*r);
      }

      inline real_type 
        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
          real_type u = 1.0/(1.0+B*r);
          dudr = A*u*u;
          d2udr2 = -B2*dudr*u;
          return A*u*r;
        }

      inline real_type f(real_type r) {
        return evaluate(r)-AoverB;
      }

      inline real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }

      bool put(xmlNodePtr cur) {return true;}
      void addOptimizables( VarRegistry<real_type>& vlist){}
    };

  /** Pade functional of \f[ u(r) = \frac{a*f(r)}{1+b*f(r)} \f] with a scale function f(r)
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct ScaledPadeFunctor:public OptimizableFunctorBase<T> {

      typedef typename OptimizableFunctorBase<T>::real_type real_type;

      ///coefficients
      real_type A, B, C; 
      real_type OneOverC, B2;

      ///constructor
      explicit ScaledPadeFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0) 
      {reset(a,b,c);}

      /** reset the internal variables.
       *
       * When RefPade is not 0, use RefPade->B to reset the values
       */
      inline void reset() {
        OneOverC=1.0/C;
        B2=2.0*B;
      }

      /** reset the internal variables.
       *@param a Pade Jastrow parameter a 
       *@param b Pade Jastrow parameter b 
       */
      void reset(real_type a, real_type b, real_type c) {
        A=a; B=b; C=c; 
        OneOverC=1.0/c;
        B2=2.0*b;
      }

      /** evaluate the value at r
       * @param r the distance
       * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
       */
      inline real_type evaluate(real_type r) {
        real_type reff((1.0-std::exp(-C*r))*OneOverC);
        return A*reff/(1.0+B*reff);
      }

      /** evaluate the value, first derivative and second derivative
       * @param r the distance
       * @param dudr return value  \f$ du/dr\f$
       * @param d2udr2 return value  \f$ d^2u/dr^2 \f$
       * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
       */
      inline real_type 
        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
          real_type reff_1(std::exp(-C*r));
          real_type reff((1.0-reff_1)*OneOverC);
          real_type u(1.0/(1.0+B*reff));
          real_type auu(A*u*u);
          dudr=reff_1*auu;
          //d2udr=auu*(-C*reff_1*reff_2-B2*reff_1*reff_1*u);
          d2udr2=-reff_1*auu*(C+B2*reff_1*u);
          return A*u*reff;
        }

      real_type f(real_type r) {
        return evaluate(r);
      }

      real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }

      bool put(xmlNodePtr cur) {return true;}
      void addOptimizables( VarRegistry<real_type>& vlist){}
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


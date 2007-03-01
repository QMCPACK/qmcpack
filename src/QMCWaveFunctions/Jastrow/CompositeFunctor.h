//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_COMPOSITEFUNCTORS_H
#define QMCPLUSPLUS_COMPOSITEFUNCTORS_H

#include "Numerics/OptimizableFunctorBase.h"
namespace qmcplusplus {

  /** Implements a truncated Pade Function \f$u[r]=A*r/(1+B*r)\f$ 
   *
   * Only valid for two-body jastrow with a "known" cusp condition.
   */
  template<class T>
    struct TruncatedPadeFunctor:public OptimizableFunctorBase<T> {

      typedef OptimizableFunctorBase<T> ThisBaseType;
      typedef typename ThisBaseType::real_type real_type;
      typedef typename ThisBaseType::OptimizableSetType OptimizableSetType;

      ///input A
      real_type A;
      ///derived parameters
      real_type B; 
      ///AB=A*B
      real_type AB;
      ///B2=2*B
      real_type B2;
      ///AoverB=A/B 
      real_type AoverB;
      ///cutoff
      real_type Rcut;
      ///maximum radius
      real_type Rmax;
      ///offset
      real_type OffSet;
      ///input function
      ThisBaseType* inFunc;

      /** constructor
       * @param a default a for e-e cusp condition
       * @param fin input function
       */
      TruncatedPadeFunctor(real_type a=-0.5, ThisBaseType* fin=0, real_type rmax=10):
        inFunc(fin),Rmax(rmax) {
        Rcut=0.0; 
        A=a;
      }

      inline real_type f(real_type r) {
        if (r>Rcut) 
          return inFunc->f(r);
        else
          return A*r/(1.0+B*r)+OffSet;
      }

      inline real_type df(real_type r) {
        if(r>Rcut)
          return inFunc->df(r);
        else 
        {
          real_type u = 1.0/(1.0+B*r);
          return A*u*u;
        }
      }

      bool put(xmlNodePtr cur) {return true;}
      void addOptimizables( VarRegistry<real_type>& vlist)
      {
        if(inFunc) inFunc->addOptimizables(vlist);
      }

      void resetParameters(OptimizableSetType& optVariables) 
      {
        if(inFunc) 
        {
          inFunc->resetParameters(optVariables);
        }
        applyCuspCondition();
      }

      void applyCuspCondition() {
        Rcut=0.0; 
        if(inFunc) {
          real_type x=0.001;
          real_type dx=0.001;
          real_type deriv0=inFunc->df(x),deriv;
          bool validInput=true;
          if(deriv0*A<0.0) {
            return;
          }
          if(A<0.0) 
          {
            bool goingdown=true;
            do {
              x+=dx;
              if(x>=Rmax) validInput=false;
              deriv=inFunc->df(x);
              if(deriv>deriv0) 
                goingdown=false;
              else 
                deriv0=deriv;
            } while(goingdown && validInput);
          }
          else
          {
            bool goingup=true;
            do {
              x+=dx;
              if(x>Rmax) validInput=false;
              deriv=inFunc->df(x);
              if(deriv<deriv0) 
                goingup=false;
              else
                deriv0=deriv;
            } while(goingup && validInput);
          }

          if(validInput)
          {
            Rcut=x+2*dx;
            deriv=inFunc->df(Rcut);
            B = (std::sqrt(A/deriv)-1.0)/Rcut;
            OffSet=inFunc->f(Rcut)-A*Rcut/(1.0+B*Rcut);
            AB = A*B; B2=2.0*B;
            AoverB=A/B;
          }
        }
      }
    };

  /** Implements \f$ u(r) = r^n*f(r) \f$ where \f$ f(r)\f$ is any OptimizableFunctorBase<T>
   * 
   * This functor is not optimized and should be used only as a temporary functor
   * for a final numerical functor.
   */
  template<class T>
    struct AnyTimesRnFunctor: public OptimizableFunctorBase<T> {
      typedef OptimizableFunctorBase<T> ThisBaseType;
      typedef typename ThisBaseType::real_type real_type;
      typedef typename ThisBaseType::OptimizableSetType OptimizableSetType;

      ///pointer to a functor
      OptimizableFunctorBase<T>* myFunc;
      ///power
      int Np;
      ///constructor
      AnyTimesRnFunctor(OptimizableFunctorBase<T>* infunc=0, int n=1): 
        myFunc(infunc),Np(n)
      {
      }

      inline real_type f(real_type r) {
        return std::pow(r,Np)*myFunc->f(r);
      }

      inline real_type df(real_type r) {
        real_type u=myFunc->f(r);
        real_type du=myFunc->df(r);
        return (Np*u+r*du)*std::pow(r,Np-1);
      }

      bool put(xmlNodePtr cur) 
      {
        return true;
      }

      void addOptimizables(VarRegistry<T>& vlist) 
      {
        if(myFunc) myFunc->addOptimizables(vlist);
      }

      void resetParameters(OptimizableSetType& optVariables) 
      { 
        if(myFunc) myFunc->resetParameters(optVariables);
      }

    };

  /** Implements \f$ u(r) = A*(-\frac{1}{B} exp(-Br))*f_c(r)\f$ with a cutoff function.
   * 
   * This functor is not optimized and should be used only as a temporary functor
   * for a final numerical functor to correct the cusp condition.
   */
  template<class T>
    struct CuspCorrectionFunctor: public OptimizableFunctorBase<T> {

      typedef OptimizableFunctorBase<T> ThisBaseType;
      typedef typename ThisBaseType::real_type real_type;
      typedef typename ThisBaseType::OptimizableSetType OptimizableSetType;

      ///scaling function or exponent
      real_type E;
      ///inverse of E
      real_type mOneOverE;
      ///maxium cutoff radius
      real_type Rmax;
      ///starting cutoff radius
      real_type Rcut;
      ///fixed maxmium cutoff provided by the user class
      real_type RmaxFixed;
      ///temporary data
      real_type R12;
      ///temporary data for derivative
      real_type dCosFac;
      ///ID for variable E 
      string ID_E;

      ///default constructor not to have nan
      CuspCorrectionFunctor():E(1.0),Rmax(1),Rcut(0), RmaxFixed(10){}
      /** constructor
       * @param c Cusp condition, i.e., the first derivative at r=0
       * @param e exponent or scaling factor
       */
      inline CuspCorrectionFunctor(real_type e, real_type rc)
      {
        E=e; 
        RmaxFixed=rc;
        resetInternals();
      }

      inline void resetInternals()
      {
        mOneOverE=-1.0/E;
        Rmax=10.0/E;
        Rmax=(Rmax>RmaxFixed)? RmaxFixed:Rmax;
        Rcut=0.9*Rmax;
        R12=1.0/(Rmax-Rcut);
        dCosFac=-0.5*M_PI*R12;
      }

      inline real_type f(real_type r)
      {
        if(r>Rmax) return 0.0;
        real_type v=mOneOverE*std::exp(-E*r);
        if(r>=Rcut)
          return v*0.5*(1.0+std::cos(M_PI*(r-Rcut)*R12));
        else
          return v;
      }

      inline real_type df(real_type r) {
        if(r>Rmax) return 0.0;
        if(r<Rcut)
        {
          return std::exp(-E*r);
        }
        else
        {
          //this may be wrong but should never be used.
          real_type rfac=M_PI*(r-Rcut)*R12;
          return std::exp(-E*r)*(mOneOverE*dCosFac*std::sin(rfac)+0.5*(1.0+std::cos(rfac)));
        }
      }

      bool put(xmlNodePtr cur) 
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(ID_E,"id"); 
        rAttrib.add(E,"exponent"); 
        rAttrib.put(cur);
        ID_E.append("_E");
        return true;
      }

      void addOptimizables(OptimizableSetType& vlist) 
      {
        vlist[ID_E]=E;
      }

      void resetParameters(OptimizableSetType& optVariables) 
      { 
        typename OptimizableSetType::iterator it_b(optVariables.find(ID_E));
        if(it_b != optVariables.end()) {
          E=(*it_b).second;
          //only change the exponent: leave the cutoff fixed
          mOneOverE=-1.0/E;
          //resetInternals();
        }
      }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: WMConstraints.h 1691 2007-02-01 21:51:50Z jnkim $ 
 ***************************************************************************/

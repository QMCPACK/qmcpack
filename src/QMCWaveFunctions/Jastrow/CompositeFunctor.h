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

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: WMConstraints.h 1691 2007-02-01 21:51:50Z jnkim $ 
 ***************************************************************************/

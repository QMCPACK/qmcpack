//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_COMPOSITEFUNCTORS_H
#define QMCPLUSPLUS_COMPOSITEFUNCTORS_H

#include <Numerics/OptimizableFunctorBase.h>
#include <config/stdlib/math.h>

namespace qmcplusplus
{

/** Implements a truncated Pade Function \f$u[r]=A*r/(1+B*r)\f$
 *
 * Only valid for two-body jastrow with a "known" cusp condition.
 */
template<class T>
struct TruncatedPadeFunctor:public OptimizableFunctorBase
{

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
  OptimizableFunctorBase* inFunc;

  /** constructor
   * @param a default a for e-e cusp condition
   * @param fin input function
   */
  TruncatedPadeFunctor(real_type a=-0.5, OptimizableFunctorBase* fin=0, real_type rmax=10):
    inFunc(fin),Rmax(rmax)
  {
    Rcut=0.0;
    A=a;
  }

  inline real_type f(real_type r)
  {
    if (r>Rcut)
      return inFunc->f(r);
    else
      return A*r/(1.0+B*r)+OffSet;
  }

  inline real_type df(real_type r)
  {
    if(r>Rcut)
      return inFunc->df(r);
    else
    {
      real_type u = 1.0/(1.0+B*r);
      return A*u*u;
    }
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void resetParameters(const opt_variables_type& active)
  {
    if(inFunc)
      inFunc->resetParameters(active);
    applyCuspCondition();
  }

  void applyCuspCondition()
  {
    Rcut=0.0;
    if(inFunc)
    {
      real_type x=0.001;
      real_type dx=0.001;
      real_type deriv0=inFunc->df(x),deriv;
      bool validInput=true;
      if(deriv0*A<0.0)
      {
        return;
      }
      if(A<0.0)
      {
        bool goingdown=true;
        do
        {
          x+=dx;
          if(x>=Rmax)
            validInput=false;
          deriv=inFunc->df(x);
          if(deriv>deriv0)
            goingdown=false;
          else
            deriv0=deriv;
        }
        while(goingdown && validInput);
      }
      else
      {
        bool goingup=true;
        do
        {
          x+=dx;
          if(x>Rmax)
            validInput=false;
          deriv=inFunc->df(x);
          if(deriv<deriv0)
            goingup=false;
          else
            deriv0=deriv;
        }
        while(goingup && validInput);
      }
      if(validInput)
      {
        Rcut=x+2*dx;
        deriv=inFunc->df(Rcut);
        B = (std::sqrt(A/deriv)-1.0)/Rcut;
        OffSet=inFunc->f(Rcut)-A*Rcut/(1.0+B*Rcut);
        AB = A*B;
        B2=2.0*B;
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
struct AnyTimesRnFunctor: public OptimizableFunctorBase
{
  ///pointer to a functor
  OptimizableFunctorBase* myFunc;
  ///power
  int Np;
  ///constructor
  AnyTimesRnFunctor(OptimizableFunctorBase* infunc=0, int n=1):
    myFunc(infunc),Np(n)
  {
  }

  OptimizableFunctorBase* makeClone() const
  {
    AnyTimesRnFunctor<T>* myclone=new AnyTimesRnFunctor<T>(*this);
    myclone->myFunc=myFunc->makeClone();
    return myclone;
  }

  inline void reset()
  {
    myFunc->reset();
  }

  inline real_type f(real_type r)
  {
    return std::pow(r,Np)*myFunc->f(r);
  }

  inline real_type df(real_type r)
  {
    real_type u=myFunc->f(r);
    real_type du=myFunc->df(r);
    return (Np*u+r*du)*std::pow(r,Np-1);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void resetParameters(const opt_variables_type& active)
  {
    if(myFunc)
      myFunc->resetParameters(active);
  }

};

/** Implements \f$ u(r) = A*(-\frac{1}{B} exp(-Br))*f_c(r)\f$ with a cutoff function.
 *
 * This functor is not optimized and should be used only as a temporary functor
 * for a final numerical functor to correct the cusp condition.
 */
template<class T>
struct CuspCorrectionFunctor: public OptimizableFunctorBase
{

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
  ///temporary data for second-derivative
  real_type d2CosFac;
  ///default constructor not to have nan
  CuspCorrectionFunctor():E(1.0),Rmax(1),Rcut(0), RmaxFixed(10) {}

  OptimizableFunctorBase* makeClone() const
  {
    return new CuspCorrectionFunctor<T>(*this);
  }

  /** constructor
   * @param c Cusp condition, i.e., the first derivative at r=0
   * @param e exponent or scaling factor
   */
  inline CuspCorrectionFunctor(real_type e, real_type rc)
  {
    real_type rin=0.5*rc;
    if(e>0)
      E=(rin)<e?4.0/rin:e;
    else
      E=5.0/rc;
    Rmax=rc;
    //E=e;
    reset();
  }

  inline void reset()
  {
    mOneOverE=-1.0/E;
    Rcut=0.7*Rmax;
    R12=M_PI/(Rmax-Rcut);
    dCosFac=0.5*R12/E;
    d2CosFac=0.5*R12*R12/E;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    if(r>Rmax)
      return 0.0;
    if(r<Rcut)
    {
      return std::exp(-E*r);
    }
    else
    {
      //this may be wrong but should never be used.
      real_type sine,cosine;
      sincos((r-Rcut)*R12,&sine,&cosine);
      return std::exp(-E*r)*(dCosFac*sine+0.5*(1.0+cosine));
    }
  }

  inline real_type evaluate(real_type r)
  {
    if(r>Rmax)
      return 0.0;
    if(r<Rcut)
      return mOneOverE*std::exp(-E*r);
    else
      return mOneOverE*std::exp(-E*r)*0.5*(1.0+std::cos((r-Rcut)*R12));
  }

  inline real_type evaluate(real_type r, real_type& du, real_type& d2u)
  {
    if(r>Rmax)
    {
      du=0.0;
      d2u=0.0;
      return 0.0;
    }
    real_type u=std::exp(-E*r);
    if(r<Rcut)
    {
      du=u;
      d2u=-E*u;
      return mOneOverE*u;
    }
    else
    {
      real_type sine,cosine;
      sincos((r-Rcut)*R12,&sine,&cosine);
      real_type fc=0.5*(1+cosine);
      du=u*(dCosFac*sine+fc);
      d2u=u*(d2CosFac*cosine+2.0*dCosFac*sine-E*fc);
      return u*fc;
    }
  }

  bool put(xmlNodePtr cur)
  {
    ///ID for variable E
    std::string ID_E;
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID_E,"id");
    rAttrib.add(E,"exponent");
    rAttrib.put(cur);
    ID_E.append("_E");
    myVars.insert(ID_E,E);
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  //cannot be optimized
  void resetParameters(const opt_variables_type& active)
  {
    int loc=myVars.where(0);
    if(loc>=0)
    {
      E=active[loc];
      reset();
    }
  }
};

/** Implements \f$ u(r) = A*(-\frac{1}{B} exp(-Br))*f_c(r)\f$ with a cutoff function.
 *
 * This functor is not optimized and should be used only as a temporary functor
 * for a final numerical functor to correct the cusp condition.
 */
template<class T>
struct DCuspCorrectionDEFunctor: public OptimizableFunctorBase
{

  ///scaling function or exponent
  real_type E;
  ///1/E
  real_type OneOverE;
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
  ///temporary data for derivative
  real_type d2CosFac;

  ///default constructor not to have nan
  DCuspCorrectionDEFunctor():E(1.0),Rmax(1),Rcut(0), RmaxFixed(10) {}
  /** constructor
   * @param c Cusp condition, i.e., the first derivative at r=0
   * @param e exponent or scaling factor
   */
  inline DCuspCorrectionDEFunctor(real_type e, real_type rc)
  {
    E=e;
    Rmax=rc;
    reset();
  }

  inline void reset()
  {
    OneOverE=1.0/E;
    Rcut=0.7*Rmax;
    R12=M_PI/(Rmax-Rcut);
    dCosFac=0.5*R12;
    d2CosFac=0.5*R12*R12;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    if(r>Rmax)
      return 0.0;
    real_type expr=std::exp(-E*r);
    if(r<Rcut)
    {
      return -r*expr;
    }
    else
    {
      real_type sine,cosine;
      sincos((r-Rcut)*R12,&sine,&cosine);
      return expr*(dCosFac*sine*OneOverE*(r+OneOverE)-r*0.5*(1+cosine));
    }
  }

  inline real_type evaluate(real_type r)
  {
    if(r>Rmax)
      return 0.0;
    real_type u=std::exp(-E*r)*OneOverE*(r+OneOverE);
    if(r<Rcut)
      return u;
    else
      return u*0.5*(1.0+std::cos((r-Rcut)*R12));
  }

  inline real_type evaluate(real_type r, real_type& du, real_type& d2u)
  {
    if(r>Rmax)
    {
      du=0.0;
      d2u=0.0;
      return 0.0;
    }
    real_type expr=std::exp(-E*r);
    real_type u=expr*OneOverE*(r+OneOverE);
    if(r<Rcut)
    {
      du=-r*expr;
      d2u=expr*(E*r-1.0);
      return u;
    }
    else
    {
      real_type sine,cosine;
      sincos((r-Rcut)*R12,&sine,&cosine);
      real_type fc=0.5*(1+cosine);
      du=u*dCosFac*sine-r*expr*fc;
      d2u=expr*((E*r-1.0)+2.0*r*sine*dCosFac)+d2CosFac*sine*u;
      return u*fc;
    }
  }

  bool put(xmlNodePtr cur)
  {
    ///ID for variable E
    std::string ID_E;
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID_E,"id");
    rAttrib.add(E,"exponent");
    rAttrib.put(cur);
    ID_E.append("_E");
    myVars.insert(ID_E,E);
    return true;
  }

  void resetParameters(const opt_variables_type& active)
  {
    int loc=myVars.where(0);
    if(loc>=0)
    {
      E=active[loc];
      reset();
    }
  }
};

}
#endif

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file LRJastrowSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LONGRANGEJASTROW_BREAKUPUTILITY_H
#define QMCPLUSPLUS_LONGRANGEJASTROW_BREAKUPUTILITY_H

#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
/** RPABreakUp
 *
 * A Func for LRHandlerTemp.  Four member functions have to be provided
 *
 * - reset(T volume) : reset the normalization factor
 * - operator() (T r, T rinv) : return a value of the original function e.g., 1.0/r
 * - Fk(T k, T rc)
 * - Xk(T k, T rc)
 *
 */
template<class T=double>
struct YukawaBreakup
{
  T Rs;
  T SqrtRs;
  T OneOverSqrtRs;
  T NormFactor;
  inline YukawaBreakup() {}

  void reset(ParticleSet& ref)
  {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    T Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    //NormFactor*=(rs*rs*rs)/(Rs*Rs*Rs);
    Rs = rs;
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
  }


  inline T operator()(T r, T rinv)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return SqrtRs-0.5*r;
    else
      return Rs*rinv*(1.0-std::exp(-r*OneOverSqrtRs));
    //if (r > 1e-10) return Rs*rinv*(1.0 - std::exp(-r*OneOverSqrtRs));
    //return 1.0 / OneOverSqrtRs - 0.5 * r;
  }

  inline T df(T r)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return -0.5+r*OneOverSqrtRs/3.0;
    else
    {
      T rinv=1.0/r;
      T exponential = std::exp(-r*OneOverSqrtRs);
      return -Rs*rinv*rinv*(1.0 - exponential) + exponential*rinv*SqrtRs;
    }
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T coskr = std::cos(k*rc);
    T sinkr = std::sin(k*rc);
    T oneOverK = 1.0/k;
    return -NormFactor * Rs *
           (coskr*oneOverK*oneOverK
            - std::exp(-rc*OneOverSqrtRs)*(coskr - OneOverSqrtRs * sinkr * oneOverK)/(k*k+1.0/Rs));
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  /** return RPA value at |k|
    * @param kk |k|^2
    */
  inline T Uk(T kk)
  {
    return NormFactor*Rs/kk;
  }

  /** return d u(k)/d rs
    *
    * Implement a correct one
    */
  inline T derivUk(T kk)
  {
    return NormFactor/kk;
  }
};

template<class T=double>
struct DerivRPABreakup
{
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline DerivRPABreakup() {}

  void reset(ParticleSet& ref)
  {
//       NormFactor= 4.0*M_PI/ref.Lattice.Volume;
//       NormFactor=4.0*M_PI/ref.getTotalNum();
    NormFactor=1.0/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = rs;
    //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv)
  {
    return 0.0;
  }

  inline T df(T r)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0)
    {
      Sy=1.0;
    }
    else
    {
      Sy = 1.5*y - 0.5*y*y*y;
    };
    //multiply by NORM?
    return  NormFactor* 3.0/( (k*k*k*k*Rs*Rs*Rs*Rs*(std::pow(1.0/(Sy*Sy)+12.0/(k*k*k*k*Rs*Rs*Rs),0.5))) );
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }
  /** return RPA value at |k|
  * @param kk |k|^2
  */
  inline T Uk(T kk)
  {
    return NormFactor*Rs/kk;
  }

  /** return d u(k)/d rs
  *
  * Implement a correct one
  */
  inline T derivUk(T kk)
  {
    return 0.0;
  }
};


template<class T=double>
struct RPABreakup
{
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline RPABreakup() {}

  void reset(ParticleSet& ref)
  {
    //       NormFactor= 4.0*M_PI/ref.Lattice.Volume;
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    NormFactor=1.0/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = rs;
    //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv)
  {
    return 0.0;
  }

  inline T df(T r)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0)
    {
      Sy=1.0;
    }
    else
    {
      Sy = 1.5*y - 0.5*y*y*y;
    };
    //multiply by NORM?
    return  NormFactor*(0.5*(-1.0/Sy + std::pow(1.0/(Sy*Sy) + 12.0/(k*k*k*k*Rs*Rs*Rs),0.5)));
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  /** return RPA value at |k|
  * @param kk |k|^2
  */
  inline T Uk(T kk)
  {
    return NormFactor*Rs/kk;
  }

  /** return d u(k)/d rs
  *
  * Implement a correct one
    */
  inline T derivUk(T kk)
  {
    return 0.0;
  }
};


template<class T=double>
struct DerivYukawaBreakup
{
  T Rc;
  T Rs;
  T SqrtRs;
  T OneOverSqrtRs;
  T NormFactor;
  T OneOverSqrtRs3;
  T OneOverRs;
  T DerivSecondTaylorTerm;
  T n2;

  inline DerivYukawaBreakup() {}

  void reset(ParticleSet& ref)
  {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    T Density=ref.getTotalNum()/ref.Lattice.Volume;
    n2 = ref.getTotalNum();
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
    OneOverSqrtRs3 = std::pow(OneOverSqrtRs,3.0);
    OneOverRs = 1.0 / Rs;
    DerivSecondTaylorTerm= (1.0/12.0)*OneOverSqrtRs3;
    Rc = ref.Lattice.LR_rc;
  }

  void reset(ParticleSet& ref, T rs)
  {
    Rs = rs;
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    T Density=ref.getTotalNum()/ref.Lattice.Volume;
    n2 = ref.getTotalNum();
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
    OneOverSqrtRs3 = std::pow(OneOverSqrtRs,3.0);
    OneOverRs = 1.0 / Rs;
    DerivSecondTaylorTerm= (1.0/12.0)*OneOverSqrtRs3;
  }

  /** need the df(r)/d(rs) */
  inline T operator()(T r, T rinv)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return 0.5*OneOverSqrtRs*(1.0-r*OneOverSqrtRs);
    else
    {
      return 0.5*OneOverSqrtRs*std::exp(-r*OneOverSqrtRs);
    }
  }

  /** need d(df(r)/d(rs))/dr */
  inline T df(T r)
  {
    if(r< std::numeric_limits<T>::epsilon())
      return -0.5*OneOverRs*(1.0-r*OneOverSqrtRs);
    else
    {
      T rinv=1.0/r;
      T exponential = std::exp(-r*OneOverSqrtRs);
      return -0.5*OneOverRs*std::exp(-r*OneOverSqrtRs);
    }
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  /** integral from kc to infinity */
  inline T Xk(T k, T rc )
  {
    T Rc=rc;
    T coskr = std::cos(k*Rc);
    T sinkr = std::sin(k*Rc);
    T Rs2 = Rs*Rs;
    T oneOverK = 1.0/k;
    T Val = (-0.5*NormFactor*OneOverSqrtRs*oneOverK)*(1.0/((1.0 + k*k* Rs)*(1.0 + k*k* Rs))*std::exp(-rc*OneOverSqrtRs)* SqrtRs* (k* SqrtRs* (rc + 2.0* SqrtRs + k*k* rc* Rs)* coskr + (rc + k*k*rc*Rs + SqrtRs* (1.0 - k*k*Rs))* sinkr));
    return Val;
  }
};

template<class T=double>
struct EPRPABreakup
{
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline EPRPABreakup() {}

  void reset(ParticleSet& ref)
  {
    NormFactor=1.0/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = rs;
    //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv)
  {
    return 0.0;
  }

  inline T df(T r)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0)
    {
      Sy=1.0;
    }
    else
    {
      Sy = 1.5*y - 0.5*y*y*y;
    };
    T val = 12.0/(k*k*k*k*Rs*Rs*Rs);
    return  -0.5*NormFactor*val*std::pow(1.0/(Sy*Sy)+val,-0.5);
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  /** return RPA value at |k|
  * @param kk |k|^2
   */
  inline T Uk(T kk)
  {
    return NormFactor*Rs/kk;
  }

  /** return d u(k)/d rs
  *
  * Implement a correct one
   */
  inline T derivUk(T kk)
  {
    return 0.0;
  }
};
template<class T=double>
struct derivEPRPABreakup
{
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline derivEPRPABreakup() {}

  void reset(ParticleSet& ref)
  {
    NormFactor=1.0/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs)
  {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = rs;
    //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv)
  {
    return 0.0;
  }

  inline T df(T r)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0)
    {
      Sy=1.0;
    }
    else
    {
      Sy = 1.5*y - 0.5*y*y*y;
    };
    T val = 12.0/(k*k*k*k*Rs*Rs*Rs);
    T uk = val*std::pow(1.0/(Sy*Sy)+val,-0.5);
    return  -0.5*NormFactor*(uk/Rs)*(1.0-0.5*val/(1.0/(Sy*Sy)+val) );
  }
  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  /** return RPA value at |k|
  * @param kk |k|^2
   */
  inline T Uk(T kk)
  {
    return NormFactor*Rs/kk;
  }

  /** return d u(k)/d rs
         *
         * Implement a correct one
   */
  inline T derivUk(T kk)
  {
    return 0.0;
  }
};

template<typename T>
struct ShortRangePartAdapter : OptimizableFunctorBase
{
  typedef LRHandlerBase HandlerType;

  explicit ShortRangePartAdapter(HandlerType* inhandler): Uconst(0), myHandler(inhandler)
  { }

  OptimizableFunctorBase* makeClone() const
  {
    return new ShortRangePartAdapter<T>(*this);
  }

  inline void reset() {}
  inline void setRmax(real_type rm)
  {
    Uconst=myHandler->evaluate(rm,1.0/rm);
  }
  inline real_type evaluate(real_type r)
  {
    return f(r);
  }
  inline real_type f(real_type r)
  {
    return myHandler->evaluate(r, 1.0/r)-Uconst;
  }
  inline real_type df(real_type r)
  {
    return myHandler->srDf(r, 1.0/r);
  }
  void checkInVariables(opt_variables_type& active) { }
  void checkOutVariables(const opt_variables_type& active) { }
  void resetParameters(const opt_variables_type& optVariables) { }
  bool put(xmlNodePtr cur)
  {
    return true;
  }
  real_type Uconst;
  HandlerType* myHandler;
};


template<class T=double>
struct RPABFeeBreakup
{
  T Rs;
  T kf;
  T kfm[2];
  T Density;
  T volume;
  T hbs2m;
  int nelec;
  int nspin;
  int nppss[2];
  inline RPABFeeBreakup() {}

  // assumes 3D here, fix
  void reset(ParticleSet& ref)
  {
    volume = ref.Lattice.Volume;
    nspin = ref.groups();
    for(int i=0; i<nspin; ++i)
      nppss[i] = ref.last(i) - ref.first(i);
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    nelec = ref.getTotalNum();
    hbs2m=0.5;
    kf = 0.0;
    kfm[0]=kfm[1]=0.0;
    for(int i=0; i<nspin; i++)
    {
      T a3=3.0*volume/(4.0*M_PI*nppss[i]);
      kfm[i]=std::pow(9.0*M_PI/(2.0*a3),1.0/OHMMS_DIM);
      kf += kfm[i]*nppss[i]/ref.getTotalNum();
    }
  }

  // assumes 3D here, fix
  void reset(ParticleSet& ref, T rs)
  {
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    volume = ref.Lattice.Volume;
    nspin = ref.groups();
    for(int i=0; i<nspin; ++i)
      nppss[i] = ref.last(i) - ref.first(i);
    Rs = rs;
    nelec = ref.getTotalNum();
    hbs2m=0.5;
    kf = 0.0;
    kfm[0]=kfm[1]=0.0;
    for(int i=0; i<nspin; i++)
    {
      T a3=3.0*volume/(4.0*M_PI*nppss[i]);
      kfm[i]=std::pow(9.0*M_PI/(2.0*a3),1.0/OHMMS_DIM);
      kf += kfm[i]*nppss[i]/ref.getTotalNum();
    }
  }


  inline T operator()(T r, T rinv)
  {
    return 0.0;
  }

  inline T df(T r)
  {
    return 0.0;
  }

  inline T Fk(T k, T rc)
  {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc)
  {
    T u = urpa(k)*volume/nelec;
    T eq=k*k*hbs2m;
    T vlr=u*u*2.0*Density*eq;
    T vq=(2.0*M_PI*(OHMMS_DIM-1.0))/(std::pow(k,OHMMS_DIM-1.0));
    T veff=vq-vlr;
    T op2=Density*2.0*vq*eq;
#if OHMMS_DIM==3
    T op2kf=Density*2.0*(2.0*M_PI*(OHMMS_DIM-1.0))*hbs2m;
    T op=(op2+1.2*( kfm[0]*kfm[0] + kfm[1]*kfm[1] )*eq*hbs2m + eq*eq );
#elif OHMMS_DIM==2
    T op2kf=Density*2.0*(2.0*M_PI*(OHMMS_DIM-1.0))*std::pow(kf,3.0-OHMMS_DIM)*hbs2m;
    T op=(op2+1.5*( kfm[0]*kfm[0] + kfm[1]*kfm[1] )*hbs2m*eq + eq*eq );
#endif
    op=std::sqrt(op);
    T denom=1.0/veff - Dlindhard(k,-eq);
    T dp=-Dlindhardp(k,-eq);
    T s0 = 0.0,ss=0.0;
    for(int i=0; i<nspin; i++)
    {
      if(nppss[i] > 0)
      {
        T y=0.5*k/kfm[i];
        if(y < 1.0)
        {
#if OHMMS_DIM==3
          ss=0.5*y*(3.0-y*y);
#elif OHMMS_DIM==2
          ss=(2.0/M_PI)*(std::asin(y)+y*std::sqrt(1.0-y*y));
#endif
        }
        else
        {
          ss=1.0;
        }
        s0 += nppss[i]*ss/nelec;
      }
    }
    T yq = s0*s0/(4.0*k*k)*(1.0/(eq*denom)+dp/(denom*denom))
           + 2.0*hbs2m*vlr/(op*(std::sqrt(op2)+eq));
    return yq/volume;
  }

  inline T integrate_r2(T rc)
  {
    return 0.0;
  }

  /** return RPA value at |k|
  * @param kk |k|^2
  */
  inline T Uk(T kk)
  {
    return 0.0;
  }

  /** return d u(k)/d rs
  *
  * Implement a correct one
    */
  inline T derivUk(T kk)
  {
    return 0.0;
  }

  T urpa(T q)
  {
    T a=0.0,vkbare;
    if(q > 0.0)
    {
      vkbare=4.0*M_PI/(q*q);
      a=2.0*Density*vkbare/(hbs2m*q*q);
    }
    return -0.5+0.5*std::sqrt(1.0+a);
  }

  // mmorales: taken from bopimc (originally by M Holzmann)
  T Dlindhard(T q, T w)
  {
    T xd1,xd2,small=0.00000001,xdummy,rq1,rq2,sd1,sd2;
    T res=0.0;
    for(int i=0; i<nspin; i++)
    {
      if(kfm[i] > 0.0)
      {
        T xq=q/kfm[i];
        if(xq < small)
          xq=small;
        T om=w/(kfm[i]*kfm[i]*hbs2m*2.0);
        xd1=om/xq-xq/2.0;
        xd2=om/xq+xq/2.0;
        if(std::abs(xd1-1.0) <= small)
        {
          if(xd1 >= 1.0)
            xd1=1.0+small;
          else
            xd1=1.0-small;
        }
        else
          if(std::abs(xd2-1.0) <= small)
          {
            if(xd2 >= 1.0)
              xd2=1.0+small;
            else
              xd2=1.0-small;
          }
#if OHMMS_DIM==3
        rq1 = std::log(std::abs(1.0+xd1)/std::abs(1.0-xd1));
        rq2 = std::log(std::abs(1.0+xd2)/std::abs(1.0-xd2));
        xdummy = -1.0+0.5*(1.0-xd1*xd1)/xq*rq1-0.50*(1.0-xd2*xd2)/xq*rq2;
        xdummy *= kfm[i]/(8.0*M_PI*M_PI*hbs2m);
#elif OHMMS_DIM==2
        if(std::abs(xd1) < small)
          sd1=0.0;
        else
          sd1=xd1/std::abs(xd1);
        if(std::abs(xd2) < small)
          sd2=0.0;
        else
          sd2=xd2/std::abs(xd2);
        if(xd1*xd1-1.0 > 1.0)
          rq1=std::sqrt(xd1*xd1-1.0)*sd1;
        else
          rq1=0.0;
        if(xd2*xd2-1.0 > 1.0)
          rq2=std::sqrt(xd2*xd2-1.0)*sd2;
        else
          rq2=0.0;
        xdummy=-(xq+rq1-rq2)/(4.0*M_PI*xq*hbs2m);
#endif
        res += xdummy;
      }
    }
    return res;
  }

  // mmorales: taken from bopimc (originally by M Holzmann)
  T Dlindhardp(T q, T w)
  {
    T xd1,xd2,small=0.00000001,xdummy,rq1,rq2,sd1,sd2;
    T res=0.0;
    for(int i=0; i<nspin; i++)
    {
      if(kfm[i] > 0.0)
      {
        T xq=q/kfm[i];
        if(xq < small)
          xq=small;
        T om=w/(kfm[i]*kfm[i]*hbs2m*2.0);
        xd1=om/xq-xq/2.0;
        xd2=om/xq+xq/2.0;
        if(std::abs(xd1-1.0) <= small)
        {
          if(xd1 >= 1.0)
            xd1=1.0+small;
          else
            xd1=1.0-small;
        }
        else
          if(std::abs(xd2-1.0) <= small)
          {
            if(xd2 >= 1.0)
              xd2=1.0+small;
            else
              xd2=1.0-small;
          }
#if OHMMS_DIM==3
        rq1 = std::log(std::abs(1.0+xd1)/std::abs(1.0-xd1));
        rq2 = std::log(std::abs(1.0+xd2)/std::abs(1.0-xd2));
        xdummy=-xd1/(xq*kfm[i]*kfm[i]*2.0*hbs2m)*rq1+ xd2/(xq*kfm[i]*kfm[i]*2.0*hbs2m)*rq2;
        xdummy *= kfm[i]/(8.0*M_PI*M_PI*hbs2m*xq);
#elif OHMMS_DIM==2
        if(std::abs(xd1) < small)
          sd1=0.0;
        else
          sd1=xd1/std::abs(xd1);
        if(std::abs(xd2) < small)
          sd2=0.0;
        else
          sd2=xd2/std::abs(xd2);
        if(xd1*xd1-1.0 > 1.0)
          rq1=sd1*xd1/std::sqrt(xd1*xd1-1.0);
        else
          rq1=0.0;
        if(xd2*xd2-1.0 > 1.0)
          rq2=sd2*xd2/std::sqrt(xd2*xd2-1.0);
        else
          rq2=0.0;
        xdummy=-(rq1-rq2)/(8.0*M_PI*(xq*kfm[i]*hbs2m)*(xq*kfm[i]*hbs2m));
#endif
        res += xdummy;
      }
    }
    return res;
  }

};

}
#endif

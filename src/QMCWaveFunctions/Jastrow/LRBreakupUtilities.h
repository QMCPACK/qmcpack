//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnum Kim
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
/** @file LRJastrowSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LONGRANGEJASTROW_BREAKUPUTILITY_H
#define QMCPLUSPLUS_LONGRANGEJASTROW_BREAKUPUTILITY_H

#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus {
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
    struct YukawaBreakup {
  T Rs;
  T SqrtRs;
  T OneOverSqrtRs;
  T NormFactor;
  inline YukawaBreakup(){}
  
  void reset(ParticleSet& ref) {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    T Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
  }
  
  void reset(ParticleSet& ref, T rs) {
    NormFactor=4.0*M_PI/ref.Lattice.Volume;
    //NormFactor*=(rs*rs*rs)/(Rs*Rs*Rs);
    Rs = rs;
    SqrtRs=std::sqrt(Rs);
    OneOverSqrtRs = 1.0 / SqrtRs;
  }

  
  inline T operator()(T r, T rinv) { 
    if(r< numeric_limits<T>::epsilon())
      return SqrtRs-0.5*r;
    else
      return Rs*rinv*(1.0-std::exp(-r*OneOverSqrtRs));
    //if (r > 1e-10) return Rs*rinv*(1.0 - std::exp(-r*OneOverSqrtRs));
    //return 1.0 / OneOverSqrtRs - 0.5 * r;
  }
  
  inline T df(T r) {
    if(r< numeric_limits<T>::epsilon())
      return -0.5+r*OneOverSqrtRs/3.0;
    else
    {
      T rinv=1.0/r;
      T exponential = std::exp(-r*OneOverSqrtRs);
      return -Rs*rinv*rinv*(1.0 - exponential) + exponential*rinv*SqrtRs;
    }
  }
  
  inline T Fk(T k, T rc) {
    return -Xk(k,rc);
  }
  
  inline T Xk(T k, T rc) {
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
    struct DerivRPABreakup {
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline DerivRPABreakup(){}

  void reset(ParticleSet& ref) {
//       NormFactor= 4.0*M_PI/ref.Lattice.Volume;
//       NormFactor=4.0*M_PI/ref.getTotalNum();
  NormFactor=1.0/ref.getTotalNum();

  Density=ref.getTotalNum()/ref.Lattice.Volume;
  Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
  //unpolarized K_f
  Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
}

void reset(ParticleSet& ref, T rs) {
  //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
  NormFactor=1.0/ref.getTotalNum();
  //       NormFactor=4.0*M_PI/ref.getTotalNum();
  Density=ref.getTotalNum()/ref.Lattice.Volume;

  Rs = rs;
  //unpolarized
  Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }
  
  
  inline T operator()(T r, T rinv) {
    return 0.0;
  }
  
  inline T df(T r) {
    return 0.0;
  }
  
  inline T Fk(T k, T rc) {
    return -Xk(k,rc);
  }
  
  inline T Xk(T k, T rc) {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0) {
      Sy=1.0;
    }
    else {
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
    struct RPABreakup {
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline RPABreakup(){}

  void reset(ParticleSet& ref) {
        //       NormFactor= 4.0*M_PI/ref.Lattice.Volume;
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
    NormFactor=1.0/ref.getTotalNum();

    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs) {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
      //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;

    Rs = rs;
    //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv) {
    return 0.0;
  }

  inline T df(T r) {
    return 0.0;
  }

  inline T Fk(T k, T rc) {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc) {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0) {
      Sy=1.0;
    }
    else {
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
      struct DerivYukawaBreakup {
    T Rc;
    T Rs;
    T SqrtRs;
    T OneOverSqrtRs;
    T NormFactor;
    T OneOverSqrtRs3;
    T OneOverRs;
    T DerivSecondTaylorTerm;
    T n2;

    inline DerivYukawaBreakup(){}

    void reset(ParticleSet& ref) {

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

    void reset(ParticleSet& ref, T rs) {
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
    inline T operator()(T r, T rinv) { 
      if(r< numeric_limits<T>::epsilon())
        return 0.5*OneOverSqrtRs*(1.0-r*OneOverSqrtRs);
      else{
        return 0.5*OneOverSqrtRs*std::exp(-r*OneOverSqrtRs);
      }
    }

    /** need d(df(r)/d(rs))/dr */
    inline T df(T r) {
      if(r< numeric_limits<T>::epsilon())
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

    inline T Fk(T k, T rc) {
      return -Xk(k,rc);
    }

    /** integral from kc to infinity */
    inline T Xk(T k, T rc ) {

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
    struct EPRPABreakup {
  T Rs;
  T Kf;
  T Density;
  T NormFactor;
  inline EPRPABreakup(){}

  void reset(ParticleSet& ref) {
    NormFactor=1.0/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;
    Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
        //unpolarized K_f
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }

  void reset(ParticleSet& ref, T rs) {
          //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
    NormFactor=1.0/ref.getTotalNum();
        //       NormFactor=4.0*M_PI/ref.getTotalNum();
    Density=ref.getTotalNum()/ref.Lattice.Volume;

    Rs = rs;
          //unpolarized
    Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
  }


  inline T operator()(T r, T rinv) {
    return 0.0;
  }

  inline T df(T r) {
    return 0.0;
  }

  inline T Fk(T k, T rc) {
    return -Xk(k,rc);
  }

  inline T Xk(T k, T rc) {
    T y = 0.5*k/Kf;
    T Sy;
    if (y >= 1.0) {
      Sy=1.0;
    }
    else {
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
        struct derivEPRPABreakup {
      T Rs;
      T Kf;
      T Density;
      T NormFactor;
      inline derivEPRPABreakup(){}

      void reset(ParticleSet& ref) {
        NormFactor=1.0/ref.getTotalNum();
        Density=ref.getTotalNum()/ref.Lattice.Volume;
        Rs = std::pow(3.0/(4.0*M_PI*Density), 1.0/3.0);
    //unpolarized K_f
        Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
      }

      void reset(ParticleSet& ref, T rs) {
    //       NormFactor=4.0*M_PI/ref.Lattice.Volume;
        NormFactor=1.0/ref.getTotalNum();
    //       NormFactor=4.0*M_PI/ref.getTotalNum();
        Density=ref.getTotalNum()/ref.Lattice.Volume;

        Rs = rs;
      //unpolarized
        Kf = std::pow(2.25*M_PI, 1.0/3.0)/Rs;
      }


      inline T operator()(T r, T rinv) {
        return 0.0;
      }

      inline T df(T r) {
        return 0.0;
      }

      inline T Fk(T k, T rc) {
        return -Xk(k,rc);
      }

      inline T Xk(T k, T rc) {
        T y = 0.5*k/Kf;
        T Sy;
        if (y >= 1.0) {
          Sy=1.0;
        }
        else {
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
struct ShortRangePartAdapter : OptimizableFunctorBase {
  typedef LRHandlerBase HandlerType;

  explicit ShortRangePartAdapter(HandlerType* inhandler): Uconst(0), myHandler(inhandler) 
  { }

  OptimizableFunctorBase* makeClone() const 
  {
    return new ShortRangePartAdapter<T>(*this);
  }

  inline void reset() {}
  inline void setRmax(real_type rm) { Uconst=myHandler->evaluate(rm,1.0/rm);}
  inline real_type evaluate(real_type r) { return f(r); }
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
  bool put(xmlNodePtr cur) {return true;}
  real_type Uconst;
  HandlerType* myHandler;
};


}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2027 $   $Date: 2007-05-31 15:37:59 -0500 (Thu, 31 May 2007) $
 * $Id: LRJastrowSingleton.h 2027 2007-05-31 20:37:59Z jnkim $
 ***************************************************************************/

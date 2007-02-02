//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_WMFUNCTOR_WITHPADECORRECTION_H
#define QMCPLUSPLUS_WMFUNCTOR_WITHPADECORRECTION_H
#include "QMCWaveFunctions/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  /** Implements a screened Function \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$
   * 
   * Short-range functor introduced by Wagner and Mitas, cond-mat/0610088
   */
  template<class T>
    struct WMFunctor: public OptimizableFunctorBase<T> {
      typedef T real_type;
      ///input B
      real_type B0;
      ///input Rcut
      real_type Rcut;
      ///1/Rcut
      real_type OneOverRc;
      ///id
      string ID;
      ///constructor
      explicit WMFunctor(real_type b, real_type rc=7.5) {
        reset(b,rc);
      }
      inline void reset() { OneOverRc=1.0/Rcut; }
      void reset(real_type b, real_type rc) { B0=b; Rcut=rc; reset(); }
      inline real_type f(real_type r) {
        real_type x=r*OneOverRc;
        real_type z=x*x*(6.0-8*x+3.0*x*x);
        return (1-z)/(1+B0*z);
      }
      inline real_type df(real_type r) {
        real_type x=r*OneOverRc;
        real_type z=x*x*(6.0-8*x+3.0*x*x);
        return -(1+B0)/(1+B0*z)/(1+B0*z)*OneOverRc*12*x*(1-2.0*x+x*x);
      }

      bool put(xmlNodePtr cur) 
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(ID,"id"); //rAttrib.add(a->B0,"b");
        rAttrib.put(cur);
        return putContent(B0,cur);
      }

      void addOptimizables(VarRegistry<T>& vlist) {
        vlist.add(ID,&B0,1);
      }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1691 $   $Date: 2007-02-01 15:51:50 -0600 (Thu, 01 Feb 2007) $
 * $Id: WMConstraints.h 1691 2007-02-01 21:51:50Z jnkim $ 
 ***************************************************************************/

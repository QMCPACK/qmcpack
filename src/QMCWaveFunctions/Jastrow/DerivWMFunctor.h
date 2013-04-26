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
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** Implements derivative of WM Functor \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$ with respect to B
 */
template<class T>
struct DWMDBFunctor: public OptimizableFunctorBase<T>
{
  ///typedef of real values
  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
  ///input B
  real_type B0;
  ///B+2
  real_type Bplus2;
  ///input Rcut
  real_type Rcut;
  ///1/Rcut
  real_type OneOverRc;
  ///id
  string ID_B;
  ///name of B-attribute
  string attribName;
  ///constructor
  WMFunctor(real_type b, real_type rc=7.5, const std::string& bname="exponent"):
    attribName(bname)
  {
    reset(b,rc);
  }
  void reset(real_type b, real_type rc)
  {
    B0=b;
    Rcut=rc;
    Bplus2=B+2.0;
    OneOverRc=1.0/Rcut;
  }

  inline real_type f(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    real_type bz=1.0/(1+B0*z);
    return z*(1-z)*bz*bz;
  }

  inline real_type df(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return (z*Bplus2-1.0)/pow(1+B0*z,3)*OneOverRc*12*x*(1-2.0*x+x*x);
  }

  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID_B,"id");
    rAttrib.add(ID_B,"ref");
    rAttrib.add(B0,"exponent");
    rAttrib.put(cur);
    ID_B.append("_E");
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    vlist[ID_B]=B0;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
    if(it_b != optVariables.end())
    {
      B0=(*it_b).second;
    }
  }
};

/** Implements a screened Function \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$
 *
 * Short-range functor introduced by Wagner and Mitas, cond-mat/0610088
 */
template<class T>
struct D2WMDB2Functor: public OptimizableFunctorBase<T>
{
  ///typedef of real values
  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
  ///input B
  real_type B0;
  ///B+2
  real_type Bplus2;
  ///B+3
  real_type Bplus3;
  ///input Rcut
  real_type Rcut;
  ///1/Rcut
  real_type OneOverRc;
  ///id
  string ID_B;
  ///name of B-attribute
  string attribName;
  ///constructor
  D2WMDB2Functor(real_type b, real_type rc=7.5, const std::string& bname="exponent"):
    attribName(bname)
  {
    reset(b,rc);
  }
  void reset(real_type b, real_type rc)
  {
    B0=b;
    Rcut=rc;
    Bplus2=B+2.0;
    Bplus3=B+3.0;
    OneOverRc=1.0/Rcut;
  }

  /** return \f$ \frac{2 z^2 (1-z)}{(1+Bz)^3}\f$
   */
  inline real_type f(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return 2.0*z*z*(1-z)/pow(1+B0*z,3);
  }

  /** return \f$ -\frac{2 z (3z+bz-2)}{(1+Bz)^4}\f$
   */
  inline real_type df(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return 2.0*z*(2.0-z*Bplus3)/pow(1+B0*z,4)*OneOverRc*12*x*(1-2.0*x+x*x);
  }

  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID_B,"id");
    rAttrib.add(ID_B,"ref");
    rAttrib.add(B0,"exponent");
    rAttrib.put(cur);
    ID_B.append("_E");
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    vlist[ID_B]=B0;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
    if(it_b != optVariables.end())
    {
      B0=(*it_b).second;
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

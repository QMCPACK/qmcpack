//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file DerivPadeFunctors.h
 * @brief Derivatives of PadeFunctor with respect to B
 *
 * DPadeDBFunctor: first derivative with respect to B
 * D2PadeDB2Functor: second derivative with respect to B
 * These are not implemented since they are the special cases for
 * existing functors.
 * DPadeDAFunctor = PadeFunctor(A=1)
 * D2PadeDABFunctor = DPadeDBFunctor(A=1)
 * D2PadeDA2Functor = 0.0
 */
#ifndef QMCPLUSPLUS_DERIV_PADEFUNCTORS_H
#define QMCPLUSPLUS_DERIV_PADEFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"

namespace qmcplusplus
{

/** Implements \f$\frac{du[r]}{dB} = -\frac{Ar^2}{(1+Br)^2}\f$
 */
template<class T>
struct DPadeDBFunctor:public OptimizableFunctorBase<T>
{

  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;

  ///A
  real_type A;
  ///B
  real_type B;
  ///A2=2*A
  real_type A2;
  ///B2=2*B
  real_type B2;
  ///id of A
  std::string ID_A;
  ///id of B
  std::string ID_B;

  ///default constructor
  DPadeDBFunctor(): ID_A("0"), ID_B("0") { }
  ///constructor
  explicit DPadeDBFunctor(real_type a, real_type b): ID_A("0"),ID_B("0")
  {
    reset(a, b);
  }
  DPadeDBFunctor(const PadeFunctor<T>& pade): ID_A(pade.ID_A), ID_B(pade.ID_B)
  {
    reset(pade.A,pade.B);
  }

  void reset(real_type a, real_type b)
  {
    A=a;
    B=b;
    A2=2.0*A;
    B2=2.0*B;
  }

  /** evaluate the value at r \f$-Ar^2/(1+Br)^2\f$
   */
  inline real_type evaluate(real_type r)
  {
    real_type br=r/(1.0+B*r);
    return -A*br*br;
  }

  /** evaluate all the values
   *
   * \f$ u=-Ar^2/(1+Br)^2\f$
   * \f$ dudr=-2*A*r/(1+Br)^3\f$
   * \f$ d2udr2=2A*(2Br-1)/(1+Br)^4\f$
   */
  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type u = 1.0/(1.0+B*r);
    real_type uu=u*u;
    dudr = -A2*r*u*uu;
    d2udr2 = A2*(B2*r-1.0)*uu*uu;
    return -A*r*r*uu;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  bool put(xmlNodePtr cur)
  {
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(cur->name));
      if(cname == "parameter" || cname == "Var")
      {
        const xmlChar* nptr=xmlGetProp(cur,(const xmlChar *)"name");
        const xmlChar* iptr=xmlGetProp(cur,(const xmlChar *)"id");
        if(nptr == NULL || iptr == NULL)
          return false;
        if(nptr[0] == 'A')
        {
          ID_A = (const char*)iptr;
          putContent(A,cur);
        }
        else
          if(nptr[0] == 'B')
          {
            ID_B = (const char*)iptr;
            putContent(B,cur);
          }
      }
      cur = cur->next;
    }
    reset(A,B);
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    if(ID_A != "0")
      vlist[ID_A]=A;
    if(ID_B != "0")
      vlist[ID_B]=B;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(const OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_a(optVariables.find(ID_A));
    if(it_a != optVariables.end())
      A=(*it_a).second;
    typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
    if(it_b != optVariables.end())
      B=(*it_b).second;
    reset(A,B);
  }
};

/** Implements \f$\frac{d^2 u[r]}{dB^2}= \frac{2Ar^3}{(1+Br)^3}\f$
 */
template<class T>
struct D2PadeDB2Functor:public OptimizableFunctorBase<T>
{

  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;

  ///A
  real_type A;
  ///B
  real_type B;
  ///A2=2*A
  real_type A2;
  ///A6=6*A
  real_type A6;
  ///A12=12*A
  real_type A12;
  ///B2=2*B
  real_type B2;
  ///id of A
  std::string ID_A;
  ///id of B
  std::string ID_B;

  ///default constructor
  D2PadeDB2Functor(): ID_A("0"), ID_B("0") { }
  ///constructor
  explicit D2PadeDB2Functor(real_type a, real_type b): ID_A("0"),ID_B("0")
  {
    reset(a, b);
  }
  D2PadeDB2Functor(const PadeFunctor<T>& pade): ID_A(pade.ID_A), ID_B(pade.ID_B)
  {
    reset(pade.A,pade.B);
  }

  void reset(real_type a, real_type b)
  {
    A=a;
    B=b;
    A2=2.0*A;
    A6=6.0*A;
    A12=12.0*A;
  }

  /* return \f$ \frac{2A r^3}{(1+Br)^3}\f$
   */
  inline real_type evaluate(real_type r)
  {
    real_type br=r/(1.0+B*r);
    return A2*br*br*br;
  }

  /* evaluate all the values
   *
   * \f$ u = \frac{2A r^3}{(1+Br)^3}\f$
   * \f$ du = \frac{6A r^2}{(1+Br)^4}\f$
   * \f$ d2udr2 = \frac{12A r(1-Br)}{(1+Br)^5}\f$
   */
  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type rr=r*r;
    real_type u = 1.0/(1.0+B*r);
    real_type uu=u*u;
    real_type uuuu=uu*uu;
    dudr = A6*rr*uuuu;
    d2udr2 = A12*r*(1-B*r)*u*uuuu;
    return A2*r*rr*u*uu;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  bool put(xmlNodePtr cur)
  {
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(cur->name));
      if(cname == "parameter")
      {
        std::string pname('0');
        std::string iname('0');
        OhmmsAttributeSet rAttrib;
        rAttrib.add(pname,"name");
        rAttrib.add(iname,"id");
        rAttrib.add(iname,"ref"); //accept either id/ref
        rAttrib.put(cur);
        if(pname[0] == 'A')
        {
          ID_A=iname;
          putContent(A,cur);
        }
        else
          if(pname[0] =='B')
          {
            ID_B=iname;
            putContent(B,cur);
          }
      }
      cur = cur->next;
    }
    reset(A,B);
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    if(ID_A != "0")
      vlist[ID_A]=A;
    if(ID_B != "0")
      vlist[ID_B]=B;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_a(optVariables.find(ID_A));
    if(it_a != optVariables.end())
      A=(*it_a).second;
    typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
    if(it_b != optVariables.end())
      B=(*it_b).second;
    reset(A,B);
  }
};
}
#endif


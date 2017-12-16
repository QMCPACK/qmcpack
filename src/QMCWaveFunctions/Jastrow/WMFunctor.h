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
    
    
#ifndef QMCPLUSPLUS_WMFUNCTOR_WITHPADECORRECTION_H
#define QMCPLUSPLUS_WMFUNCTOR_WITHPADECORRECTION_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** Implements a screened Function \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$
 *
 * Short-range functor introduced by Wagner and Mitas, cond-mat/0610088
 */
template<class T>
struct WMFunctor: public OptimizableFunctorBase
{
  ///input B
  real_type B0;
  ///input Rcut
  real_type Rcut;
  ///1/Rcut
  real_type OneOverRc;
  ///12/Rcut
  real_type DxDr;
  ///12/Rcut/Rcut
  real_type D2xDr2;
  ///id
  std::string ID_B;
  ///name of B-attribute
  std::string attribName;
  ///constructor
  WMFunctor(real_type b, real_type rc=7.5, const std::string& bname="exponent"):
    B0(b),Rcut(rc),attribName(bname)
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new WMFunctor<T>(*this);
  }

  void reset()
  {
    //B0=b; Rcut=rc;
    OneOverRc=1.0/Rcut;
    DxDr=12*OneOverRc;
    D2xDr2=OneOverRc*DxDr;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return -(1+B0)/(1+B0*z)/(1+B0*z)*OneOverRc*12*x*(1-2.0*x+x*x);
  }

  inline real_type evaluate(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type z=x*x*(6.0-8*x+3.0*x*x);
    return (1-z)/(1+B0*z);
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if(r>Rcut)
    {
      dudr=0.0;
      d2udr2=0.0;
      return 0.0;
    }
    real_type x=r*OneOverRc;
    real_type xx=x*x;
    real_type z=xx*(6.0-8*x+3.0*xx);
    real_type u=1.0/(1+B0*z);
    real_type du=(1.0+B0)*u*u;
    real_type dzdr=DxDr*x*(1.0-2.0*x+xx);
    dudr=-du*dzdr;
    d2udr2=du*(2.0*B0*u*dzdr*dzdr-D2xDr2*(1.0-4.0*x+3.0*xx));
    return (1.0-z)*u;
  }

  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet rAttrib;
    rAttrib.add(ID_B,"id");
    rAttrib.add(ID_B,"ref");
    rAttrib.add(B0,"exponent");
    rAttrib.put(cur);
    ID_B.append("_E");
    myVars.insert(ID_B,B0,true);
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    //disable optimization of E
    //active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    //disable optimization of E
    //myVars.getIndex(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    //disable optimization of E
    //int loc=myVars.where(0);
    //if(loc>=0) {myVars[0]=B0=active[loc];}
  }
};

template<class T>
struct WMFunctorSet: public OptimizableFunctorBase
{
  //contraction and exponent pair
  typedef TinyVector<real_type,2>  variable_type;
  ///input Rcut
  real_type Rcut;
  ///1/Rcut
  real_type OneOverRc;
  ///12/Rcut
  real_type DxDr;
  ///12/Rcut/Rcut
  real_type D2xDr2;
  ///input  Param[i][0]=B; Param[i][1]=C for i-th component
  std::vector<variable_type> Params;

  WMFunctorSet(real_type rc=7.5):Rcut(rc)
  {
    OneOverRc=1.0/Rcut;
    DxDr=12*OneOverRc;
    D2xDr2=OneOverRc*DxDr;
    Params.reserve(8);//limit to 8
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new WMFunctorSet<T>(*this);
  }

  inline void reset() {}

  inline void addComponent(real_type c, real_type b, const std::string& aname)
  {
    int i=Params.size();
    Params.push_back(variable_type(c,b));
    std::string cname=aname+"_C";
    std::string ename=aname+"_E";
    myVars.insert(cname,Params[i][0]);
    myVars.insert(ename,Params[i][1]);
  }

  inline real_type evaluate(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type xx=x*x;
    real_type z=xx*(6.0-8*x+3.0*xx);
    real_type value=0.0;
    for(int i=0; i<Params.size(); ++i)
      value += Params[i][0]*(1.0-z)/(1+Params[i][1]*z);
    return value;
  }

  inline real_type evaluate(real_type r, real_type& dv, real_type& d2v)
  {
    dv=0.0;
    d2v=0.0;
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type xx=x*x;
    real_type z=xx*(6.0-8*x+3.0*xx);
    real_type dzdr=DxDr*x*(1.0-2.0*x+xx);
    real_type dzdrsq=dzdr*dzdr;
    real_type d2zdr2=D2xDr2*(1.0-4.0*x+3.0*xx);
    real_type v=0.0;
    for(int i=0; i<Params.size(); ++i)
    {
      real_type c=Params[i][0];
      real_type b=Params[i][1];
      real_type u=1.0/(1+b*z);
      real_type du=c*(1.0+b)*u*u;
      v  += c*(1.0-z)*u;
      dv -= du*dzdr;
      d2v+= du*(2.0*b*u*dzdrsq-d2zdr2);
    }
    return v;
  }

  inline bool
  evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    if(r>=Rcut)
      return false;
    real_type x=r*OneOverRc;
    real_type xx=x*x;
    real_type z=xx*(6.0-8*x+3.0*xx);
    real_type dzdr=DxDr*x*(1.0-2.0*x+xx);
    real_type dzdrsq=dzdr*dzdr;
    real_type d2zdr2=D2xDr2*(1.0-4.0*x+3.0*xx);
    typedef TinyVector<real_type,3> e_type;
    typename std::vector<e_type>::iterator dit(derivs.begin());
    for(int i=0,ii=0; i<Params.size(); ++i)
    {
      real_type c=Params[i][0];
      real_type b=Params[i][1];
      real_type u=1.0/(1+b*z);
      real_type uu=u*u;
      real_type du=(1.0+b)*uu;
      (*dit++)=e_type((1.0-z)*u,-du*dzdr,du*(2.0*b*u*dzdrsq-d2zdr2));
      //Disable derivatives with respect to exponents: incorrect and unstable
      //real_type dudz=(z*(b+2)-1.0)*uu*u;
      //(*dit++)=e_type(c*z*(z-1)*uu,c*dudz*dzdr,c*uu*uu*((b+2)*(1.0-2*b*z)+3*b)*dzdrsq+dudz*d2zdr2);
    }
    return true;
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }
  inline real_type df(real_type r)
  {
    if(r>Rcut)
      return 0.0;
    real_type x=r*OneOverRc;
    real_type xx=x*x;
    real_type z=xx*(6.0-8*x+3.0*xx);
    real_type dzdr=DxDr*x*(1.0-2.0*x+xx);
    real_type dv=0.0;
    for(int i=0; i<Params.size(); ++i)
    {
      real_type c=Params[i][0];
      real_type b=Params[i][1];
      real_type u=1.0/(1+b*z);
      dv -= c*(1.0+b)*u*u*dzdr;
    }
    return dv;
  }

  bool put(xmlNodePtr cur)
  {
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

  void resetParameters(const opt_variables_type& active)
  {
    int ii=0;
    for(int i=0; i<Params.size(); ++i)
    {
      int loc_c=myVars.where(ii++);
      int loc_e=myVars.where(ii++);
      if(loc_c>=0)
        Params[i][0]=active[loc_c];
      if(loc_e>=0)
        Params[i][1]=active[loc_e];
    }
  }
};
}
#endif

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file NewPadeFunctors.h
 * @brief Functors which implement Pade functions
 */
#ifndef QMCPLUSPLUS_NEW_PADEFUNCTORS_H
#define QMCPLUSPLUS_NEW_PADEFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus
{

/** Implements a Pade Function \f$u[r]=A*r/(1+B*r) - A/B\f$
 *  cusp = A, so only one optimizable parameter
 */
template<class T>
struct pade1functor:public OptimizableFunctorBase
{
  ///true, if B is optimizable
  bool Opt_B;
  ///input A
  real_type A;
  ///input B
  real_type B;
  ///AoverB=A/B
  real_type AoverB;
  ///id of B
  std::string ID_B;

  ///default constructor
  pade1functor(): ID_B("0") { }

  ///constructor
  explicit pade1functor(real_type cusp, real_type b):
  A(cusp),B(b),AoverB(cusp/b),Opt_B(true),ID_B("0"))
  {
    reset();
  }


  OptimizableFunctorBase* makeClone() const
  {
    return new pade1functor(*this);
  }

  void reset()
  {
    AoverB=A/B;
  }

  inline real_type evaluate(real_type r) const
  {
    return A*r/(1.0+B*r)-AoverB;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    const real_type u = 1.0/(1.0+B*r);
    dudr = A*u*u;
    d2udr2 = -2*B*u*dudr
    return A*u*r-AoverB;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    real_type u = 1.0/(1.0+B*r);
    dudr = A*u*u;
    d2udr2 = -2*B*u*dudr;
    d3udr3 = -3.0*B*u*d2udr2;
    return A*u*r-AoverB;
  }

  inline real_type evaluateV(const int iat, const int iStart, const int iEnd,
    const T* restrict _distArray, T* restrict distArrayCompressed ) const
  {
    real_type sum(0);
    for(int idx=iStart; idx<iEnd; idx++)
      if (idx!=iat) sum += evaluate(_distArray[idx]);
    return sum;
  }

  inline void evaluateVGL(const int iat, const int iStart, const int iEnd,
    const T* distArray,  T* restrict valArray,
    T* restrict gradArray, T* restrict laplArray,
    T* restrict distArrayCompressed, int* restrict distIndices ) const
  {
    for(int idx=iStart; idx<iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if ( iat>=iStart && iat<iEnd )
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
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

  /// compute derivatives with respect to B
  inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    int i=0;
    real_type u = 1.0/(1.0+B*r);
    if(Opt_B)
    {
      const real_type bru = B*r*u;
      derivs[i][0]= -A*r*r*u*u; //du/db
      derivs[i][1]= -2*A*r*u*u*u; //d(du/db)/dr
      derivs[i][2]= 2.0*A*u*u*(-3.0*bru*bru + 4.0*bru - 1.0); //d^2(du/db)/dr^2
      ++i;
    }
    return true;
  }

  /// compute derivatives with respect to B
  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    int i=0;
    real_type u = 1.0/(1.0+B*r);
    if(Opt_B)
    {
      derivs[i]= -A*r*r*u*u; //du/db
      ++i;
    }
    return true;
  }

  bool put(xmlNodePtr cur)
  {
    real_type Btemp(B0);
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "var" || cname == "parameter") 
      {
        std::string id_in("0");
        std::string p_name("B");
	std::string doopt("no");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
	rAttrib.add(doopt, "optimize");
        rAttrib.put(cur);
        if(p_name=="B")
        else if(p_name == "B")
        {
          ID_B = id_in;
          putContent(Btemp,cur);
          Opt_B=(doopt != "no");
        }
      }
      cur = cur->next;
    }
    B=Btemp;
    reset();
    myVars.clear();
    if(Opt_B)
      myVars.insert(ID_B,B0,Opt_B,optimize::OTHER_P);
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if(myVars.size())
    {
      int ia=myVars.where(0);
      if(ia>-1)
      {
        int i=0;
        if(Opt_B)
          B0=myVars[i++]=active[ia];
      }
      reset();
    }
  }
};


/** 2nd order Pade function of \f[ u(r) = \frac{a*r+b*r^2}{1+c*r+d*r^2} - B/D\f]
 *  cusp = A, so three optimizable parameters
 */
template<class T>
struct pade2functor:public OptimizableFunctorBase
{

  ///coefficients
  real_type A, B, C, D, BoverD;
  bool Opt_B, Opt_C, Opt_D;
  ///id for A
  std::string ID_B;
  ///id for B
  std::string ID_C;
  ///id for C
  std::string ID_D;
  
  //default constructor
 pade2functor(): ID_B("0"), ID_C("1"), ID_D("2") { }

  ///constructor
 pade2functor(real_type cusp=0.0, real_type b=1.0, real_type c=1.0, real_type d=1.0) :
  A(cusp),B(b),C(c),D(d),ID_B("0"),ID_C("1"),ID_C("2"), Opt_B(true), Opt_C(true), Opt_D(true)
    {
      reset();
    }
  
  OptimizableFunctorBase* makeClone() const
  {
    return new pade2functor(*this);
  }

  /** reset the internal variables.
   */
  void reset()
  {
    BoverD=B/D;
  }

  inline real_type evaluate(real_type r)
  {
    const real_type denom = 1.0/(1.0+C*r+D*r*r);
    const real_type abr = A+B*r;
    return abr*denom;
  }

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {

    const real_type denom = 1.0/(1.0+C*r+D*r*r);
    const real_type abr = A+B*r;
    const real_type a2br = A+2*B*r;
    const real_type c2dr = C+2*D*r;
    dudr = denom*(a2br - abr*c2dr*denom*r);
    d2udr2 = 2*denom*(-a2br*c2dr*denom + abr*c2dr*c2dr*denom*denom*r - abr*D*denom*r + B)
    return abr*denom;
  }

  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type& d3udr3)
  {
    const real_type denom = 1.0/(1.0+C*r+D*r*r);
    const real_type abr = A+B*r;
    const real_type a2br = A+2*B*r;
    const real_type c2dr = C+2*D*r;
    dudr = denom*(a2br - abr*c2dr*denom*r);
    d2udr2 = 2*denom*(-a2br*c2dr*denom + abr*c2dr*c2dr*denom*denom*r - abr*D*denom*r + B);
    d3udr3 = 6*denom*denom*(a2br*c2dr*c2dr*denom - a2br*D - abr*c2dr*c2dr*c2dr*denom*denom*r + 2*abr*c2dr*D*denom*r - B*c2dr);
    return abr*denom;
  }

  real_type f(real_type r)
  {
    return evaluate(r);
  }

  real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    const real_type denom = 1.0/(1.0+C*r+D*r*r);
    const real_type abr = A+B*r;
    const real_type a2br = A+2.0*B*r;
    const real_type c2dr = C+2.0*D*r;
    const real_type combo = c2dr*denom*r;
    int i=0;
    if(Opt_B)
    {
      derivs[i][0]= denom*r*r;
      derivs[i][1]= denom*r*(-combo + 2.0);
      derivs[i][2]= 2*denom*(combo*combo - 2.0*combo - D*denom*r*r + 1.0);
      i++;
    }

    if(Opt_C)
    {
      derivs[i][0]= -abr*denom*denom*r*r;
      derivs[i][1]= denom*denom*r*(2.0*abr*combo - (2.0*A+3.0*B*r));
      derivs[i][2]= 2.0*denom*denom*((2.0*combo - 1.0)*(2.0*A+3.0*B*r) - 3.0*abr*combo*combo + 2.0*abr*d*denom*r*r);
      i++;
    }

    if(Opt_D)
    {
      derivs[i][0]= -denom*denom*r*r*r*r*abr;
      derivs[i][1]= denom*denom*r*r*(2.0*abr*combo - (3.0*A+4*B*r));
      derivs[i][2]= -2.0*denom*denom*r*(3.0*a2br - 2.0*a2br*combo + 3.0*abr*combo*combo - 4.0*abr*combo - 2.0*abr*D*denom*r*r);
;
      i++;
    }
    return true;
  }



  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    const real_type denom = 1.0/(1.0+C*r+D*r*r);
    const real_type abr = A+B*r;
    int i=0;
    if(Opt_B)
    {
      derivs[i]= denom*r*r;
      i++;
    }
    if(Opt_C)
    {
      derivs[i]= -abr*denom*denom*r*r;
      i++;
    }
    if(Opt_D)
    {
      derivs[i]= -denom*denom*r*r*r*r*abr;
      i++;
    }
    return true;
  }

  inline real_type evaluateV(const int iat, const int iStart, const int iEnd,
    const T* restrict _distArray, T* restrict distArrayCompressed ) const
  {
    real_type sum(0);
    for(int idx=iStart; idx<iEnd; idx++)
      if (idx!=iat) sum += evaluate(_distArray[idx]);
    return sum;
  }

  inline void evaluateVGL(const int iat, const int iStart, const int iEnd,
    const T* distArray,  T* restrict valArray,
    T* restrict gradArray, T* restrict laplArray,
    T* restrict distArrayCompressed, int* restrict distIndices ) const
  {
    for(int idx=iStart; idx<iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if ( iat>=iStart && iat<iEnd )
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
  }

  /** process input xml node
   * @param cur current xmlNode from which the data members are reset
   *
   * T1 is the type of VarRegistry, typically double.
   * Read in the Pade parameters from the xml input file.
   */
  bool put(xmlNodePtr cur)
  {
    real_type Btemp,Ctemp, Dtemp;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    bool renewed=false;
    while(tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "var")
      {
        std::string doopt("yes");
        std::string id_in("0");
        std::string p_name("B");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.add(doopt, "optimize");
        rAttrib.put(tcur);
        if(p_name=="B")
        {
          ID_B=id_in;
          Opt_B=(doopt != "no");
          putContent(Btemp,tcur);
          renewed=true;
        }
        else if(p_name == "C")
        {
          ID_C = id_in;
          Opt_C=(doopt != "no");
          putContent(Ctemp,tcur);
          renewed=true;
        }
        else if(p_name == "D")
        {
          ID_D = id_in;
          Opt_D=(doopt != "no");
          putContent(Dtemp,tcur);
          renewed=true;
        }
      }
      tcur = tcur->next;
    }
    if (renewed)
    {
      B=Btemp;
      C=Ctemp;
      D=Dtemp;
      reset();
      //these are always active
      myVars.clear();
      if(Opt_B)
        myVars.insert(ID_B,B,true,optimize::OTHER_P);
      if(Opt_C)
        myVars.insert(ID_C,C,true,optimize::OTHER_P);
      if(Opt_D)
        myVars.insert(ID_D,D,true,optimize::OTHER_P);
    }
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
    int i=0;
    if (ID_B!="0")
    {
      int ib=myVars.where(i);
      if(ib>-1)
        B=myVars[i]=active[ib];
      i++;
    }
    if (ID_C!="0")
    {
      int ic=myVars.where(i);
      if(ib>-1)
        C=myVars[i]=active[ic];
      i++;
    }
    if (ID_D!="0")
    {
      int id=myVars.where(i);
      if(id>-1)
        D=myVars[i]=active[id];
      i++;
    }
  }
};

#endif


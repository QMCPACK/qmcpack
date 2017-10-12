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
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file PadeFunctors.h
 * @brief Functors which implement Pade functions
 */
#ifndef QMCPLUSPLUS_PADEFUNCTORS_H
#define QMCPLUSPLUS_PADEFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus
{

/** Implements a Pade Function \f$u[r]=A*r/(1+B*r)\f$
 *
 * Similar to PadeJastrow with a scale.
 */
template<class T>
struct PadeFunctor:public OptimizableFunctorBase
{
  ///true, if A is optimizable
  bool Opt_A;
  ///true, if B is optimizable
  bool Opt_B;
  ///input A
  real_type A;
  ///input B
  real_type B0;
  ///input scaling, default=1.0
  real_type Scale;
  ///B=B0*Scale
  real_type B;
  ///AB=A*B
  real_type AB;
  ///B2=2*B
  real_type B2;
  ///AoverB=A/B
  real_type AoverB;
  ///id of A
  std::string ID_A;
  ///id of B
  std::string ID_B;

  ///default constructor
  PadeFunctor(): Scale(1.0),ID_A("0"),ID_B("0") { }

  ///constructor
  explicit PadeFunctor(real_type a, real_type b, real_type s=1.0):
    A(a),B0(b),Scale(s),Opt_A(true),Opt_B(true)
  {
    reset();
  }

  /** constructor with A
   * @param a value of A
   * @param ida id of A
   *
   * Special constructor for two-body Jastrow for Coulomb interactions
   * Automatically fix the cusp conditions
   */
  explicit PadeFunctor(real_type a, const std::string& ida)
    :A(a),B0(1.0),Scale(1.0),ID_A(ida),Opt_A(false),Opt_B(true)
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new PadeFunctor(*this);
  }

  void reset()
  {
    cutoff_radius=1.0e4; //some big range
    //A=a; B0=b; Scale=s;
    B = B0*Scale;
    AB = A*B;
    B2=2.0*B;
    AoverB=A/B;
  }

  inline real_type evaluate(real_type r) const
  {
    return A*r/(1.0+B*r)-AoverB;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    real_type u = 1.0/(1.0+B*r);
    dudr = A*u*u;
    d2udr2 = -B2*dudr*u;
    return A*u*r-AoverB;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    real_type u = 1.0/(1.0+B*r);
    dudr = A*u*u;
    d2udr2 = -B2*dudr*u;
    d3udr3 = -3.0*B*d2udr2*u;
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
    return evaluate(r)-AoverB;
  }

  inline real_type df(real_type r)
  {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  /// compute derivatives with respect to A and B
  inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    int i=0;
    real_type u = 1.0/(1.0+B*r);
    if(Opt_A)
    {
      derivs[i][0]= r*u-1/B; //du/da
      derivs[i][1]= u*u; //d(du/da)/dr
      derivs[i][2]= -B2*u*u*u; //d^2 (du/da)/dr
      ++i;
    }
    if(Opt_B)
    {
      derivs[i][0]= -A*r*r*u*u+AoverB/B; //du/db
      derivs[i][1]= -2.0*A*r*u*u*u; //d(du/db)/dr
      derivs[i][2]=  2.0*A*(B2*r-1)*u*u*u*u; //d^2(du/db)/dr^2
    }
    return true;
  }

  /// compute derivatives with respect to A and B
  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    int i=0;
    real_type u = 1.0/(1.0+B*r);
    if(Opt_A)
    {
      derivs[i]= r*u-1/B; //du/da
      ++i;
    }
    if(Opt_B)
    {
      derivs[i]= -A*r*r*u*u+AoverB/B; //du/db
    }
    return true;
  }

  bool put(xmlNodePtr cur)
  {
    real_type Atemp(A),Btemp(B0);
    cur = cur->xmlChildrenNode;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "var") //only accept var
      {
        std::string id_in("0");
        std::string p_name("B");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.put(cur);
        if(p_name=="A")
        {
          ID_A = id_in;
          putContent(Atemp,cur);
          Opt_A=true;
        }
        else if(p_name == "B")
        {
          ID_B = id_in;
          putContent(Btemp,cur);
          Opt_B=true;
        }
      }
      cur = cur->next;
    }
    A=Atemp;
    B0=Btemp;
    reset();
    myVars.clear();
    if(Opt_A)
      myVars.insert(ID_A,A, Opt_A,optimize::LOGLINEAR_P);
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
        if(Opt_A)
          A=myVars[i++]=active[ia++];
        if(Opt_B)
          B0=myVars[i]=active[ia];
      }
      reset();
    }
  }
};


/** Pade function of \f[ u(r) = \frac{a*r+c*r^2}{1+b*r} \f]
 *
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct Pade2ndOrderFunctor:public OptimizableFunctorBase
{

  ///coefficients
  real_type A, B, C, C2;
  bool Opt_A, Opt_B, Opt_C;
  ///id for A
  std::string ID_A;
  ///id for B
  std::string ID_B;
  ///id for C
  std::string ID_C;

  ///constructor
  Pade2ndOrderFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0)
    : A(a),B(b),C(c),ID_A("0"),ID_B("0"),ID_C("0")
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new Pade2ndOrderFunctor(*this);
  }

  /** reset the internal variables.
   */
  void reset()
  {
    // A = a; B=b; C = c;
    C2 = 2.0*C;
  }

  /**@param r the distance
    @return \f$ u(r) = a*r/(1+b*r) \f$
    */
  inline real_type evaluate(real_type r)
  {
    real_type br(B*r);
    return (A+br)*r/(1.0+br);
  }

  /** evaluate the value at r
   * @param r the distance
   @param dudr return value  \f$ du/dr = a/(1+br)^2 \f$
   @param d2udr2 return value  \f$ d^2u/dr^2 = -2ab/(1+br)^3 \f$
   @return \f$ u(r) = a*r/(1+b*r) \f$
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type u = 1.0/(1.0+B*r);
    real_type v = A*r+C*r*r;
    real_type w = A+C2*r;
    dudr = u*(w-B*u*v);
    d2udr2 = 2.0*u*u*u*(C-B*A);
    return u*v;
  }

  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type& d3udr3)
  {
    real_type u = 1.0/(1.0+B*r);
    real_type v = A*r+C*r*r;
    real_type w = A+C2*r;
    dudr = u*(w-B*u*v);
    d2udr2 = 2.0*u*u*u*(C-B*A);
    std::cerr << "Third derivative not imlemented for Pade functor.\n";
    return u*v;
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
    real_type u = 1.0/(1.0+B*r);
    real_type u2 = u*u;
    real_type u3 = u*u*u;
    real_type u4 = u*u*u*u;
    int i=0;
    if(Opt_A)
    {
      derivs[i][0]= r*u;
      derivs[i][1]= u2;
      derivs[i][2]= -2.0*B*u3;
      i++;
    }

    if(Opt_B)
    {
      derivs[i][0]= -r*r*(A+C*r)*u2;
      derivs[i][1]= -r*(2.0*A+C*r*(3.0+B*r))*u3;
      derivs[i][2]= -2.0*u4*(A-2.0*A*B*r+3.0*C*r);
      i++;
    }

    if(Opt_C)
    {
      derivs[i][0]= r*r*u;
      derivs[i][1]= r*(2.0+B*r)*u2;
      derivs[i][2]= 2.0*u3;
      i++;
    }
    return true;
  }



  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    real_type u = 1.0/(1.0+B*r);
    int i=0;
    if(Opt_A)
    {
      derivs[i]= r*u;
      i++;
    }
    if(Opt_B)
    {
      derivs[i]= -r*r*(A+C*r)*u*u;
      i++;
    }
    if(Opt_C)
    {
      derivs[i]= r*r*u;
      i++;
    }
    return true;
  }


  /** process input xml node
   * @param cur current xmlNode from which the data members are reset
   *
   * T1 is the type of VarRegistry, typically double.
   * Read in the Pade parameters from the xml input file.
   */
  bool put(xmlNodePtr cur)
  {
    real_type Atemp,Btemp, Ctemp;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    bool renewed=false;
    while(tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string doopt("yes");
        std::string id_in("0");
        std::string p_name("B");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.add(doopt, "optimize");
        rAttrib.put(tcur);
        if(p_name=="A")
        {
          ID_A=id_in;
          Opt_A=(doopt != "no");
          putContent(Atemp,tcur);
          renewed=true;
        }
        else if(p_name == "B")
        {
          ID_B = id_in;
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
      }
      tcur = tcur->next;
    }
    if (renewed)
    {
      A=Atemp;
      B=Btemp;
      C=Ctemp;
      reset();
      //these are always active
      myVars.clear();
      if(Opt_A)
        myVars.insert(ID_A,A,true,optimize::LOGLINEAR_P);
      if(Opt_B)
        myVars.insert(ID_B,B,true,optimize::OTHER_P);
      if(Opt_C)
        myVars.insert(ID_C,C,true,optimize::LOGLINEAR_P);
    }
    //LOGMSG("Jastrow (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")")
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
    if (ID_A!="0")
    {
      int ia=myVars.where(i);
      if(ia>-1)
        A=myVars[i]=active[ia];
      i++;
    }
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
      if(ic>-1)
        C=myVars[i]=active[ic];
      i++;
    }
    C2 = 2.0*C;
  }
};

/** Pade function of \f[ u(r) = \frac{a*r+b*r^2}{1+c*r+d*r^2} \f]
*
* Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
*/
template<class T>
struct PadeTwo2ndOrderFunctor:public OptimizableFunctorBase
{

  ///coefficients
  real_type A, B, C, D;
  bool Opt_A, Opt_B, Opt_C, Opt_D;
  ///id for A
  std::string ID_A;
  ///id for B
  std::string ID_B;
  ///id for C
  std::string ID_C;
  ///id for D
  std::string ID_D;

  ///constructor
  PadeTwo2ndOrderFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0, real_type d=1.0)
    : A(a),B(b),C(c),D(d), Opt_A(false),Opt_B(false),Opt_C(false),Opt_D(false),
    ID_A("0"),ID_B("0"),ID_C("0"),ID_D("0")
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new PadeTwo2ndOrderFunctor(*this);
  }

  /** reset the internal variables.
   */
  void reset()
  {
    // A = a; B=b; C = c;
  }

  inline real_type evaluate(real_type r)
  {
    real_type br(B*r);
    real_type dr(D*r);
    return (A*r+br*r)/(1.0+C*C*r+dr*dr);
  }

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type ar(A*r);
    real_type br(B*r);
    real_type cr(C*r);
    real_type dr(D*r);
    real_type bttm( 1.0/(1.0+C*cr+dr*dr) );
    dudr = (A - A*dr*dr + br*(2.0 + C*cr))*bttm*bttm;
    d2udr2 = -2.0*(A*(C*C + 3.0*dr*D - dr*dr*dr*D) + B*(-1.0 + dr*dr*(3.0 + C*cr)))*bttm*bttm*bttm;
    return (A*r+br*r)*bttm;
  }

  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type& d3udr3)
  {
    d3udr3 = 0;
    return evaluate(r,dudr,d2udr2);
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
    real_type ar(A*r);
    real_type br(B*r);
    real_type cr(C*r);
    real_type dr(D*r);
    real_type r2(r*r);
    real_type dr2(D*r2);
    real_type d2r2(dr*dr);
    real_type c2(C*C);
    real_type c2r(c2*r);
    real_type bttm(1.0/(1.0+c2r+d2r2));
    real_type tp( A*r+br*r );
    real_type c2r2(cr*cr);
    real_type d4r4(d2r2*d2r2);
    real_type bttm2(bttm*bttm);
    real_type bttm3(bttm*bttm2);
    real_type bttm4(bttm2*bttm2);
    int i=0;
    if(Opt_A)
    {
      derivs[i][0]= r2*bttm;
      derivs[i][1]= r*(2.0+c2r)*bttm2;
      derivs[i][2]= (2.0 - 2.0*d2r2*(3.0+c2r))*bttm3;
      i++;
    }

    if(Opt_B)
    {
      derivs[i][0]= -2.0*cr*tp*bttm2;
      derivs[i][1]= -2*cr*(A*(2.0 - 2.0*d2r2) + br*(3.0 + c2r - d2r2))*bttm3;
      derivs[i][2]= 4.0*C*(A*(-1.0 + 2.0*c2r + 8.0*d2r2 - 3.0*d4r4) + br*(-3.0 - d4r4 + 2.0*d2r2 *(4.0 + c2r)))*bttm4;
      i++;
    }

    if(Opt_C)
    {
      derivs[i][0]= -2.0*dr2*tp*bttm2;
      derivs[i][1]= -2.0*dr2*(2.0*br*(2.0 + c2r) + A*(3.0 + c2r - d2r2))*bttm3;
      derivs[i][2]= -4.0*dr*(br*(6.0 + 4.0* c2r + c2*c2r2 - 6.0*d2r2 - 2.0*c2r*d2r2) + A*(3.0 + d4r4 - 2.0 *d2r2 *(4.0 + c2r))) *bttm4;
      i++;
    }
    return true;
  }


  inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
  {
    real_type ar(A*r);
    real_type br(B*r);
    real_type cr(C*r);
    real_type dr(D*r);
    real_type r2(r*r);
    real_type dr2(D*r2);
    real_type d2r2(dr*dr);
    real_type c2(C*C);
    real_type c2r(c2*r);
    real_type bttm(1.0/(1.0+c2r+d2r2));
    real_type tp( A*r+br*r );
    real_type c2r2(cr*cr);
    real_type d4r4(d2r2*d2r2);
    real_type bttm2(bttm*bttm);
    int i=0;
    if(Opt_A)
    {
      derivs[i]= r2*bttm;
      i++;
    }

    if(Opt_B)
    {
      derivs[i]= -2.0*cr*tp*bttm2;
      i++;
    }

    if(Opt_C)
    {
      derivs[i]= -2.0*dr2*tp*bttm2;
      i++;
    }
    return true;
  }




  /** process input xml node
   * @param cur current xmlNode from which the data members are reset
   *
   * T1 is the type of VarRegistry, typically double.
   * Read in the Pade parameters from the xml input file.
   */
  bool put(xmlNodePtr cur)
  {
    std::string fcup("yes");
    OhmmsAttributeSet p;
    p.add(fcup,"fixcusp");
    p.put(cur);
    if (fcup=="true")
      fcup="yes";
    // if (fcup=="yes") app_log()<<" fixing cusp conditions"<< std::endl;
    real_type Atemp=A,Btemp=B, Ctemp=C, Dtemp=D;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    bool renewed=false;
    while(tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string id_in("0");
        std::string p_name("B");
        std::string doopt("yes");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.add(doopt, "optimize");
        rAttrib.put(tcur);
        if(p_name=="A")
        {
          ID_A=id_in;
          Opt_A=(doopt!="no");
          putContent(Atemp,tcur);
          renewed=true;
        }
        else if(p_name == "B")
        {
          ID_B=id_in;
          Opt_B=(doopt!="no");
          putContent(Btemp,tcur);
          renewed=true;
        }
        else if(p_name == "C")
        {
          ID_C=id_in;
          Opt_C=(doopt!="no");
          putContent(Ctemp,tcur);
          renewed=true;
        }
        else if(p_name == "D")
        {
          ID_D=id_in;
          Opt_D=(doopt!="no");
          putContent(Dtemp,tcur);
          renewed=true;
        }
      }
      tcur = tcur->next;
    }
    if (renewed)
    {
      A=Atemp;
      B=Btemp;
      C=Ctemp;
      D=Dtemp;
      reset();
      //these are always active
      myVars.clear();
      if(Opt_A)
        myVars.insert(ID_A,B,true,optimize::LOGLINEAR_P);
      if(Opt_B)
        myVars.insert(ID_B,B,true,optimize::LOGLINEAR_P);
      if(Opt_C)
        myVars.insert(ID_C,C,true,optimize::OTHER_P);
      if(Opt_D)
        myVars.insert(ID_D,D,true,optimize::OTHER_P);
      //myVars.insert(ID_A,A,fcup!="yes");
    }
    //LOGMSG("Jastrow (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")")
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
    if(myVars.size()==0) return;

    int ia=myVars.where(0);
    if(ia<0) return;

    int i=0;
    if(Opt_A)
      A=myVars[i++]=active[ia++];
    if(Opt_B)
      B=myVars[i++]=active[ia++];
    if(Opt_C)
      C=myVars[i++]=active[ia++];
    if(Opt_D)
      D=myVars[i++]=active[ia++];

    reset();
  }
};

/** Pade functional of \f[ u(r) = \frac{a*f(r)}{1+b*f(r)} \f] with a scale function f(r)
 *
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct ScaledPadeFunctor:public OptimizableFunctorBase
{

  ///coefficients
  real_type A, B, C;
  real_type OneOverC, B2;

  ///constructor
  explicit ScaledPadeFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0)
    : A(a),B(b),C(c)
  {
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new ScaledPadeFunctor(*this);
  }

  /** reset the internal variables.
   */
  void reset()
  {
    OneOverC=1.0/C;
    B2=2.0*B;
  }

  /** evaluate the value at r
   * @param r the distance
   * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
   */
  inline real_type evaluate(real_type r)
  {
    real_type reff((1.0-std::exp(-C*r))*OneOverC);
    return A*reff/(1.0+B*reff);
  }

  /** evaluate the value, first derivative and second derivative
   * @param r the distance
   * @param dudr return value  \f$ du/dr\f$
   * @param d2udr2 return value  \f$ d^2u/dr^2 \f$
   * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
   */
  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type reff_1(std::exp(-C*r));
    real_type reff((1.0-reff_1)*OneOverC);
    real_type u(1.0/(1.0+B*reff));
    real_type auu(A*u*u);
    dudr=reff_1*auu;
    //d2udr=auu*(-C*reff_1*reff_2-B2*reff_1*reff_1*u);
    d2udr2=-reff_1*auu*(C+B2*reff_1*u);
    return A*u*reff;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type d3udr3)
  {
    real_type reff_1(std::exp(-C*r));
    real_type reff((1.0-reff_1)*OneOverC);
    real_type u(1.0/(1.0+B*reff));
    real_type auu(A*u*u);
    dudr=reff_1*auu;
    //d2udr=auu*(-C*reff_1*reff_2-B2*reff_1*reff_1*u);
    d2udr2=-reff_1*auu*(C+B2*reff_1*u);
    std::cerr << "Third derivative not imlemented for Pade functor.\n";
    return A*u*reff;
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

  inline void resetParameters(const opt_variables_type& active)
  {
    OneOverC=1.0/C;
    B2=2.0*B;
  }
};
}
#endif


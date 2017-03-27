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
    
    
#ifndef QMCPLUSPLUS_GAUSSIAN_H
#define QMCPLUSPLUS_GAUSSIAN_H
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
template<class T>
struct GaussianFunctor: public OptimizableFunctorBase
{
/// f(r) = A exp(-(r)^2/C^2)
///for smooth truncation f(r)->f(r)+f(2rc-r)-2f(rc)
  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A;
  real_type C;
  real_type RC;
  real_type c0,c1,c2,c3,c4;

  std::string ID_A;
  std::string ID_C;
  std::string ID_RC;
  bool optimizable;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  GaussianFunctor(real_type a=1.0, real_type c=1.0,real_type rc=1.0):ID_A("G_A"), ID_C("G_C"),ID_RC("G_RC"), optimizable(true)
  {
    A=a;
    C=c;
    RC=rc;
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new GaussianFunctor<T>(*this);
  }

  inline void reset()
  {
    c0 = -1.0/(C*C);
    c1 = 2*RC;
    c2 = 2*A*std::exp(-(RC*RC)/(C*C));
    c3 = -2.0*A/(C*C);
    c4 = 4.0*A/(C*C*C*C);
  }

  /** reset the internal variables.
  *
  * USE_resetParameters
   */
  void resetParameters(const opt_variables_type& active)
  {
    if (!optimizable)
      return;
    int ia=myVars.where(0);
    if (ia>-1)
      A=myVars[0]=active[ia];
    int ic=myVars.where(1);
    if (ic>-1)
      C=myVars[1]=active[ic];
//         int id=myVars.where(2);
//         if (id>-1) RC=myVars[2]=active[id];
    reset();
  }

  void checkInVariables(opt_variables_type& active)
  {
    if (optimizable)
      active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    if (optimizable)
      myVars.getIndex(active);
  }

  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    if (r<RC)
    {
      real_type expart = std::exp(c0*r*r);
      real_type r1=c1-r;
      real_type expartRC = std::exp(c0*r1*r1);
      return A*(expart+expartRC)+c2;
    }
    else
    {
      return 0.0;
    }
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @param d3udr3 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type& d3udr3)
  {
    std::cerr << "Third derivative not implemented for GaussianFunctor.\n";
    if (r<RC)
    {
      real_type expart = std::exp(c0*r*r);
      real_type r1=c1-r;
      real_type expartRC = std::exp(c0*r1*r1);
      dudr =  c3*(r*expart - r1*expartRC);
      d2udr2 = (c4*r*r+c3)*expart + (c4*r1*r1+c3)*expartRC;
      return A*(expart+expartRC)+c2;
    }
    else
    {
      dudr = 0.0;
      d2udr2 = 0.0;
      return 0.0;
    }
  }


  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if (r<RC)
    {
      real_type expart = std::exp(c0*r*r);
      real_type r1=c1-r;
      real_type expartRC = std::exp(c0*r1*r1);
      dudr =  c3*(r*expart - r1*expartRC);
      d2udr2 = (c4*r*r+c3)*expart + (c4*r1*r1+c3)*expartRC;
      return A*(expart+expartRC)+c2;
    }
    else
    {
      dudr = 0.0;
      d2udr2 = 0.0;
      return 0.0;
    }
  }

  /** return a value at r
   */
  real_type f(real_type r)
  {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  real_type df(real_type r)
  {
    if (r<RC)
    {
      real_type r1=(c1-r);
      return c3*(r*std::exp(c0*r*r) - r1*std::exp(c0*r1*r1));
    }
    else
      return 0.0;
  }
  inline bool
  evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    if ((r >= RC)||(!optimizable))
      return false;
    else
    {
      real_type rcmr =(c1-r);
      real_type rcmr2 =rcmr*rcmr;
      real_type r2   =r*r;
      real_type expart = std::exp(c0*r2);
      real_type expartRC = std::exp(c0*rcmr2);
      real_type Am1 = -2.0/(C*C);
      real_type cm2 = 1.0/(C*C);
      real_type cm3 = 1.0/(C*C*C);
      real_type cm4 = cm2*cm2;
      real_type cm5 = cm2*cm3;
      //Derivs of function
//
//  df_dA
      derivs[0][0]= (expart+expartRC)+c2*Am1 ;
//  df_dC
      derivs[1][0]= 2*A*( expart*r2 - 2.0*expartRC*rcmr2 + expartRC*rcmr2)*cm3 ;
//  df_dRC
//             derivs[2][0]= (2.0*c2*RC - 2.0*A*expartRC*rcmr)*cm2 ;
      //grad derivs
//  dr_df_dA
      derivs[0][1]= Am1*(r*expart - rcmr*expartRC) ;
//  dr_df_dC
      derivs[1][1]= 4.0*A*(expart*r*cm3*(1.0-r2*cm2)-expartRC*rcmr*cm3*(1.0-rcmr2*cm2)) ;
//  dr_df_dRC
//             derivs[2][1]= 4.0*A*expartRC*cm2*(1.0-2.0*rcmr2*cm2) ;
      //lap derivs
//  dr2_df_dA
      derivs[0][2]= Am1*((Am1*r2+1)*expart + (Am1*rcmr2 + 1)*expartRC);
//  dr2_df_dC
      derivs[1][2]= 4.0*A*cm3*(expart*(1.0-5.0*r2*cm2+2.0*r2*r2*cm4)+expartRC*(1.0-5.0*rcmr2*cm2+2.0*rcmr2*rcmr2*cm4)) ;
//  dr2_df_dRC
//             derivs[2][2]= 24*A*expartRC*rcmr*cm4 - 16.0*A*rcmr*rcmr2*expartRC*cm3*cm3;
      return true;
    }
  }


  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    RC = cutoff_radius;
    std::string opt("yes");
    OhmmsAttributeSet Tattrib;
    Tattrib.add(opt,"optimize");
    Tattrib.put(cur);
    optimizable = (opt=="yes");
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while (tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if (cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
//            std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if (aname == "A")
        {
          ID_A = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
          putContent(A,tcur);
        }
        else
          if (aname == "C")
          {
            ID_C = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
            putContent(C,tcur);
          }
          else
            if (aname == "RC")
            {
              ID_RC = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
              putContent(RC,tcur);
            }
      }
      tcur = tcur->next;
    }
    if (optimizable)
    {
      myVars.insert(ID_A,A);
      myVars.insert(ID_C,C);
//          myVars.insert(ID_RC,RC);
    }
    reset();
    return true;
  }
};

template<class T>
struct TruncatedShiftedGaussianFunctor: public OptimizableFunctorBase
{

  typedef typename OptimizableFunctorBase::real_type real_type;
  ///coefficients
  real_type A, B, C, RC;
  real_type c0,c1,c2,c3,c4,c5,c6;
  ///id of A
  std::string ID_A;
  ///id of B
  std::string ID_B;
  ///id of rc
  std::string ID_RC;
  std::string ID_C;
  ///constructor
  TruncatedShiftedGaussianFunctor(real_type a=1.0 , real_type b=1.0 , real_type c=1.0):ID_A("SG_A"), ID_B("SG_B"),ID_C("SG_C")
  {
    A=a;
    B=b;
    C=c;
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new TruncatedShiftedGaussianFunctor(*this);
  }

  /**
  *@brief reset the internal variables.
  */
  void reset()
  {
    c0 = -1.0/(C*C);
    c1 = std::exp(c0*B*B);
    c2 = A*c1/(C*C);
    c3 = -B*B*c2;
    c4 = 2.0*A*c0;
    c5 = 2.0*c2;
    c6 = 4.0*A*c0*c0;
  }

  /**@param r the distance
  @return \f$ u(r) = a/(1+br^2) \f$
   */
  inline real_type evaluate(real_type r)
  {
    if (r<RC)
    {
      r-=B;
      real_type r2=r*r;
      return A*(std::exp(c0*r2)-c1)+c2*r2+c3;
    }
    else
      return 0.0;
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @param d3udr3 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr,
                            real_type& d2udr2, real_type& d3udr3)
  {
    std::cerr << "Third derivative not implemented for GaussianFunctor.\n";
    return 0.0;
  }


  /**@param r the distance
  @param dudr return value  \f$ du/dr = -2abr/(1+br^2)^2 \f$
  @param d2udr2 return value  \f$ d^2u/dr^2 =
  -2ab(1-3br^2)/(1+br^2)^3 \f$
  @return \f$ u(r) = a/(1+br^2) \f$
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if (r<RC)
    {
      r-=B;
      real_type r2=r*r;
      real_type expart = std::exp(c0*r2);
      dudr = c4*r*expart+c5*r;
      d2udr2 = (c4+c6*r2)*expart+c5;
      return  A*(expart-c1)+c2*r2+c3;
    }
    else
    {
      dudr=0.0;
      d2udr2 = 0.0;
      return 0.0;
    };
  }

  real_type f(real_type r)
  {
    return evaluate(r);
  }

  real_type df(real_type r)
  {
    if (r<RC)
    {
      r-=B;
      return c4*r*std::exp(c0*r*r)+c5*r;
    }
    else
      return 0.0;
  }

  /** implements virtual function
  * @param cur xml node
   */
  bool put(xmlNodePtr cur)
  {
    RC = cutoff_radius;
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while (tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if (cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if (aname == "A")
        {
          ID_A = idname;
          putContent(A,tcur);
        }
        else
          if (aname == "B")
          {
            ID_B = idname;
            putContent(B,tcur);
          }
          else
            if (aname == "C")
            {
              ID_C = idname;
              putContent(C,tcur);
            }
      }
      tcur = tcur->next;
    }
    ///Sim cell too small
    if (B<0.5*cutoff_radius)
      RC=2.0*B;
    else
      return false;
    reset();
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    myVars.insert(ID_C,C);
    LOGMSG("  TruncatedShiftedGaussianFunctor Parameters ")
    LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B<<"  C (" << ID_C << ") =  " << C  << "  R_c (" << ID_RC << ") =  " << RC)
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

  /** reset the internal variables.
  *
  * USE_resetParameters
   */
  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if (ia>-1)
      A=active[ia];
    int ib=myVars.where(1);
    if (ib>-1)
      B=active[ib];
    int ic=myVars.where(2);
    if (ic>-1)
      C=active[ic];
    if (B>0.5*cutoff_radius)
    {
      B=0.5*cutoff_radius;
      RC=cutoff_radius;
    }
    else
      RC=2.0*B;
    reset();
  }
};

}
#endif

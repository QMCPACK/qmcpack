//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_MCMILLANJ2G_H
#define QMCPLUSPLUS_MCMILLANJ2G_H
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
/**
 *  McMillan Jastrow functor with cusp condition that matches HFD-BHE
 *  N = 64, density equilibrium
 *  Will be tuned up for general N, rho values
 */
template<class T>
struct McMillanJ2GFunctor: public OptimizableFunctorBase
{

  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A, B, R0, c1, c2;
  // setting R0 = 2.5 by default, in the constructor

  std::string ID_A;
  std::string ID_B;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  McMillanJ2GFunctor(real_type a=5.0, real_type b=5.7448):ID_A("McA"),ID_B("McB"),R0(2.5)
  {
    reset(a,b);
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new McMillanJ2GFunctor<T>(*this);
  }
  inline void reset()
  {
    reset(A,B);
    real_type Y;
    c1 = evaluate(R0, Y, c2);
    c2 = -Y/(2.0*R0*c1);
    c1 *= std::exp(-c2*R0*R0);
  }

  /** reset the internal variables.
  *
  * USE_resetParameters
   */
  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if(ia>-1)
      A=active[ia];
    int ib=myVars.where(1);
    if(ib>-1)
      B=active[ib];
    reset(A,B);
    real_type Y;
    c1 = evaluate(R0, Y, c2);
    c2 = -Y/(2.0*R0*c1);
    c1 *= std::exp(-c2*R0*R0);
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
  *@param a New Jastrow parameter a
   */
  inline void reset(real_type a, real_type b)
  {
    A = a;
    B = b;
  }

  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    // 64 particles, equil density for now
    return ( r < 2.5 ? c1*std::exp(-c2*r*r) :
             ( r < cutoff_radius ? std::pow(B,A)*(std::pow(r,-A) + std::pow(2.0*cutoff_radius-r,-A) - 2.0*std::pow(cutoff_radius,-A)) :
               0.0 ) );
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type wert;
    if (r < 2.5)
    {
      wert = c1*std::exp(-c2*r*r);
      dudr = -2.0*c1*c2*r*std::exp(-c2*r*r);
      d2udr2 = -2.0*c1*c2*std::exp(-c2*r*r) + 4.0*c1*c2*c2*r*r*std::exp(-c2*r*r);
      /*
      wert = 17318.0912945689*std::exp(-1.00041061270143*r*r);
      dudr = -34650.4046456380*r*std::exp(-1.00041061270143*r*r);
      d2udr2 = -34650.4046456380*std::exp(-1.00041061270143*r*r) + 69329.2650837902*r*r*std::exp(-1.00041061270143*r*r);
      */
    }
    else
      if (r < cutoff_radius)
      {
        real_type rd = 2.0*cutoff_radius - r;
        wert = std::pow(B,A)*(std::pow(r,-A) + std::pow(rd,-A) - 2.0*std::pow(cutoff_radius,-A));
        dudr = A*std::pow(B,A)*(std::pow(rd,-(A+1.0))-std::pow(r,-(A+1.0)));
        d2udr2 = (A+1.0)*A*std::pow(B,A)*(std::pow(rd,-(A+2.0))+std::pow(r,-(A+2.0)));
      }
      else
      {
        wert = 0.0;
        dudr = 0.0;
        d2udr2 = 0.0;
      }
    return wert;
    /*
    real_type rinv = 1.0/r, rdinv = 1.0/(2.0*cutoff_radius - r), wert1 = std::pow(B*rinv,A), wert2 = std::pow(B*rdinv,A);
    dudr = -A*(wert1*rinv - wert2*rdinv);
    d2udr2= (A+1.)*A*(wert1*rinv*rinv + wert2*rdinv*rdinv);
    return (wert1 + wert2 - 2.0*std::pow(B/cutoff_radius,A));
    */
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
    return ( r < 2.5 ? -2.0*c1*c2*r*std::exp(-c2*r*r) :
             ( r < cutoff_radius ? A*std::pow(B,A)*(std::pow(2.0*cutoff_radius-r,-(A+1.0))-std::pow(r,-(A+1.0))) :
               0.0 ) );
  }

  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
//            std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "a")
        {
          putContent(A,tcur);
        }
        else
          if(aname == "b")
          {
            ID_B = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
            putContent(B,tcur);
          }
      }
      tcur = tcur->next;
    }
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    reset(A,B);
    real_type Y;
    c1 = evaluate(R0, Y, c2);
    c2 = -Y/(2.0*R0*c1);
    c1 *= std::exp(c2*R0*R0);
//      std::cout << "ChangMo's test: " << R0 << ", " << cutoff_radius << ", " << c1 << ", " << c2 << std::endl << evaluate(R0-0.0001) << ", " << evaluate(R0+0.0001) << std::endl;
    return true;
  }


};
}
#endif

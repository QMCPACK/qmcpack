//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GAUSSIANSLATER_H
#define QMCPLUSPLUS_GAUSSIANSLATER_H
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
template<class T>
struct OpenGaussianSlaterFunctor: public OptimizableFunctorBase
{
/// f(r) = A exp(-(r)^2/C^2)
///for smooth truncation f(r)->f(r)+f(2rc-r)-2f(rc)
  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A;
  real_type B;
  std::string ID_A;
  std::string ID_B;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  OpenGaussianSlaterFunctor(real_type a=1.0, real_type b=1.0 ):ID_A("G_A"), ID_B("G_B")
  {
    A=a;
    B=b;
    reset();
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new OpenGaussianSlaterFunctor<T>(*this);
  }

  inline void reset()
  {
  }

  /** reset the internal variables.
  *
  * USE_resetParameters
   */
  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if (ia>-1)
      A=myVars[0]=active[ia];
    int ib=myVars.where(1);
    if (ib>-1)
      B=myVars[1]=active[ib];
    reset();
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    return A*r*r/(1.0+B*r);
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
    real_type part1 = 1.0/(1.0+B*r);
    dudr = A*r*(2.0+B*r)*part1*part1;
    d2udr2 = 2.0*A*part1*part1*part1;
    d3udr3 = -6.0*A*B*part1*part1*part1*part1;
    return A*r*r*part1;
  }


  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type part1 = 1.0/(1.0+B*r);
    dudr = A*r*(2.0+B*r)*part1*part1;
    d2udr2 = 2.0*A*part1*part1*part1;
    return A*r*r*part1;
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
    real_type part1 = 1.0/(1.0+B*r);
    return A*r*(2.0+B*r)*part1*part1;
  }

  inline bool
  evaluateDerivatives(real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    real_type u = 1.0/(1.0+B*r);
    derivs[0][0]= r*r*u;
    derivs[1][0]= -A*r*r*r*u*u;
    derivs[0][1]= r*(2.0+B*r)*u*u;
    derivs[1][1]= -A*r*r*(3.0+B*r)*u*u*u;
    derivs[0][2]= 2.0*u*u*u;
    derivs[1][2]= -6.0*A*r*u*u*u*u;
    return true;
  }


  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
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
          if (aname == "B")
          {
            ID_B = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
            putContent(B,tcur);
          }
      }
      tcur = tcur->next;
    }
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    reset();
    return true;
  }
};
}
#endif

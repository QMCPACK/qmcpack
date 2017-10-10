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
    
    
#ifndef QMCPLUSPLUS_RPA_FUNCTOR_H
#define QMCPLUSPLUS_RPA_FUNCTOR_H

#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
/** RPA Jastrow functor
 *
 * \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right], \f$
 * where \f$F\f$ satisfies the cusp conditions: \f$F_{\uparrow \uparrow} = \sqrt{2A}\f$
 * and \f$F_{\uparrow \downarrow} = \sqrt{A} \f$.
 * Prototype of the template parameter of TwoBodyJastrow
 * Reformulated by using B=1/F.
 */
template<class T>
struct RPAJastrowFunctor: public OptimizableFunctorBase<T>
{

  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;

  bool SameSpin;

  ///coefficients
  real_type A, B, AB, ABB;
  ///id of A for optimization
  std::string ID_A;
  /** constructor
   * @param a A coefficient
   * @param samespin boolean to indicate if this function is for parallel spins
   */
  RPAJastrowFunctor(bool samespin=true): SameSpin(samespin)
  {
    reset(1.0);
  }

  /** dummy constructor to handle referenced case
  */
  RPAJastrowFunctor(RPAJastrow<T>* func): SameSpin(true)
  {
    reset(1.0);
  }

  /** reset the internal variables.
   *@param a New Jastrow parameter a
   */
  void reset(real_type a)
  {
    A = a;
    resetInternals();
  }

  /** reset internal variables for the current value of A
   */
  void resetInternals()
  {
    real_type F=std::sqrt(std::abs(A));
    if(SameSpin)
      F*=std::sqrt(2.0);
    B=1.0/F;
    AB=A*B;
    ABB=AB*B;
  }

  /** evaluate the value at r
   * @param r the distance
   * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    return A/r*(1.0-std::exp(-B*r));
  }

  /**@param r the distance
    @param dudr first derivative
    @param d2udr second derivative
    @return the value
    */
  inline real_type evaluate(real_type r, T& dudr, T& d2udr2)
  {
    real_type rinv=1.0/r;
    real_type expbr=std::exp(-B*r);
    real_type u = A*rinv*(1.0-expbr);
    dudr=-rinv*(u-AB*expbr);
    d2udr2=-rinv*(2.0*dudr+ABB*expbr);
    return u;
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
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  void setDensity(real_type n)
  {
    A = std::pow(4.0*M_PI*n/3.0,-1.0/3.0);
    resetInternals();
  }

  /** Read in the parameter from the xml input file.
    @param cur current xmlNode from which the data members are reset
    @param vlist VarRegistry<T1> to which the variable A will be added for optimization
    */
  bool put(xmlNodePtr cur)
  {
    real_type Atemp=-100;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "A")
        {
          ID_A = idname;
          putContent(Atemp,tcur);
        }
      }
      tcur = tcur->next;
    }
    //only if Atemp is given
    if(Atemp>0)
      reset(Atemp);
    app_log() << "  RPA Jastrow A/r[1-exp(-r/F)] = (" << A << "," << 1.0/B << ")" << std::endl;
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    vlist[ID_A]=A;
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    typename OptimizableSetType::iterator it_a(optVariables.find(ID_A));
    if(it_a != optVariables.end())
    {
      A=(*it_a).second;
      resetInternals();
    }
  }
};





}

#endif


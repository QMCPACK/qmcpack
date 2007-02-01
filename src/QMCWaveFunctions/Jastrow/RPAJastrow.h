//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_RPAFUNCTION_H
#define QMCPLUSPLUS_RPAFUNCTION_H

#include "QMCWaveFunctions/OptimizableFunctorBase.h"

namespace qmcplusplus {
/** RPA Jastrow functional
 *
 * \f[ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right], \f]
 * where $F$ satisfies the cusp conditions: \f$F_{\uparrow \uparrow} = \sqrt{2A}\f$ 
 * and \f$F_{\uparrow \downarrow} = \sqrt{A} \f$.
 * Prototype of the template parameter of TwoBodyJastrow 
 * Reformulated by using B=1/F.
 */
template<class T>
struct RPAJastrow: public OptimizableFunctorBase<T> {

  typedef typename OptimizableFunctorBase<T>::real_type real_type;

  bool SameSpin;

  ///coefficients
  T A, B, AB, ABB;
  ///id of A for optimization
  string ID_A;
  /** constructor
   * @param a A coefficient
   * @param samespin boolean to indicate if this function is for parallel spins
   */
  RPAJastrow(bool samespin=true): SameSpin(samespin) {reset(1.0);}

  /** dummy constructor to handle referenced case
   */
  RPAJastrow(RPAJastrow<T>* func): SameSpin(true) {
    reset(1.0);
  }

  /** reset the internal variables.
   */
  inline void reset() {
    T F=std::sqrt(std::abs(A));
    if(SameSpin) F*=std::sqrt(2.0);
    B=1.0/F;
    AB=A*B;
    ABB=AB*B;
  }

  /** reset the internal variables.
   *@param a New Jastrow parameter a 
   */
  void reset(T a) {
    A = a;
    reset();
  }

  /** evaluate the value at r
   * @param r the distance
   * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline T evaluate(T r) {
    return A/r*(1.0-std::exp(-B*r));
  }

  /**@param r the distance
     @param dudr first derivative
     @param d2udr second derivative
     @return the value
  */
  inline T evaluate(T r, T& dudr, T& d2udr2) {
    T rinv=1.0/r;
    T expbr=std::exp(-B*r);
    T u = A*rinv*(1.0-expbr);
    dudr=-rinv*(u-AB*expbr);
    d2udr2=-rinv*(2.0*dudr+ABB*expbr);
    return u;
  }

  /** return a value at r
   */
  real_type f(real_type r) {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  real_type df(real_type r) {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  void setDensity(real_type n) {
    A = std::pow(4.0*M_PI*n/3.0,-1.0/3.0);
    reset();
  }

  /** Read in the parameter from the xml input file.
   @param cur current xmlNode from which the data members are reset
   @param vlist VarRegistry<T1> to which the variable A will be added for optimization
  */
  bool put(xmlNodePtr cur){
    T Atemp=-100;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL) {
      //@todo Var -> <param(eter) role="opt"/>
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var") {
	string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	if(aname == "A") {
	  ID_A = idname;
	  putContent(Atemp,tcur);
	} 
      }
      tcur = tcur->next;
    }
    //only if Atemp is given
    if(Atemp>0) reset(Atemp);
    app_log() << "  RPA Jastrow A/r[1-exp(-r/F)] = (" << A << "," << 1.0/B << ")" << std::endl;
    return true;
  }

  void addOptimizables( VarRegistry<real_type>& vlist)
  {
    vlist.add(ID_A,&A,1);
  }
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


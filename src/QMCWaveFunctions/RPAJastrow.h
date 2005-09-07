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
#ifndef OHMMS_QMC_RPAFUNCTION_H
#define OHMMS_QMC_RPAFUNCTION_H
#include "OhmmsData/libxmldefs.h"
#include "Optimize/VarList.h"

/**RPA Jastrow functional
 *
 * \f[ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right], \f]
 * where $F$ satisfies the cusp conditions: \f$F_{\uparrow \uparrow} = \sqrt{2A}\f$ 
 * and \f$F_{\uparrow \downarrow} = \sqrt{A} \f$.
 * Prototype of the template parameter of TwoBodyJastrow 
 * Reformulated by using B=1/F.
 */
template<class T>
struct RPAJastrow {

  bool SameSpin;

  ///coefficients
  T A, B, AB, ABB;

  /** constructor
   * @param a A coefficient
   * @param samespin boolean to indicate if this function is for parallel spins
   */
  RPAJastrow(T a=1.0, bool samespin=true): SameSpin(samespin) {reset(a);}

  /**
   *@brief reset the internal variables.
   */
  inline void reset() {
    T F=sqrt(abs(a));
    if(SameSpin) F*=sqrt(2.0);
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

  /**@param r the distance
     @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
  */
  inline T evaluate(T r) {
    return A/r*(1.0-exp(-B*r));
  }

  /**@param r the distance
     @param dudr first derivative
     @param d2udr second derivative
     @return the value
  */
  inline T evaluate(T r, T& dudr, T& d2udr2) {
    T rinv=1.0/r;
    T expbr=exp(-B*r);
    T u = A*rinv*(1.0-expbr);
    dudr=-rinv*(u-AB*expbr);
    d2udr2=-rinv*(2.0*dudr+ABB*expbr);
    return u;
  }

  /**@param cur current xmlNode from which the data members are reset
   @param vlist VarRegistry<T1> to which the variable A will be
   added for optimization
   @brief T1 is the type of VarRegistry, typically double.  Read 
   in the parameter from the xml input file.
  */
  template<class T1>
  void put(xmlNodePtr cur, VarRegistry<T1>& vlist){
    T Atemp;
    string ida;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL) {
      //@todo Var -> <param(eter) role="opt"/>
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var") {
	string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	if(aname == "A") {
	  ida = idname;
	  putContent(Atemp,tcur);
	} 
      }
      tcur = tcur->next;
    }
    reset(Atemp);
    vlist.add(ida,&A,1);
    XMLReport("Jastrow A/r[1-exp(-r/F)] = (" << A << "," << 1.0/B << ")") 
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


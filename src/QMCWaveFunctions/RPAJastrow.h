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

/**class RPA Jastrow functional
 *@brief \f[ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right], \f]
 * where $F$ satisfies the cusp conditions: \f$F_{\uparrow \uparrow} = 
  \sqrt{2A}\f$ and \f$F_{\uparrow \downarrow} = \sqrt{A} \f$.
 * Prototype of the template parameter of TwoBodyJastrow 
 */
template<class T>
struct RPAJastrow {

  ///coefficients
  T A, A2, B, Finv;

  ///constructor
  RPAJastrow(T a=1.0, T b=1.0) {reset(a,b);}

  /**
   *@brief reset the internal variables.
   */
  inline void reset() {
    A2=2.0*A; Finv=1.0/sqrt(A*B);
  }

  /**
   *@param a New Jastrow parameter a 
   *@brief reset the internal variables.
   */
  void reset(T a, T b) {
    A=a; B=b; A2=2.0*A; Finv=1.0/sqrt(A*B);
  }
  /**@param r the distance
     @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
  */
  inline T evaluate(T r) {
    return A/r*(1.0-exp(-r*Finv));
  }

  /**@param r the distance
     @param dudr return value 
     \f$ du/dr = -\frac{A}{r^2}\left[1-\exp(-\frac{r}{F}) -
     \frac{r}{F}\exp(-\frac{r}{F})\right]\f$
     @param d2udr2 return value  \f$ d^2u/dr^2 =
     \frac{2A}{r^3}\left[1-\exp(-\frac{r}{F}) -
     \frac{r}{F}\exp(-\frac{r}{F}) -
     \frac{r^2}{2F^2}\exp(-\frac{r}{F})\right]\f$
     @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right] \f$ 
  */
  inline T evaluate(T r, T& dudr, T& d2udr2) {
    T rinv = 1.0/r;
    T rinv2 = rinv*rinv;
    T rinv3 = rinv*rinv2;
    T rFinv = r*Finv;
    T e_minus_rFinv = exp(-rFinv);
    T u = 1.0-e_minus_rFinv;//1-exp(-r/F)
    dudr = u-rFinv*e_minus_rFinv;//1-exp(-r/F)-r/F*exp(-r/F)
    //2A/r^3*[1-exp(-r/F)-r/F*exp(-r/F)-r^2/(2F^2)*exp(-r/F)]
    d2udr2 = A2*rinv3*(dudr-0.5*rFinv*rFinv*e_minus_rFinv);
    dudr *= -A*rinv2;//multiply by -A/r^2
    return A*rinv*u;//multiply by A/r
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
    T Btemp = 1.0;
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
	} else if(aname == "B") {
	  putContent(Btemp,tcur);
	}
      }
      tcur = tcur->next;
    }
    reset(Atemp,Btemp);
    vlist.add(ida,&A,1);
    XMLReport("Jastrow A/r[1-exp(-r/F)] = (" << A << "," << sqrt(A*B) << ")") 
      }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


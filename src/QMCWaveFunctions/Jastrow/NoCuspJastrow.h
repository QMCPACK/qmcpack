//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_NOCUSPJASTROW_H
#define QMCPLUSPLUS_NOCUSPJASTROW_H
#include "QMCWaveFunctions/OptimizableFunctorBase.h"

namespace qmcplusplus {
/**Class NoCusp functional
 *@brief \f[ u(r) = \frac{a}{1+br^2} \f]
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct NoCuspJastrow: public OptimizableFunctorBase<T> {

  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  ///coefficients
  real_type A, B, AB2;
  string ID_A;
  string ID_B;

  ///constructor
  NoCuspJastrow(real_type a=1.0, real_type b=1.0) {reset(a,b);}

  /**
   *@brief reset the internal variables.
   */
  inline void reset() {
    AB2 = 2.0*A*B;
  }

  /**
   *@param a NoCusp Jastrow parameter a 
   *@param b NoCusp Jastrow parameter b 
   *@brief reset the internal variables.
   */
  void reset(real_type a, real_type b) {
    A=a; B=b; AB2=2.0*a*b;
  }
  /**@param r the distance
     @return \f$ u(r) = a/(1+br^2) \f$
  */
  inline real_type evaluate(real_type r) {
    return A/(1.0+B*r*r);
  }

  /**@param r the distance
     @param dudr return value  \f$ du/dr = -2abr/(1+br^2)^2 \f$
     @param d2udr2 return value  \f$ d^2u/dr^2 =
     -2ab(1-3br^2)/(1+br^2)^3 \f$
     @return \f$ u(r) = a/(1+br^2) \f$
  */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
    real_type rsq = r*r;
    real_type u = 1.0/(1.0+B*rsq);
    real_type usq = u*u;
    dudr = -AB2*r*usq;
    d2udr2 = -AB2*usq*u*(1.0-3.0*B*rsq);
    return A*u;
  }

  real_type f(real_type r) {
    return evaluate(r);
  }

  real_type df(real_type r) {
    real_type dudr,d2udr2;
    real_type res=evaluate(r,dudr,d2udr2);
    return dudr;
  }

  /** implements virtual function
   * @param cur xml node
   */
  bool put(xmlNodePtr cur)
  {
    real_type Atemp,Btemp;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL) {
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var") {
	string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	if(aname == "A") {
	  ID_A = idname;
	  putContent(Atemp,tcur);
	} else if(aname == "B"){
	  ID_B = idname;
	  putContent(Btemp,tcur);
	}
      }
      tcur = tcur->next;
    }
    reset(Atemp,Btemp);
    LOGMSG("  NoCuspJastrow Parameters ")
    LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B)
    return true;
  }

  void addOptimizables( VarRegistry<real_type>& vlist){
    vlist.add(ID_A,&A,1);
    vlist.add(ID_B,&B,1);
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


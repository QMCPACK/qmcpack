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
#ifndef OHMMS_QMC_NOCUSPJASTROW_H
#define OHMMS_QMC_NOCUSPJASTROW_H
#include "OhmmsData/libxmldefs.h"
#include "Optimize/VarList.h"

/**Class NoCusp functional
 *@brief \f[ u(r) = \frac{a}{1+br^2} \f]
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct NoCuspJastrow {

  ///coefficients
  T A, B, AB2;

  ///constructor
  NoCuspJastrow(T a=1.0, T b=1.0) {reset(a,b);}

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
  void reset(T a, T b) {
    A=a; B=b; AB2=2.0*a*b;
  }
  /**@param r the distance
     @return \f$ u(r) = a/(1+br^2) \f$
  */
  inline T evaluate(T r) {
    return A/(1.0+B*r*r);
  }

  /**@param r the distance
     @param dudr return value  \f$ du/dr = -2abr/(1+br^2)^2 \f$
     @param d2udr2 return value  \f$ d^2u/dr^2 =
     -2ab(1-3br^2)/(1+br^2)^3 \f$
     @return \f$ u(r) = a/(1+br^2) \f$
  */
  inline T evaluate(T r, T& dudr, T& d2udr2) {
    T rsq = r*r;
    T u = 1.0/(1.0+B*rsq);
    T usq = u*u;
    dudr = -AB2*r*usq;
    d2udr2 = -AB2*usq*u*(1.0-3.0*B*rsq);
    return A*u;
  }

  /**@param cur current xmlNode from which the data members are reset
   @param vlist VarRegistry<T1> to which the Pade variables A and B
   are added for optimization
   @brief T1 is the type of VarRegistry, typically double.  Read 
   in the Pade parameters from the xml input file.
  */
  template<class T1>
  void put(xmlNodePtr cur, VarRegistry<T1>& vlist){
    T Atemp,Btemp;
    string ida, idb;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL) {
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var") {
	string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	if(aname == "A") {
	  ida = idname;
	  putContent(Atemp,tcur);
	} else if(aname == "B"){
	  idb = idname;
	  putContent(Btemp,tcur);
	}
      }
      tcur = tcur->next;
    }
    reset(Atemp,Btemp);
    vlist.add(ida,&A,1);
    vlist.add(idb,&B,1);
    XMLReport("Jastrow Parameters = (" << A << "," << B << ")") 
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


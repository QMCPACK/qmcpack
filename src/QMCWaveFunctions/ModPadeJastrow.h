//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_MODIFIED_PADEFUNCTION_H
#define QMCPLUSPLUS_MODIFIED_PADEFUNCTION_H

#include "QMCWaveFunctions/JastrowFunctorBase.h"

namespace qmcplusplus {
/** ModPade Jastrow functional
 *
 * \f[ u(r) = \frac{1}{2*A}\left[1-\exp(-A*r)\right], \f]
 */
template<class T>
struct ModPadeJastrow: public JastrowFunctorBase<T> {

  typedef typename JastrowFunctorBase<T>::real_type real_type;
  typedef typename JastrowFunctorBase<T>::value_type value_type;

  ///coefficients
  T A;
  T B;
  T Zeff;

  T Coeff;
  T mAB;

  /** constructor
   * @param a A coefficient
   * @param samespin boolean to indicate if this function is for parallel spins
   */
  ModPadeJastrow(T a=-0.5, T b=1): Zeff(1.0) {reset(a,b);}

  /** reset the internal variables.
   */
  inline void reset() {
    Coeff=A/B;
    mAB = -A*B;
  }

  /** reset the internal variables.
   *@param a New Jastrow parameter a 
   */
  inline void reset(T a, T b) {
    A = a;
    B = b;
    reset();
  }

  /** evaluate the value at r
   * @param r the distance
   * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline T evaluate(T r) {
    return Coeff*(1.0-exp(-B*r));
  }

  /**@param r the distance
     @param dudr first derivative
     @param d2udr second derivative
     @return the value
  */
  inline T evaluate(T r, T& dudr, T& d2udr2) {
    T expmar=exp(-B*r);
    dudr=A*expmar;
    d2udr2=mAB*expmar;
    return Coeff*(1.0-expmar);
  }

  /** return a value at r
   */
  value_type f(real_type r) {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  value_type df(real_type r) {
    return A*exp(-B*r);
  }

  /** Read in the parameter from the xml input file.
   * @param cur current xmlNode from which the data members are reset
   * @param vlist VarRegistry<T1> to which the variable A will be added for optimization
   */
  void put(xmlNodePtr cur, VarRegistry<real_type>& vlist){
    Zeff=1.0;
    string idb;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL) {
      //@todo Var -> <param(eter) role="opt"/>
      string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var") {
	string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	if(aname == "A") {
	  putContent(A,tcur);
        } else if(aname == "B") {
	  idb = idname;
	  putContent(B,tcur);
	} else if(aname == "Z") {
	  putContent(Zeff,tcur);
        }
      }
      tcur = tcur->next;
    }
    reset(A,B);
    vlist.add(idb,&A,1);
    XMLReport("Jastrow A/B [1-exp(-B*r)] = (" << A << "," << B << ")") 
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


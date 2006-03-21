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
#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_SLATERBASISSET_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_SLATERBASISSET_H

#include <cmath>
#include "Numerics/SlaterTypeOrbital.h"
#include "OhmmsData/AttributeSet.h"

template<class T>
struct SlaterCombo: public RadialOrbitalBase<T> {

  typedef T value_type;
  typedef GenericSTO<T> Component_t;

  int L;
  bool Normalized;

  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::vector<xmlNodePtr> InParam;
  std::vector<Component_t> sset;
  T Y, dY, d2Y;

  explicit 
    SlaterCombo(int l=0, 
                bool normalized=true,
                const char* node_name="radfunc",
                const char* exp_name="exponent", 
                const char* c_name="contraction");

  ~SlaterCombo(){ }

  void reset();

  inline value_type f(value_type r) const {
    value_type res=0;
    typename std::vector<Component_t>::const_iterator it(sset.begin());
    typename std::vector<Component_t>::const_iterator it_end(sset.end());
    while(it != it_end) {
      res += (*it).f(r); ++it;
    }
    return res;
  }

  inline value_type df(value_type r) const {
    value_type res=0;
    typename std::vector<Component_t>::const_iterator it(sset.begin());
    typename std::vector<Component_t>::const_iterator it_end(sset.end());
    while(it != it_end) {
      res += (*it).df(r); ++it;
    }
    return res;
  }

  inline value_type evaluate(T r, T rinv) {
    Y=0.0;dY=0.0;d2Y=0.0;
    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
    while(it != it_end) {
      Y+=(*it).evaluate(r,rinv); ++it;
    }
    return Y;
  }

  inline void evaluateAll(T r, T rinv) {
    Y=0.0;dY=0.0;d2Y=0.0;
    T du, d2u;
    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
    while(it != it_end) {
      Y+=(*it).evaluate(r,rinv,du,d2u); dY+=du; d2Y+=d2u;
      ++it;
    }
  }

//  inline value_type evaluate(T r, T rinv, T& drnl, T& d2rnl) {
//    Y=0.0;drnl=0.0;d2rnl=0.0;
//    T du, d2u;
//    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
//    while(it != it_end) {
//      Y+=(*it).evaluate(r,rinv,du,d2u); dY+=du; d2Y+=d2u;
//      ++it;
//    }
//    return Y;
//  }

  bool putBasisGroup(xmlNodePtr cur);
};

template<class T>
SlaterCombo<T>::SlaterCombo(int l, bool normalized,
    const char* node_name, const char* exp_name, const char* c_name): 
  L(l), Normalized(normalized), 
  nodeName(node_name), expName(exp_name), coeffName(c_name)
{
}

template<class T>
void SlaterCombo<T>::reset() {
 //Add reset function
}

template<class T>
bool SlaterCombo<T>::putBasisGroup(xmlNodePtr cur) {
  cur = cur->children;
  while(cur != NULL) {
    string cname((const char*)cur->name);
    if(cname == "radfunc" || cname == "Rnl") {
      T zeta(1.0),c(1.0);
      int qN=1;
      OhmmsAttributeSet radAttrib;
      radAttrib.add(zeta,expName); radAttrib.add(zeta,"alpha");
      radAttrib.add(c,coeffName); radAttrib.add(zeta,"c");
      radAttrib.add(qN,"node"); radAttrib.add(qN,"n");
      radAttrib.put(cur);
      if(Normalized) {
        //z is not right
        sset.push_back(Component_t(qN-1,zeta,c));
      } else {
        STONorm<T> anorm(qN);
        //multiply a normalization factor to the contraction factor
        //anorm(n,\zeta) = 1/\sqrt((2n+2)!/(2*\zeta)^{2*n+3))
        c *= anorm(qN-1,zeta);
        sset.push_back(Component_t(qN-L-1,zeta,c));
      }
      LOGMSG(" Slater Component (n,zeta,c)= " << qN << " " << zeta << " " << c)
    }
    cur=cur->next;
  }
  reset();
  return true;
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

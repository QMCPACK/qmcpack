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
#ifndef OHMMS_QMC_RADIALGRIDFUNCTOR_SLATERBASISSET_H
#define OHMMS_QMC_RADIALGRIDFUNCTOR_SLATERBASISSET_H

#include "Numerics/SlaterTypeOrbital.h"
#include "OhmmsData/AttributeSet.h"

template<class T>
struct SlaterCombo: public RadialOrbitalBase<T> {

  typedef T value_type;
  T L;
  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::vector<xmlNodePtr> InParam;
  std::vector<GenericSTO> sset;

  explicit 
    SlaterCombo(int l=0, 
                const char* node_name="radfunc",
                const char* exp_name="exponent", 
                const char* c_name="contraction");

  ~SlaterCombo(){ }

  void reset();

  inline value_type f(value_type r) const {
    value_type res=0;
    typename std::vector<GenericSTO>::const_iterator it(gset.begin());
    typename std::vector<GenericSTO>::const_iterator it_end(gset.end());
    while(it != it_end) {
      res += (*it).f(r); ++it;
    }
    return res;
  }

  inline value_type df(value_type r) const {
    value_type res=0;
    typename std::vector<GenericSTO>::const_iterator it(gset.begin());
    typename std::vector<GenericSTO>::const_iterator it_end(gset.end());
    while(it != it_end) {
      res += (*it).df(r); ++it;
    }
    return res;
  }

  bool put(xmlNodePtr cur);
};

template<class T>
SlaterCombo<T>::SlaterCombo(int l,
    const char* node_name, const char* exp_name, const char* c_name): 
  Normalized(normalized), nodeName(node_name),
  expName(exp_name), coeffName(c_name)
{
  L = static_cast<T>(l);
}

template<class T>
void SlaterCombo<T>::reset() {
  int n=gset.size();
  while(n<InParam.size()) {
    gset.push_back(BasicGaussian());
    n++;
  }

  OhmmsAttributeSet radAttrib;
  T zeta(1.0),c(1.0);
  int qN=0;

  OhmmsAttributeSet radAttrib;
  radAttrib.add(zeta,expName); radAttrib.add(zeta,"alpha");
  radAttrib.add(c,coeffName); radAttrib.add(zeta,"c");
  radAttrib.add(qN,"node"); radAttrib.add(qN,"n");

  for(int i=0; i<InParam.size(); i++) {
    radAttrib.put(InParam[i]);
    LOGMSG(" Slater Component (n,zeta,c)= " << qN << " " << zeta << " " << c)
    STONorm<T> anorm(qN);
    InFunc.add(new GenericSTO<T>(qN-L-1,zeta,c*anorm(n-1,zeta)));
    gset[i].reset(zeta,c);
  }
}

template<class T>
bool SlaterCombo<T>::putBasisGroup(xmlNodePtr cur) {
  cur = cur->children;
  while(cur != NULL) {
    string cname((const char*)cur->name);
    if(cname == "radfunc" || cname == "Rnl") {
      put(cur);
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

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
#ifndef OHMMS_QMC_RADIALGRIDFUNCTOR_GAUSSIANBASISSET_H
#define OHMMS_QMC_RADIALGRIDFUNCTOR_GAUSSIANBASISSET_H
#include "Numerics/RadialOrbitalBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

template<class T>
struct GaussianCombo: public RadialOrbitalBase<T> {

  typedef T value_type;

  struct BasicGaussian {
    T Sigma;
    T Coeff;
    T CoeffP;
    BasicGaussian(): Sigma(1.0), Coeff(1.0) { } 
    inline BasicGaussian(T sig, T c) { 
      reset(sig,c);
    } 
    void reset(T sig, T c){
      Sigma = sig; Coeff = c;
      CoeffP = -2.0*Sigma*Coeff;
    }
    inline void setgrid(T r) { }
    inline T f(T r2) const {
      return Coeff*exp(-Sigma*r2);
    }
    inline T df(T r, T r2) const {
      return CoeffP*r*exp(-Sigma*r2);
    }
  };

  ///Boolean
  bool Normalized;
  T L;
  T NormL;
  T NormPow;
  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::vector<xmlNodePtr> InParam;
  std::vector<BasicGaussian> gset;

  explicit 
    GaussianCombo(int l=0, bool normalized=false, 
                  const char* node_name="radfunc",
                  const char* exp_name="exponent", 
                  const char* c_name="contraction");

  ~GaussianCombo(){ }

  void reset();

  inline value_type f(value_type r) const {
    value_type res=0;
    value_type r2 = r*r;
    typename std::vector<BasicGaussian>::const_iterator it(gset.begin());
    typename std::vector<BasicGaussian>::const_iterator it_end(gset.end());
    while(it != it_end) {
      res += (*it).f(r2); ++it;
    }
    return res;
  }
  inline value_type df(value_type r) const {
    value_type res=0;
    value_type r2 = r*r;
    typename std::vector<BasicGaussian>::const_iterator it(gset.begin());
    typename std::vector<BasicGaussian>::const_iterator it_end(gset.end());
    while(it != it_end) {
      res += (*it).df(r,r2); ++it;
    }
    return res;
  }

  bool put(xmlNodePtr cur);

  bool putBasisGroup(xmlNodePtr cur);

  int DFactorial(int l) {
    return (l==1)? 1: DFactorial(l-2);
  }
};

template<class T>
GaussianCombo<T>::GaussianCombo(int l, bool normalized,
    const char* node_name, const char* exp_name, const char* c_name): 
  Normalized(normalized), nodeName(node_name),
  expName(exp_name), coeffName(c_name)
{
  L = static_cast<T>(l);
  //Everything related to L goes to NormL and NormPow
  const T pi = 4.0*atan(1.0);
  NormL = pow(2,L+1)*sqrt(2.0/static_cast<T>(DFactorial(2*l+1)))*pow(2.0/pi,0.25);
  NormPow = 0.5*(L+1.0)+0.25;
}

template<class T>
bool GaussianCombo<T>::put(xmlNodePtr cur) {
  InParam.push_back(cur);
  return true;
}

template<class T>
void GaussianCombo<T>::reset() {
  int n=gset.size();
  while(n<InParam.size()) {
    gset.push_back(BasicGaussian());
    n++;
  }

  T alpha(1.0),c(1.0);
  OhmmsAttributeSet radAttrib;
  radAttrib.add(alpha,expName); 
  radAttrib.add(c,coeffName);

  for(int i=0; i<InParam.size(); i++) {
    radAttrib.put(InParam[i]);
    //const xmlChar* aptr = xmlGetProp(InParam[i],(const xmlChar *)expName.c_str());
    //const xmlChar* cptr = xmlGetProp(InParam[i],(const xmlChar *)coeffName.c_str());
    //if(aptr == 0) {
    //  cout << "Exponent is not found. Skip this" << endl;
    //  continue;
    //}
    //T alpha = atof((const char*)aptr);
    //T c(1.0);
    //if(cptr) c = atof((const char*)cptr);
    //get the normalization factor
    if(!Normalized) c *= NormL*pow(alpha,NormPow); 
    LOGMSG(" Gaussian exponent = " << alpha << " contraction=" << c)
    gset[i].reset(alpha,c);
  }
}

template<class T>
bool GaussianCombo<T>::putBasisGroup(xmlNodePtr cur) {
  cur = cur->children;
  while(cur != NULL) {
    string cname((const char*)cur->name);
    if(cname == "radfunc") {
      put(cur);
    }
    cur=cur->next;
  }
  reset();
//  if(nptr) nprim = atoi((const char*)nptr);
//
//  int nprim=1;
//  const xmlChar* nptr = xmlGetProp(cur,(const xmlChar*)"nprim");
//  if(nptr) nprim = atoi((const char*)nptr);
//  T minL(L);
//  const xmlChar* minLptr = xmlGetProp(cur,(const xmlChar*)"minL");
//  if(minLptr) minL = atof((const char*)minLptr);
//  vector<T> alpha(nprim),c(nprim,1.0);
//  cur = cur->children;
//  while(cur != NULL) {
//    string cname((const char*)cur->name);
//    if(cname == "exponents") {
//      putContent(alpha,cur);
//    } else if(cname == "contraction") {
//      if(minL == L) {
//        putContent(c,cur);
//        if(!Normalized) 
//          for(int i=0; i<nprim; i++) c[i] *= NormL*pow(alpha[i],NormPow); 
//      }
//      ++minL;
//    }
//    cur = cur->next;
//  }
//  int i=0;
//  while(i<nprim) {
//    gset.push_back(BasicGaussian());
//    gset[i].reset(alpha[i],c[i]);
//    ++i;
//  }
  return true;
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

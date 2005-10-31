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
#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANBASISSET_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANBASISSET_H
#include "Numerics/RadialOrbitalBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

template<class T>
struct GaussianCombo: public RadialOrbitalBase<T> {

  typedef T value_type;
  T Y, dY, d2Y;

  struct BasicGaussian {
    T Sigma;
    T Coeff;

    T MinusSigma;
    T CoeffP;
    T CoeffPP;
    BasicGaussian(): Sigma(1.0), Coeff(1.0) { } 

    inline BasicGaussian(T sig, T c) { 
      reset(sig,c);
    } 

    inline void reset(T sig, T c){
      Sigma = sig; 
      MinusSigma=-sig;
      Coeff = c;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
    }

    inline void reset() {
      MinusSigma=-Sigma;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
    }

    inline void setgrid(T r) { }

    inline T f(T rr) const {
      return Coeff*exp(MinusSigma*rr);
    }
    inline T df(T r, T rr) const {
      return CoeffP*r*exp(MinusSigma*rr);
    }
    inline T evaluate(T r, T rr, T& du, T& d2u) {
      T v=exp(MinusSigma*rr);
      du += CoeffP*r*v;
      d2u += (CoeffP+CoeffPP*rr)*v;
      return Coeff*v;
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

  /** return the number Gaussians
   */
  inline int size() const {
    return gset.size();
  }

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

  inline value_type evaluate(T r, T rinv) {
    Y=0.0;
    value_type rr = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin()),it_end(gset.end());
    while(it != it_end) {
      Y+=(*it).f(rr); ++it;
    }
    return Y;
  }

  inline void evaluateAll(T r, T rinv) {
    Y=0.0;dY=0.0;d2Y=0.0;
    value_type rr = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin()),it_end(gset.end());
    while(it != it_end) {
      Y+=(*it).evaluate(r,rr,dY,d2Y); ++it;
    }
  }

  //inline value_type evaluate(T r, T rinv, T& drnl, T& d2rnl) {
  //  Y=0.0;drnl=0.0;d2rnl=0.0;
  //  T du, d2u;
  //  typename std::vector<BasicGaussian>::iterator 
  //    it(sset.begin()),it_end(sset.end());
  //  while(it != it_end) {
  //    Y+=(*it).evaluate(r,rinv,du,d2u); dY+=du; d2Y+=d2u;
  //    ++it;
  //  }
  //  return Y;
  //}


  bool put(xmlNodePtr cur);

  bool putBasisGroup(xmlNodePtr cur);

  /**  double factorial of num
   * @param num integer to be factored
   * @return num!!
   *
   * \if num == odd,
   * \f$ num!! = 1\cdot 3\cdot ... \cdot num-2 \cdot num\f$
   * \else num == even,
   * \f$ num!! = 2\cdot 4\cdot ... \cdot num-2 \cdot num\f$
   */
  int DFactorial(int num) {
    return (num<2)? 1: num*DFactorial(num-2);
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
  if(InParam.empty()) {
    for(int i=0; i<gset.size(); i++) gset[i].reset();
  } else {
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
      if(!Normalized) c *= NormL*pow(alpha,NormPow); 
      LOGMSG(" Gaussian exponent = " << alpha << " contraction=" << c)
      gset[i].reset(alpha,c);
    }
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
  return true;
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

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
#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANTIMESRN_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANTIMESRN_H
#include "Numerics/RadialOrbitalBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

template<class T>
struct GaussianTimesRN: public RadialOrbitalBase<T> {

  typedef T value_type;
  T Y, dY, d2Y;

  struct BasicGaussian {
    T Sigma;
    T Coeff;
    int Power;

    T MinusSigma;
    T CoeffP;
    T CoeffPP;
    T PowerC;
    BasicGaussian(): Sigma(1.0), Coeff(1.0), Power(0) { } 

    inline BasicGaussian(T sig, T c, int p=0) { 
      reset(sig,c,p);
    } 

    inline void reset(T sig, T c, int p=0){
      Sigma = sig; 
      Coeff = c;
      Power = p;
      LOGMSG(" Gaussian exponent = " << sig << " contraction=" << c << " power = " << Power)
      reset();
    }

    inline void reset() {
      MinusSigma=-Sigma;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
      PowerC=static_cast<T>(Power)*Coeff;
    }

    inline void setgrid(T r) { }

    inline T f(T r, T rr) const {
      if(Power==0) 
        return Coeff*exp(MinusSigma*rr);
      else if(Power==1)
        return r*Coeff*exp(MinusSigma*rr);
      else 
        return std::pow(r,Power)*Coeff*exp(MinusSigma*rr);
    }

    inline T df(T r, T rr) const {
      if(Power==0) 
        return CoeffP*r*exp(MinusSigma*rr);
      else if(Power==1) 
        return (Coeff+CoeffP*r)*exp(MinusSigma*rr);
      else  {
        return 
          std::pow(r,Power-1)*(PowerC+CoeffP*rr)*exp(MinusSigma*rr);
      }
    }

    inline T evaluate(T r, T rr, T& du, T& d2u) {
      T v=exp(MinusSigma*rr);
      if(Power==0) {
        du += CoeffP*r*v;
        d2u += (CoeffP+CoeffPP*rr)*v;
        return Coeff*v;
      } else {
        return std::pow(r,Power)*Coeff*v;
      }
    }
  };

  int basePower;
  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::string  powerName;
  std::vector<xmlNodePtr> InParam;
  std::vector<BasicGaussian> gset;

  explicit 
    GaussianTimesRN(const char* node_name="radfunc",
                  const char* exp_name="exponent", 
                  const char* c_name="contraction",
                  const char* p_name="power"): 
      basePower(0),
      nodeName(node_name), expName(exp_name), 
      coeffName(c_name), powerName(p_name)
    {}

  ~GaussianTimesRN(){ }

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
      res += (*it).f(r,r2); ++it;
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
      Y+=(*it).f(r,rr); ++it;
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

  bool put(xmlNodePtr cur);

  bool putBasisGroup(xmlNodePtr cur);

};

template<class T>
bool GaussianTimesRN<T>::put(xmlNodePtr cur) {
  InParam.push_back(cur);
  return true;
}

template<class T>
void GaussianTimesRN<T>::reset() {
  if(InParam.empty()) {
    for(int i=0; i<gset.size(); i++) gset[i].reset();
  } else {
    int n=gset.size();
    while(n<InParam.size()) {
      gset.push_back(BasicGaussian());
      n++;
    }
    T alpha(1.0),c(1.0); 
    int np(0);
    OhmmsAttributeSet radAttrib;
    radAttrib.add(alpha,expName); 
    radAttrib.add(c,coeffName);
    radAttrib.add(np,powerName);

    for(int i=0; i<InParam.size(); i++) {
      radAttrib.put(InParam[i]);
      gset[i].reset(alpha,c,np+basePower);
    }
  }
}

template<class T>
bool GaussianTimesRN<T>::putBasisGroup(xmlNodePtr cur) {
  const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"basePower");
  if(t!= NULL) basePower = atoi((const char*)t);
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

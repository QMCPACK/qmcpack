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
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

template<class T>
struct GaussianCombo: public OptimizableFunctorBase
{

  typedef T value_type;
  real_type Y, dY, d2Y, d3Y;

  struct BasicGaussian
  {
    real_type Sigma;
    real_type Coeff;

    real_type MinusSigma;
    real_type CoeffP;
    real_type CoeffPP;
    real_type CoeffPPP1;
    real_type CoeffPPP2;
    BasicGaussian(): Sigma(1.0), Coeff(1.0) { }

    inline BasicGaussian(real_type sig, real_type c)
    {
      reset(sig,c);
    }

    inline void reset(real_type sig, real_type c)
    {
      Sigma = sig;
      MinusSigma=-sig;
      Coeff = c;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
      CoeffPPP1= 12.0*Sigma*Sigma*Coeff;
      CoeffPPP2= -8.0*Sigma*Sigma*Sigma*Coeff;
    }

    inline void reset()
    {
      MinusSigma=-Sigma;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
      CoeffPPP1= 12.0*Sigma*Sigma*Coeff;
      CoeffPPP2= -8.0*Sigma*Sigma*Sigma*Coeff;
    }

    void resetParameters(const opt_variables_type& active)
    {
      //DO NOTHING
    }

    inline void setgrid(real_type r) { }

    inline real_type f(real_type rr) const
    {
      return Coeff*exp(MinusSigma*rr);
    }
    inline real_type df(real_type r, real_type rr) const
    {
      return CoeffP*r*exp(MinusSigma*rr);
    }
    inline real_type evaluate(real_type r, real_type rr, real_type& du, real_type& d2u)
    {
      real_type v=exp(MinusSigma*rr);
      du += CoeffP*r*v;
      d2u += (CoeffP+CoeffPP*rr)*v;
      return Coeff*v;
    }
    inline real_type evaluate(real_type r, real_type rr, real_type& du, real_type& d2u, real_type& d3u)
    {
      real_type v=exp(MinusSigma*rr);
      du += CoeffP*r*v;
      d2u += (CoeffP+CoeffPP*rr)*v;
      d3u += (CoeffPPP1*r+CoeffPPP2*r*rr)*v;
      return Coeff*v;
    }
  };

  ///Boolean
  bool Normalized;
  real_type L;
  real_type NormL;
  real_type NormPow;
  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::vector<BasicGaussian> gset;

  explicit
  GaussianCombo(int l=0, bool normalized=false,
                const char* node_name="radfunc",
                const char* exp_name="exponent",
                const char* c_name="contraction");

  ~GaussianCombo() { }

  OptimizableFunctorBase* makeClone() const
  {
    return new GaussianCombo<T>(*this);
  }

  void reset();

  /** return the number Gaussians
   */
  inline int size() const
  {
    return gset.size();
  }

  inline real_type f(real_type r)
  {
    real_type res=0;
    real_type r2 = r*r;
    typename std::vector<BasicGaussian>::const_iterator it(gset.begin());
    typename std::vector<BasicGaussian>::const_iterator it_end(gset.end());
    while(it != it_end)
    {
      res += (*it).f(r2);
      ++it;
    }
    return res;
  }
  inline real_type df(real_type r)
  {
    real_type res=0;
    real_type r2 = r*r;
    typename std::vector<BasicGaussian>::const_iterator it(gset.begin());
    typename std::vector<BasicGaussian>::const_iterator it_end(gset.end());
    while(it != it_end)
    {
      res += (*it).df(r,r2);
      ++it;
    }
    return res;
  }

  inline real_type evaluate(real_type r, real_type rinv)
  {
    Y=0.0;
    real_type rr = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin()),it_end(gset.end());
    while(it != it_end)
    {
      Y+=(*it).f(rr);
      ++it;
    }
    return Y;
  }

  inline void evaluateAll(real_type r, real_type rinv)
  {
    Y=0.0;
    dY=0.0;
    d2Y=0.0;
    real_type rr = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin()),it_end(gset.end());
    while(it != it_end)
    {
      Y+=(*it).evaluate(r,rr,dY,d2Y);
      ++it;
    }
  }

  inline void evaluateWithThirdDeriv(real_type r, real_type rinv)
  {
    Y=0.0;
    dY=0.0;
    d2Y=0.0,d3Y=0.0;
    real_type rr = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin()),it_end(gset.end());
    while(it != it_end)
    {
      Y+=(*it).evaluate(r,rr,dY,d2Y,d3Y);
      ++it;
    }
  }

  //inline real_type evaluate(real_type r, real_type rinv, real_type& drnl, real_type& d2rnl) {
  //  Y=0.0;drnl=0.0;d2rnl=0.0;
  //  real_type du, d2u;
  //  typename std::vector<BasicGaussian>::iterator
  //    it(sset.begin()),it_end(sset.end());
  //  while(it != it_end) {
  //    Y+=(*it).evaluate(r,rinv,du,d2u); dY+=du; d2Y+=d2u;
  //    ++it;
  //  }
  //  return Y;
  //}


  bool put(xmlNodePtr cur);

  void checkInVariables(opt_variables_type& active) { }
  void checkOutVariables(const opt_variables_type& active) { }
  void resetParameters(const opt_variables_type& active)
  {
    //DO NOTHING FOR NOW
  }

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
  int DFactorial(int num)
  {
    return (num<2)? 1: num*DFactorial(num-2);
  }
};

template<class T>
GaussianCombo<T>::GaussianCombo(int l, bool normalized,
                                const char* node_name, const char* exp_name, const char* c_name):
  Normalized(normalized), nodeName(node_name),
  expName(exp_name), coeffName(c_name)
{
  L = static_cast<real_type>(l);
  //Everything related to L goes to NormL and NormPow
  const real_type pi = 4.0*atan(1.0);
  NormL = pow(2,L+1)*sqrt(2.0/static_cast<real_type>(DFactorial(2*l+1)))*pow(2.0/pi,0.25);
  NormPow = 0.5*(L+1.0)+0.25;
}

template<class T>
bool GaussianCombo<T>::put(xmlNodePtr cur)
{
  real_type alpha(1.0),c(1.0);
  OhmmsAttributeSet radAttrib;
  radAttrib.add(alpha,expName);
  radAttrib.add(c,coeffName);
  radAttrib.put(cur);
  real_type c0=c;
  if(!Normalized)
    c *= NormL*pow(alpha,NormPow);
  LOGMSG("    Gaussian exponent = " << alpha
         << "\n              contraction=" << c0 <<  " nomralized contraction = " << c)
  gset.push_back(BasicGaussian(alpha,c));
  return true;
}

template<class T>
void GaussianCombo<T>::reset()
{
  for(int i=0; i<gset.size(); i++)
    gset[i].reset();
}

template<class T>
bool GaussianCombo<T>::putBasisGroup(xmlNodePtr cur)
{
  cur = cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "radfunc")
    {
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

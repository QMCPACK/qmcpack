//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANTIMESRN_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_GAUSSIANTIMESRN_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>

template<class T>
struct GaussianTimesRN: public OptimizableFunctorBase
{

  typedef T value_type;
  real_type Y, dY, d2Y;

  struct BasicGaussian
  {
    real_type Sigma;
    real_type Coeff;
    int Power;

    real_type MinusSigma;
    real_type CoeffP;
    real_type CoeffPP;
    real_type PowerC;

    //BasicGaussian(): Sigma(1.0), Coeff(1.0), Power(0) { }
    inline BasicGaussian(real_type sig=1.0, real_type c=1.0, int p=0):
      Sigma(sig),Coeff(c),Power(p)
    {
      reset();
    }

    inline void resetGaussian(real_type sig, real_type c, int p)
    {
      Sigma=sig;
      Coeff=c;
      Power=p;
      reset();
    }

    inline void reset()
    {
      MinusSigma=-Sigma;
      CoeffP = -2.0*Sigma*Coeff;
      CoeffPP = 4.0*Sigma*Sigma*Coeff;
      PowerC=static_cast<T>(Power)*Coeff;
    }

    inline void setgrid(real_type r) { }

    inline real_type f(real_type r, real_type rr)
    {
      if(Power==0)
        return Coeff*exp(MinusSigma*rr);
      else
        if(Power==1)
          return r*Coeff*exp(MinusSigma*rr);
        else
          return std::pow(r,Power)*Coeff*exp(MinusSigma*rr);
    }

    inline real_type df(real_type r, real_type rr)
    {
      if(Power==0)
        return CoeffP*r*exp(MinusSigma*rr);
      else
        if(Power==1)
          return (Coeff+CoeffP*r)*exp(MinusSigma*rr);
        else
        {
          return
            std::pow(r,Power-1)*(PowerC+CoeffP*rr)*exp(MinusSigma*rr);
        }
    }

    inline real_type evaluate(real_type r, real_type rr, real_type& du, real_type& d2u)
    {
      T v=exp(MinusSigma*rr);
      if(Power==0)
      {
        du += CoeffP*r*v;
        d2u += (CoeffP+CoeffPP*rr)*v;
        return Coeff*v;
      }
      else
      {
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

  ~GaussianTimesRN() { }

  OptimizableFunctorBase* makeClone() const
  {
    return new GaussianTimesRN<T>(*this);
  }

  void reset();

  void checkInVariables(opt_variables_type& active) { }
  void checkOutVariables(const opt_variables_type& active) { }
  void resetParameters(const opt_variables_type& active)
  {
    ///DO NOTHING FOR NOW
  }

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
    typename std::vector<BasicGaussian>::iterator it(gset.begin());
    typename std::vector<BasicGaussian>::iterator it_end(gset.end());
    while(it != it_end)
    {
      res += (*it).f(r,r2);
      ++it;
    }
    return res;
  }

  inline real_type df(real_type r)
  {
    real_type res=0;
    real_type r2 = r*r;
    typename std::vector<BasicGaussian>::iterator it(gset.begin());
    typename std::vector<BasicGaussian>::iterator it_end(gset.end());
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
      Y+=(*it).f(r,rr);
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

  bool put(xmlNodePtr cur);

  void addOptimizables( VarRegistry<real_type>& vlist) {}

  /** process cur xmlnode
   * @param cur root node
   * @param baseOff offset to the basePower
   */
  bool putBasisGroup(xmlNodePtr cur, int baseOff=0);

};

template<class T>
bool GaussianTimesRN<T>::put(xmlNodePtr cur)
{
  InParam.push_back(cur);
  return true;
}

template<class T>
void GaussianTimesRN<T>::reset()
{
  if(InParam.empty())
  {
    for(int i=0; i<gset.size(); i++)
      gset[i].reset();
  }
  else
  {
    int n=gset.size();
    while(n<InParam.size())
    {
      gset.push_back(BasicGaussian());
      n++;
    }
    for(int i=0; i<InParam.size(); i++)
    {
      T alpha(1.0),c(1.0);
      int np(0);
      OhmmsAttributeSet radAttrib;
      radAttrib.add(alpha,expName);
      radAttrib.add(c,coeffName);
      radAttrib.add(np,powerName);
      radAttrib.put(InParam[i]);
      gset[i].resetGaussian(alpha,c,np+basePower);
    }
  }
}

template<class T>
bool GaussianTimesRN<T>::putBasisGroup(xmlNodePtr cur, int baseOff)
{
  const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"basePower");
  if(t!= NULL)
    basePower = atoi((const char*)t);
  basePower += baseOff;
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

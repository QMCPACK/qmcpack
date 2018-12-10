//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_SLATERBASISSET_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_SLATERBASISSET_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/SlaterTypeOrbital.h"
#include "OhmmsData/AttributeSet.h"
#include "io/hdf_archive.h"

namespace qmcplusplus
{
template<class T>
struct SlaterCombo: public OptimizableFunctorBase
{

  typedef T value_type;
  typedef GenericSTO<T> Component_t;

  int L;
  bool Normalized;

  std::string  nodeName;
  std::string  expName;
  std::string  coeffName;
  std::vector<xmlNodePtr> InParam;
  std::vector<Component_t> sset;
  real_type Y, dY, d2Y, d3Y;

  explicit
  SlaterCombo(int l=0,
              bool normalized=true,
              const char* node_name="radfunc",
              const char* exp_name="exponent",
              const char* c_name="contraction");

  ~SlaterCombo() { }

  OptimizableFunctorBase* makeClone() const
  {
    return new SlaterCombo<value_type>(*this);
  }

  inline real_type f(real_type r)
  {
    real_type res=0;
    typename std::vector<Component_t>::iterator it(sset.begin());
    typename std::vector<Component_t>::iterator it_end(sset.end());
    while (it != it_end)
    {
      res += (*it).f(r);
      ++it;
    }
    return res;
  }

  inline real_type df(real_type r)
  {
    real_type res=0;
    typename std::vector<Component_t>::iterator it(sset.begin());
    typename std::vector<Component_t>::iterator it_end(sset.end());
    while (it != it_end)
    {
      res += (*it).df(r);
      ++it;
    }
    return res;
  }

  inline real_type evaluate(real_type r, real_type rinv)
  {
    Y=0.0;
    dY=0.0;
    d2Y=0.0;
    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
    while (it != it_end)
    {
      Y+=(*it).evaluate(r,rinv);
      ++it;
    }
    return Y;
  }

  inline void evaluateAll(real_type r, real_type rinv)
  {
    Y=0.0;
    dY=0.0;
    d2Y=0.0;
    real_type du, d2u;
    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
    while (it != it_end)
    {
      Y+=(*it).evaluate(r,rinv,du,d2u);
      dY+=du;
      d2Y+=d2u;
      ++it;
    }
  }

  inline void  evaluateWithThirdDeriv(real_type r, real_type rinv)
  {
    Y=0.0;
    dY=0.0;
    d2Y=0.0;
    d3Y=0.0;
    real_type du, d2u, d3u;
    typename std::vector<Component_t>::iterator it(sset.begin()),it_end(sset.end());
    while (it != it_end)
    {
      Y+=(*it).evaluate(r,rinv,du,d2u,d3u);
      dY+=du;
      d2Y+=d2u;
      d3Y+=d3u;
      ++it;
    }
  }

//  inline real_type evaluate(T r, T rinv, T& drnl, T& d2rnl) {
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
  bool putBasisGroupH5(hdf_archive &hin)
  {
    APP_ABORT(" Error: Slater Orbitals with HDF5 not implemented. Please contact developers. Aborting.\n");
    return true;
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  void reset()
  {
    //DO NOTHING FOR NOW
  }
  void checkInVariables(opt_variables_type& active) { }
  void checkOutVariables(const opt_variables_type& active) { }
  void resetParameters(const opt_variables_type& active)
  {
    //DO NOTHING FOR NOW
  }

};

template<class T>
SlaterCombo<T>::SlaterCombo(int l, bool normalized,
                            const char* node_name, const char* exp_name, const char* c_name):
  L(l), Normalized(normalized),
  nodeName(node_name), expName(exp_name), coeffName(c_name)
{
}

template<class T>
bool SlaterCombo<T>::putBasisGroup(xmlNodePtr cur)
{
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "radfunc" || cname == "Rnl")
    {
      real_type zeta(1.0),c(1.0);
      int qN=1;
      OhmmsAttributeSet radAttrib;
      radAttrib.add(zeta,expName);
      radAttrib.add(zeta,"alpha");
      radAttrib.add(c,coeffName);
      radAttrib.add(zeta,"c");
      radAttrib.add(qN,"node");
      radAttrib.add(qN,"n");
      radAttrib.put(cur);
      if (Normalized)
      {
        //z is not right
        sset.push_back(Component_t(qN-1,zeta,c));
        LOGMSG(" Slater Component (n,zeta,c)= " << qN-1 << " " << zeta << " " << c)
      }
      else
      {
        STONorm<T> anorm(qN);
        //multiply a normalization factor to the contraction factor
        //anorm(n,\zeta) = 1/\sqrt((2n+2)!/(2*\zeta)^{2*n+3))
        c *= anorm(qN-1,zeta);
        sset.push_back(Component_t(qN-L-1,zeta,c));
        LOGMSG(" Slater Component (n,zeta,c)= " << qN << " " << zeta << " " << c)
      }
    }
    cur=cur->next;
  }
  //reset();
  return true;
}

}
#endif

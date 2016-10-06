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
    
    



/** @file TestFunc.h
 *
 * Test functor in a PW basis set
 */
#ifndef QCMPLUSPLUS_PW_COMPONENT_BENCH_H
#define QCMPLUSPLUS_PW_COMPONENT_BENCH_H
#include <type_traits/scalar_traits.h>
#include <config/stdlib/math.h>

namespace qmcplusplus
{
/** a PW function at K
 */
template<typename T>
struct PWComponent
{

  TinyVector<T,3> K;
  T kk;

  /** default constructor
   *
   * Multiply 2pi
   */
  PWComponent(int nk0=1, int nk1=1, int nk2=1)
  {
    K[0]=TWOPI*static_cast<T>(nk0);
    K[1]=TWOPI*static_cast<T>(nk1);
    K[2]=TWOPI*static_cast<T>(nk2);
    kk=-(K[0]*K[0]+K[1]*K[1]+K[2]*K[2]);
  }

  template<typename PV, typename VT>
  inline void v(const PV& pos, VT& val)
  {
    T kdotp=dot(K,pos);
    std::complex<T> res(std::cos(kdotp),std::sin(kdotp));
    convert(res,val);
  }

  template<class PV, typename VT, typename GT, typename LT>
  inline void vgl(const PV& pos, VT& v, GT& g, LT& l)
  {
    T kdotp=dot(K,pos);
    T sinx=std::sin(kdotp);
    T cosx=std::cos(kdotp);
    std::complex<T> eikr(cosx,sinx);
    std::complex<T> ieikr(-sinx,cosx);
    TinyVector<std::complex<T>,3> g1(K[0]*ieikr,K[1]*ieikr,K[2]*ieikr);
    convert(eikr,v);
    convert(g1,g);
    convert(kk*eikr,l);
  }

  template<class PV, typename VT, typename GT, typename HT>
  inline void vgh(const PV& pos, VT& v, GT& g, HT& h)
  {
    T kdotp=dot(K,pos);
    T sinx=std::sin(kdotp);
    T cosx=std::cos(kdotp);
    std::complex<T> eikr(cosx,sinx);
    std::complex<T> ieikr(-sinx,cosx);
    TinyVector<std::complex<T>,3> g1(K[0]*ieikr,K[1]*ieikr,K[2]*ieikr);
    Tensor<std::complex<T>,3> h1(-K[0]*K[0]*eikr,-K[0]*K[1]*eikr,-K[0]*K[2]*eikr
                            ,-K[1]*K[0]*eikr,-K[1]*K[1]*eikr,-K[1]*K[2]*eikr
                            ,-K[2]*K[0]*eikr,-K[2]*K[1]*eikr,-K[2]*K[2]*eikr
                           );
    convert(eikr,v);
    convert(g1,g);
    convert(h1,h);
  }
};

template<typename T>
struct PW
{

  typedef typename scalar_traits<T>::real_type real_type;
  typedef typename scalar_traits<T>::value_type value_type;
  typedef PWComponent<real_type> basis_type;

  std::vector<value_type> C;
  std::vector<basis_type*> F;

  PW() {}
  ~PW()
  {
    for(int i=0; i<F.size(); i++)
      delete F[i];
  }

  void push_back(value_type c, basis_type* fn)
  {
    C.push_back(c);
    F.push_back(fn);
  }

  template<typename PV>
  inline void v(const PV& pos,value_type& res)
  {
    value_type t;
    res=value_type();
    for(int i=0; i<C.size(); i++)
    {
      F[i]->v(pos,t);
      res+=C[i]*t;
    }
  }

  template<typename PV>
  inline void vgl(const PV& pos, value_type& res, TinyVector<value_type,3>& grad, value_type& lap)
  {
    res=0.0;
    grad=0.0;
    lap=0.0;
    value_type v,l;
    TinyVector<value_type,3> g;
    for(int i=0; i<C.size(); i++)
    {
      F[i]->vgl(pos,v,g,l);
      res+=C[i]*v;
      grad+=C[i]*g;
      lap+=C[i]*l;
    }
  }

  template<typename PV>
  inline void vgh(const PV& pos,value_type& res,TinyVector<value_type,3>& grad, Tensor<value_type,3>& hess)
  {
    res=0.0;
    grad=0.0;
    hess=0.0;
    value_type v;
    TinyVector<value_type,3> g;
    Tensor<value_type,3> h;
    for(int i=0; i<C.size(); i++)
    {
      F[i]->vgh(pos,v,g,h);
      res+=C[i]*v;
      grad+=C[i]*g;
      hess+=C[i]*h;
    }
  }

};

}
#endif


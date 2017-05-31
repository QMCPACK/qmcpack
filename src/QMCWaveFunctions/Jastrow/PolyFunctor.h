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
    
    
#ifndef QMCPLUSPLUS_POLYFUNCTOR_H
#define QMCPLUSPLUS_POLYFUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** Implements a screened Function \f$u[r]=\left(\frac{r-L}{L}\right)^k\sum_i C_i r^i\f$
 *
 * Short-range functor using polynomial functions. Reference to casino
 */
template<class T>
struct PolyFunctor: public OptimizableFunctorBase<T>
{
  ///typedef of real values
  typedef typename OptimizableFunctorBase<T>::real_type real_type;
  typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
  ///power for \f$\left(\frac{r-L}{L}\right)^k\f$
  int K;
  ///maximum power for \f$\sum_i C_i r^i\f$
  int N;
  ///input Rcut
  real_type L;
  ///1/L
  real_type OneOverL;
  ///K/L
  real_type KOverL;
  ///K*(K-1)/L/L
  real_type KKOverL2;
  ///coefficient
  std::vector<T> C;
  ///first derivative of coefficient
  std::vector<T> dC;
  ///second derivative of coefficient
  std::vector<T> d2C;
  ///id
  std::vector<std::string> ID_C;
  ///default constructor
  inline PolyFunctor(): L(10.0), K(0), N(0), OneOverL(0),KOverL(0), KKOverL2(0)
  {}
  ///constructor
  inline PolyFunctor(real_type rc, int k): L(rc), K(k), N(0)
  {
    OneOverL=1.0/L;
    KOverL=static_cast<real_type>(K)/L;
    KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
  }

  inline real_type evaluate(real_type r)
  {
    if(r>L)
      return 0.0;
    real_type rp=r;
    real_type s=C[0];
    for(int i=1; i<N; i++)
    {
      s += C[i]*rp;
      rp*=r;
    }
    return std::pow((r-L)*OneOverL,K)*s;
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if(r<L)
    {
      real_type rp=r, rp2=r*r;
      real_type s=C[0]+C[1]*r+C[2]*rp2;
      real_type ds=dC[1]+dC[2]*r;
      real_type d2s=d2C[2];
      for(int i=3; i<N; i++)
      {
        ds += dC[i]*rp2;
        d2s += d2C[i]*rp;
        real_type rp3=rp2*r;
        s+=C[i]*rp3;
        rp=rp2;
        rp2=rp3;
      }
      real_type x=(r-L)*OneOverL;
      real_type xpk=std::pow(x,K-2);
      dudr=xpk*x*(KOverL*s+x*ds);
      d2udr2=xpk*(KKOverL2*s+x*(2.0*KOverL*ds+x*d2s));
      return xpk*x*x*s;
    }
    else
    {
      dudr=0.0;
      d2udr2=0.0;
      return 0.0;
    }
  }

  inline real_type f(real_type r)
  {
    return evaluate(r);
  }

  inline real_type df(real_type r)
  {
    if(r>L)
      return 0.0;
    real_type rp=r;
    real_type s=C[0]+C[1]*rp;
    real_type ds=dC[1];
    for(int i=2; i<N; i++)
    {
      real_type rp2=rp*r;
      s+=C[i]*rp2;
      ds += dC[i]*rp;
      rp=rp2;
    }
    real_type x=(r-L)*OneOverL;
    return std::pow(x,K-1)*(KOverL*s+x*ds);
  }

  /** basisGroup
   */
  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet aAttrib;
    aAttrib.add(K,"n");
    aAttrib.put(cur);
    std::map<int,std::pair<real_type,std::string> > ctemp;
    cur=cur->children;
    N=0;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "radfunc")
      {
        int i=0;
        real_type c=0.0;
        std::string aname("pn");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(c,"contraction");
        rAttrib.add(i,"node");
        rAttrib.add(aname,"id");
        rAttrib.put(cur);
        aname.append("_C");
        ctemp[i]= std::pair<real_type,std::string>(c,aname);
        N=std::max(i,N);
      }
      cur=cur->next;
    }
    N++;
    C.resize(N,0.0);
    ID_C.resize(N,"0");
    typename std::map<int,std::pair<real_type,std::string> >::iterator it(ctemp.begin());
    typename std::map<int,std::pair<real_type,std::string> >::iterator it_end(ctemp.end());
    app_log() << "    PolyFunctor [(r-L)/L]^K  K=" << K << " L=" << L << std::endl;
    while(it != it_end)
    {
      int i=(*it).first;
      C[i]=(*it).second.first;
      ID_C[i]=(*it).second.second;
      app_log() << "      i="<< i << " "  << ID_C[i] << "=" << C[i] << std::endl;
      ++it;
    }
    resetInternals();
    return true;
  }

  void addOptimizables(OptimizableSetType& vlist)
  {
    for(int i=0; i<N; i++)
      vlist[ID_C[i]]=C[i];
  }

  /** reset the internal variables.
   *
   * USE_resetParameters
   */
  void resetParameters(OptimizableSetType& optVariables)
  {
    for(int i=0; i<N; i++)
    {
      typename OptimizableSetType::iterator it_b(optVariables.find(ID_C[i]));
      if(it_b != optVariables.end())
      {
        C[i]=(*it_b).second;
      }
    }
    resetInternals();
  }

  inline void setL(real_type rc)
  {
    L=rc;
  }


  inline void resetInternals()
  {
    OneOverL=1.0/L;
    KOverL=static_cast<real_type>(K)/L;
    KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
    dC.resize(N,0.0);
    d2C.resize(N,0.0);
    for(int i=1; i<N; i++)
      dC[i]=static_cast<T>(i)*C[i];
    for(int i=2; i<N; i++)
      dC[i]=static_cast<T>(i*(i-1))*C[i];
  }

  inline void add(real_type c, int i)
  {
    C.push_back(c);
    dC.push_back(i*c);
    d2C.push_back(i*(i-1)*c);
    N++;
  }

};

//  //consider full recursive template
//  template<class T, unsigned K>
//    struct RminusLOverLFunctor: public OptimizableFunctorBase<T> {
//
//      typedef typename OptimizableFunctorBase<T>::real_type real_type;
//      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
//
//      ///input Rcut
//      real_type L;
//      ///1/L
//      real_type OneOverL;
//      ///K/L
//      real_type KOverL;
//      ///K*(K-1)/L/L
//      real_type KKOverL2;
//      ///constructor
//      inline RminusLOverLFunctor(real_type rc): L(rc)
//      {
//        OneOverL=1.0/L;
//        KOverL=static_cast<real_type>(K)/L;
//        KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
//      }
//
//      inline real_type evaluate(real_type r)
//      {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K);
//      }
//
//      inline real_type
//        evaluate(real_type r, real_type& dudr, real_type& d2udr2)
//        {
//          if(r<L)
//          {
//            real_type x=(r-L)*OneOverL;
//            real_type xpk=std::pow(x,K-2);
//            d2udr2=KKOverL2*xpk;
//            dudr=KOverL*x*xpk;
//            return xpk*x*x;
//          }
//          else
//          {dudr=0.0; d2udr2=0.0; return 0.0;}
//        }
//
//      inline real_type f(real_type r)
//      {
//        return evaluate(r);
//      }
//
//      inline real_type df(real_type r)
//      {
//        if(r>L) return 0.0;
//        return KOverL*std::pow(x,K-1);
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        return true;
//      }
//
//      void addOptimizables(OptimizableSetType& vlist)
//      {
//      }
//
//      /** reset the internal variables.
//       *
//       * USE_resetParameters
//       */
//      void resetParameters(OptimizableSetType& optVariables)
//      {
//      }
//    };

//  template<class T>
//    struct PowerFunctor<T,0>: public OptimizableFunctorBase<T> {
//
//      typedef typename OptimizableFunctorBase<T>::real_type real_type;
//      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
//
//      ///power for \f$\left(\frac{r-L}{L}\right)^k\f$
//      int K;
//      ///input Rcut
//      real_type L;
//      ///1/L
//      real_type OneOverL;
//      ///K/L
//      real_type KOverL;
//      ///K*(K-1)/L/L
//      real_type KKOverL2;
//
//      ///constructor
//      inline PolyFunctor(real_type rc, int k): L(rc), K(k)
//      {
//        OneOverL=1.0/L;
//        KOverL=static_cast<real_type>(K)/L;
//        KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
//      }
//
//      inline real_type evaluate(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K);
//      }
//
//      inline real_type
//        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//          if(r<L)
//          {
//            real_type x=(r-L)*OneOverL;
//            real_type xpk=std::pow(x,K-2);
//            d2udr2=xpk*KKOverL2;
//            dudr=xpk*x*KOverL;
//            return xpk*x*x;
//          }
//          else
//          {dudr=0.0; d2udr2=0.0; return 0.0;}
//        }
//
//      inline real_type f(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K);
//      }
//
//      inline real_type df(real_type r) {
//        if(r>L) return 0.0;
//        return KOverL*std::pow((r-L)*OneOverL,K-1);
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        return true;
//      }
//
//      void addOptimizables(OptimizableSetType& vlist)
//      {
//      }
//
//      /** reset the internal variables.
//       *
//       * USE_resetParameters
//       */
//      void resetParameters(OptimizableSetType& optVariables)
//      {
//      }
//    };
//
//  template<class T>
//    struct PowerFunctor<T,1>: public OptimizableFunctorBase<T> {
//
//      typedef typename OptimizableFunctorBase<T>::real_type real_type;
//      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
//
//      ///power for \f$\left(\frac{r-L}{L}\right)^k\f$
//      int K;
//      ///input Rcut
//      real_type L;
//      ///1/L
//      real_type OneOverL;
//      ///K/L
//      real_type KOverL;
//      ///K*(K-1)/L/L
//      real_type KKOverL2;
//
//      ///constructor
//      inline PolyFunctor(real_type rc, int k): L(rc), K(k)
//      {
//        OneOverL=1.0/L;
//        KOverL=static_cast<real_type>(K)/L;
//        KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
//      }
//
//      inline real_type evaluate(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K)*r;
//      }
//
//      inline real_type
//        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//          if(r<L)
//          {
//            real_type x=(r-L)*OneOverL;
//            real_type xpk=std::pow(x,K-2);
//            d2udr2=xpk*(KKOverL2+2.0*KOverL*x);
//            dudr=xpk*x*(KOverL*r+x);
//            return xpk*x*x*r;
//          }
//          else
//          {dudr=0.0; d2udr2=0.0; return 0.0;}
//        }
//
//      inline real_type f(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K);
//      }
//
//      inline real_type df(real_type r) {
//        if(r>L) return 0.0;
//        real_type x=(r-L)*OneOverL;
//        real_type xpk=std::pow(x,K-1);
//        return xpk*(KOverL*r+x);
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        return true;
//      }
//
//      void addOptimizables(OptimizableSetType& vlist)
//      {
//      }
//
//      /** reset the internal variables.
//       *
//       * USE_resetParameters
//       */
//      void resetParameters(OptimizableSetType& optVariables)
//      {
//      }
//    };
//
//  template<class T>
//    struct PowerFunctor<T,2>: public OptimizableFunctorBase<T> {
//
//      typedef typename OptimizableFunctorBase<T>::real_type real_type;
//      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
//
//      ///power for \f$\left(\frac{r-L}{L}\right)^k\f$
//      int K;
//      ///input Rcut
//      real_type L;
//      ///1/L
//      real_type OneOverL;
//      ///K/L
//      real_type KOverL;
//      ///K*(K-1)/L/L
//      real_type KKOverL2;
//
//      ///constructor
//      inline PolyFunctor(real_type rc, int k): L(rc), K(k)
//      {
//        OneOverL=1.0/L;
//        KOverL=static_cast<real_type>(K)/L;
//        KKOverL2=static_cast<real_type>(K*(K-1))/L/L;
//      }
//
//      inline real_type evaluate(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K)*r*r;
//      }
//
//      inline real_type
//        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//          if(r<L)
//          {
//            real_type x=(r-L)*OneOverL;
//            real_type xpk=std::pow(x,K-2);
//            d2udr2=xpk*(KKOverL2+2.0*KOverL*x);
//            dudr=xpk*x*(KOverL*r+x);
//            return xpk*x*x*r;
//          }
//          else
//          {dudr=0.0; d2udr2=0.0; return 0.0;}
//        }
//
//      inline real_type f(real_type r) {
//        if(r>L) return 0.0;
//        return std::pow((r-L)*OneOverL,K);
//      }
//
//      inline real_type df(real_type r) {
//        if(r>L) return 0.0;
//        real_type x=(r-L)*OneOverL;
//        real_type xpk=std::pow(x,K-1);
//        return xpk*(KOverL*r+x);
//      }
//
//      bool put(xmlNodePtr cur)
//      {
//        return true;
//      }
//
//      void addOptimizables(OptimizableSetType& vlist)
//      {
//      }
//
//      /** reset the internal variables.
//       *
//       * USE_resetParameters
//       */
//      void resetParameters(OptimizableSetType& optVariables)
//      {
//      }
//    };
}
#endif

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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#define QMCPLUSPLUS_WMFUNCTOR_SM_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/NumericalJastrowFunctor.h"

namespace qmcplusplus {

  /** Implements a screened Function \f$u[r]=(1-z(r/rcut))/(1+B*z(r/rcut)\f$
   * 
   * Short-range functor introduced by Wagner and Mitas, cond-mat/0610088
   */
  template<class T>
    struct WMFunctor {
      typedef T real_type;
      ///input B
      real_type B0;
      ///input Rcut
      real_type Rcut;
      ///1/Rcut
      real_type OneOverRc;
      ///id
      string ID;
      ///constructor
      explicit WMFunctor(real_type b, real_type rc=7.5) {
        reset(b,rc);
      }
      inline void reset() { OneOverRc=1.0/Rcut; }
      void reset(real_type b, real_type rc) { B0=b; Rcut=rc; reset(); }
      inline real_type f(real_type r) {
        real_type x=r*OneOverRc;
        real_type z=x*x*(6.0-8*x+3.0*x*x);
        return (1-z)/(1+B0*z);
      }
      inline real_type df(real_type r) {
        real_type x=r*OneOverRc;
        real_type z=x*x*(6.0-8*x+3.0*x*x);
        return -(1+B0)/(1+B0*z)*(1+B0*z)*OneOverRc*12*x*(1-2.0*x+x*x);
      }
    };

  /** Implements a linear combination of WMFunctors for Two-body and One-body Jastrow orbitals
   */
  template<class T>
    struct WMComboFunctor: public JastrowFunctorBase<T> {
      typedef typename JastrowFunctorBase<T>::real_type real_type;
      vector<T> C;
      vector<WMFunctor*> Phi;
      vector<string> ID;

      WMComboFunctor() { 
        C.reserve(8);
        Phi.reserve(8);
        ID.reserve(8);
      }

      int size() const { return Phi.size();}

      void add(WMFunctor* func, T c,  const string& id) {
        C.push_back(c);
        Phi.push_back(func);
        ID.push_back(id);
      }

      inline void reset() {
        for(int i=0; i<Phi.size(); i++) Phi.[i]->reset();
      }

      inline real_type f(real_type r) {
        real_type res=0;
        for(int i=0; i<Phi.size(); i++) { res += C[i]*Phi[i]->f(r);}
        return res;
      }

      inline real_type df(real_type r) {
        real_type res(0);
        for(int i=0; i<Phi.size(); i++) { res += C[i]*Phi[i]->df(r);}
        return res;
      }

      void put(xmlNodePtr cur, VarRegistry<real_type>& vlist) {}
    };


  struct WMConstraints: public OrbitalConstraintsBase {
    ///analytic functor
    typedef WMComboFunctor<RealType> InFuncType;
    ///numerical functor
    typedef SplineJastrow<RealType> FuncType;
    bool IgnoreSpin;
    RealType Rcut;

    typedef vector<WMFunctor*> BasisSetType;
    map<string,BasisSetType*> myBasisSet;
    vector<InFuncType*> InFuncList;
    vector<FuncType*> FuncList;

    ~WMConstraints();
    WMConstraints(bool nospin=true):IgnoreSpin(nospin) {}
    void apply();
    void addOptimizables(VarRegistry<RealType>& outVars);
    OrbitalBase* createTwoBody(ParticleSet& target);
    OrbitalBase* createOneBody(ParticleSet& target, ParticleSet& source);
    bool put(xmlNodePtr cur);
    void addBasisGroup(xmlNodePtr cur);
    InFuncType* createCorrelation(xmlNodePtr cur, BasisSetType* basis);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

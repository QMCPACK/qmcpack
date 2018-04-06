//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SOA_MULTIANALYTICFUNCTOR_BUILDER_H
#define QMCPLUSPLUS_SOA_MULTIANALYTICFUNCTOR_BUILDER_H

#include "Configuration.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/GaussianBasisSet.h"

namespace qmcplusplus
{

  /** generic functor that computes a set of 1D functors
   * @tparam FN analytic functor, SlaterCombo<T>, GaussianCombo<T>
   *
   * Analytic functors are light and no state but not efficient.
   * Only for benchmarking.
   */
  template<typename FN>
    struct MultiFunctorAdapter
    {
      typedef typename FN::real_type value_type;
      typedef LogGridLight<value_type> grid_type;
      typedef FN single_type;
      aligned_vector<single_type*> Rnl;

      MultiFunctorAdapter<FN>* makeClone() const
      {
        MultiFunctorAdapter<FN>* clone=new MultiFunctorAdapter<FN>(*this);
        for(size_t i=0; i<Rnl.size(); ++i)
          clone->Rnl[i]=new single_type(*Rnl[i]);
        return clone;
      }

      ~MultiFunctorAdapter()
      {
        for(size_t i=0; i<Rnl.size(); ++i) delete Rnl[i];
      }

      inline value_type rmax() const
      {
        CONSTEXPR value_type r0(100);
        return r0;
      }

      inline void evaluate(value_type r, value_type* restrict u) 
      {
        for(size_t i=0, n=Rnl.size(); i<n; ++i)
          u[i]=Rnl[i]->f(r);
      }

      inline void evaluate(value_type r, value_type* restrict u, value_type* restrict du, value_type* restrict d2u) 
      {
        const value_type rinv=value_type(1)/r;
        for(size_t i=0, n=Rnl.size(); i<n; ++i)
        {
          Rnl[i]->evaluateAll(r,rinv);
          u[i]=Rnl[i]->Y;
          du[i]=Rnl[i]->dY;
          d2u[i]=Rnl[i]->d2Y;
        }
      }
    };

  template<typename FN, typename SH>
    struct RadialOrbitalSetBuilder<SoaAtomicBasisSet<MultiFunctorAdapter<FN>, SH> >
    {
      typedef SoaAtomicBasisSet<MultiFunctorAdapter<FN>,SH> COT;
      typedef MultiFunctorAdapter<FN>                       RadialOrbital_t;
      typedef typename RadialOrbital_t::single_type    single_type;
   
      ///true, if the RadialOrbitalType is normalized
      bool Normalized;
      ///orbitals to build
      COT* m_orbitals;
      ///temporary
      RadialOrbital_t* m_multiset;

      ///constructor
      RadialOrbitalSetBuilder(xmlNodePtr cur=NULL):m_multiset(nullptr)
      {
        if(cur != NULL) putCommon(cur);
        Normalized=true;
      }

      ///implement functions used by AOBasisBuilder
      void setOrbitalSet(COT* oset, const std::string& acenter) { m_orbitals = oset; }
      bool addGrid(xmlNodePtr cur) { return true;}
      bool addGridH5(hdf_archive &hin) { return true;}
      inline bool putCommon(xmlNodePtr cur) 
      { 
        const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
        if(a)
        {
          if(xmlStrEqual(a,(const xmlChar*)"no"))
            Normalized=false;
        }
        return true;
      }
      bool put(xmlNodePtr cur)
      {
        const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
        if(a)
        {
          if(xmlStrEqual(a,(const xmlChar*)"no"))
            Normalized=false;
        }
        return true;
      }

      bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
      {

        if(m_multiset==nullptr)
          m_multiset=new RadialOrbital_t;

        single_type* radorb= new single_type(nlms[q_l],Normalized);
        radorb->putBasisGroup(cur);

        m_orbitals->RnlID.push_back(nlms);
        m_multiset->Rnl.push_back(radorb);
        return true;
      }

      bool addRadialOrbitalH5(hdf_archive & hin, const QuantumNumberType& nlms)
      {
        if(m_multiset==nullptr)
          m_multiset=new RadialOrbital_t;

        single_type* radorb= new single_type(nlms[q_l],Normalized);
        radorb->putBasisGroupH5(hin);

        m_orbitals->RnlID.push_back(nlms);
        m_multiset->Rnl.push_back(radorb);

        return true;
      }

      void finalize() 
      {
        m_orbitals->MultiRnl=m_multiset;
        m_orbitals->setRmax(m_multiset->rmax()); //set Rmax
      }

    };
}
#endif

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
    
    
#ifndef QMCPLUSPLUS_SOA_GAUSSIANTYPEORBITAL_BUILDER_H
#define QMCPLUSPLUS_SOA_GAUSSIANTYPEORBITAL_BUILDER_H

#include "Configuration.h"
#include "Numerics/GaussianBasisSet.h"

namespace qmcplusplus
{

  /** multi version of SlaterCombo for SoaAtomicBasisSet 
   */
  template<typename T>
    struct MultiGTO
    {
      typedef GaussianCombo<T> gto_type;
      aligned_vector<gto_type*> Rnl;
      inline void evaluate(T r, T* restrict u) 
      {
        for(size_t i=0, n=Rnl.size(); i<n; ++i)
          u[i]=Rnl[i]->f(r);
      }
      inline void evaluate(T r, T* restrict u, T* restrict du, T* restrict d2u) 
      {
        const size_t n=Rnl.size();
        for(size_t i=0; i<n; ++i)
          u[i]=Rnl[i]->evaluate(r,du[i],d2u[i]);
      }
    };

  /** scpecialization with STO 
   */
  template<>
    class RadialOrbitalSetBuilder<MultiGTO<QMCTraits::RealType> >
    {
      typedef OHMMS_PRECISION RealType;
      typedef MultiGTO<RealType> RadialOrbital_t;
      typedef MultiGTO<RealType>::sto_type sto_type;
      typedef SoaAtomicBasisSet<RadialOrbital_t> COT;
      typedef LogGridLight<RealType> GridType;

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

        sto_type* radorb= new sto_type(nlms[q_l],Normalized);
        radorb->putBasisGroup(cur);

        m_orbitals->RnlID.push_back(nlms);
        m_multiset->push_back(radorb);
        return true;
      }

      void finalize() 
      {
        m_orbitals->MultiRnl=m_multiset;
        m_orbitals->setRmax(50.); //set Rmax
      }

    };
}
#endif

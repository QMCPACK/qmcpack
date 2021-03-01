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
  using RealType = typename FN::real_type;
  using GridType = LogGridLight<RealType>;
  typedef FN single_type;
  aligned_vector<single_type*> Rnl;

  MultiFunctorAdapter<FN>* makeClone() const
  {
    MultiFunctorAdapter<FN>* clone = new MultiFunctorAdapter<FN>(*this);
    for (size_t i = 0; i < Rnl.size(); ++i)
      clone->Rnl[i] = new single_type(*Rnl[i]);
    return clone;
  }

  ~MultiFunctorAdapter()
  {
    for (size_t i = 0; i < Rnl.size(); ++i)
      delete Rnl[i];
  }

  inline RealType rmax() const
  {
    //Another magic r_max
    constexpr RealType r0(100);
    return r0;
  }

  inline void evaluate(RealType r, RealType* restrict u)
  {
    for (size_t i = 0, n = Rnl.size(); i < n; ++i)
      u[i] = Rnl[i]->f(r);
  }

  inline void evaluate(RealType r, RealType* restrict u, RealType* restrict du, RealType* restrict d2u)
  {
    const RealType rinv = RealType(1) / r;
    for (size_t i = 0, n = Rnl.size(); i < n; ++i)
    {
      Rnl[i]->evaluateAll(r, rinv);
      u[i]   = Rnl[i]->Y;
      du[i]  = Rnl[i]->dY;
      d2u[i] = Rnl[i]->d2Y;
    }
  }
  inline void evaluate(RealType r,
                       RealType* restrict u,
                       RealType* restrict du,
                       RealType* restrict d2u,
                       RealType* restrict d3u)
  {
    const RealType rinv = RealType(1) / r;
    for (size_t i = 0, n = Rnl.size(); i < n; ++i)
    {
      Rnl[i]->evaluateWithThirdDeriv(r, rinv);
      u[i]   = Rnl[i]->Y;
      du[i]  = Rnl[i]->dY;
      d2u[i] = Rnl[i]->d2Y;
      d3u[i] = Rnl[i]->d3Y;
    }
  }
};

template<typename FN, typename SH>
struct RadialOrbitalSetBuilder<SoaAtomicBasisSet<MultiFunctorAdapter<FN>, SH>> : public MPIObjectBase
{
  typedef SoaAtomicBasisSet<MultiFunctorAdapter<FN>, SH> COT;
  typedef MultiFunctorAdapter<FN> RadialOrbital_t;
  typedef typename RadialOrbital_t::single_type single_type;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///orbitals to build
  COT* m_orbitals;
  ///temporary
  RadialOrbital_t* m_multiset;

  ///constructor
  RadialOrbitalSetBuilder(Communicate* comm) : MPIObjectBase(comm), Normalized(true), m_multiset(nullptr) {}

  ///implement functions used by AOBasisBuilder
  void setOrbitalSet(COT* oset, const std::string& acenter) { m_orbitals = oset; }
  bool addGrid(xmlNodePtr cur, const std::string& rad_type) { return true; }
  bool addGridH5(hdf_archive& hin) { return true; }
  bool openNumericalBasisH5(xmlNodePtr cur) { return true; }
  bool put(xmlNodePtr cur)
  {
    const XMLAttrString a(cur, "normalized");
    if (a == "no")
      Normalized = false;
    return true;
  }

  bool addRadialOrbital(xmlNodePtr cur, const std::string& rad_type, const QuantumNumberType& nlms)
  {
    if (m_multiset == nullptr)
      m_multiset = new RadialOrbital_t;

    single_type* radorb = new single_type(nlms[q_l], Normalized);
    radorb->putBasisGroup(cur);

    m_orbitals->RnlID.push_back(nlms);
    m_multiset->Rnl.push_back(radorb);
    return true;
  }

  bool addRadialOrbitalH5(hdf_archive& hin, const std::string& rad_type, const QuantumNumberType& nlms)
  {
    if (m_multiset == nullptr)
      m_multiset = new RadialOrbital_t;

    single_type* radorb = new single_type(nlms[q_l], Normalized);
    radorb->putBasisGroupH5(hin);

    m_orbitals->RnlID.push_back(nlms);
    m_multiset->Rnl.push_back(radorb);

    return true;
  }

  void finalize()
  {
    m_orbitals->MultiRnl = m_multiset;
    m_orbitals->setRmax(m_multiset->rmax()); //set Rmax
  }
};
} // namespace qmcplusplus
#endif

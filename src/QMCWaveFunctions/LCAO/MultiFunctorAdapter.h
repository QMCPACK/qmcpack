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
#include "Message/MPIObjectBase.h"
#include "ModernStringUtils.hpp"
#include "hdf/hdf_archive.h"

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
  using RealType       = typename FN::real_type;
  using GridType       = LogGridLight<RealType>;
  using single_type    = FN;
  using OffloadArray2D = Array<RealType, 2, OffloadPinnedAllocator<RealType>>;
  using OffloadArray3D = Array<RealType, 3, OffloadPinnedAllocator<RealType>>;
  aligned_vector<std::unique_ptr<single_type>> Rnl;

  MultiFunctorAdapter() = default;
  MultiFunctorAdapter(const MultiFunctorAdapter& other)
  {
    Rnl.reserve(other.Rnl.size());
    for (size_t i = 0; i < other.Rnl.size(); ++i)
      Rnl.push_back(std::make_unique<single_type>(*other.Rnl[i]));
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

  inline void batched_evaluate(OffloadArray2D& r, OffloadArray3D& u, RealType Rmax) const
  {
    r.updateFrom(); // TODO: remove after offload
    const size_t nElec = r.size(0);
    const size_t Nxyz  = r.size(1); // number of PBC images
    assert(nElec == u.size(0));
    assert(Nxyz == u.size(1));
    const size_t nRnl = u.size(2);    // number of splines
    const size_t nR   = nElec * Nxyz; // total number of positions to evaluate

    auto* r_ptr = r.data();
    auto* u_ptr = u.data();


    for (size_t ir = 0; ir < nR; ir++)
    {
      if (r_ptr[ir] >= Rmax)
      {
        for (size_t i = 0, n = Rnl.size(); i < n; ++i)
          u_ptr[ir * nRnl + i] = 0.0;
      }
      else
      {
        for (size_t i = 0, n = Rnl.size(); i < n; ++i)
          u_ptr[ir * nRnl + i] = Rnl[i]->f(r_ptr[ir]);
      }
    }

    u.updateTo(); // TODO: remove after offload
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

template<typename COT>
class RadialOrbitalSetBuilder;

template<typename FN, typename SH>
class RadialOrbitalSetBuilder<SoaAtomicBasisSet<MultiFunctorAdapter<FN>, SH>> : public MPIObjectBase
{
public:
  using COT             = SoaAtomicBasisSet<MultiFunctorAdapter<FN>, SH>;
  using RadialOrbital_t = MultiFunctorAdapter<FN>;
  using single_type     = typename RadialOrbital_t::single_type;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///orbitals to build
  COT& m_orbitals;

  ///constructor
  RadialOrbitalSetBuilder(Communicate* comm, COT& aos) : MPIObjectBase(comm), Normalized(true), m_orbitals(aos) {}

  ///implement functions used by AOBasisBuilder
  bool addGrid(xmlNodePtr cur, const std::string& rad_type) { return true; }
  bool addGridH5(hdf_archive& hin) { return true; }
  bool openNumericalBasisH5(xmlNodePtr cur) { return true; }
  bool put(xmlNodePtr cur)
  {
    const std::string a(lowerCase(getXMLAttributeValue(cur, "normalized")));
    if (a == "no")
      Normalized = false;
    return true;
  }

  bool addRadialOrbital(xmlNodePtr cur, const std::string& rad_type, const QuantumNumberType& nlms)
  {
    auto radorb = std::make_unique<single_type>(nlms[q_l], Normalized);
    radorb->putBasisGroup(cur);

    m_orbitals.RnlID.push_back(nlms);
    m_orbitals.MultiRnl.Rnl.push_back(std::move(radorb));
    return true;
  }

  bool addRadialOrbitalH5(hdf_archive& hin, const std::string& rad_type, const QuantumNumberType& nlms)
  {
    auto radorb = std::make_unique<single_type>(nlms[q_l], Normalized);
    radorb->putBasisGroupH5(hin, *myComm);

    m_orbitals.RnlID.push_back(nlms);
    m_orbitals.MultiRnl.Rnl.push_back(std::move(radorb));

    return true;
  }

  void finalize()
  {
    m_orbitals.setRmax(0); //set Rmax
  }
};
} // namespace qmcplusplus
#endif

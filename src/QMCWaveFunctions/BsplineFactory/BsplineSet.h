//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineSet.h
 *
 * BsplineSet is a SPOSet derived class and serves as a base class for B-spline SPO C2C/C2R/R2R implementation
 */
#ifndef QMCPLUSPLUS_BSPLINESET_H
#define QMCPLUSPLUS_BSPLINESET_H

#include "QMCWaveFunctions/SPOSet.h"
#include "spline/einspline_engine.hpp"
#include "spline/einspline_util.hpp"

namespace qmcplusplus
{
/** BsplineSet is the base class for SplineC2C, SplineC2R, SplineR2R.
 * Its derived template classes manage the storage and evaluation at given precision.
 * BsplineSet also implements a few fallback routines in case optimized implementation is not necessary in the derived class.
 */
class BsplineSet : public SPOSet
{
protected:
  static const int D = DIM;
  ///Index of this adoptor, when multiple adoptors are used for NUMA or distributed cases
  size_t MyIndex;
  ///first index of the SPOs this Spline handles
  size_t first_spo;
  ///last index of the SPOs this Spline handles
  size_t last_spo;
  ///sign bits at the G/2 boundaries
  TinyVector<int, D> HalfG;
  ///flags to unpack sin/cos
  std::vector<bool> MakeTwoCopies;
  /** kpoints for each unique orbitals.
   * Note: for historic reason, this sign is opposite to what was used in DFT when orbitals were generated.
   * Changing the sign requires updating all the evaluation code.
   */
  std::vector<SPOSet::PosType> kPoints;
  ///remap splines to orbitals
  aligned_vector<int> BandIndexMap;
  ///band offsets used for communication
  std::vector<int> offset;

public:
  BsplineSet(const std::string& my_name) : SPOSet(my_name), MyIndex(0), first_spo(0), last_spo(0) {}

  virtual bool isComplex() const         = 0;
  virtual std::string getKeyword() const = 0;

  auto& getHalfG() const { return HalfG; }

  inline void init_base(int n)
  {
    kPoints.resize(n);
    MakeTwoCopies.resize(n);
    BandIndexMap.resize(n);
    for (int i = 0; i < n; i++)
      BandIndexMap[i] = i;
  }

  ///remap kpoints to group general kpoints & special kpoints
  int remap_kpoints()
  {
    std::vector<SPOSet::PosType> k_copy(kPoints);
    const int nk = kPoints.size();
    int nCB      = 0;
    //two pass
    for (int i = 0; i < nk; ++i)
    {
      if (MakeTwoCopies[i])
      {
        kPoints[nCB]        = k_copy[i];
        BandIndexMap[nCB++] = i;
      }
    }
    int nRealBands = nCB;
    for (int i = 0; i < nk; ++i)
    {
      if (!MakeTwoCopies[i])
      {
        kPoints[nRealBands]        = k_copy[i];
        BandIndexMap[nRealBands++] = i;
      }
    }
    return nCB; //return the number of complex bands
  }

  // propagate SPOSet virtual functions
  using SPOSet::evaluateDetRatios;
  using SPOSet::evaluateValue;
  using SPOSet::evaluateVGH;
  using SPOSet::evaluateVGHGH;
  using SPOSet::evaluateVGL;
  using SPOSet::finalizeConstruction;
  using SPOSet::mw_evaluateDetRatios;
  using SPOSet::mw_evaluateVGL;
  using SPOSet::mw_evaluateVGLandDetRatioGrads;

  using SPOSet::acquireResource;
  using SPOSet::createResource;
  using SPOSet::releaseResource;

  std::unique_ptr<SPOSet> makeClone() const override = 0;

  void setOrbitalSetSize(int norbs) override { OrbitalSetSize = norbs; }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override
  {
    using value_type = ValueMatrix::value_type;
    using grad_type  = GradMatrix::value_type;
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector v(logdet[i], logdet.cols());
      GradVector g(dlogdet[i], dlogdet.cols());
      ValueVector l(d2logdet[i], d2logdet.cols());
      evaluateVGL(P, iat, v, g, l);
    }
  }

  void mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                               const RefVectorWithLeader<ParticleSet>& P_list,
                               int first,
                               int last,
                               const RefVector<ValueMatrix>& logdet_list,
                               const RefVector<GradMatrix>& dlogdet_list,
                               const RefVector<ValueMatrix>& d2logdet_list) const override
  {
    assert(this == &spo_list.getLeader());
    using value_type = ValueMatrix::value_type;
    using grad_type  = GradMatrix::value_type;

    const size_t nw = spo_list.size();
    std::vector<ValueVector> mw_psi_v;
    std::vector<GradVector> mw_dpsi_v;
    std::vector<ValueVector> mw_d2psi_v;
    RefVector<ValueVector> psi_v_list;
    RefVector<GradVector> dpsi_v_list;
    RefVector<ValueVector> d2psi_v_list;
    mw_psi_v.reserve(nw);
    mw_dpsi_v.reserve(nw);
    mw_d2psi_v.reserve(nw);
    psi_v_list.reserve(nw);
    dpsi_v_list.reserve(nw);
    d2psi_v_list.reserve(nw);

    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      mw_psi_v.clear();
      mw_dpsi_v.clear();
      mw_d2psi_v.clear();
      psi_v_list.clear();
      dpsi_v_list.clear();
      d2psi_v_list.clear();

      for (int iw = 0; iw < nw; iw++)
      {
        mw_psi_v.emplace_back(logdet_list[iw].get()[i], logdet_list[iw].get().cols());
        mw_dpsi_v.emplace_back(dlogdet_list[iw].get()[i], dlogdet_list[iw].get().cols());
        mw_d2psi_v.emplace_back(d2logdet_list[iw].get()[i], d2logdet_list[iw].get().cols());
        psi_v_list.push_back(mw_psi_v.back());
        dpsi_v_list.push_back(mw_dpsi_v.back());
        d2psi_v_list.push_back(mw_d2psi_v.back());
      }

      mw_evaluateVGL(spo_list, P_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list);
    }
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector v(logdet[i], logdet.cols());
      GradVector g(dlogdet[i], dlogdet.cols());
      HessVector h(grad_grad_logdet[i], grad_grad_logdet.cols());
      evaluateVGH(P, iat, v, g, h);
    }
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector v(logdet[i], logdet.cols());
      GradVector g(dlogdet[i], dlogdet.cols());
      HessVector h(grad_grad_logdet[i], grad_grad_logdet.cols());
      GGGVector gh(grad_grad_grad_logdet[i], grad_grad_grad_logdet.cols());
      evaluateVGHGH(P, iat, v, g, h, gh);
    }
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& gradphi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi,
                          HessMatrix& grad_grad_phi,
                          GradMatrix& grad_lapl_phi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  template<class BSPLINESPO>
  friend struct SplineSetReader;
  friend struct BsplineReaderBase;
};

} // namespace qmcplusplus
#endif

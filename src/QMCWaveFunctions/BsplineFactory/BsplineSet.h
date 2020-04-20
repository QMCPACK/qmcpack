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
  ///true if the computed values are complex
  bool is_complex;
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
  ///kpoints for each unique orbitals
  std::vector<SPOSet::PosType> kPoints;
  ///remap splines to orbitals
  aligned_vector<int> BandIndexMap;
  ///band offsets used for communication
  std::vector<int> offset;
  ///keyword used to match hdf5
  std::string KeyWord;

public:
  BsplineSet(bool use_OMP_offload = false, bool ion_deriv = false, bool optimizable = false)
      : SPOSet(use_OMP_offload, ion_deriv, optimizable), is_complex(false), MyIndex(0), first_spo(0), last_spo(0) {}

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
  using SPOSet::mw_evaluateDetRatios;
  using SPOSet::evaluateValue;
  using SPOSet::evaluateVGH;
  using SPOSet::evaluateVGHGH;
  using SPOSet::evaluateVGL;
  using SPOSet::finalizeConstruction;
  using SPOSet::mw_evaluateVGL;

  virtual SPOSet* makeClone() const override = 0;

  void resetParameters(const opt_variables_type& active) override {}

  void resetTargetParticleSet(ParticleSet& e) override {}

  void setOrbitalSetSize(int norbs) override { OrbitalSetSize = norbs; }

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet) override
  {
    typedef ValueMatrix_t::value_type value_type;
    typedef GradMatrix_t::value_type grad_type;
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      ValueVector_t l(d2logdet[i], OrbitalSetSize);
      evaluateVGL(P, iat, v, g, l);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      HessVector_t h(grad_grad_logdet[i], OrbitalSetSize);
      evaluateVGH(P, iat, v, g, h);
    }
  }

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    HessMatrix_t& grad_grad_logdet,
                                    GGGMatrix_t& grad_grad_grad_logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector_t v(logdet[i], OrbitalSetSize);
      GradVector_t g(dlogdet[i], OrbitalSetSize);
      HessVector_t h(grad_grad_logdet[i], OrbitalSetSize);
      GGGVector_t gh(grad_grad_grad_logdet[i], OrbitalSetSize);
      evaluateVGHGH(P, iat, v, g, h, gh);
    }
  }

  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& gradphi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  virtual void evaluateGradSource(const ParticleSet& P,
                                  int first,
                                  int last,
                                  const ParticleSet& source,
                                  int iat_src,
                                  GradMatrix_t& grad_phi,
                                  HessMatrix_t& grad_grad_phi,
                                  GradMatrix_t& grad_lapl_phi) override
  {
    //Do nothing, since Einsplines don't explicitly depend on ion positions.
  }

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend class BsplineReaderBase;
};

} // namespace qmcplusplus
#endif

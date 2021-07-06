//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_TWOBODYJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_TWOBODYJASTROW_OPTIMIZED_SOA_H

#include <map>
#include <numeric>
#include "Configuration.h"
#if !defined(QMC_BUILD_SANDBOX_ONLY)
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#endif
#include "Particle/DistanceTableData.h"
#include "LongRange/StructFact.h"
#include "CPU/SIMD/aligned_allocator.hpp"

namespace qmcplusplus
{
// helper class to activate KEcorr during optimizing Jastrow
template<typename RT, class FT>
class J2KECorrection
{
  size_t num_groups_;
  std::vector<size_t> num_elec_in_groups_;
  RT num_elecs_;
  RT vol;
  RT G0mag;
  const std::vector<FT*>& F_;
  bool SK_enabled;

public:
  J2KECorrection(const ParticleSet& targetPtcl, const std::vector<FT*>& F)
      : num_groups_(targetPtcl.groups()),
        num_elecs_(targetPtcl.getTotalNum()),
        vol(targetPtcl.Lattice.Volume),
        F_(F),
        SK_enabled(targetPtcl.SK != nullptr)
  {
    // compute num_elec_in_groups_
    num_elec_in_groups_.reserve(3);
    for (int i = 0; i < num_groups_; i++)
      num_elec_in_groups_.push_back(targetPtcl.last(i) - targetPtcl.first(i));

    if (SK_enabled)
      G0mag = std::sqrt(targetPtcl.SK->KLists.ksq[0]);
  }

  RT computeKEcorr()
  {
    if (!SK_enabled)
      return 0;

    const int numPoints = 1000;
    RT uk               = 0.0;
    RT a                = 1.0;

    for (int i = 0; i < num_groups_; i++)
    {
      int Ni = num_elec_in_groups_[i];
      for (int j = 0; j < num_groups_; j++)
      {
        int Nj = num_elec_in_groups_[j];
        if (F_[i * num_groups_ + j])
        {
          FT& ufunc = *(F_[i * num_groups_ + j]);
          RT radius = ufunc.cutoff_radius;
          RT k      = G0mag;
          RT dr     = radius / (RT)(numPoints - 1);
          for (int ir = 0; ir < numPoints; ir++)
          {
            RT r = dr * (RT)ir;
            RT u = ufunc.evaluate(r);
            uk += 0.5 * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr * (RT)Nj / (RT)(Ni + Nj);
          }
        }
      }
    }
    for (int iter = 0; iter < 20; iter++)
      a = uk / (4.0 * M_PI * (1.0 / (G0mag * G0mag) - 1.0 / (G0mag * G0mag + 1.0 / a)));
    return 4.0 * M_PI * a / (4.0 * vol) * num_elecs_;
  }
};

/** @ingroup WaveFunctionComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 *
 * Each pair-type can have distinct function \f$u(r_{ij})\f$.
 * For electrons, distinct pair correlation functions are used
 * for spins up-up/down-down and up-down/down-up.
 *
 * Based on J2OrbitalSoA.h with these considerations
 * - DistanceTableData using SoA containers
 * - support mixed precision: FT::real_type != OHMMS_PRECISION
 * - loops over the groups: elminated PairID
 * - support simd function
 * - double the loop counts
 * - Memory use is O(N). 
 */
template<class FT>
class J2OrbitalSoA : public WaveFunctionComponent
{
public:
  ///alias FuncType
  using FuncType = FT;
  ///type of each component U, dU, d2U;
  using valT = typename FT::real_type;
  ///element position type
  using posT = TinyVector<valT, OHMMS_DIM>;
  ///use the same container
  using DistRow         = DistanceTableData::DistRow;
  using DisplRow        = DistanceTableData::DisplRow;
  using gContainer_type = VectorSoaContainer<valT, OHMMS_DIM>;

protected:
  ///number of particles
  size_t N;
  ///number of particles + padded
  size_t N_padded;
  ///number of groups of the target particleset
  size_t NumGroups;
  ///diff value
  RealType DiffVal;
  ///Correction
  RealType KEcorr;
  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Uat;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  gContainer_type dUat;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  Vector<valT> d2Uat;
  valT cur_Uat;
  aligned_vector<valT> cur_u, cur_du, cur_d2u;
  aligned_vector<valT> old_u, old_du, old_d2u;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  ///Uniquue J2 set for cleanup
  std::map<std::string, std::unique_ptr<FT>> J2Unique;
  ///Container for \f$F[ig*NumGroups+jg]\f$. treat every pointer as a reference.
  std::vector<FT*> F;
  /// e-e table ID
  const int my_table_ID_;
  // helper for compute J2 Chiesa KE correction
  J2KECorrection<RealType, FT> j2_ke_corr_helper;

public:
  J2OrbitalSoA(const std::string& obj_name, ParticleSet& p);
  J2OrbitalSoA(const J2OrbitalSoA& rhs) = delete;
  ~J2OrbitalSoA() override;

  /* initialize storage */
  void init(ParticleSet& p);

  /** add functor for (ia,ib) pair */
  void addFunc(int ia, int ib, std::unique_ptr<FT> j);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active) override;

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active) override;

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active) override;

  inline void finalizeOptimization() override { KEcorr = j2_ke_corr_helper.computeKEcorr(); }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi) override;

  /** recompute internal data assuming distance table is fully ready */
  void recompute(const ParticleSet& P) override;

  PsiValueType ratio(ParticleSet& P, int iat) override;
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override;
  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  GradType evalGrad(ParticleSet& P, int iat) override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;
  inline void restore(int iat) override {}

  /** compute G and L after the sweep
   */
  LogValueType evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch = false) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  /*@{ internal compute engines*/
  valT computeU(const ParticleSet& P, int iat, const DistRow& dist);

  void computeU3(const ParticleSet& P,
                 int iat,
                 const DistRow& dist,
                 RealType* restrict u,
                 RealType* restrict du,
                 RealType* restrict d2u,
                 bool triangle = false);

  /** compute gradient
   */
  posT accumulateG(const valT* restrict du, const DisplRow& displ) const;
  /**@} */

  inline RealType ChiesaKEcorrection() { return KEcorr = j2_ke_corr_helper.computeKEcorr(); }

  inline RealType KECorrection() override { return KEcorr; }

  const std::vector<FT*>& getPairFunctions() const { return F; }
};

} // namespace qmcplusplus
#endif

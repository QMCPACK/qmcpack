//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EEIJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_EEIJASTROW_OPTIMIZED_SOA_H
#include "Configuration.h"
#if !defined(QMC_BUILD_SANDBOX_ONLY)
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#endif
#include "Particle/DistanceTable.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include <map>
#include <numeric>
#include <memory>

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief Specialization for three-body Jastrow function using multiple functors
 *
 *Each pair-type can have distinct function \f$u(r_{ij})\f$.
 *For electrons, distinct pair correlation functions are used
 *for spins up-up/down-down and up-down/down-up.
 */
template<class FT>
class JeeIOrbitalSoA : public WaveFunctionComponent
{
  ///type of each component U, dU, d2U;
  using valT = typename FT::real_type;
  ///element position type
  using posT = TinyVector<valT, OHMMS_DIM>;
  ///use the same container
  using DistRow  = DistanceTable::DistRow;
  using DisplRow = DistanceTable::DisplRow;
  ///table index for el-el
  const int ee_Table_ID_;
  ///table index for i-el
  const int ei_Table_ID_;
  //number of particles
  int Nelec, Nion;
  ///number of particles + padded
  size_t Nelec_padded;
  //number of groups of the target particleset
  int eGroups, iGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;
  ///diff value
  RealType DiffVal;

  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Uat, oldUk, newUk;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  using gContainer_type = VectorSoaContainer<valT, OHMMS_DIM>;
  gContainer_type dUat, olddUk, newdUk;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  Vector<valT> d2Uat, oldd2Uk, newd2Uk;
  /// current values during PbyP
  valT cur_Uat, cur_d2Uat;
  posT cur_dUat, dUat_temp;
  ///container for the Jastrow functions
  Array<FT*, 3> F;

  std::map<std::string, std::unique_ptr<FT>> J3Unique;
  //YYYY
  std::map<FT*, int> J3UniqueIndex;

  /// the cutoff for e-I pairs
  std::vector<valT> Ion_cutoff;
  /// the electrons around ions within the cutoff radius, grouped by species
  Array<std::vector<int>, 2> elecs_inside;
  Array<std::vector<valT>, 2> elecs_inside_dist;
  Array<std::vector<posT>, 2> elecs_inside_displ;
  /// the ids of ions within the cutoff radius of an electron on which a move is proposed
  std::vector<int> ions_nearby_old, ions_nearby_new;

  /// work buffer size
  size_t Nbuffer;
  /// compressed distances
  aligned_vector<valT> Distjk_Compressed, DistkI_Compressed, DistjI_Compressed;
  std::vector<int> DistIndice_k;
  /// compressed displacements
  gContainer_type Disp_jk_Compressed, Disp_jI_Compressed, Disp_kI_Compressed;
  /// work result buffer
  VectorSoaContainer<valT, 9> mVGL;

  // Used for evaluating derivatives with respect to the parameters
  Array<std::pair<int, int>, 3> VarOffset;
  Vector<RealType> dLogPsi;
  Array<PosType, 2> gradLogPsi;
  Array<RealType, 2> lapLogPsi;

  // Temporary store for parameter derivatives of functor
  // The first index is the functor index in J3Unique.  The second is the parameter index w.r.t. to that
  // functor
  std::vector<std::vector<RealType>> du_dalpha;
  std::vector<std::vector<PosType>> dgrad_dalpha;
  std::vector<std::vector<Tensor<RealType, 3>>> dhess_dalpha;

  void resizeWFOptVectors()
  {
    dLogPsi.resize(myVars.size());
    gradLogPsi.resize(myVars.size(), Nelec);
    lapLogPsi.resize(myVars.size(), Nelec);

    du_dalpha.resize(J3Unique.size());
    dgrad_dalpha.resize(J3Unique.size());
    dhess_dalpha.resize(J3Unique.size());

    int ifunc = 0;
    for (auto& j3UniquePair : J3Unique)
    {
      auto functorPtr           = j3UniquePair.second.get();
      J3UniqueIndex[functorPtr] = ifunc;
      const int numParams       = functorPtr->getNumParameters();
      du_dalpha[ifunc].resize(numParams);
      dgrad_dalpha[ifunc].resize(numParams);
      dhess_dalpha[ifunc].resize(numParams);
      ifunc++;
    }
  }

  /// compute G and L from internally stored data
  QTFull::RealType computeGL(ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L) const
  {
    for (int iat = 0; iat < Nelec; ++iat)
    {
      G[iat] += dUat[iat];
      L[iat] += d2Uat[iat];
    }
    return -0.5 * simd::accumulate_n(Uat.data(), Nelec, QTFull::RealType());
  }


public:
  ///alias FuncType
  using FuncType = FT;

  JeeIOrbitalSoA(const std::string& obj_name, const ParticleSet& ions, ParticleSet& elecs, bool is_master = false)
      : WaveFunctionComponent("JeeIOrbitalSoA", obj_name),
        ee_Table_ID_(elecs.addTable(elecs, DTModes::NEED_TEMP_DATA_ON_HOST | DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
        ei_Table_ID_(elecs.addTable(ions, DTModes::NEED_FULL_TABLE_ANYTIME | DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
        Ions(ions)
  {
    if (myName.empty())
      throw std::runtime_error("JeeIOrbitalSoA object name cannot be empty!");
    init(elecs);
  }

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& elecs) const override
  {
    auto eeIcopy = std::make_unique<JeeIOrbitalSoA<FT>>(myName, Ions, elecs, false);
    std::map<const FT*, FT*> fcmap;
    for (int iG = 0; iG < iGroups; iG++)
      for (int eG1 = 0; eG1 < eGroups; eG1++)
        for (int eG2 = 0; eG2 < eGroups; eG2++)
        {
          if (F(iG, eG1, eG2) == nullptr)
            continue;
          auto fit = fcmap.find(F(iG, eG1, eG2));
          if (fit == fcmap.end())
          {
            auto fc                = std::make_unique<FT>(*F(iG, eG1, eG2));
            fcmap[F(iG, eG1, eG2)] = fc.get();
            eeIcopy->addFunc(iG, eG1, eG2, std::move(fc));
          }
        }
    // Ye: I don't like the following memory allocated by default.
    eeIcopy->myVars.clear();
    eeIcopy->myVars.insertFrom(myVars);
    eeIcopy->VarOffset   = VarOffset;
    eeIcopy->Optimizable = Optimizable;
    return eeIcopy;
  }

  void init(ParticleSet& p)
  {
    Nelec        = p.getTotalNum();
    Nelec_padded = getAlignedSize<valT>(Nelec);
    Nion         = Ions.getTotalNum();
    iGroups      = Ions.getSpeciesSet().getTotalNum();
    eGroups      = p.groups();

    Uat.resize(Nelec);
    dUat.resize(Nelec);
    d2Uat.resize(Nelec);

    oldUk.resize(Nelec);
    olddUk.resize(Nelec);
    oldd2Uk.resize(Nelec);
    newUk.resize(Nelec);
    newdUk.resize(Nelec);
    newd2Uk.resize(Nelec);

    F.resize(iGroups, eGroups, eGroups);
    F = nullptr;
    elecs_inside.resize(eGroups, Nion);
    elecs_inside_dist.resize(eGroups, Nion);
    elecs_inside_displ.resize(eGroups, Nion);
    ions_nearby_old.resize(Nion);
    ions_nearby_new.resize(Nion);
    Ion_cutoff.resize(Nion, 0.0);

    //initialize buffers
    Nbuffer = Nelec;
    mVGL.resize(Nbuffer);
    Distjk_Compressed.resize(Nbuffer);
    DistjI_Compressed.resize(Nbuffer);
    DistkI_Compressed.resize(Nbuffer);
    Disp_jk_Compressed.resize(Nbuffer);
    Disp_jI_Compressed.resize(Nbuffer);
    Disp_kI_Compressed.resize(Nbuffer);
    DistIndice_k.resize(Nbuffer);
  }

  void addFunc(int iSpecies, int eSpecies1, int eSpecies2, std::unique_ptr<FT> j)
  {
    if (eSpecies1 == eSpecies2)
    {
      //if only up-up is specified, assume spin-unpolarized correlations
      if (eSpecies1 == 0)
        for (int eG1 = 0; eG1 < eGroups; eG1++)
          for (int eG2 = 0; eG2 < eGroups; eG2++)
          {
            if (F(iSpecies, eG1, eG2) == 0)
              F(iSpecies, eG1, eG2) = j.get();
          }
    }
    else
    {
      F(iSpecies, eSpecies1, eSpecies2) = j.get();
      F(iSpecies, eSpecies2, eSpecies1) = j.get();
    }
    if (j)
    {
      RealType rcut = 0.5 * j->cutoff_radius;
      for (int i = 0; i < Nion; i++)
        if (Ions.GroupID[i] == iSpecies)
          Ion_cutoff[i] = rcut;
    }
    else
    {
      APP_ABORT("JeeIOrbitalSoA::addFunc  Jastrow function pointer is NULL");
    }
    std::stringstream aname;
    aname << iSpecies << "_" << eSpecies1 << "_" << eSpecies2;
    J3Unique.emplace(aname.str(), std::move(j));
  }


  /** check that correlation information is complete
   */
  void check_complete()
  {
    //check that correlation pointers are either all 0 or all assigned
    bool complete = true;
    for (int i = 0; i < iGroups; ++i)
    {
      int nfilled = 0;
      bool partial;
      for (int e1 = 0; e1 < eGroups; ++e1)
        for (int e2 = 0; e2 < eGroups; ++e2)
          if (F(i, e1, e2) != 0)
            nfilled++;
      partial = nfilled > 0 && nfilled < eGroups * eGroups;
      if (partial)
        app_log() << "J3 eeI is missing correlation for ion " << i << std::endl;
      complete = complete && !partial;
    }
    if (!complete)
    {
      APP_ABORT("JeeIOrbitalSoA::check_complete  J3 eeI is missing correlation components\n  see preceding messages "
                "for details");
    }
    //first set radii
    for (int i = 0; i < Nion; ++i)
    {
      FT* f = F(Ions.GroupID[i], 0, 0);
      if (f != 0)
        Ion_cutoff[i] = .5 * f->cutoff_radius;
    }
    //then check radii
    bool all_radii_match = true;
    for (int i = 0; i < iGroups; ++i)
    {
      if (F(i, 0, 0) != 0)
      {
        bool radii_match = true;
        RealType rcut    = F(i, 0, 0)->cutoff_radius;
        for (int e1 = 0; e1 < eGroups; ++e1)
          for (int e2 = 0; e2 < eGroups; ++e2)
            radii_match = radii_match && F(i, e1, e2)->cutoff_radius == rcut;
        if (!radii_match)
          app_log() << "eeI functors for ion species " << i << " have different radii" << std::endl;
        all_radii_match = all_radii_match && radii_match;
      }
    }
    if (!all_radii_match)
    {
      APP_ABORT("JeeIOrbitalSoA::check_radii  J3 eeI are inconsistent for some ion species\n  see preceding messages "
                "for details");
    }
  }

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active) override
  {
    myVars.clear();

    for (auto& ftPair : J3Unique)
    {
      ftPair.second->checkInVariables(active);
      ftPair.second->checkInVariables(myVars);
    }
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.clear();

    for (auto& ftPair : J3Unique)
    {
      ftPair.second->myVars.getIndex(active);
      myVars.insertFrom(ftPair.second->myVars);
    }

    myVars.getIndex(active);
    const size_t NumVars = myVars.size();
    if (NumVars)
    {
      VarOffset.resize(iGroups, eGroups, eGroups);
      int varoffset = myVars.Index[0];
      for (int ig = 0; ig < iGroups; ig++)
        for (int jg = 0; jg < eGroups; jg++)
          for (int kg = 0; kg < eGroups; kg++)
          {
            FT* func_ijk = F(ig, jg, kg);
            if (func_ijk == nullptr)
              continue;
            VarOffset(ig, jg, kg).first  = func_ijk->myVars.Index.front() - varoffset;
            VarOffset(ig, jg, kg).second = func_ijk->myVars.Index.size() + VarOffset(ig, jg, kg).first;
          }
    }
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active) override
  {
    if (!Optimizable)
      return;

    for (auto& ftPair : J3Unique)
      ftPair.second->resetParameters(active);

    for (int i = 0; i < myVars.size(); ++i)
    {
      int ii = myVars.Index[i];
      if (ii >= 0)
        myVars[i] = active[ii];
    }
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os) override
  {
    for (auto& ftPair : J3Unique)
      ftPair.second->myVars.print(os);
  }

  void build_compact_list(const ParticleSet& P)
  {
    const auto& eI_dists  = P.getDistTableAB(ei_Table_ID_).getDistances();
    const auto& eI_displs = P.getDistTableAB(ei_Table_ID_).getDisplacements();

    for (int iat = 0; iat < Nion; ++iat)
      for (int jg = 0; jg < eGroups; ++jg)
      {
        elecs_inside(jg, iat).clear();
        elecs_inside_dist(jg, iat).clear();
        elecs_inside_displ(jg, iat).clear();
      }

    for (int jg = 0; jg < eGroups; ++jg)
      for (int jel = P.first(jg); jel < P.last(jg); jel++)
        for (int iat = 0; iat < Nion; ++iat)
          if (eI_dists[jel][iat] < Ion_cutoff[iat])
          {
            elecs_inside(jg, iat).push_back(jel);
            elecs_inside_dist(jg, iat).push_back(eI_dists[jel][iat]);
            elecs_inside_displ(jg, iat).push_back(eI_displs[jel][iat]);
          }
  }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    recompute(P);
    return log_value_ = computeGL(G, L);
  }

  PsiValueType ratio(ParticleSet& P, int iat) override
  {
    UpdateMode = ORB_PBYP_RATIO;

    const auto& eI_table = P.getDistTableAB(ei_Table_ID_);
    const auto& ee_table = P.getDistTableAA(ee_Table_ID_);
    cur_Uat = computeU(P, iat, P.GroupID[iat], eI_table.getTempDists(), ee_table.getTempDists(), ions_nearby_new);
    DiffVal = Uat[iat] - cur_Uat;
    return std::exp(static_cast<PsiValueType>(DiffVal));
  }

  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios) override
  {
    for (int k = 0; k < ratios.size(); ++k)
      ratios[k] = std::exp(Uat[VP.refPtcl] -
                           computeU(VP.refPS, VP.refPtcl, VP.refPS.GroupID[VP.refPtcl],
                                    VP.getDistTableAB(ei_Table_ID_).getDistRow(k),
                                    VP.getDistTableAB(ee_Table_ID_).getDistRow(k), ions_nearby_old));
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override
  {
    const auto& eI_table = P.getDistTableAB(ei_Table_ID_);
    const auto& eI_dists = eI_table.getDistances();
    const auto& ee_table = P.getDistTableAA(ee_Table_ID_);

    for (int jg = 0; jg < eGroups; ++jg)
    {
      const valT sumU = computeU(P, -1, jg, eI_table.getTempDists(), ee_table.getTempDists(), ions_nearby_new);

      for (int j = P.first(jg); j < P.last(jg); ++j)
      {
        // remove self-interaction
        valT Uself(0);
        for (int iat = 0; iat < Nion; ++iat)
        {
          const valT& r_Ij = eI_table.getTempDists()[iat];
          const valT& r_Ik = eI_dists[j][iat];
          if (r_Ij < Ion_cutoff[iat] && r_Ik < Ion_cutoff[iat])
          {
            const int ig = Ions.GroupID[iat];
            Uself += F(ig, jg, jg)->evaluate(ee_table.getTempDists()[j], r_Ij, r_Ik);
          }
        }
        ratios[j] = std::exp(Uat[j] + Uself - sumU);
      }
    }
  }

  GradType evalGrad(ParticleSet& P, int iat) override { return GradType(dUat[iat]); }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    UpdateMode = ORB_PBYP_PARTIAL;

    const auto& eI_table = P.getDistTableAB(ei_Table_ID_);
    const auto& ee_table = P.getDistTableAA(ee_Table_ID_);
    computeU3(P, iat, eI_table.getTempDists(), eI_table.getTempDispls(), ee_table.getTempDists(),
              ee_table.getTempDispls(), cur_Uat, cur_dUat, cur_d2Uat, newUk, newdUk, newd2Uk, ions_nearby_new);
    DiffVal = Uat[iat] - cur_Uat;
    grad_iat += cur_dUat;
    return std::exp(static_cast<PsiValueType>(DiffVal));
  }

  inline void restore(int iat) override {}

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override
  {
    const auto& eI_table = P.getDistTableAB(ei_Table_ID_);
    const auto& ee_table = P.getDistTableAA(ee_Table_ID_);
    // get the old value, grad, lapl
    computeU3(P, iat, eI_table.getDistRow(iat), eI_table.getDisplRow(iat), ee_table.getOldDists(),
              ee_table.getOldDispls(), Uat[iat], dUat_temp, d2Uat[iat], oldUk, olddUk, oldd2Uk, ions_nearby_old);
    if (UpdateMode == ORB_PBYP_RATIO)
    { //ratio-only during the move; need to compute derivatives
      computeU3(P, iat, eI_table.getTempDists(), eI_table.getTempDispls(), ee_table.getTempDists(),
                ee_table.getTempDispls(), cur_Uat, cur_dUat, cur_d2Uat, newUk, newdUk, newd2Uk, ions_nearby_new);
    }

#pragma omp simd
    for (int jel = 0; jel < Nelec; jel++)
    {
      Uat[jel] += newUk[jel] - oldUk[jel];
      d2Uat[jel] += newd2Uk[jel] - oldd2Uk[jel];
    }
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
    {
      valT* restrict save_g      = dUat.data(idim);
      const valT* restrict new_g = newdUk.data(idim);
      const valT* restrict old_g = olddUk.data(idim);
#pragma omp simd aligned(save_g, new_g, old_g : QMC_SIMD_ALIGNMENT)
      for (int jel = 0; jel < Nelec; jel++)
        save_g[jel] += new_g[jel] - old_g[jel];
    }

    log_value_ += Uat[iat] - cur_Uat;
    Uat[iat]   = cur_Uat;
    dUat(iat)  = cur_dUat;
    d2Uat[iat] = cur_d2Uat;

    const int ig = P.GroupID[iat];
    // update compact list elecs_inside
    // if the old position exists in elecs_inside
    for (int iind = 0; iind < ions_nearby_old.size(); iind++)
    {
      int jat         = ions_nearby_old[iind];
      auto iter       = std::find(elecs_inside(ig, jat).begin(), elecs_inside(ig, jat).end(), iat);
      auto iter_dist  = elecs_inside_dist(ig, jat).begin() + std::distance(elecs_inside(ig, jat).begin(), iter);
      auto iter_displ = elecs_inside_displ(ig, jat).begin() + std::distance(elecs_inside(ig, jat).begin(), iter);
// sentinel code
#ifndef NDEBUG
      if (iter == elecs_inside(ig, jat).end())
      {
        std::cerr << std::setprecision(std::numeric_limits<valT>::digits10 + 1) << "updating electron iat = " << iat
                  << " near ion " << jat << " dist " << eI_table.getDistRow(iat)[jat] << std::endl;
        throw std::runtime_error("BUG electron not found in elecs_inside");
      }
      else if (std::abs(eI_table.getDistRow(iat)[jat] - *iter_dist) >= std::numeric_limits<valT>::epsilon())
      {
        std::cerr << std::setprecision(std::numeric_limits<valT>::digits10 + 1) << "inconsistent electron iat = " << iat
                  << " near ion " << jat << " dist " << eI_table.getDistRow(iat)[jat]
                  << " stored value = " << *iter_dist << std::endl;
        throw std::runtime_error("BUG eI distance stored value elecs_inside_dist not matching distance table");
      }
#endif

      if (eI_table.getTempDists()[jat] < Ion_cutoff[jat]) // the new position is still inside
      {
        *iter_dist                                                      = eI_table.getTempDists()[jat];
        *iter_displ                                                     = eI_table.getTempDispls()[jat];
        *std::find(ions_nearby_new.begin(), ions_nearby_new.end(), jat) = -1;
      }
      else
      {
        *iter = elecs_inside(ig, jat).back();
        elecs_inside(ig, jat).pop_back();
        *iter_dist = elecs_inside_dist(ig, jat).back();
        elecs_inside_dist(ig, jat).pop_back();
        *iter_displ = elecs_inside_displ(ig, jat).back();
        elecs_inside_displ(ig, jat).pop_back();
      }
    }

    // if the old position doesn't exist in elecs_inside but the new position do
    for (int iind = 0; iind < ions_nearby_new.size(); iind++)
    {
      int jat = ions_nearby_new[iind];
      if (jat >= 0)
      {
        elecs_inside(ig, jat).push_back(iat);
        elecs_inside_dist(ig, jat).push_back(eI_table.getTempDists()[jat]);
        elecs_inside_displ(ig, jat).push_back(eI_table.getTempDispls()[jat]);
      }
    }
  }

  inline void recompute(const ParticleSet& P) override
  {
    const auto& eI_table = P.getDistTableAB(ei_Table_ID_);
    const auto& ee_table = P.getDistTableAA(ee_Table_ID_);

    build_compact_list(P);

    for (int jel = 0; jel < Nelec; ++jel)
    {
      computeU3(P, jel, eI_table.getDistRow(jel), eI_table.getDisplRow(jel), ee_table.getDistRow(jel),
                ee_table.getDisplRow(jel), Uat[jel], dUat_temp, d2Uat[jel], newUk, newdUk, newd2Uk, ions_nearby_new,
                true);
      dUat(jel) = dUat_temp;
// add the contribution from the upper triangle
#pragma omp simd
      for (int kel = 0; kel < jel; kel++)
      {
        Uat[kel] += newUk[kel];
        d2Uat[kel] += newd2Uk[kel];
      }
      for (int idim = 0; idim < OHMMS_DIM; ++idim)
      {
        valT* restrict save_g      = dUat.data(idim);
        const valT* restrict new_g = newdUk.data(idim);
#pragma omp simd aligned(save_g, new_g : QMC_SIMD_ALIGNMENT)
        for (int kel = 0; kel < jel; kel++)
          save_g[kel] += new_g[kel];
      }
    }
  }

  inline valT computeU(const ParticleSet& P,
                       int jel,
                       int jg,
                       const DistRow& distjI,
                       const DistRow& distjk,
                       std::vector<int>& ions_nearby)
  {
    ions_nearby.clear();
    for (int iat = 0; iat < Nion; ++iat)
      if (distjI[iat] < Ion_cutoff[iat])
        ions_nearby.push_back(iat);

    valT Uj = valT(0);
    for (int kg = 0; kg < eGroups; ++kg)
    {
      int kel_counter = 0;
      for (int iind = 0; iind < ions_nearby.size(); ++iind)
      {
        const int iat   = ions_nearby[iind];
        const int ig    = Ions.GroupID[iat];
        const valT r_jI = distjI[iat];
        for (int kind = 0; kind < elecs_inside(kg, iat).size(); kind++)
        {
          const int kel = elecs_inside(kg, iat)[kind];
          if (kel != jel)
          {
            DistkI_Compressed[kel_counter] = elecs_inside_dist(kg, iat)[kind];
            Distjk_Compressed[kel_counter] = distjk[kel];
            DistjI_Compressed[kel_counter] = r_jI;
            kel_counter++;
            if (kel_counter == Nbuffer)
            {
              const FT& feeI(*F(ig, jg, kg));
              Uj += feeI.evaluateV(kel_counter, Distjk_Compressed.data(), DistjI_Compressed.data(),
                                   DistkI_Compressed.data());
              kel_counter = 0;
            }
          }
        }
        if ((iind + 1 == ions_nearby.size() || ig != Ions.GroupID[ions_nearby[iind + 1]]) && kel_counter > 0)
        {
          const FT& feeI(*F(ig, jg, kg));
          Uj +=
              feeI.evaluateV(kel_counter, Distjk_Compressed.data(), DistjI_Compressed.data(), DistkI_Compressed.data());
          kel_counter = 0;
        }
      }
    }
    return Uj;
  }

  inline void computeU3_engine(const ParticleSet& P,
                               const FT& feeI,
                               int kel_counter,
                               valT& Uj,
                               posT& dUj,
                               valT& d2Uj,
                               Vector<valT>& Uk,
                               gContainer_type& dUk,
                               Vector<valT>& d2Uk)
  {
    constexpr valT czero(0);
    constexpr valT cone(1);
    constexpr valT ctwo(2);
    constexpr valT lapfac = OHMMS_DIM - cone;

    valT* restrict val     = mVGL.data(0);
    valT* restrict gradF0  = mVGL.data(1);
    valT* restrict gradF1  = mVGL.data(2);
    valT* restrict gradF2  = mVGL.data(3);
    valT* restrict hessF00 = mVGL.data(4);
    valT* restrict hessF11 = mVGL.data(5);
    valT* restrict hessF22 = mVGL.data(6);
    valT* restrict hessF01 = mVGL.data(7);
    valT* restrict hessF02 = mVGL.data(8);

    feeI.evaluateVGL(kel_counter, Distjk_Compressed.data(), DistjI_Compressed.data(), DistkI_Compressed.data(), val,
                     gradF0, gradF1, gradF2, hessF00, hessF11, hessF22, hessF01, hessF02);

    // compute the contribution to jel, kel
    Uj               = simd::accumulate_n(val, kel_counter, Uj);
    valT gradF0_sum  = simd::accumulate_n(gradF0, kel_counter, czero);
    valT gradF1_sum  = simd::accumulate_n(gradF1, kel_counter, czero);
    valT hessF00_sum = simd::accumulate_n(hessF00, kel_counter, czero);
    valT hessF11_sum = simd::accumulate_n(hessF11, kel_counter, czero);
    d2Uj -= hessF00_sum + hessF11_sum + lapfac * (gradF0_sum + gradF1_sum);
    std::fill_n(hessF11, kel_counter, czero);
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
    {
      valT* restrict jk = Disp_jk_Compressed.data(idim);
      valT* restrict jI = Disp_jI_Compressed.data(idim);
      valT* restrict kI = Disp_kI_Compressed.data(idim);
      valT dUj_x(0);
#pragma omp simd aligned(gradF0, gradF1, gradF2, hessF11, jk, jI, kI : QMC_SIMD_ALIGNMENT) reduction(+ : dUj_x)
      for (int kel_index = 0; kel_index < kel_counter; kel_index++)
      {
        // recycle hessF11
        hessF11[kel_index] += kI[kel_index] * jk[kel_index];
        dUj_x += gradF1[kel_index] * jI[kel_index];
        // destroy jk, kI
        const valT temp = jk[kel_index] * gradF0[kel_index];
        dUj_x += temp;
        jk[kel_index] *= jI[kel_index];
        kI[kel_index] = kI[kel_index] * gradF2[kel_index] - temp;
      }
      dUj[idim] += dUj_x;

      valT* restrict jk0 = Disp_jk_Compressed.data(0);
      if (idim > 0)
      {
#pragma omp simd aligned(jk, jk0 : QMC_SIMD_ALIGNMENT)
        for (int kel_index = 0; kel_index < kel_counter; kel_index++)
          jk0[kel_index] += jk[kel_index];
      }

      valT* restrict dUk_x = dUk.data(idim);
      for (int kel_index = 0; kel_index < kel_counter; kel_index++)
        dUk_x[DistIndice_k[kel_index]] += kI[kel_index];
    }
    valT sum(0);
    valT* restrict jk0 = Disp_jk_Compressed.data(0);
#pragma omp simd aligned(jk0, hessF01 : QMC_SIMD_ALIGNMENT) reduction(+ : sum)
    for (int kel_index = 0; kel_index < kel_counter; kel_index++)
      sum += hessF01[kel_index] * jk0[kel_index];
    d2Uj -= ctwo * sum;

#pragma omp simd aligned(hessF00, hessF22, gradF0, gradF2, hessF02, hessF11 : QMC_SIMD_ALIGNMENT)
    for (int kel_index = 0; kel_index < kel_counter; kel_index++)
      hessF00[kel_index] = hessF00[kel_index] + hessF22[kel_index] + lapfac * (gradF0[kel_index] + gradF2[kel_index]) -
          ctwo * hessF02[kel_index] * hessF11[kel_index];

    for (int kel_index = 0; kel_index < kel_counter; kel_index++)
    {
      const int kel = DistIndice_k[kel_index];
      Uk[kel] += val[kel_index];
      d2Uk[kel] -= hessF00[kel_index];
    }
  }

  inline void computeU3(const ParticleSet& P,
                        int jel,
                        const DistRow& distjI,
                        const DisplRow& displjI,
                        const DistRow& distjk,
                        const DisplRow& displjk,
                        valT& Uj,
                        posT& dUj,
                        valT& d2Uj,
                        Vector<valT>& Uk,
                        gContainer_type& dUk,
                        Vector<valT>& d2Uk,
                        std::vector<int>& ions_nearby,
                        bool triangle = false)
  {
    constexpr valT czero(0);

    Uj   = czero;
    dUj  = posT();
    d2Uj = czero;

    const int jg = P.GroupID[jel];

    const int kelmax = triangle ? jel : Nelec;
    std::fill_n(Uk.data(), kelmax, czero);
    std::fill_n(d2Uk.data(), kelmax, czero);
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
      std::fill_n(dUk.data(idim), kelmax, czero);

    ions_nearby.clear();
    for (int iat = 0; iat < Nion; ++iat)
      if (distjI[iat] < Ion_cutoff[iat])
        ions_nearby.push_back(iat);

    for (int kg = 0; kg < eGroups; ++kg)
    {
      int kel_counter = 0;
      for (int iind = 0; iind < ions_nearby.size(); ++iind)
      {
        const int iat      = ions_nearby[iind];
        const int ig       = Ions.GroupID[iat];
        const valT r_jI    = distjI[iat];
        const posT disp_Ij = displjI[iat];
        for (int kind = 0; kind < elecs_inside(kg, iat).size(); kind++)
        {
          const int kel = elecs_inside(kg, iat)[kind];
          if (kel < kelmax && kel != jel)
          {
            DistkI_Compressed[kel_counter]  = elecs_inside_dist(kg, iat)[kind];
            DistjI_Compressed[kel_counter]  = r_jI;
            Distjk_Compressed[kel_counter]  = distjk[kel];
            Disp_kI_Compressed(kel_counter) = elecs_inside_displ(kg, iat)[kind];
            Disp_jI_Compressed(kel_counter) = disp_Ij;
            Disp_jk_Compressed(kel_counter) = displjk[kel];
            DistIndice_k[kel_counter]       = kel;
            kel_counter++;
            if (kel_counter == Nbuffer)
            {
              const FT& feeI(*F(ig, jg, kg));
              computeU3_engine(P, feeI, kel_counter, Uj, dUj, d2Uj, Uk, dUk, d2Uk);
              kel_counter = 0;
            }
          }
        }
        if ((iind + 1 == ions_nearby.size() || ig != Ions.GroupID[ions_nearby[iind + 1]]) && kel_counter > 0)
        {
          const FT& feeI(*F(ig, jg, kg));
          computeU3_engine(P, feeI, kel_counter, Uj, dUj, d2Uj, Uk, dUk, d2Uk);
          kel_counter = 0;
        }
      }
    }
  }

  inline void registerData(ParticleSet& P, WFBufferType& buf) override
  {
    if (Bytes_in_WFBuffer == 0)
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Uat.begin(), Uat.end());
      buf.add(dUat.data(), dUat.end());
      buf.add(d2Uat.begin(), d2Uat.end());
      Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
      // free local space
      Uat.free();
      dUat.free();
      d2Uat.free();
    }
    else
    {
      buf.forward(Bytes_in_WFBuffer);
    }
  }

  inline LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override
  {
    log_value_ = computeGL(P.G, P.L);
    buf.forward(Bytes_in_WFBuffer);
    return log_value_;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override
  {
    Uat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    dUat.attachReference(Nelec, Nelec_padded, buf.lendReference<valT>(Nelec_padded * OHMMS_DIM));
    d2Uat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    build_compact_list(P);
  }

  LogValueType evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient& G,
                          ParticleSet::ParticleLaplacian& L,
                          bool fromscratch = false) override
  {
    return log_value_ = computeGL(G, L);
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
    resizeWFOptVectors();

    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (optvars.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }

    if (recalculate)
    {
      constexpr valT czero(0);
      constexpr valT cone(1);
      constexpr valT cminus(-1);
      constexpr valT ctwo(2);
      constexpr valT lapfac = OHMMS_DIM - cone;

      const auto& ee_table  = P.getDistTableAA(ee_Table_ID_);
      const auto& ee_dists  = ee_table.getDistances();
      const auto& ee_displs = ee_table.getDisplacements();

      build_compact_list(P);

      dLogPsi    = czero;
      gradLogPsi = PosType();
      lapLogPsi  = czero;

      for (int iat = 0; iat < Nion; ++iat)
      {
        const int ig = Ions.GroupID[iat];
        for (int jg = 0; jg < eGroups; ++jg)
          for (int jind = 0; jind < elecs_inside(jg, iat).size(); jind++)
          {
            const int jel       = elecs_inside(jg, iat)[jind];
            const valT r_Ij     = elecs_inside_dist(jg, iat)[jind];
            const posT disp_Ij  = cminus * elecs_inside_displ(jg, iat)[jind];
            const valT r_Ij_inv = cone / r_Ij;

            for (int kg = 0; kg < eGroups; ++kg)
              for (int kind = 0; kind < elecs_inside(kg, iat).size(); kind++)
              {
                const int kel = elecs_inside(kg, iat)[kind];
                if (kel < jel)
                {
                  const valT r_Ik     = elecs_inside_dist(kg, iat)[kind];
                  const posT disp_Ik  = cminus * elecs_inside_displ(kg, iat)[kind];
                  const valT r_Ik_inv = cone / r_Ik;

                  const valT r_jk     = ee_dists[jel][kel];
                  const posT disp_jk  = ee_displs[jel][kel];
                  const valT r_jk_inv = cone / r_jk;

                  FT& func = *F(ig, jg, kg);
                  int idx  = J3UniqueIndex[F(ig, jg, kg)];
                  func.evaluateDerivatives(r_jk, r_Ij, r_Ik, du_dalpha[idx], dgrad_dalpha[idx], dhess_dalpha[idx]);
                  int first                               = VarOffset(ig, jg, kg).first;
                  int last                                = VarOffset(ig, jg, kg).second;
                  std::vector<RealType>& dlog             = du_dalpha[idx];
                  std::vector<PosType>& dgrad             = dgrad_dalpha[idx];
                  std::vector<Tensor<RealType, 3>>& dhess = dhess_dalpha[idx];

                  for (int p = first, ip = 0; p < last; p++, ip++)
                  {
                    RealType& dval          = dlog[ip];
                    PosType& dg             = dgrad[ip];
                    Tensor<RealType, 3>& dh = dhess[ip];

                    dg[0] *= r_jk_inv;
                    dg[1] *= r_Ij_inv;
                    dg[2] *= r_Ik_inv;

                    PosType gr_ee = dg[0] * disp_jk;

                    gradLogPsi(p, jel) -= dg[1] * disp_Ij - gr_ee;
                    lapLogPsi(p, jel) -=
                        (dh(0, 0) + lapfac * dg[0] - ctwo * dh(0, 1) * dot(disp_jk, disp_Ij) * r_jk_inv * r_Ij_inv +
                         dh(1, 1) + lapfac * dg[1]);

                    gradLogPsi(p, kel) -= dg[2] * disp_Ik + gr_ee;
                    lapLogPsi(p, kel) -=
                        (dh(0, 0) + lapfac * dg[0] + ctwo * dh(0, 2) * dot(disp_jk, disp_Ik) * r_jk_inv * r_Ik_inv +
                         dh(2, 2) + lapfac * dg[2]);

                    dLogPsi[p] -= dval;
                  }
                }
              }
          }
      }

      for (int k = 0; k < myVars.size(); ++k)
      {
        int kk = myVars.where(k);
        if (kk < 0)
          continue;
        dlogpsi[kk]  = (ValueType)dLogPsi[k];
        RealType sum = 0.0;
        for (int i = 0; i < Nelec; i++)
        {
#if defined(QMC_COMPLEX)
          sum -= 0.5 * lapLogPsi(k, i);
          for (int jdim = 0; jdim < OHMMS_DIM; ++jdim)
            sum -= P.G[i][jdim].real() * gradLogPsi(k, i)[jdim];
#else
          sum -= 0.5 * lapLogPsi(k, i) + dot(P.G[i], gradLogPsi(k, i));
#endif
        }
        dhpsioverpsi[kk] = (ValueType)sum;
      }
    }
  }

  inline GradType evalGradSource(ParticleSet& P, ParticleSet& source, int isrc) override
  {
    ParticleSet::ParticleGradient tempG;
    ParticleSet::ParticleLaplacian tempL;
    tempG.resize(P.getTotalNum());
    tempL.resize(P.getTotalNum());
    QTFull::RealType delta = 0.00001;
    QTFull::RealType c1    = 1.0 / delta / 2.0;
    QTFull::RealType c2    = 1.0 / delta / delta;

    GradType g_return(0.0);
    // GRAD TEST COMPUTATION
    PosType rI = source.R[isrc];
    for (int iondim = 0; iondim < 3; iondim++)
    {
      source.R[isrc][iondim] = rI[iondim] + delta;
      source.update();
      P.update();

      LogValueType log_p = evaluateLog(P, tempG, tempL);

      source.R[isrc][iondim] = rI[iondim] - delta;
      source.update();
      P.update();
      LogValueType log_m = evaluateLog(P, tempG, tempL);

      QTFull::RealType log_p_r(0.0), log_m_r(0.0);

      log_p_r = log_p.real();
      log_m_r = log_m.real();
      //symmetric finite difference formula for gradient.
      g_return[iondim] = c1 * (log_p_r - log_m_r);

      //reset everything to how it was.
      source.R[isrc][iondim] = rI[iondim];
    }
    // this last one makes sure the distance tables and internal neighbourlist correspond to unperturbed source.
    source.update();
    P.update();
    build_compact_list(P);
    return g_return;
  }

  inline GradType evalGradSource(ParticleSet& P,
                                 ParticleSet& source,
                                 int isrc,
                                 TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
                                 TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) override
  {
    ParticleSet::ParticleGradient Gp, Gm, dG;
    ParticleSet::ParticleLaplacian Lp, Lm, dL;
    Gp.resize(P.getTotalNum());
    Gm.resize(P.getTotalNum());
    dG.resize(P.getTotalNum());
    Lp.resize(P.getTotalNum());
    Lm.resize(P.getTotalNum());
    dL.resize(P.getTotalNum());

    QTFull::RealType delta = 0.00001;
    QTFull::RealType c1    = 1.0 / delta / 2.0;
    QTFull::RealType c2    = 1.0 / delta / delta;
    GradType g_return(0.0);
    // GRAD TEST COMPUTATION
    PosType rI = source.R[isrc];
    for (int iondim = 0; iondim < 3; iondim++)
    {
      Lp                     = 0;
      Gp                     = 0;
      Lm                     = 0;
      Gm                     = 0;
      source.R[isrc][iondim] = rI[iondim] + delta;
      source.update();
      P.update();

      LogValueType log_p = evaluateLog(P, Gp, Lp);

      source.R[isrc][iondim] = rI[iondim] - delta;
      source.update();
      P.update();
      LogValueType log_m = evaluateLog(P, Gm, Lm);
      QTFull::RealType log_p_r(0.0), log_m_r(0.0);

      log_p_r = log_p.real();
      log_m_r = log_m.real();
      dG      = Gp - Gm;
      dL      = Lp - Lm;
      //symmetric finite difference formula for gradient.
      g_return[iondim] = c1 * (log_p_r - log_m_r);
      grad_grad[iondim] += c1 * dG;
      lapl_grad[iondim] += c1 * dL;
      //reset everything to how it was.
      source.R[isrc][iondim] = rI[iondim];
    }
    // this last one makes sure the distance tables and internal neighbourlist correspond to unperturbed source.
    source.update();
    P.update();
    build_compact_list(P);
    return g_return;
  }
};

} // namespace qmcplusplus
#endif

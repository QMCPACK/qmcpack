//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 */
template<class FT>
class DiffTwoBodyJastrowOrbital : public DiffWaveFunctionComponent
{
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///number of groups, e.g., for the up/down electrons
  int NumGroups;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  std::vector<FT*> F;
  /// e-e table ID
  const int my_table_ID_;
  /// Map indices from subcomponent variables to component variables
  std::vector<std::pair<int, int>> OffSet;
  Vector<RealType> dLogPsi;
  std::vector<GradVectorType*> gradLogPsi;
  std::vector<ValueVectorType*> lapLogPsi;
  std::map<std::string, std::unique_ptr<FT>> J2Unique;

public:
  // return for testing
  const std::vector<FT*>& getPairFunctions() const { return F; }

  ///constructor
  DiffTwoBodyJastrowOrbital(ParticleSet& p) : NumVars(0), my_table_ID_(p.addTable(p))
  {
    NumPtcls  = p.getTotalNum();
    NumGroups = p.groups();
    F.resize(NumGroups * NumGroups, 0);
  }

  ~DiffTwoBodyJastrowOrbital() override
  {
    delete_iter(gradLogPsi.begin(), gradLogPsi.end());
    delete_iter(lapLogPsi.begin(), lapLogPsi.end());
  }

  // Accessors for unit testing
  std::pair<int, int> getComponentOffset(int index) { return OffSet.at(index); }

  opt_variables_type& getComponentVars() { return myVars; }


  void addFunc(int ia, int ib, std::unique_ptr<FT> j)
  {
    // make all pair terms equal to uu initially
    //   in case some terms are not provided explicitly
    if (ia == ib)
    {
      if (ia == 0) //first time, assign everything
      {
        int ij = 0;
        for (int ig = 0; ig < NumGroups; ++ig)
          for (int jg = 0; jg < NumGroups; ++jg, ++ij)
            if (F[ij] == nullptr)
              F[ij] = j.get();
      }
      else
        F[ia * NumGroups + ib] = j.get();
    }
    else
    {
      // a very special case, 1 particle of each type (e.g. 1 up + 1 down)
      // uu/dd/etc. was prevented by the builder
      if (NumPtcls == NumGroups)
        for (int ig = 0; ig < NumGroups; ++ig)
          F[ig * NumGroups + ig] = j.get();
      // generic case
      F[ia * NumGroups + ib] = j.get();
      F[ib * NumGroups + ia] = j.get();
    }
    std::stringstream aname;
    aname << ia << ib;
    J2Unique[aname.str()] = std::move(j);
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active) override
  {
    auto it     = J2Unique.begin();
    auto it_end = J2Unique.end();
    while (it != it_end)
    {
      (*it++).second->resetParameters(active);
    }
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.clear();
    auto it     = J2Unique.begin();
    auto it_end = J2Unique.end();
    while (it != it_end)
    {
      (*it).second->myVars.getIndex(active);
      myVars.insertFrom((*it).second->myVars);
      ++it;
    }
    // Remove inactive variables so the mappings are correct
    myVars.removeInactive();

    myVars.getIndex(active);
    NumVars = myVars.size();

    //myVars.print(std::cout);

    if (NumVars && dLogPsi.size() == 0)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars, 0);
      lapLogPsi.resize(NumVars, 0);
      for (int i = 0; i < NumVars; ++i)
      {
        gradLogPsi[i] = new GradVectorType(NumPtcls);
        lapLogPsi[i]  = new ValueVectorType(NumPtcls);
      }
      OffSet.resize(F.size());

      // Find first active variable for the starting offset
      int varoffset = -1;
      for (int i = 0; i < myVars.size(); i++)
      {
        varoffset = myVars.Index[i];
        if (varoffset != -1)
          break;
      }

      for (int i = 0; i < F.size(); ++i)
      {
        if (F[i] && F[i]->myVars.Index.size())
        {
          OffSet[i].first  = F[i]->myVars.Index.front() - varoffset;
          OffSet[i].second = F[i]->myVars.Index.size() + OffSet[i].first;
        }
        else
        {
          OffSet[i].first = OffSet[i].second = -1;
        }
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
    if (myVars.size() == 0)
      return;
    evaluateDerivativesWF(P, active, dlogpsi);
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (active.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate)
    {
      for (int k = 0; k < myVars.size(); ++k)
      {
        int kk = myVars.where(k);
        if (kk < 0)
          continue;
        if (rcsingles[k])
        {
          dhpsioverpsi[kk] = -RealType(0.5) * ValueType(Sum(*lapLogPsi[k])) - ValueType(Dot(P.G, *gradLogPsi[k]));
        }
      }
    }
  }

  void evaluateDerivativesWF(ParticleSet& P, const opt_variables_type& active, std::vector<ValueType>& dlogpsi) override
  {
    if (myVars.size() == 0)
      return;
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(), false);
    for (int k = 0; k < myVars.size(); ++k)
    {
      int kk = myVars.where(k);
      if (kk < 0)
        continue;
      if (active.recompute(kk))
        recalculate = true;
      rcsingles[k] = true;
    }
    if (recalculate)
    {
      ///precomputed recalculation switch
      std::vector<bool> RecalcSwitch(F.size(), false);
      for (int i = 0; i < F.size(); ++i)
      {
        if (OffSet[i].first < 0)
        {
          // nothing to optimize
          RecalcSwitch[i] = false;
        }
        else
        {
          bool recalcFunc(false);
          for (int rcs = OffSet[i].first; rcs < OffSet[i].second; rcs++)
            if (rcsingles[rcs] == true)
              recalcFunc = true;
          RecalcSwitch[i] = recalcFunc;
        }
      }
      dLogPsi = 0.0;
      for (int p = 0; p < NumVars; ++p)
        (*gradLogPsi[p]) = 0.0;
      for (int p = 0; p < NumVars; ++p)
        (*lapLogPsi[p]) = 0.0;
      std::vector<TinyVector<RealType, 3>> derivs(NumVars);
      const auto& d_table = P.getDistTableAA(my_table_ID_);
      constexpr RealType cone(1);
      constexpr RealType lapfac(OHMMS_DIM - cone);
      const size_t n  = d_table.sources();
      const size_t ng = P.groups();
      for (size_t i = 1; i < n; ++i)
      {
        const size_t ig   = P.GroupID[i] * ng;
        const auto& dist  = d_table.getDistRow(i);
        const auto& displ = d_table.getDisplRow(i);
        for (size_t j = 0; j < i; ++j)
        {
          const size_t ptype = ig + P.GroupID[j];
          if (RecalcSwitch[ptype])
          {
            std::fill(derivs.begin(), derivs.end(), 0.0);
            if (!F[ptype]->evaluateDerivatives(dist[j], derivs))
              continue;
            RealType rinv(cone / dist[j]);
            PosType dr(displ[j]);
            for (int p = OffSet[ptype].first, ip = 0; p < OffSet[ptype].second; ++p, ++ip)
            {
              RealType dudr(rinv * derivs[ip][1]);
              RealType lap(derivs[ip][2] + lapfac * dudr);
              //RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
              PosType gr(dudr * dr);
              dLogPsi[p] -= derivs[ip][0];
              (*gradLogPsi[p])[i] += gr;
              (*gradLogPsi[p])[j] -= gr;
              (*lapLogPsi[p])[i] -= lap;
              (*lapLogPsi[p])[j] -= lap;
            }
          }
        }
      }
      for (int k = 0; k < myVars.size(); ++k)
      {
        int kk = myVars.where(k);
        if (kk < 0)
          continue;
        if (rcsingles[k])
        {
          dlogpsi[kk] = dLogPsi[k];
        }
        //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
      }
    }
  }

  std::unique_ptr<DiffWaveFunctionComponent> makeClone(ParticleSet& tqp) const override
  {
    auto j2copy = std::make_unique<DiffTwoBodyJastrowOrbital<FT>>(tqp);
    std::map<const FT*, FT*> fcmap;
    for (int ig = 0; ig < NumGroups; ++ig)
      for (int jg = ig; jg < NumGroups; ++jg)
      {
        int ij = ig * NumGroups + jg;
        if (F[ij] == nullptr)
          continue;
        auto fit = fcmap.find(F[ij]);
        if (fit == fcmap.end())
        {
          auto fc      = std::make_unique<FT>(*F[ij]);
          fcmap[F[ij]] = fc.get();
          j2copy->addFunc(ig, jg, std::move(fc));
        }
      }
    j2copy->myVars.clear();
    j2copy->myVars.insertFrom(myVars);
    j2copy->NumVars   = NumVars;
    j2copy->NumPtcls  = NumPtcls;
    j2copy->NumGroups = NumGroups;
    j2copy->dLogPsi.resize(NumVars);
    j2copy->gradLogPsi.resize(NumVars, 0);
    j2copy->lapLogPsi.resize(NumVars, 0);
    for (int i = 0; i < NumVars; ++i)
    {
      j2copy->gradLogPsi[i] = new GradVectorType(NumPtcls);
      j2copy->lapLogPsi[i]  = new ValueVectorType(NumPtcls);
    }
    j2copy->OffSet = OffSet;
    return j2copy;
  }
};
} // namespace qmcplusplus
#endif

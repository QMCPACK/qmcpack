//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DIFFERENTIAL_ONEBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_ONEBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"
#include "Particle/DistanceTableData.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"


namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 */
template<class FT>
class DiffOneBodyJastrowOrbital : public DiffWaveFunctionComponent
{
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///index of the table
  const int myTableIndex;
  ///reference to the ions
  const ParticleSet& CenterRef;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  std::vector<FT*> Fs;
  ///container for the unique Jastrow functions
  std::vector<FT*> Funique;
  std::vector<std::pair<int, int>> OffSet;
  Vector<RealType> dLogPsi;
  std::vector<GradVectorType*> gradLogPsi;
  std::vector<ValueVectorType*> lapLogPsi;

public:
  ///constructor
  DiffOneBodyJastrowOrbital(const ParticleSet& centers, ParticleSet& els)
      : NumVars(0), myTableIndex(els.addTable(centers)), CenterRef(centers)
  {
    NumPtcls = els.getTotalNum();
  }

  ~DiffOneBodyJastrowOrbital() override
  {
    delete_iter(gradLogPsi.begin(), gradLogPsi.end());
    delete_iter(lapLogPsi.begin(), lapLogPsi.end());
  }

  /** Add a radial functor for a group
   * @param source_type group index of the center species
   * @param afunc radial functor
   */
  void addFunc(int source_type, FT* afunc, int target_type = -1)
  {
    if (Fs.empty())
    {
      Fs.resize(CenterRef.getTotalNum(), 0);
      Funique.resize(CenterRef.getSpeciesSet().size(), 0);
    }
    for (int i = 0; i < Fs.size(); i++)
      if (CenterRef.GroupID[i] == source_type)
        Fs[i] = afunc;
    Funique[source_type] = afunc;
  }


  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active) override
  {
    for (int i = 0; i < Funique.size(); ++i)
      if (Funique[i])
        Funique[i]->resetParameters(active);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.clear();
    for (int i = 0; i < Funique.size(); ++i)
    {
      if (Funique[i])
      {
        Funique[i]->myVars.getIndex(active);
        myVars.insertFrom(Funique[i]->myVars);
      }
    }
    myVars.getIndex(active);
    NumVars = myVars.size();
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
      OffSet.resize(Fs.size());
      int varoffset = myVars.Index[0];
      for (int i = 0; i < Fs.size(); ++i)
      {
        if (Fs[i])
        {
          OffSet[i].first  = Fs[i]->myVars.Index.front() - varoffset;
          OffSet[i].second = Fs[i]->myVars.Index.size() + OffSet[i].first;
        }
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override
  {
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
      const auto& d_table = P.getDistTable(myTableIndex);
      dLogPsi             = 0.0;
      for (int p = 0; p < NumVars; ++p)
        (*gradLogPsi[p]) = 0.0;
      for (int p = 0; p < NumVars; ++p)
        (*lapLogPsi[p]) = 0.0;
      std::vector<TinyVector<RealType, 3>> derivs(NumVars);

      constexpr RealType cone(1);
      constexpr RealType lapfac(OHMMS_DIM - cone);
      const size_t ns = d_table.sources();
      const size_t nt = P.getTotalNum();

      aligned_vector<int> iadj(nt);
      aligned_vector<RealType> dist(nt);
      std::vector<PosType> displ(nt);

      for (size_t i = 0; i < ns; ++i)
      {
        FT* func = Fs[i];
        if (func == 0)
          continue;
        int first(OffSet[i].first);
        int last(OffSet[i].second);
        bool recalcFunc(false);
        for (int rcs = first; rcs < last; rcs++)
          if (rcsingles[rcs] == true)
            recalcFunc = true;
        if (recalcFunc)
        {
          size_t nn = d_table.get_neighbors(i, func->cutoff_radius, iadj.data(), dist.data(), displ.data());
          for (size_t nj = 0; nj < nn; ++nj)
          {
            std::fill(derivs.begin(), derivs.end(), 0);
            if (!func->evaluateDerivatives(dist[nj], derivs))
              continue;
            int j = iadj[nj];
            RealType rinv(cone / dist[nj]);
            PosType& dr = displ[nj];
            for (int p = first, ip = 0; p < last; ++p, ++ip)
            {
              dLogPsi[p] -= derivs[ip][0];
              RealType dudr(rinv * derivs[ip][1]);
              (*gradLogPsi[p])[j] -= dudr * dr;
              (*lapLogPsi[p])[j] -= derivs[ip][2] + lapfac * dudr;
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
          dlogpsi[kk] = ValueType(dLogPsi[k]);
        }
      }
    }
  }

  inline void setVars(const opt_variables_type& vars)
  {
    NumVars = vars.size();
    if (NumVars == 0)
      return;
    myVars = vars;
    dLogPsi.resize(NumVars);
    gradLogPsi.resize(NumVars, 0);
    lapLogPsi.resize(NumVars, 0);
    for (int i = 0; i < NumVars; ++i)
    {
      gradLogPsi[i] = new GradVectorType(NumPtcls);
      lapLogPsi[i]  = new ValueVectorType(NumPtcls);
    }
  }

  std::unique_ptr<DiffWaveFunctionComponent> makeClone(ParticleSet& tqp) const override
  {
    auto j1copy = std::make_unique<DiffOneBodyJastrowOrbital<FT>>(CenterRef, tqp);
    for (int i = 0; i < Funique.size(); ++i)
    {
      if (Funique[i])
        j1copy->addFunc(i, new FT(*Funique[i]), -1);
    }
    j1copy->setVars(myVars);
    j1copy->OffSet = OffSet;
    return j1copy;
  }
};


} // namespace qmcplusplus
#endif

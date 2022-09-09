#ifndef QMCPLUSPLUS_AFQMC_MIXEDRDMESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_MIXEDRDMESTIMATOR_H

#include "Message/MPIObjectBase.h"
#include "AFQMC/config.h"
#include <vector>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/libxmldefs.h"

#include "AFQMC/Wavefunctions/Wavefunction.hpp"
#include "AFQMC/Walkers/WalkerSet.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class MixedRDMEstimator : public EstimatorBase
{
  using CMatrix_ref = boost::multi::array_ref<ComplexType, 2>;
  using CVector     = boost::multi::array<ComplexType, 1>;
  using stdCVector  = boost::multi::array<ComplexType, 1>;

public:
  MixedRDMEstimator(afqmc::TaskGroup_& tg_,
                    AFQMCInfo& info,
                    std::string title,
                    xmlNodePtr cur,
                    WALKER_TYPES wlk,
                    Wavefunction& wfn,
                    bool impsamp_ = true)
      : EstimatorBase(info), TG(tg_), writer(false), wfn0(wfn), block_size(1), importanceSampling(impsamp_)
  {
    if (cur != NULL)
    {
      ParameterSet m_param;
      std::string restore_paths;
      m_param.add(block_size, "block_size");
      m_param.put(cur);
    }

    ncores_per_TG = TG.getNCoresPerTG();
    core_rank     = TG.getLocalTGRank();
    writer        = (TG.getGlobalRank() == 0);
    if (wlk == CLOSED)
    {
      dm_size = NMO * NMO;
      dm_dims = {NMO, NMO};
    }
    else if (wlk == COLLINEAR)
    {
      dm_size = 2 * NMO * NMO;
      dm_dims = {2 * NMO, NMO};
    }
    else if (wlk == NONCOLLINEAR)
    {
      dm_size = 4 * NMO * NMO;
      dm_dims = {2 * NMO, 2 * NMO};
    }
    if (DMBuffer.size() < dm_size)
    {
      DMBuffer.reextent(iextensions<1u>{dm_size});
    }
    if (DMAverage.size() < dm_size)
    {
      DMAverage.reextent(iextensions<1u>{dm_size});
    }
    std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0, 0.0));
    std::fill(DMAverage.begin(), DMAverage.end(), ComplexType(0.0, 0.0));
    denom.reextent({1});
    denom_average.reextent({1});
    denom_average[0] = ComplexType(0.0);
  }

  ~MixedRDMEstimator() {}

  void accumulate_step(WalkerSet& wset, std::vector<ComplexType>& curData) override {}

  void accumulate_block(WalkerSet& wset) override
  {
    // check to see whether we should be accumulating estimates.
    CMatrix_ref OneRDM(DMBuffer.data(), {dm_dims.first, dm_dims.second});
    denom[0] = ComplexType(0.0, 0.0);
    std::fill(DMBuffer.begin(), DMBuffer.end(), ComplexType(0.0, 0.0));
    stdCVector wgt(iextensions<1u>{wset.size()});
    wset.getProperty(WEIGHT, wgt);

    int nx((wset.getWalkerType() == COLLINEAR) ? 2 : 1);
    if (std::get<0>(wDMsum.sizes()) != wset.size() || std::get<1>(wDMsum.sizes()) != nx)
      wDMsum.reextent({wset.size(), nx});
    if (std::get<0>(wOvlp.sizes()) != wset.size() || std::get<1>(wOvlp.sizes()) != nx)
      wOvlp.reextent({wset.size(), nx});

    if (!importanceSampling)
    {
      stdCVector phase(iextensions<1u>{wset.size()});
      wset.getProperty(PHASE, phase);
      for (int i = 0; i < wgt.size(); i++)
        wgt[i] *= phase[i];
    }
    wfn0.WalkerAveragedDensityMatrix(wset, wgt, OneRDM, denom, wOvlp, wDMsum, !importanceSampling);
    iblock++;
  }

  void tags(std::ofstream& out) override {}

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset) override
  {
    // I doubt we will ever collect a billion blocks of data.
    int n_zero = 9;
    if (writer)
    {
      for (int i = 0; i < DMBuffer.size(); i++)
        DMAverage[i] += DMBuffer[i];
      denom_average[0] += denom[0];
      if (iblock % block_size == 0)
      {
        for (int i = 0; i < DMAverage.size(); i++)
          DMAverage[i] /= block_size;
        denom_average[0] /= block_size;
        dump.push("Mixed");
        std::string padded_iblock = std::string(n_zero - std::to_string(iblock).length(), '0') + std::to_string(iblock);
        boost::multi::array_ref<ComplexType, 1> wOvlp_(wOvlp.origin(), {std::get<0>(wOvlp.sizes()) * std::get<1>(wOvlp.sizes())});
        boost::multi::array_ref<ComplexType, 1> wDMsum_(wDMsum.origin(), {std::get<0>(wDMsum.sizes()) * std::get<1>(wDMsum.sizes())});
        dump.write(DMAverage, "one_rdm_" + padded_iblock);
        dump.write(denom_average, "one_rdm_denom_" + padded_iblock);
        dump.write(wOvlp_, "one_rdm_walker_overlaps_" + padded_iblock);
        dump.write(wDMsum_, "one_rdm_walker_dm_sums_" + padded_iblock);
        dump.pop();
        std::fill(DMAverage.begin(), DMAverage.end(), ComplexType(0.0, 0.0));
        std::fill(denom_average.begin(), denom_average.end(), ComplexType(0.0, 0.0));
      }
    }
  }

private:
  TaskGroup_& TG;

  bool writer;

  Wavefunction& wfn0;

  // The first element of data stores the denominator of the estimator (i.e., the total
  // walker weight including rescaling factors etc.). The rest of the elements store the
  // averages of the various elements of the green's function.
  CVector DMBuffer, DMAverage;

  int core_rank;
  int ncores_per_TG;
  int iblock       = 0;
  ComplexType zero = ComplexType(0.0, 0.0);
  ComplexType one  = ComplexType(1.0, 0.0);

  // Block size over which RDM will be averaged.
  int block_size;
  bool importanceSampling;
  std::vector<ComplexType> weights;
  int dm_size;
  std::pair<int, int> dm_dims;
  CVector denom, denom_average;

  boost::multi::array<ComplexType, 2> wDMsum;
  boost::multi::array<ComplexType, 2> wOvlp;
};
} // namespace afqmc
} // namespace qmcplusplus

#endif

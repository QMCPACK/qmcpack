#ifndef QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H
#define QMCPLUSPLUS_AFQMC_ENERGYESTIMATOR_H

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
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{
namespace afqmc
{
class EnergyEstimator : public EstimatorBase
{
public:
  EnergyEstimator(afqmc::TaskGroup_& tg_,
                  AFQMCInfo info,
                  xmlNodePtr cur,
                  Wavefunction& wfn,
                  bool impsamp_ = true,
                  bool timer    = true)
      : EstimatorBase(info), TG(tg_), wfn0(wfn), importanceSampling(impsamp_), energy_components(false)
  {
    if (cur != NULL)
    {
      ParameterSet m_param;
      std::string print_components;
      m_param.add(print_components, "print_components");
      m_param.put(cur);
      if (print_components == "true" || print_components == "yes")
      {
        energy_components = true;
      }
      else
      {
        energy_components = false;
      }
    }
    data.resize(5);
  }

  ~EnergyEstimator() {}

  void accumulate_step(WalkerSet& wlks, std::vector<ComplexType>& curData) override {}

  void accumulate_block(WalkerSet& wset) override
  {
    ScopedTimer local_timer(AFQMCTimers[energy_timer]);
    size_t nwalk = wset.size();
    if (eloc.size(0) != nwalk || eloc.size(1) != 3)
      eloc.reextent({static_cast<boost::multi::size_t>(nwalk), 3});
    if (ovlp.size(0) != nwalk)
      ovlp.reextent(iextensions<1u>(nwalk));
    if (wprop.size(0) != 4 || wprop.size(1) != nwalk)
      wprop.reextent({4, static_cast<boost::multi::size_t>(nwalk)});

    ComplexType dum, et;
    wfn0.Energy(wset, eloc, ovlp);
    // in case GPU
    ComplexMatrix<std::allocator<ComplexType>> eloc_(eloc);
    ComplexVector<std::allocator<ComplexType>> ovlp_(ovlp);
    if (TG.TG_local().root())
    {
      wset.getProperty(WEIGHT, wprop[0]);
      wset.getProperty(OVLP, wprop[1]);
      wset.getProperty(PHASE, wprop[2]);
      std::fill_n(data.begin(), data.size(), ComplexType(0.0));
      for (int i = 0; i < nwalk; i++)
      {
        if (std::isnan(real(wprop[0][i])))
          continue;
        if (importanceSampling)
        {
          dum = (wprop[0][i]) * ovlp_[i] / (wprop[1][i]);
        }
        else
        {
          dum = (wprop[0][i]) * ovlp_[i] * (wprop[2][i]);
        }
        et = eloc_[i][0] + eloc_[i][1] + eloc_[i][2];
        if ((!std::isfinite(real(dum))) || (!std::isfinite(real(et * dum))))
          continue;
        data[1] += dum;
        data[0] += et * dum;
        data[2] += eloc_[i][0] * dum;
        data[3] += eloc_[i][1] * dum;
        data[4] += eloc_[i][2] * dum;
      }
      TG.TG_heads().all_reduce_in_place_n(data.begin(), data.size(), std::plus<>());
    }
  }

  void tags(std::ofstream& out) override
  {
    if (TG.Global().root())
    {
      out << "EnergyEstim_" << name << "_nume_real  EnergyEstim_" << name << "_nume_imag "
          << "EnergyEstim_" << name << "_deno_real  EnergyEstim_" << name << "_deno_imag "
          << "EnergyEstim_" << name << "_timer ";
      if (energy_components)
      {
        out << "OneBodyEnergyEstim__nume_real "
            << "ECoulEnergyEstim__nume_real "
            << "EXXEnergyEstim__nume_real ";
      }
    }
  }

  void print(std::ofstream& out, hdf_archive& dump, WalkerSet& wset) override
  {
    if (TG.Global().root())
    {
      int n = wset.get_global_target_population();
      out << data[0].real() / n << " " << data[0].imag() / n << " " << data[1].real() / n << " " << data[1].imag() / n
          << " " << AFQMCTimers[energy_timer].get().get_total() << " ";
      if (energy_components)
      {
        out << data[2].real() / n << " " << data[3].real() / n << " " << data[4].real() / n << " ";
      }
      AFQMCTimers[energy_timer].get().reset();
    }
  }

private:
  std::string name;

  TaskGroup_& TG;

  Wavefunction& wfn0;

  ComplexMatrix<device_allocator<ComplexType>> eloc;
  ComplexVector<device_allocator<ComplexType>> ovlp;
  ComplexMatrix<std::allocator<ComplexType>> wprop;

  std::vector<std::complex<double>> data;

  bool importanceSampling;
  bool energy_components;
};
} // namespace afqmc
} // namespace qmcplusplus

#endif

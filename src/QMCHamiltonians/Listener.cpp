//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Listener.hpp"
#include <unordered_map>
#include <OhmmsPETE/OhmmsVector.h>
#include <string>

namespace qmcplusplus
{

template void combinePerParticleEnergies<double>(const CrowdEnergyValues<double>& cev_in,
                                                 std::vector<Vector<double>>& values_out);

template void combinePerParticleEnergies<float>(const CrowdEnergyValues<float>& cev_in,
                                                std::vector<Vector<float>>& values_out);

template void combinePerParticleEnergies<std::complex<double>>(const CrowdEnergyValues<std::complex<double>>& cev_in,
                                                               std::vector<Vector<std::complex<double>>>& values_out);

template void combinePerParticleEnergies<std::complex<float>>(const CrowdEnergyValues<std::complex<float>>& cev_in,
                                                              std::vector<Vector<std::complex<float>>>& values_out);

template<typename T>
void combinePerParticleEnergies(const CrowdEnergyValues<T>& cev_in, std::vector<Vector<T>>& values_out)
{
  const auto num_walkers = cev_in.begin()->second.size();
  const auto num_particles = cev_in.begin()->second[0].size();

  values_out.resize(num_walkers);
  for (int iw = 0; iw < num_walkers; ++iw)
    values_out[iw].resize(num_particles);

  assert(cev_in.begin()->second.size() == values_out.size());
  assert(cev_in.begin()->second[0].size() == values_out[0].size());

  for (auto& [component, values] : cev_in)
  {
    for (int iw = 0; iw < values.size(); ++iw)
    {
      // using Vector operator+=
      for (int ip = 0; ip < values[iw].size(); ++ip)
        values_out[iw][ip] += values[iw][ip];
    }
  }
}

} // namespace qmcplusplus

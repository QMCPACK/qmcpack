//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: DensityMatrices1b.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PER_PARTICLE_HAMILTONIAN_LOGGER_INPUT_H
#define QMCPLUSPLUS_PER_PARTICLE_HAMILTONIAN_LOGGER_INPUT_H

#include "Configuration.h"
#include "InputSection.h"
#include "OperatorEstBase.h"

namespace qmcplusplus
{
class PerParticleHamiltonianLogger;
class PerParticleHamiltonianLoggerInput
{
public:
  using Consumer = PerParticleHamiltonianLogger;
  using Real     = QMCTraits::RealType;

  class PerParticleHamiltonianLoggerInputSection : public InputSection
  {
  public:
    PerParticleHamiltonianLoggerInputSection()
    {
      section_name = "PerParticleHamiltonianLogger";
      attributes   = {"to_stdout", "validate_per_particle_sum", "type", "name"};
      bools        = {"to_stdout", "validate_per_particle_sum"};
      strings      = {"type", "name"};
    }
    PerParticleHamiltonianLoggerInputSection(const PerParticleHamiltonianLoggerInputSection& other) = default;
  };
  PerParticleHamiltonianLoggerInput(const PerParticleHamiltonianLoggerInput& other) = default;
  PerParticleHamiltonianLoggerInput(xmlNodePtr cur);

  /** For this input class its valid with just its defaults
   */
  PerParticleHamiltonianLoggerInput() = default;
  
  const std::string& get_name() const { return name_; }
  bool get_to_stdout() const { return to_stdout_; }

private:
  PerParticleHamiltonianLoggerInputSection input_section_;
  std::string name_               = "per_particle_log";
  bool to_stdout_                 = false;
};
} // namespace qmcplusplus
#endif

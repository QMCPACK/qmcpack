//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/DeepQMC/DeepQMCWaveFunctionBuilder.h"

#include <stdexcept>

#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCBridge.h"
#include "QMCWaveFunctions/DeepQMC/DeepQMCWF.h"

namespace qmcplusplus
{

DeepQMCWaveFunctionBuilder::DeepQMCWaveFunctionBuilder(Communicate* comm, ParticleSet& target, const PSetMap& psets)
    : WaveFunctionComponentBuilder(comm, target), ptcl_pool_(psets)
{}

std::unique_ptr<WaveFunctionComponent> DeepQMCWaveFunctionBuilder::buildComponent(xmlNodePtr cur)
{
  std::string name("deepqmc");
  std::string source_name("ion0");
  std::string model_path;
  std::string python_module_path;
  int mol_idx = 0;

  OhmmsAttributeSet attribs;
  attribs.add(name, "name");
  attribs.add(source_name, "source");
  attribs.add(model_path, "model");
  attribs.add(python_module_path, "python_module_path");
  attribs.add(mol_idx, "mol_idx");
  attribs.put(cur);

  auto ion_it = ptcl_pool_.find(source_name);
  if (ion_it == ptcl_pool_.end())
    throw std::runtime_error("DeepQMC wavefunction source particle set not found: " + source_name);

  auto bridge = makePythonDeepQMCBridge(model_path, python_module_path);
  return std::make_unique<DeepQMCWF>(name, *ion_it->second, std::move(bridge), mol_idx);
}

} // namespace qmcplusplus

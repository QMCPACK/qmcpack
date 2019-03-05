//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

/**@file WaveFunctionComponentBuilder.cpp
  *@brief Initialization of static data members for wavefunction-related tags.
  *
  *The main input files for qmcplusplus applications should use the matching
  *tags defined here.
  */
namespace qmcplusplus
{
//int WaveFunctionComponentBuilder::print_level=1;

std::string WaveFunctionComponentBuilder::wfs_tag = "wavefunction";

std::string WaveFunctionComponentBuilder::param_tag = "parameter";

std::string WaveFunctionComponentBuilder::dtable_tag = "distancetable";

std::string WaveFunctionComponentBuilder::jastrow_tag = "jastrow";

std::string WaveFunctionComponentBuilder::fdlrwfn_tag = "fdlrwfn";

std::string WaveFunctionComponentBuilder::detset_tag = "determinantset";

std::string WaveFunctionComponentBuilder::sd_tag = "slaterdeterminant";

std::string WaveFunctionComponentBuilder::det_tag = "determinant";

std::string WaveFunctionComponentBuilder::rn_tag = "determinant_rn";

std::string WaveFunctionComponentBuilder::spo_tag = "psi";

std::string WaveFunctionComponentBuilder::sposet_tag = "sposet";

std::string WaveFunctionComponentBuilder::basisset_tag = "basisset";

std::string WaveFunctionComponentBuilder::basis_tag = "basis";

std::string WaveFunctionComponentBuilder::basisfunc_tag = "phi";

std::string WaveFunctionComponentBuilder::ionorb_tag = "ionwf";

std::string WaveFunctionComponentBuilder::backflow_tag = "backflow";

std::string WaveFunctionComponentBuilder::multisd_tag = "multideterminant";

std::string WaveFunctionComponentBuilder::sposcanner_tag = "spo_scanner";

WaveFunctionComponentBuilder::WaveFunctionComponentBuilder(ParticleSet& p, TrialWaveFunction& psi)
    : MPIObjectBase(psi.getCommunicator()), targetPtcl(p), targetPsi(psi), myNode(NULL)
{}

WaveFunctionComponentBuilder::~WaveFunctionComponentBuilder() {}
} // namespace qmcplusplus

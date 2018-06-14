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
    
    
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

/**@file OrbitalBuilderBase.cpp
  *@brief Initialization of static data members for wavefunction-related tags.
  *
  *The main input files for qmcplusplus applications should use the matching
  *tags defined here.
  */
namespace qmcplusplus
{

//int OrbitalBuilderBase::print_level=1;

std::string OrbitalBuilderBase::wfs_tag="wavefunction";

std::string OrbitalBuilderBase::param_tag="parameter";

std::string OrbitalBuilderBase::dtable_tag="distancetable";

std::string OrbitalBuilderBase::jastrow_tag="jastrow";

std::string OrbitalBuilderBase::fdlrwfn_tag="fdlrwfn";

std::string OrbitalBuilderBase::detset_tag="determinantset";

std::string OrbitalBuilderBase::sd_tag="slaterdeterminant";

std::string OrbitalBuilderBase::det_tag="determinant";

std::string OrbitalBuilderBase::rn_tag="determinant_rn";

std::string OrbitalBuilderBase::spo_tag="psi";

std::string OrbitalBuilderBase::sposet_tag="sposet";

std::string OrbitalBuilderBase::basisset_tag="basisset";

std::string OrbitalBuilderBase::basis_tag="basis";

std::string OrbitalBuilderBase::basisfunc_tag="phi";

std::string OrbitalBuilderBase::ionorb_tag="ionwf";

std::string OrbitalBuilderBase::backflow_tag="backflow";

std::string OrbitalBuilderBase::multisd_tag="multideterminant";

std::string OrbitalBuilderBase::sposcanner_tag="spo_scanner";

OrbitalBuilderBase::OrbitalBuilderBase(ParticleSet& p, TrialWaveFunction& psi):
  MPIObjectBase(psi.getCommunicator()),
  targetPtcl(p), targetPsi(psi), myNode(NULL)
{
}

OrbitalBuilderBase::~OrbitalBuilderBase()
{
}
}

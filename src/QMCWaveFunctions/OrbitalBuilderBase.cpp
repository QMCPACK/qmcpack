//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

string OrbitalBuilderBase::wfs_tag="wavefunction";

string OrbitalBuilderBase::param_tag="parameter";

string OrbitalBuilderBase::dtable_tag="distancetable";

string OrbitalBuilderBase::jastrow_tag="jastrow";

string OrbitalBuilderBase::detset_tag="determinantset";

string OrbitalBuilderBase::sd_tag="slaterdeterminant";

string OrbitalBuilderBase::det_tag="determinant";

string OrbitalBuilderBase::rn_tag="determinant_rn";

string OrbitalBuilderBase::spo_tag="psi";

string OrbitalBuilderBase::sposet_tag="sposet";

string OrbitalBuilderBase::basisset_tag="basisset";

string OrbitalBuilderBase::basis_tag="basis";

string OrbitalBuilderBase::basisfunc_tag="phi";

string OrbitalBuilderBase::ionorb_tag="ionwf";

string OrbitalBuilderBase::backflow_tag="backflow";

string OrbitalBuilderBase::multisd_tag="multideterminant";

OrbitalBuilderBase::OrbitalBuilderBase(ParticleSet& p, TrialWaveFunction& psi):
  MPIObjectBase(psi.getCommunicator()),
  targetPtcl(p), targetPsi(psi), myNode(NULL)
{
}

OrbitalBuilderBase::~OrbitalBuilderBase()
{
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

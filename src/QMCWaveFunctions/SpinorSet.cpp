//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:  Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/SpinorSet.h"

namespace qmcplusplus
{
SpinorSet::SpinorSet():SPOSet(),className("SpinorSet")
{

}

SpinorSet::~SpinorSet()
{

}

void SpinorSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{

}
}

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

void SpinorSet::resetParameters(const opt_variables_type& optVariables){};

void SpinorSet::resetTargetParticleSet(ParticleSet& P){};

void SpinorSet::setOrbitalSetSize(int norbs){};


void SpinorSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{

}

void SpinorSet::evaluate(const ParticleSet& P,
                      int iat,
                      ValueVector_t& psi,
                      GradVector_t& dpsi,
                      ValueVector_t& d2psi)
{

}

void SpinorSet::evaluate_notranspose(const ParticleSet& P,
                                  int first,
                                  int last,
                                  ValueMatrix_t& logdet,
                                  GradMatrix_t& dlogdet,
                                  ValueMatrix_t& d2logdet)
{
}

}

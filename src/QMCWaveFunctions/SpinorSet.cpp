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
SpinorSet::SpinorSet():SPOSet(),className("SpinorSet"),spo_up(nullptr),spo_dn(nullptr)
{

}

SpinorSet::~SpinorSet()
{

}

void SpinorSet::set_spos(std::shared_ptr<SPOSet> up, std::shared_ptr<SPOSet> dn)
{
  //Sanity check for input SPO's.  They need to be the same size or 
  IndexType spo_size_up=up->getOrbitalSetSize();
  IndexType spo_size_down=dn->getOrbitalSetSize();

  if(spo_size_up != spo_size_down)
    APP_ABORT("SpinorSet::set_spos(...):  up and down SPO components have different sizes.");

  setOrbitalSetSize(spo_size_up);

  spo_up=up;
  spo_dn=dn;

  psi_work_up.resize(OrbitalSetSize);
  psi_work_down.resize(OrbitalSetSize);
  
  dpsi_work_up.resize(OrbitalSetSize);
  dpsi_work_down.resize(OrbitalSetSize);

  d2psi_work_up.resize(OrbitalSetSize);
  d2psi_work_down.resize(OrbitalSetSize);
}

void SpinorSet::resetParameters(const opt_variables_type& optVariables){};

void SpinorSet::resetTargetParticleSet(ParticleSet& P){};

void SpinorSet::setOrbitalSetSize(int norbs){ OrbitalSetSize=norbs; };


void SpinorSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  psi_work_up=0.0;
  psi_work_down=0.0;
 
  spo_up->evaluate(P,iat,psi_work_up);
  spo_dn->evaluate(P,iat,psi_work_down);
  
  psi = psi_work_up+psi_work_down;
  
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
  IndexType nelec=P.getTotalNum();

  logpsi_work_up.resize(nelec, OrbitalSetSize);
  logpsi_work_down.resize(nelec, OrbitalSetSize);

  dlogpsi_work_up.resize(nelec, OrbitalSetSize);
  dlogpsi_work_down.resize(nelec, OrbitalSetSize);

  d2logpsi_work_up.resize(nelec, OrbitalSetSize);
  d2logpsi_work_down.resize(nelec, OrbitalSetSize);

  spo_up->evaluate_notranspose(P, first, last, logpsi_work_up, dlogpsi_work_up, d2logpsi_work_up);
  spo_dn->evaluate_notranspose(P, first, last, logpsi_work_down, dlogpsi_work_down, d2logpsi_work_down);

  logdet=logpsi_work_up+logpsi_work_down;
  dlogdet=dlogpsi_work_up+dlogpsi_work_down;
  d2logdet=d2logpsi_work_up+d2logpsi_work_down;
   
}

}

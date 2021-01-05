//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorEstBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

/**@file
 */
#include "Message/Communicate.h"
#include "OperatorEstBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
OperatorEstBase::OperatorEstBase() : Value(0.0), walkers_weight_(0)
{
  quantum_domain = no_quantum_domain;
  energy_domain  = no_energy_domain;
  data_locality_ = DataLocality::crowd;
  UpdateMode.set(PRIMARY, 1);
}

OperatorEstBase::OperatorEstBase(const OperatorEstBase& oth) : walkers_weight_(0)
{
  quantum_domain = oth.quantum_domain;
  energy_domain = oth.energy_domain;
  data_locality_ = oth.data_locality_;     
  UpdateMode.set(PRIMARY, 1);
}

// I suspect this can be a pure function outside of the class.
// In this case at least we don't care to copy the data_ as we are going to reduce these later and don't want
// to end up with a multiplicative factor if we already have data.
OperatorEstBase::Data OperatorEstBase::createLocalData(size_t size, DataLocality data_locality)
{
  Data new_data;
  if (data_locality == DataLocality::crowd)
  {
    new_data = std::make_unique<std::vector<QMCT::RealType>>(size, 0);
  }
  else
  {
    throw std::runtime_error("currently SpinDensityNew only supports crowd level datalocality");
  }
  return new_data;
}

void OperatorEstBase::set_energy_domain(energy_domains edomain)
{
  if (energy_domain_valid(edomain))
    energy_domain = edomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_energy_domain\n  input energy domain is invalid");
}

void OperatorEstBase::set_quantum_domain(quantum_domains qdomain)
{
  if (quantum_domain_valid(qdomain))
    quantum_domain = qdomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_quantum_domain\n  input quantum domain is invalid");
}

void OperatorEstBase::one_body_quantum_domain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain = classical;
  else if (P.is_quantum())
    quantum_domain = quantum;
  else
    APP_ABORT("OperatorEstBase::one_body_quantum_domain\n  quantum domain of input particles is invalid");
}

void OperatorEstBase::two_body_quantum_domain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain = classical_classical;
  else if (P.is_quantum())
    quantum_domain = quantum_quantum;
  else
    APP_ABORT("OperatorEstBase::two_body_quantum_domain(P)\n  quantum domain of input particles is invalid");
}

void OperatorEstBase::two_body_quantum_domain(const ParticleSet& P1, const ParticleSet& P2)
{
  bool c1 = P1.is_classical();
  bool c2 = P2.is_classical();
  bool q1 = P1.is_quantum();
  bool q2 = P2.is_quantum();
  if (c1 && c2)
    quantum_domain = classical_classical;
  else if ((q1 && c2) || (c1 && q2))
    quantum_domain = quantum_classical;
  else if (q1 && q2)
    quantum_domain = quantum_quantum;
  else
    APP_ABORT("OperatorEstBase::two_body_quantum_domain(P1,P2)\n  quantum domain of input particles is invalid");
}

bool OperatorEstBase::quantum_domain_valid(quantum_domains qdomain) { return qdomain != no_quantum_domain; }

} // namespace qmcplusplus

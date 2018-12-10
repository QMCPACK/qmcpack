//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file QMCHamiltonianBase.cpp
 *@brief Definition of QMCHamiltonianBase
 */
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>

namespace qmcplusplus
{

QMCHamiltonianBase::QMCHamiltonianBase()
  :myIndex(-1),Value(0.0),Dependants(0),tWalker(0)
{
  quantum_domain       = no_quantum_domain;
  energy_domain        = no_energy_domain;

#if !defined(REMOVE_TRACEMANAGER)
  streaming_scalars    = false;
  streaming_particles  = false;
  have_required_traces = false;
#endif
  UpdateMode.set(PRIMARY,1);
}

void QMCHamiltonianBase::set_energy_domain(energy_domains edomain)
{
  if(energy_domain_valid(edomain))
    energy_domain = edomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_energy_domain\n  input energy domain is invalid");
}

void QMCHamiltonianBase::set_quantum_domain(quantum_domains qdomain)
{
  if(quantum_domain_valid(qdomain))
    quantum_domain = qdomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_quantum_domain\n  input quantum domain is invalid");
}

void QMCHamiltonianBase::one_body_quantum_domain(const ParticleSet& P)
{
  if(P.is_classical())
    quantum_domain = classical;
  else if(P.is_quantum())
    quantum_domain = quantum;
  else
    APP_ABORT("QMCHamiltonianBase::one_body_quantum_domain\n  quantum domain of input particles is invalid");
}

void QMCHamiltonianBase::two_body_quantum_domain(const ParticleSet& P)
{
  if(P.is_classical())
    quantum_domain = classical_classical;
  else if(P.is_quantum())
    quantum_domain = quantum_quantum;
  else
    APP_ABORT("QMCHamiltonianBase::two_body_quantum_domain(P)\n  quantum domain of input particles is invalid");
}

void QMCHamiltonianBase::two_body_quantum_domain(const ParticleSet& P1,const ParticleSet& P2)
{
  bool c1 = P1.is_classical();
  bool c2 = P2.is_classical();
  bool q1 = P1.is_quantum();
  bool q2 = P2.is_quantum();
  if(c1 && c2)
    quantum_domain = classical_classical;
  else if((q1 && c2) || (c1 && q2))
    quantum_domain = quantum_classical;
  else if(q1 && q2)
    quantum_domain = quantum_quantum;
  else
    APP_ABORT("QMCHamiltonianBase::two_body_quantum_domain(P1,P2)\n  quantum domain of input particles is invalid");
}

bool QMCHamiltonianBase::quantum_domain_valid(quantum_domains qdomain)
{
  return qdomain!=no_quantum_domain;
}

void QMCHamiltonianBase::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi
    ,QMCHamiltonian& targetH)
{
  QMCHamiltonianBase* myclone=makeClone(qp,psi);
  if(myclone)
    targetH.addOperator(myclone,myName,UpdateMode[PHYSICAL]);
}

void QMCHamiltonianBase::registerObservables(std::vector<observable_helper*>& h5desc
    , hid_t gid)  const
{
  bool collect=UpdateMode.test(COLLECTABLE);
  //exclude collectables
  if(!collect)
  {
    int loc=h5desc.size();
    h5desc.push_back(new observable_helper(myName));
    std::vector<int> onedim(1,1);
    h5desc[loc]->set_dimensions(onedim,myIndex);
    h5desc[loc]->open(gid);
  }
}

void QMCHamiltonianBase::addEnergy(MCWalkerConfiguration &W,
    std::vector<RealType> &LocalEnergy)
{
  app_error() << "Need specialization for " << myName
    << "::addEnergy(MCWalkerConfiguration &W).\n";
}

}


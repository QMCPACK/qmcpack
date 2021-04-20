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


/**@file OperatorBase.cpp
 *@brief Definition of OperatorBase
 */
#include "Message/Communicate.h"
#include "OperatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
OperatorBase::OperatorBase() : myIndex(-1), Dependants(0), Value(0.0), tWalker(0)
{
  quantum_domain = no_quantum_domain;
  energy_domain  = no_energy_domain;

#if !defined(REMOVE_TRACEMANAGER)
  streaming_scalars    = false;
  streaming_particles  = false;
  have_required_traces = false;
#endif
  UpdateMode.set(PRIMARY, 1);
}

/** The correct behavior of this routine requires estimators with non-deterministic components
 * in their evaluate() function to override this function.
 */
OperatorBase::Return_t OperatorBase::evaluateDeterministic(ParticleSet& P) { return evaluate(P); }
/** Take o_list and p_list update evaluation result variables in o_list?
 *
 * really should reduce vector of local_energies. matching the ordering and size of o list
 * the this can be call for 1 or more QMCHamiltonians
 */
void OperatorBase::mw_evaluate(const RefVectorWithLeader<OperatorBase>& O_list,
                               const RefVectorWithLeader<ParticleSet>& P_list) const
{
  assert(this == &O_list.getLeader());
/**  Temporary raw omp pragma for simple thread parallelism
   *   ignoring the driver level concurrency
   *   
   *  \todo replace this with a proper abstraction. It should adequately describe the behavior
   *  and strictly limit the activation of this level concurrency to when it is intended.
   *  It is unlikely to belong in this function.
   *  
   *  This implicitly depends on openmp work division logic. Essentially adhoc runtime
   *  crowds over which we have given up control of thread/global scope.
   *  How many walkers per thread? How to handle their data movement if any of these
   *  hamiltonians should be accelerated? We can neither reason about or describe it in C++
   *
   *  As I understand it it should only be required for as long as the AMD openmp offload 
   *  compliler is incapable of running multiple threads. They should/must fix their compiler
   *  before delivery of frontier and it should be removed at that point at latest
   *
   *  If you want 16 threads of 1 walker that should be 16 crowds of 1
   *  not one crowd of 16 with openmp thrown in at hamiltonian level.
   *  If this must be different from the other crowd batching. Make this a reasoned about
   *  and controlled level of concurency blocking at the driver level.
   *
   *  This is only thread safe only if each walker has a complete
   *  set of anything involved in an Operator.evaluate.
   */
#pragma omp parallel for
  for (int iw = 0; iw < O_list.size(); iw++)
    O_list[iw].evaluate(P_list[iw]);
}

void OperatorBase::mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& O_list,
                                                       const RefVectorWithLeader<ParticleSet>& P_list,
                                                       const opt_variables_type& optvars,
                                                       RecordArray<ValueType>& dlogpsi,
                                                       RecordArray<ValueType>& dhpsioverpsi) const
{
  const int nparam = dlogpsi.nparam();
  std::vector<ValueType> tmp_dlogpsi(nparam);
  std::vector<ValueType> tmp_dhpsioverpsi(nparam);
  for (int iw = 0; iw < O_list.size(); iw++)
  {
    for (int j = 0; j < nparam; j++)
    {
      tmp_dlogpsi[j] = dlogpsi.getValue(j, iw);
    }

    O_list[iw].evaluateValueAndDerivatives(P_list[iw], optvars, tmp_dlogpsi, tmp_dhpsioverpsi);

    for (int j = 0; j < nparam; j++)
    {
      dhpsioverpsi.setValue(j, iw, dhpsioverpsi.getValue(j, iw) + tmp_dhpsioverpsi[j]);
    }
  }
}


void OperatorBase::set_energy_domain(energy_domains edomain)
{
  if (energy_domain_valid(edomain))
    energy_domain = edomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_energy_domain\n  input energy domain is invalid");
}

void OperatorBase::set_quantum_domain(quantum_domains qdomain)
{
  if (quantum_domain_valid(qdomain))
    quantum_domain = qdomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_quantum_domain\n  input quantum domain is invalid");
}

void OperatorBase::one_body_quantum_domain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain = classical;
  else if (P.is_quantum())
    quantum_domain = quantum;
  else
    APP_ABORT("OperatorBase::one_body_quantum_domain\n  quantum domain of input particles is invalid");
}

void OperatorBase::two_body_quantum_domain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain = classical_classical;
  else if (P.is_quantum())
    quantum_domain = quantum_quantum;
  else
    APP_ABORT("OperatorBase::two_body_quantum_domain(P)\n  quantum domain of input particles is invalid");
}

void OperatorBase::two_body_quantum_domain(const ParticleSet& P1, const ParticleSet& P2)
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
    APP_ABORT("OperatorBase::two_body_quantum_domain(P1,P2)\n  quantum domain of input particles is invalid");
}

bool OperatorBase::quantum_domain_valid(quantum_domains qdomain) { return qdomain != no_quantum_domain; }

void OperatorBase::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  std::shared_ptr<OperatorBase> myclone = makeClone(qp, psi);
  if (myclone)
    targetH.addOperator(myclone, myName, UpdateMode[PHYSICAL]);
}

void OperatorBase::registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  bool collect = UpdateMode.test(COLLECTABLE);
  //exclude collectables
  if (!collect)
  {
    h5desc.emplace_back(myName);
    auto& oh = h5desc.back();
    std::vector<int> onedim(1, 1);
    oh.set_dimensions(onedim, myIndex);
    oh.open(gid);
  }
}

void OperatorBase::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  APP_ABORT("Need specialization for " + myName +
            "::addEnergy(MCWalkerConfiguration &W).\n Required functionality not implemented\n");
}

} // namespace qmcplusplus

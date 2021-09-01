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
//NonLocalData
NonLocalData::NonLocalData() : pid(-1), weight(1.0) {}

NonLocalData::NonLocalData(IndexType id, RealType w, const PosType& d) : pid(id), weight(w), delta(d) {}

//OperatorBase

//PUBLIC
OperatorBase::OperatorBase()
    : updateMode(std::move(std::bitset<8>{}.set(PRIMARY, 1))),
      value(0.0),
      myIndex(-1),
      newValue(0.0),
      tWalker(0),
#if !defined(REMOVE_TRACEMANAGER)
      haveRequiredTraces(false),
      streamingParticles(false),
#endif
      quantumDomain(quantum_domains::no_quantum_domain),
      energyDomain(energy_domains::no_energy_domain)
#if !defined(REMOVE_TRACEMANAGER)
      ,
      streamingScalars(false),
      valueSample(nullptr)
#endif
{}

OperatorBase::Return_t OperatorBase::getEnsembleAverage() { return 0.0; }

void OperatorBase::updateSource(ParticleSet& s) {}

void OperatorBase::setObservables(PropertySetType& plist) { plist[myIndex] = value; }

void OperatorBase::addObservables(PropertySetType& plist, BufferType& collectables) { addValue(plist); }

void OperatorBase::registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  bool collect = updateMode.test(COLLECTABLE);
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

void OperatorBase::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const {}

void OperatorBase::getRequiredTraces(TraceManager& tm) {}

void OperatorBase::setParticlePropertyList(PropertySetType& plist, int offset) { plist[myIndex + offset] = value; }

void OperatorBase::setRandomGenerator(RandomGenerator_t* rng) {}

OperatorBase::Return_t OperatorBase::evaluateDeterministic(ParticleSet& P) { return evaluate(P); }

OperatorBase::Return_t OperatorBase::evaluateWithToperator(ParticleSet& P) { return evaluate(P); }

// FIXME: unused variables commented off
OperatorBase::Return_t OperatorBase::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& /*ions*/,
                                                           TrialWaveFunction& /*psi*/,
                                                           ParticleSet::ParticlePos_t& /*hf_term*/,
                                                           ParticleSet::ParticlePos_t& /*pulay_term*/)
{
  return evaluate(P);
}

OperatorBase::Return_t OperatorBase::evaluateWithIonDerivsDeterministic(ParticleSet& P,
                                                                        ParticleSet& ions,
                                                                        TrialWaveFunction& psi,
                                                                        ParticleSet::ParticlePos_t& hf_term,
                                                                        ParticleSet::ParticlePos_t& pulay_term)
{
  //If there's no stochastic component, defaults to above defined evaluateWithIonDerivs.
  //If not otherwise specified, this defaults to evaluate().
  return evaluateWithIonDerivs(P, ions, psi, hf_term, pulay_term);
}

// FIXME: unused variables commented off
OperatorBase::Return_t OperatorBase::evaluateValueAndDerivatives(ParticleSet& P,
                                                                 const opt_variables_type& /*optvars*/,
                                                                 const std::vector<ValueType>& /*dlogpsi*/,
                                                                 std::vector<ValueType>& /*dhpsioverpsi*/)
{
  return evaluate(P);
}


void OperatorBase::mwEvaluate(const RefVectorWithLeader<OperatorBase>& O_list,
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

void OperatorBase::mwEvaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& O_list,
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

void OperatorBase::mwEvaluateWithToperator(const RefVectorWithLeader<OperatorBase>& O_list,
                                            const RefVectorWithLeader<ParticleSet>& P_list) const
{
  mwEvaluate(O_list, P_list);
}

void OperatorBase::setHistories(Walker_t& ThisWalker) { tWalker = &(ThisWalker); }

OperatorBase::Return_t OperatorBase::rejectedMove(ParticleSet& P) { return 0; }

void OperatorBase::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  std::unique_ptr<OperatorBase> myclone = makeClone(qp, psi);
  if (myclone)
  {
    targetH.addOperator(std::move(myclone), myName, updateMode[PHYSICAL]);
  }
}

void OperatorBase::createResource(ResourceCollection& collection) const {}

void OperatorBase::acquireResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& O_list) const
{}

void OperatorBase::releaseResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& O_list) const
{}

bool OperatorBase::isClassical() const noexcept { return quantumDomain == quantum_domains::classical; }

bool OperatorBase::isQuantum() const noexcept { return quantumDomain == quantum_domains::quantum; }

bool OperatorBase::isClassicalClassical() const noexcept
{
  return quantumDomain == quantum_domains::classical_classical;
}

bool OperatorBase::isQuantumClassical() const noexcept
{
  return quantumDomain == quantum_domains::quantum_classical;
}

bool OperatorBase::isQuantumQuantum() const noexcept { return quantumDomain == quantum_domains::quantum_quantum; }

bool OperatorBase::isNonLocal() const noexcept { return updateMode[NONLOCAL]; }

bool OperatorBase::getMode(int i) const noexcept { return updateMode[i]; }

#if !defined(REMOVE_TRACEMANAGER)
void OperatorBase::contributeTraceQuantities()
{
  contributeScalarQuantities();
  contributeParticleQuantities();
}

void OperatorBase::collectScalarTraces() { collectScalarQuantities(); }
#endif

void OperatorBase::checkoutTraceQuantities(TraceManager& tm)
{
  //derived classes must guard individual checkouts using request info
  checkoutScalarQuantities(tm);
  checkoutParticleQuantities(tm);
}

void OperatorBase::deleteTraceQuantities()
{
  deleteScalarQuantities();
  deleteParticleQuantities();
  streamingScalars    = false;
  streamingParticles  = false;
  haveRequiredTraces = false;
  request.reset();
}

// PROTECTED
#if !defined(REMOVE_TRACEMANAGER)
void OperatorBase::contributeScalarQuantities() { request.contribute_scalar(myName); }

void OperatorBase::checkoutScalarQuantities(TraceManager& tm)
{
  streamingScalars = request.streaming_scalar(myName);
  if (streamingScalars)
    valueSample = tm.checkout_real<1>(myName);
}

void OperatorBase::collectScalarQuantities()
{
  if (streamingScalars)
    (*valueSample)(0) = value;
}

void OperatorBase::deleteScalarQuantities()
{
  if (streamingScalars)
    delete valueSample;
}

void OperatorBase::contributeParticleQuantities(){};

void OperatorBase::checkoutParticleQuantities(TraceManager& /*tm*/){};

void OperatorBase::deleteParticleQuantities(){};

#endif

void OperatorBase::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  APP_ABORT("Need specialization for " + myName +
            "::addEnergy(MCWalkerConfiguration &W).\n Required functionality not implemented\n");
}

void OperatorBase::addEnergy(MCWalkerConfiguration& W,
                             std::vector<RealType>& LocalEnergy,
                             std::vector<std::vector<NonLocalData>>& /*Txy*/)
{
  addEnergy(W, LocalEnergy);
}

void OperatorBase::setComputeForces(bool compute) {}

void OperatorBase::setEnergyDomain(energy_domains edomain)
{
  if (energyDomainValid(edomain))
    energyDomain = edomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_energy_domain\n  input energy domain is invalid");
}

void OperatorBase::setQuantumDomain(quantum_domains qdomain)
{
  if (quantumDomainValid(qdomain))
    quantumDomain = qdomain;
  else
    APP_ABORT("QMCHamiltonainBase::set_quantum_domain\n  input quantum domain is invalid");
}

void OperatorBase::oneBodyQuantumDomain(const ParticleSet& P)
{
  if (P.is_classical())
    quantumDomain = quantum_domains::classical;
  else if (P.is_quantum())
    quantumDomain = quantum_domains::quantum;
  else
    APP_ABORT("OperatorBase::one_body_quantum_domain\n  quantum domain of input particles is invalid");
}

void OperatorBase::twoBodyQuantumDomain(const ParticleSet& P)
{
  if (P.is_classical())
    quantumDomain = quantum_domains::classical_classical;
  else if (P.is_quantum())
    quantumDomain = quantum_domains::quantum_quantum;
  else
    APP_ABORT("OperatorBase::two_body_quantum_domain(P)\n  quantum domain of input particles is invalid");
}

void OperatorBase::twoBodyQuantumDomain(const ParticleSet& P1, const ParticleSet& P2)
{
  bool c1 = P1.is_classical();
  bool c2 = P2.is_classical();
  bool q1 = P1.is_quantum();
  bool q2 = P2.is_quantum();
  if (c1 && c2)
    quantumDomain = quantum_domains::classical_classical;
  else if ((q1 && c2) || (c1 && q2))
    quantumDomain = quantum_domains::quantum_classical;
  else if (q1 && q2)
    quantumDomain = quantum_domains::quantum_quantum;
  else
    APP_ABORT("OperatorBase::two_body_quantum_domain(P1,P2)\n  quantum domain of input particles is invalid");
}

void OperatorBase::addValue(PropertySetType& plist)
{
  if (!updateMode[COLLECTABLE])
    myIndex = plist.add(myName.c_str());
}


// PRIVATE
bool OperatorBase::energyDomainValid(const energy_domains edomain) const noexcept
{
  return edomain != energy_domains::no_energy_domain;
}

bool OperatorBase::energyDomainValid() const noexcept { return energyDomainValid(energyDomain); }

bool OperatorBase::quantumDomainValid(const quantum_domains qdomain) const noexcept
{
  return qdomain != quantum_domains::no_quantum_domain;
}

bool OperatorBase::quantumDomainValid() const noexcept { return quantumDomainValid(quantumDomain); }


} // namespace qmcplusplus

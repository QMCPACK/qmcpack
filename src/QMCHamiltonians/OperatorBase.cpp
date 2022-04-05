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
// PUBLIC

OperatorBase::OperatorBase()
    : value_(0.0),
      my_index_(-1),
      t_walker_(0),
#if !defined(REMOVE_TRACEMANAGER)
      streaming_particles_(false),
      have_required_traces_(false),
      streaming_scalars_(false),
#endif
      quantum_domain_(NO_QUANTUM_DOMAIN),
      energy_domain_(NO_ENERGY_DOMAIN)
{
  update_mode_.set(PRIMARY, 1);
}

std::bitset<8>& OperatorBase::getUpdateMode() noexcept { return update_mode_; }

OperatorBase::Return_t OperatorBase::getValue() const noexcept { return value_; }

std::string OperatorBase::getName() const noexcept { return name_; }

void OperatorBase::setName(const std::string name) noexcept { name_ = name; }

#if !defined(REMOVE_TRACEMANAGER)
TraceRequest& OperatorBase::getRequest() noexcept { return request_; }
#endif

////////  FUNCTIONS ////////////////
void OperatorBase::addObservables(PropertySetType& plist, BufferType& collectables) { addValue(plist); }

void OperatorBase::registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  const bool collect = update_mode_.test(COLLECTABLE);
  //exclude collectables
  if (!collect)
  {
    h5desc.emplace_back(name_);
    auto& oh = h5desc.back();
    std::vector<int> onedim(1, 1);
    oh.set_dimensions(onedim, my_index_);
    oh.open(gid);
  }
}

void OperatorBase::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const {}

void OperatorBase::setObservables(PropertySetType& plist) { plist[my_index_] = value_; }

void OperatorBase::setParticlePropertyList(PropertySetType& plist, int offset) { plist[my_index_ + offset] = value_; }

void OperatorBase::setHistories(Walker_t& ThisWalker) { t_walker_ = &(ThisWalker); }

OperatorBase::Return_t OperatorBase::evaluateDeterministic(ParticleSet& P) { return evaluate(P); }

void OperatorBase::mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                               const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                               const RefVectorWithLeader<ParticleSet>& p_list) const
{
  assert(this == &o_list.getLeader());
/**  Temporary raw omp pragma for simple thread parallelism
   *   ignoring the driver level concurrency
   *   
   *  TODO: replace this with a proper abstraction. It should adequately describe the behavior
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
  for (int iw = 0; iw < o_list.size(); iw++)
    o_list[iw].evaluate(p_list[iw]);
}

void OperatorBase::mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& o_list,
                                                       const RefVectorWithLeader<ParticleSet>& p_list,
                                                       const opt_variables_type& optvars,
                                                       RecordArray<ValueType>& dlogpsi,
                                                       RecordArray<ValueType>& dhpsioverpsi) const
{
  const int nparam = dlogpsi.nparam();
  std::vector<ValueType> tmp_dlogpsi(nparam);
  std::vector<ValueType> tmp_dhpsioverpsi(nparam);
  for (int iw = 0; iw < o_list.size(); iw++)
  {
    for (int j = 0; j < nparam; j++)
    {
      tmp_dlogpsi[j] = dlogpsi.getValue(j, iw);
    }

    o_list[iw].evaluateValueAndDerivatives(p_list[iw], optvars, tmp_dlogpsi, tmp_dhpsioverpsi);

    for (int j = 0; j < nparam; j++)
    {
      dhpsioverpsi.setValue(j, iw, dhpsioverpsi.getValue(j, iw) + tmp_dhpsioverpsi[j]);
    }
  }
}

OperatorBase::Return_t OperatorBase::rejectedMove(ParticleSet& P) { return 0; }

OperatorBase::Return_t OperatorBase::evaluateWithToperator(ParticleSet& P) { return evaluate(P); }

void OperatorBase::mw_evaluateWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                            const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                            const RefVectorWithLeader<ParticleSet>& p_list) const
{
  mw_evaluate(o_list, wf_list, p_list);
}

OperatorBase::Return_t OperatorBase::evaluateValueAndDerivatives(ParticleSet& P,
                                                                 const opt_variables_type& optvars,
                                                                 const std::vector<ValueType>& dlogpsi,
                                                                 std::vector<ValueType>& dhpsioverpsi)
{
  return evaluate(P);
}

OperatorBase::Return_t OperatorBase::evaluateWithIonDerivs(ParticleSet& P,
                                                           ParticleSet& ions,
                                                           TrialWaveFunction& psi,
                                                           ParticleSet::ParticlePos& hf_term,
                                                           ParticleSet::ParticlePos& pulay_term)
{
  return evaluate(P);
}

OperatorBase::Return_t OperatorBase::evaluateWithIonDerivsDeterministic(ParticleSet& P,
                                                                        ParticleSet& ions,
                                                                        TrialWaveFunction& psi,
                                                                        ParticleSet::ParticlePos& hf_term,
                                                                        ParticleSet::ParticlePos& pulay_term)
{
  return evaluateWithIonDerivs(P, ions, psi, hf_term, pulay_term);
}

void OperatorBase::updateSource(ParticleSet& s) {}

OperatorBase::Return_t OperatorBase::getEnsembleAverage() { return 0.0; }

void OperatorBase::createResource(ResourceCollection& collection) const {}

void OperatorBase::acquireResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{}

void OperatorBase::releaseResource(ResourceCollection& collection,
                                   const RefVectorWithLeader<OperatorBase>& o_list) const
{}

void OperatorBase::setRandomGenerator(RandomGenerator* rng) {}

void OperatorBase::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH)
{
  std::unique_ptr<OperatorBase> myclone = makeClone(qp, psi);
  if (myclone)
  {
    targetH.addOperator(std::move(myclone), name_, update_mode_[PHYSICAL]);
  }
}

#if !defined(REMOVE_TRACEMANAGER)
void OperatorBase::getRequiredTraces(TraceManager& tm){};
#endif

void OperatorBase::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  APP_ABORT("Need specialization for " + name_ +
            "::addEnergy(MCWalkerConfiguration &W).\n Required functionality not implemented\n");
}

void OperatorBase::addEnergy(MCWalkerConfiguration& W,
                             std::vector<RealType>& LocalEnergy,
                             std::vector<std::vector<NonLocalData>>& Txy)
{
  addEnergy(W, LocalEnergy);
}

// END  FUNCTIONS //

bool OperatorBase::isClassical() const noexcept { return quantum_domain_ == CLASSICAL; }

bool OperatorBase::isQuantum() const noexcept { return quantum_domain_ == QUANTUM; }

bool OperatorBase::isClassicalClassical() const noexcept { return quantum_domain_ == CLASSICAL_CLASSICAL; }

bool OperatorBase::isQuantumClassical() const noexcept { return quantum_domain_ == QUANTUM_CLASSICAL; }

bool OperatorBase::isQuantumQuantum() const noexcept { return quantum_domain_ == QUANTUM_QUANTUM; }

bool OperatorBase::getMode(const int i) const noexcept { return update_mode_[i]; }

bool OperatorBase::isNonLocal() const noexcept { return update_mode_[NONLOCAL]; }


#if !defined(REMOVE_TRACEMANAGER)

void OperatorBase::contributeTraceQuantities()
{
  contributeScalarQuantities();
  contributeParticleQuantities();
}

void OperatorBase::checkoutTraceQuantities(TraceManager& tm)
{
  checkoutScalarQuantities(tm);
  checkoutParticleQuantities(tm);
}

void OperatorBase::collectScalarTraces() { collectScalarQuantities(); }

void OperatorBase::deleteTraceQuantities()
{
  deleteScalarQuantities();
  deleteParticleQuantities();
  streaming_scalars_    = false;
  streaming_particles_  = false;
  have_required_traces_ = false;
  request_.reset();
}

#endif

////// PROTECTED FUNCTIONS
#if !defined(REMOVE_TRACEMANAGER)
void OperatorBase::contributeScalarQuantities() { request_.contribute_scalar(name_); }

void OperatorBase::checkoutScalarQuantities(TraceManager& tm)
{
  streaming_scalars_ = request_.streaming_scalar(name_);
  if (streaming_scalars_)
    value_sample_ = tm.checkout_real<1>(name_);
}

void OperatorBase::collectScalarQuantities()
{
  if (streaming_scalars_)
    (*value_sample_)(0) = value_;
}

void OperatorBase::deleteScalarQuantities()
{
  if (streaming_scalars_)
    delete value_sample_;
}

void OperatorBase::contributeParticleQuantities(){};
void OperatorBase::checkoutParticleQuantities(TraceManager& tm){};
void OperatorBase::deleteParticleQuantities(){};
#endif

void OperatorBase::setComputeForces(bool compute) {}

void OperatorBase::setEnergyDomain(EnergyDomains edomain)
{
  if (energyDomainValid(edomain))
    energy_domain_ = edomain;
  else
    APP_ABORT("QMCHamiltonainBase::setEnergyDomain\n  input energy domain is invalid");
}

void OperatorBase::setQuantumDomain(QuantumDomains qdomain)
{
  if (quantumDomainValid(qdomain))
    quantum_domain_ = qdomain;
  else
    APP_ABORT("QMCHamiltonainBase::setQuantumDomain\n  input quantum domain is invalid");
}

void OperatorBase::oneBodyQuantumDomain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain_ = CLASSICAL;
  else if (P.is_quantum())
    quantum_domain_ = QUANTUM;
  else
    APP_ABORT("OperatorBase::oneBodyQuantumDomain\n  quantum domain of input particles is invalid");
}

void OperatorBase::twoBodyQuantumDomain(const ParticleSet& P)
{
  if (P.is_classical())
    quantum_domain_ = CLASSICAL_CLASSICAL;
  else if (P.is_quantum())
    quantum_domain_ = QUANTUM_QUANTUM;
  else
    APP_ABORT("OperatorBase::twoBodyQuantumDomain(P)\n  quantum domain of input particles is invalid");
}

void OperatorBase::twoBodyQuantumDomain(const ParticleSet& P1, const ParticleSet& P2)
{
  const bool c1 = P1.is_classical();
  const bool c2 = P2.is_classical();
  const bool q1 = P1.is_quantum();
  const bool q2 = P2.is_quantum();
  if (c1 && c2)
    quantum_domain_ = CLASSICAL_CLASSICAL;
  else if ((q1 && c2) || (c1 && q2))
    quantum_domain_ = QUANTUM_CLASSICAL;
  else if (q1 && q2)
    quantum_domain_ = QUANTUM_QUANTUM;
  else
    APP_ABORT("OperatorBase::twoBodyQuantumDomain(P1,P2)\n  quantum domain of input particles is invalid");
}

void OperatorBase::addValue(PropertySetType& plist)
{
  if (!update_mode_[COLLECTABLE])
    my_index_ = plist.add(name_.c_str());
}

////// PRIVATE FUNCTIONS
bool OperatorBase::energyDomainValid(EnergyDomains edomain) const noexcept { return edomain != NO_ENERGY_DOMAIN; }

bool OperatorBase::quantumDomainValid(QuantumDomains qdomain) const noexcept { return qdomain != NO_QUANTUM_DOMAIN; }

} // namespace qmcplusplus

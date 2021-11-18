//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCHamiltonian.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "Utilities/TimerManager.h"
#include "Containers/MinimalContainers/RecordArray.hpp"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif

namespace qmcplusplus
{
/** constructor
*/
QMCHamiltonian::QMCHamiltonian(const std::string& aname)
    : myIndex(0),
      numCollectables(0),
      myName(aname),
      nlpp_ptr(nullptr),
      l2_ptr(nullptr),
      ham_timer_(*timer_manager.createTimer("Hamiltonian:" + aname, timer_level_medium))
#if !defined(REMOVE_TRACEMANAGER)
      ,
      streaming_position(false),
      id_sample(nullptr),
      pid_sample(nullptr),
      step_sample(nullptr),
      gen_sample(nullptr),
      age_sample(nullptr),
      mult_sample(nullptr),
      weight_sample(nullptr),
      position_sample(nullptr)
#endif
{}

///// copy constructor is distable by declaring it as private
//QMCHamiltonian::QMCHamiltonian(const QMCHamiltonian& qh) {}

/** destructor
 */
QMCHamiltonian::~QMCHamiltonian()
{
  //@todo clean up H and auxH
}

bool QMCHamiltonian::get(std::ostream& os) const
{
  for (int i = 0; i < H.size(); i++)
  {
    os.setf(std::ios::left);
    os << "  " << std::setw(16) << H[i]->getName();
    H[i]->get(os);
    os << "\n";
  }
  return true;
}

/** add a new Hamiltonian the the list of Hamiltonians.
 * @param h an operator
 * @param aname name of h
 * @param physical if true, a physical operator
 */
void QMCHamiltonian::addOperator(std::unique_ptr<OperatorBase>&& h, const std::string& aname, bool physical)
{
  //change UpdateMode[PHYSICAL] of h so that cloning can be done correctly
  h->getUpdateMode()[OperatorBase::PHYSICAL] = physical;
  if (physical)
  {
    for (int i = 0; i < H.size(); ++i)
    {
      if (H[i]->getName() == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << std::endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to H, physical Hamiltonian " << std::endl;
    h->setName(aname);
    H.push_back(std::move(h));
    std::string tname = "Hamiltonian:" + aname;
    my_timers_.push_back(*timer_manager.createTimer(tname, timer_level_fine));
  }
  else
  {
    //ignore timers for now
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (auxH[i]->getName() == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << std::endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to auxH " << std::endl;
    h->setName(aname);
    auxH.push_back(std::move(h));
  }

  //assign save NLPP if found
  //  name is fixed in ECPotentialBuilder::put()
  if (aname == "NonLocalECP")
  {
    if (nlpp_ptr == nullptr)
    {
      // original h arguments moved to either H or auxH
      nlpp_ptr = physical ? dynamic_cast<NonLocalECPotential*>(H.back().get())
                          : dynamic_cast<NonLocalECPotential*>(auxH.back().get());
    }
    else
    {
      APP_ABORT("QMCHamiltonian::addOperator nlpp_ptr is supposed to be null. Something went wrong!");
    }
  }

  //save L2 potential if found
  //  name is fixed in ECPotentialBuilder::put()
  if (aname == "L2")
  {
    if (l2_ptr == nullptr)
    {
      l2_ptr = physical ? dynamic_cast<L2Potential*>(H.back().get()) : dynamic_cast<L2Potential*>(auxH.back().get());
    }
    else
    {
      APP_ABORT("QMCHamiltonian::addOperator l2_ptr is supposed to be null. Something went wrong!");
    }
  }
}


void QMCHamiltonian::addOperatorType(const std::string& name, const std::string& type)
{
  app_log() << "QMCHamiltonian::addOperatorType added type " << type << " named " << name << std::endl;
  operator_types[name] = type;
}


const std::string& QMCHamiltonian::getOperatorType(const std::string& name)
{
  std::map<std::string, std::string>::iterator type = operator_types.find(name);
  if (type == operator_types.end())
    APP_ABORT("QMCHamiltonian::getOperatorType\n  operator type not found for name " + name);
  return type->second;
}

///** remove a named Hamiltonian from the list
// *@param aname the name of the Hamiltonian
// *@return true, if the request hamiltonian exists and is removed.
// */
//bool
//QMCHamiltonian::remove(const std::string& aname)
//{
//  return false;
//}

void QMCHamiltonian::updateSource(ParticleSet& s)
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->updateSource(s);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->updateSource(s);
}

/** add a number of properties to the ParticleSet
 * @param P ParticleSet to which multiple columns to be added
 *
 * QMCHamiltonian can add any number of properties to a ParticleSet.
 * Hindex contains the index map to the ParticleSet::PropertyList.
 * This enables assigning the properties evaluated by each OperatorBase
 * object to the correct property column.
 */
//void
//QMCHamiltonian::addObservables(PropertySetType& plist)
//{
//  //first add properties to Observables
//  Observables.clear();
//  for(int i=0; i<H.size(); ++i) H[i]->addObservables(Observables);
//  for(int i=0; i<auxH.size(); ++i) auxH[i]->addObservables(Observables);
//
//  myIndex=plist.add(Observables.Names[0]);
//  for(int i=1; i<Observables.size(); ++i) plist.add(Observables.Names[i]);
//
//  app_log() << "  QMCHamiltonian::add2WalkerProperty added " << Observables.size() << " data to PropertyList" << std::endl;
//  app_log() << "    starting Index = " << myIndex << std::endl;
//}
//
int QMCHamiltonian::addObservables(ParticleSet& P)
{
  //first add properties to Observables
  Observables.clear();
  //ParticleSet::mcObservables (large data, e.g. density) are accumulated while evaluations
  P.Collectables.clear();
  P.Collectables.rewind();
  for (int i = 0; i < H.size(); ++i)
    H[i]->addObservables(Observables, P.Collectables);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->addObservables(Observables, P.Collectables);
  myIndex = P.PropertyList.add(Observables.Names[0]);
  for (int i = 1; i < Observables.size(); ++i)
    P.PropertyList.add(Observables.Names[i]);
  numCollectables = P.Collectables.size();
  app_log() << "\n  QMCHamiltonian::add2WalkerProperty added"
            << "\n    " << Observables.size() << " to P::PropertyList "
            << "\n    " << P.Collectables.size() << " to P::Collectables "
            << "\n    starting Index of the observables in P::PropertyList = " << myIndex << std::endl;
  return Observables.size();
}

void QMCHamiltonian::resetObservables(int start, int ncollects)
{
  Observables.clear();
  BufferType collectables;
  collectables.rewind();
  for (int i = 0; i < H.size(); ++i)
    H[i]->addObservables(Observables, collectables);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->addObservables(Observables, collectables);
  if (collectables.size() != ncollects)
  {
    APP_ABORT("  QMCHamiltonian::resetObservables numCollectables != ncollects");
  }
  myIndex         = start;
  numCollectables = ncollects;
}

void QMCHamiltonian::registerObservables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->registerObservables(h5desc, gid);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->registerObservables(h5desc, gid);
}

void QMCHamiltonian::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  //The physical operators cannot add to collectables
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->registerCollectables(h5desc, gid);
}


#if !defined(REMOVE_TRACEMANAGER)
void QMCHamiltonian::initialize_traces(TraceManager& tm, ParticleSet& P)
{
  static bool first_init = true;
  bool trace_log         = first_init && tm.verbose && omp_get_thread_num() == 0;
  if (trace_log)
    app_log() << "\n  Hamiltonian is initializing traces" << std::endl;

  //fill std::string vectors for combined trace quantities
  std::vector<std::string> Eloc;
  std::vector<std::string> Vloc;
  std::vector<std::string> Vq, Vc, Vqq, Vqc, Vcc;
  for (int i = 0; i < H.size(); ++i)
    Eloc.push_back(H[i]->getName());
  for (int i = 1; i < H.size(); ++i)
    Vloc.push_back(H[i]->getName());

  // These contributions are based on the potential energy components.
  // Loop starts at one to skip the kinetic energy component.
  for (int i = 1; i < H.size(); ++i)
  {
    OperatorBase& h = *H[i];
    if (h.isQuantum())
      Vq.push_back(h.getName());
    else if (h.isClassical())
      Vc.push_back(h.getName());
    else if (h.isQuantumQuantum())
      Vqq.push_back(h.getName());
    else if (h.isQuantumClassical())
      Vqc.push_back(h.getName());
    else if (h.isClassicalClassical())
      Vcc.push_back(h.getName());
    else if (omp_get_thread_num() == 0)
      app_log() << "  warning: potential named " << h.getName()
                << " has not been classified according to its quantum domain (q,c,qq,qc,cc)\n    estimators depending "
                   "on this classification may not function properly"
                << std::endl;
  }


  //make trace quantities available
  request.contribute_scalar("id", true);           //default trace quantity
  request.contribute_scalar("parent_id", true);    //default trace quantity
  request.contribute_scalar("step", true);         //default trace quantity
  request.contribute_scalar("generation", true);   //default trace quantity
  request.contribute_scalar("age", true);          //default trace quantity
  request.contribute_scalar("multiplicity", true); //default trace quantity
  request.contribute_scalar("weight", true);       //default trace quantity
  request.contribute_array("position");
  for (int i = 0; i < H.size(); ++i)
    H[i]->contributeTraceQuantities();
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->contributeTraceQuantities();


  //note availability of combined quantities
  request.contribute_combined("LocalEnergy", Eloc, true);
  request.contribute_combined("LocalPotential", Vloc, true, true);
  if (Vq.size() > 0)
    request.contribute_combined("Vq", Vq, true, true);
  if (Vc.size() > 0)
    request.contribute_combined("Vc", Vc, true, true);
  if (Vqq.size() > 0)
    request.contribute_combined("Vqq", Vqq, true, true);
  if (Vqc.size() > 0)
    request.contribute_combined("Vqc", Vqc, true, true);
  if (Vcc.size() > 0)
    request.contribute_combined("Vcc", Vcc, true, true);


  //collect trace requests
  std::vector<TraceRequest*> requests;
  //  Hamiltonian request (id, step, weight, positions)
  requests.push_back(&request);
  //  requests from Hamiltonian components
  for (int i = 0; i < H.size(); ++i)
    requests.push_back(&H[i]->getRequest());
  //  requests from other observables
  for (int i = 0; i < auxH.size(); ++i)
    requests.push_back(&auxH[i]->getRequest());

  //collect trace quantity availability/requests from contributors/requestors
  for (int i = 0; i < requests.size(); ++i)
    tm.request.incorporate(*requests[i]);

  //balance requests with availability, mark quantities as streaming/writing
  tm.request.determine_stream_write();

  //relay updated streaming information to all contributors/requestors
  for (int i = 0; i < requests.size(); ++i)
    tm.request.relay_stream_info(*requests[i]);

  //set streaming/writing traces in general
  tm.update_status();

  // setup traces, if any quantities should be streaming

  // tracing
  bool tracing = request.streaming();
  if (tracing != tm.streaming_traces)
    APP_ABORT("QMCHamiltonian::initialize_traces  trace request failed to initialize properly");
  if (!tracing)
  {
    // Empty. Do not log if nothing will be done

    if (trace_log)
      app_log() << "    no traces streaming" << std::endl;
  }
  else
  {
    if (trace_log)
      app_log() << "    traces streaming" << std::endl;
    //checkout trace quantities
    //(requested sources checkout arrays to place samples in for streaming)
    //  checkout walker trace quantities
    streaming_position = request.streaming_array("position");
    if (request.streaming_default_scalars)
    {
      id_sample     = tm.checkout_int<1>("id");
      pid_sample    = tm.checkout_int<1>("parent_id");
      step_sample   = tm.checkout_int<1>("step");
      gen_sample    = tm.checkout_int<1>("generation");
      age_sample    = tm.checkout_int<1>("age");
      mult_sample   = tm.checkout_int<1>("multiplicity");
      weight_sample = tm.checkout_real<1>("weight");
    }
    if (streaming_position)
      position_sample = tm.checkout_real<2>("position", P, DIM);
    //  checkout observable trace quantities
    for (int i = 0; i < H.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::checkoutTraceQuantities  " << H[i]->getName() << std::endl;
      H[i]->checkoutTraceQuantities(tm);
    }
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::checkoutTraceQuantities  " << auxH[i]->getName() << std::endl;
      auxH[i]->checkoutTraceQuantities(tm);
    }
    //setup combined traces that depend on H information
    //  LocalEnergy, LocalPotential, Vq, Vc, Vqq, Vqc, Vcc
    if (Vloc.size() > 0 && request.streaming("LocalPotential"))
      tm.make_combined_trace("LocalPotential", Vloc);
    if (Eloc.size() > 0 && request.streaming("LocalEnergy"))
      tm.make_combined_trace("LocalEnergy", Eloc);
    if (Vq.size() > 0 && request.streaming("Vq"))
      tm.make_combined_trace("Vq", Eloc);
    if (Vc.size() > 0 && request.streaming("Vc"))
      tm.make_combined_trace("Vc", Eloc);
    if (Vqq.size() > 0 && request.streaming("Vqq"))
      tm.make_combined_trace("Vqq", Eloc);
    if (Vqc.size() > 0 && request.streaming("Vqc"))
      tm.make_combined_trace("Vqc", Eloc);
    if (Vcc.size() > 0 && request.streaming("Vcc"))
      tm.make_combined_trace("Vcc", Eloc);

    //all trace samples have been created ( streaming instances)
    //  mark the ones that will be writing also
    tm.screen_writes();

    //observables that depend on traces check them out
    if (trace_log)
      app_log() << "\n  Hamiltonian is fulfilling trace requests from observables" << std::endl;
    for (int i = 0; i < auxH.size(); ++i)
    {
      if (trace_log)
        app_log() << "    OperatorBase::getRequiredTraces  " << auxH[i]->getName() << std::endl;
      auxH[i]->getRequiredTraces(tm);
    }
    //report

    //write traces status to the log
    if (trace_log)
      tm.user_report();

    first_init = false;
  }
}


void QMCHamiltonian::collect_walker_traces(Walker_t& walker, int step)
{
  if (request.streaming_default_scalars)
  {
    (*id_sample)(0)     = walker.ID;
    (*pid_sample)(0)    = walker.ParentID;
    (*step_sample)(0)   = step;
    (*gen_sample)(0)    = walker.Generation;
    (*age_sample)(0)    = walker.Age;
    (*mult_sample)(0)   = walker.Multiplicity;
    (*weight_sample)(0) = walker.Weight;
  }
  if (streaming_position)
    for (int i = 0; i < walker.R.size(); ++i)
      for (int d = 0; d < DIM; ++d)
        (*position_sample)(i, d) = walker.R[i][d];
}


void QMCHamiltonian::finalize_traces()
{
  if (request.streaming_default_scalars)
  {
    delete id_sample;
    delete pid_sample;
    delete step_sample;
    delete gen_sample;
    delete age_sample;
    delete mult_sample;
    delete weight_sample;
  }
  if (streaming_position)
    delete position_sample;
  if (request.streaming())
  {
    for (int i = 0; i < H.size(); ++i)
      H[i]->deleteTraceQuantities();
    for (int i = 0; i < auxH.size(); ++i)
      auxH[i]->deleteTraceQuantities();
  }
  streaming_position = false;
  request.reset();
}
#endif

/** Evaluate all the Hamiltonians for the N-particle  configuration
 *@param P input configuration containing N particles
 *@return the local energy
 */
QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluate(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    const auto LocalEnergyComponent = H[i]->evaluate(P);
    if (std::isnan(LocalEnergyComponent))
      APP_ABORT("QMCHamiltonian::evaluate component " + H[i]->getName() + " returns NaN\n");
    LocalEnergy += LocalEnergyComponent;
    H[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
    H[i]->setParticlePropertyList(P.PropertyList, myIndex);
  }
  KineticEnergy                      = H[0]->getValue();
  P.PropertyList[WP::LOCALENERGY]    = LocalEnergy;
  P.PropertyList[WP::LOCALPOTENTIAL] = LocalEnergy - KineticEnergy;
  // auxHevaluate(P);
  return LocalEnergy;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateDeterministic(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    const auto LocalEnergyComponent = H[i]->evaluateDeterministic(P);
    if (std::isnan(LocalEnergyComponent))
      APP_ABORT("QMCHamiltonian::evaluate component " + H[i]->getName() + " returns NaN\n");
    LocalEnergy += LocalEnergyComponent;
    H[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
    H[i]->setParticlePropertyList(P.PropertyList, myIndex);
  }
  KineticEnergy                      = H[0]->getValue();
  P.PropertyList[WP::LOCALENERGY]    = LocalEnergy;
  P.PropertyList[WP::LOCALPOTENTIAL] = LocalEnergy - KineticEnergy;
  // auxHevaluate(P);
  return LocalEnergy;
}
void QMCHamiltonian::updateNonKinetic(OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset)
{
  if (std::isnan(op.getValue()))
    APP_ABORT("QMCHamiltonian::evaluate component " + op.getName() + " returns NaN\n");
  // The following is a ridiculous breach of encapsulation.
  ham.LocalEnergy += op.getValue();
  op.setObservables(ham.Observables);
  op.setParticlePropertyList(pset.PropertyList, ham.myIndex);
}

void QMCHamiltonian::updateKinetic(OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset)
{
  ham.KineticEnergy                     = op.getValue();
  pset.PropertyList[WP::LOCALENERGY]    = ham.LocalEnergy;
  pset.PropertyList[WP::LOCALPOTENTIAL] = ham.LocalEnergy - ham.KineticEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluate(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& ham_leader = ham_list.getLeader();
  ScopedTimer local_timer(ham_leader.ham_timer_);
  for (QMCHamiltonian& ham : ham_list)
    ham.LocalEnergy = 0.0;

  const int num_ham_operators = ham_leader.H.size();
  for (int i_ham_op = 0; i_ham_op < num_ham_operators; ++i_ham_op)
  {
    ScopedTimer h_timer(ham_leader.my_timers_[i_ham_op]);
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));

    // // This lambda accomplishes two things
    // // 1. It makes clear T& and not std::reference_wrapper<T> is desired removing need for gets.
    // // 2. [] captures nothing insuring that we know these updates only depend on the three object involved.
    // auto updateNonKinetic = [](OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset) {
    //   // both hamiltonian and operatorbase should have operator<< overides
    //   if (std::isnan(op.Value))
    //     APP_ABORT("QMCHamiltonian::evaluate component " + op.myName + " returns NaN\n");

    //   // The following is a ridiculous breach of encapsulation.
    //   ham.LocalEnergy += op.Value;
    //   op.setObservables(ham.Observables);
    //   op.setParticlePropertyList(pset.PropertyList, ham.myIndex);
    // };
    ham_leader.H[i_ham_op]->mw_evaluate(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); iw++)
      updateNonKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  // auto updateKinetic = [](OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset) {
  //   ham.KineticEnergy                 = op.Value;
  //   pset.PropertyList[WP::LOCALENERGY]    = ham.LocalEnergy;
  //   pset.PropertyList[LOCALPOTENTIAL] = ham.LocalEnergy - ham.KineticEnergy;
  // };

  for (int iw = 0; iw < ham_list.size(); iw++)
  {
    const auto HC_list(extract_HC_list(ham_list, 0));
    updateKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  std::vector<FullPrecRealType> local_energies(ham_list.size(), 0.0);
  for (int iw = 0; iw < ham_list.size(); ++iw)
    local_energies[iw] = ham_list[iw].get_LocalEnergy();

  return local_energies;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateValueAndDerivatives(ParticleSet& P,
                                                                             const opt_variables_type& optvars,
                                                                             std::vector<ValueType>& dlogpsi,
                                                                             std::vector<ValueType>& dhpsioverpsi,
                                                                             bool compute_deriv)
{
  LocalEnergy = KineticEnergy = H[0]->evaluate(P);
  if (compute_deriv)
    for (int i = 1; i < H.size(); ++i)
      LocalEnergy += H[i]->evaluateValueAndDerivatives(P, optvars, dlogpsi, dhpsioverpsi);
  else
    for (int i = 1; i < H.size(); ++i)
      LocalEnergy += H[i]->evaluate(P);
  return LocalEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluateValueAndDerivativesInner(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const opt_variables_type& optvars,
    RecordArray<ValueType>& dlogpsi,
    RecordArray<ValueType>& dhpsioverpsi)
{
  std::vector<FullPrecRealType> local_energies(ham_list.size(), 0.0);
  for (int iw = 0; iw < ham_list.size(); iw++)
    ham_list[iw].LocalEnergy = 0.0;

  if (ham_list.size() > 0)
  {
    auto& ham_leader            = ham_list.getLeader();
    const int num_ham_operators = ham_leader.H.size();
    for (int i_ham_op = 0; i_ham_op < num_ham_operators; ++i_ham_op)
    {
      ScopedTimer local_timer(ham_leader.my_timers_[i_ham_op]);
      const auto HC_list(extract_HC_list(ham_list, i_ham_op));

      ham_leader.H[i_ham_op]->mw_evaluateWithParameterDerivatives(HC_list, p_list, optvars, dlogpsi, dhpsioverpsi);

      for (int iw = 0; iw < ham_list.size(); iw++)
        updateNonKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
    }

    for (int iw = 0; iw < ham_list.size(); iw++)
    {
      const auto HC_list(extract_HC_list(ham_list, 0));
      updateKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
    }

    for (int iw = 0; iw < ham_list.size(); ++iw)
      local_energies[iw] = ham_list[iw].get_LocalEnergy();
  }

  return local_energies;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluateValueAndDerivatives(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const opt_variables_type& optvars,
    RecordArray<ValueType>& dlogpsi,
    RecordArray<ValueType>& dhpsioverpsi,
    bool compute_deriv)
{
  std::vector<FullPrecRealType> local_energies(ham_list.size(), 0.0);
  if (compute_deriv)
    local_energies =
        QMCHamiltonian::mw_evaluateValueAndDerivativesInner(ham_list, wf_list, p_list, optvars, dlogpsi, dhpsioverpsi);
  else
    local_energies = QMCHamiltonian::mw_evaluate(ham_list, wf_list, p_list);

  return local_energies;
}


QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateVariableEnergy(ParticleSet& P, bool free_nlpp)
{
  RealType nlpp = 0.0;
  RealType ke   = H[0]->evaluate(P);
  if (free_nlpp)
    for (int i = 1; i < H.size(); ++i)
    {
      if (H[i]->isNonLocal())
        nlpp += H[i]->evaluate(P);
    }
  return ke + nlpp;
}

void QMCHamiltonian::auxHevaluate(ParticleSet& P)
{
  for (int i = 0; i < auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    auxH[i]->collectScalarTraces();
#endif
    auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
    //H[i]->setParticlePropertyList(P.PropertyList,myIndex);
  }
}

///This is more efficient. Only calculate auxH elements if moves are accepted.
void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker)
{
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  for (int i = 0; i < auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    auxH[i]->collectScalarTraces();
#endif
    auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
  }
}
///Evaluate properties only.
void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker, bool do_properties, bool do_collectables)
{
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  for (int i = 0; i < auxH.size(); ++i)
  {
    bool is_property    = !(auxH[i]->getMode(OperatorBase::COLLECTABLE));
    bool is_collectable = (auxH[i]->getMode(OperatorBase::COLLECTABLE));
    if ((is_property && do_properties) || (is_collectable && do_collectables))
    {
      auxH[i]->setHistories(ThisWalker);
      RealType sink = auxH[i]->evaluate(P);
      auxH[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
      auxH[i]->collectScalarTraces();
#endif
      auxH[i]->setParticlePropertyList(P.PropertyList, myIndex);
    }
  }
}

/** Looks like a hack see DMCBatched.cpp and DMC.cpp weight is used like temporary flag
 *  from DMC.
 */
void QMCHamiltonian::rejectedMove(ParticleSet& P, Walker_t& ThisWalker)
{
  // weight should be 0 from DMC
  //   since other traced properties will be invalid
  //   (they will be from the walker moved before this one)
#if !defined(REMOVE_TRACEMANAGER)
  collect_walker_traces(ThisWalker, P.current_step);
#endif
  //   ThisWalker.rejectedMove();
  for (int i = 0; i < auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->rejectedMove(P);
  }
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateWithToperator(ParticleSet& P)
{
  ScopedTimer local_timer(ham_timer_);
  LocalEnergy = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    LocalEnergy += H[i]->evaluateWithToperator(P);
    H[i]->setObservables(Observables);
#if !defined(REMOVE_TRACEMANAGER)
    H[i]->collectScalarTraces();
#endif
  }
  KineticEnergy                      = H[0]->getValue();
  P.PropertyList[WP::LOCALENERGY]    = LocalEnergy;
  P.PropertyList[WP::LOCALPOTENTIAL] = LocalEnergy - KineticEnergy;
  //   auxHevaluate(P);
  return LocalEnergy;
}

std::vector<QMCHamiltonian::FullPrecRealType> QMCHamiltonian::mw_evaluateWithToperator(
    const RefVectorWithLeader<QMCHamiltonian>& ham_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list)
{
  for (QMCHamiltonian& ham : ham_list)
    ham.LocalEnergy = 0.0;

  auto& ham_leader            = ham_list.getLeader();
  const int num_ham_operators = ham_leader.H.size();
  for (int i_ham_op = 0; i_ham_op < num_ham_operators; ++i_ham_op)
  {
    ScopedTimer local_timer(ham_leader.my_timers_[i_ham_op]);
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));

    ham_leader.H[i_ham_op]->mw_evaluateWithToperator(HC_list, wf_list, p_list);
    for (int iw = 0; iw < ham_list.size(); ++iw)
      updateNonKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  for (int iw = 0; iw < ham_list.size(); iw++)
  {
    const auto HC_list(extract_HC_list(ham_list, 0));
    updateKinetic(HC_list[iw], ham_list[iw], p_list[iw]);
  }

  std::vector<FullPrecRealType> local_energies(ham_list.size());
  for (int iw = 0; iw < ham_list.size(); ++iw)
    local_energies[iw] = ham_list[iw].get_LocalEnergy();

  return local_energies;
}
void QMCHamiltonian::evaluateElecGrad(ParticleSet& P,
                                      TrialWaveFunction& psi,
                                      ParticleSet::ParticlePos_t& Egrad,
                                      RealType delta)
{
  int nelec = P.getTotalNum();
  RealType ep(0.0);
  RealType em(0.0);
  RealType e0(0.0);
  for (int iel = 0; iel < nelec; iel++)
  {
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      RealType r0 = P.R[iel][dim];
      ep          = 0;
      em          = 0;
      //Plus
      RealType rp   = r0 + delta;
      P.R[iel][dim] = rp;
      P.update();
      psi.evaluateLog(P);
      ep = evaluateDeterministic(P);

      //minus
      RealType rm   = r0 - delta;
      P.R[iel][dim] = rm;
      P.update();
      psi.evaluateLog(P);
      em = evaluateDeterministic(P);

      Egrad[iel][dim] = (ep - em) / (2.0 * delta);
      P.R[iel][dim]   = r0;
      P.update();
      psi.evaluateLog(P);
    }
  }
}
QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateIonDerivs(ParticleSet& P,
                                                                   ParticleSet& ions,
                                                                   TrialWaveFunction& psi,
                                                                   ParticleSet::ParticlePos_t& hf_term,
                                                                   ParticleSet::ParticlePos_t& pulay_terms,
                                                                   ParticleSet::ParticlePos_t& wf_grad)
{
  ParticleSet::ParticleGradient_t wfgradraw_(ions.getTotalNum());
  wfgradraw_           = 0.0;
  RealType localEnergy = 0.0;

  for (int i = 0; i < H.size(); ++i)
    localEnergy += H[i]->evaluateWithIonDerivs(P, ions, psi, hf_term, pulay_terms);

  for (int iat = 0; iat < ions.getTotalNum(); iat++)
  {
    wfgradraw_[iat] = psi.evalGradSource(P, ions, iat);
    convert(wfgradraw_[iat], wf_grad[iat]);
  }
  return localEnergy;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::evaluateIonDerivsDeterministic(ParticleSet& P,
                                                                                ParticleSet& ions,
                                                                                TrialWaveFunction& psi,
                                                                                ParticleSet::ParticlePos_t& hf_term,
                                                                                ParticleSet::ParticlePos_t& pulay_terms,
                                                                                ParticleSet::ParticlePos_t& wf_grad)
{
  ParticleSet::ParticleGradient_t wfgradraw_(ions.getTotalNum());
  wfgradraw_           = 0.0;
  RealType localEnergy = 0.0;

  for (int i = 0; i < H.size(); ++i)
    localEnergy += H[i]->evaluateWithIonDerivsDeterministic(P, ions, psi, hf_term, pulay_terms);

  for (int iat = 0; iat < ions.getTotalNum(); iat++)
  {
    wfgradraw_[iat] = psi.evalGradSource(P, ions, iat);
    convert(wfgradraw_[iat], wf_grad[iat]);
  }
  return localEnergy;
}

QMCHamiltonian::FullPrecRealType QMCHamiltonian::getEnsembleAverage()
{
  FullPrecRealType sum = 0.0;
  for (int i = 0; i < H.size(); i++)
    sum += H[i]->getEnsembleAverage();
  return sum;
}

/** return pointer to the QMCHamtiltonian with the name
 *@param aname the name of Hamiltonian
 *@return the pointer to the named term.
 *
 * If not found, return 0
 */
OperatorBase* QMCHamiltonian::getHamiltonian(const std::string& aname)
{
  for (int i = 0; i < H.size(); ++i)
    if (H[i]->getName() == aname)
      return H[i].get();
  for (int i = 0; i < auxH.size(); ++i)
    if (auxH[i]->getName() == aname)
      return auxH[i].get();
  return nullptr;
}

void QMCHamiltonian::resetTargetParticleSet(ParticleSet& P)
{
  for (int i = 0; i < H.size(); i++)
    H[i]->resetTargetParticleSet(P);
  for (int i = 0; i < auxH.size(); i++)
    auxH[i]->resetTargetParticleSet(P);
}

void QMCHamiltonian::setRandomGenerator(RandomGenerator_t* rng)
{
  for (int i = 0; i < H.size(); i++)
    H[i]->setRandomGenerator(rng);
  for (int i = 0; i < auxH.size(); i++)
    auxH[i]->setRandomGenerator(rng);
  if (nlpp_ptr)
    nlpp_ptr->setRandomGenerator(rng);
}

void QMCHamiltonian::setNonLocalMoves(xmlNodePtr cur)
{
  if (nlpp_ptr != nullptr)
    nlpp_ptr->setNonLocalMoves(cur);
}

void QMCHamiltonian::setNonLocalMoves(const std::string& non_local_move_option,
                                      const double tau,
                                      const double alpha,
                                      const double gamma)
{
  if (nlpp_ptr != nullptr)
    nlpp_ptr->setNonLocalMoves(non_local_move_option, tau, alpha, gamma);
}

int QMCHamiltonian::makeNonLocalMoves(ParticleSet& P)
{
  if (nlpp_ptr == nullptr)
    return 0;
  else
    return nlpp_ptr->makeNonLocalMovesPbyP(P);
}


std::vector<int> QMCHamiltonian::mw_makeNonLocalMoves(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                      const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                      const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& ham_leader = ham_list.getLeader();

  std::vector<int> num_accepts(ham_list.size(), 0);
  if (ham_list.getLeader().nlpp_ptr)
  {
    for (int iw = 0; iw < ham_list.size(); ++iw)
      num_accepts[iw] = ham_list[iw].nlpp_ptr->makeNonLocalMovesPbyP(p_list[iw]);
  }
  return num_accepts;
}

void QMCHamiltonian::createResource(ResourceCollection& collection) const
{
  for (int i = 0; i < H.size(); ++i)
    H[i]->createResource(collection);
}

void QMCHamiltonian::acquireResource(ResourceCollection& collection,
                                     const RefVectorWithLeader<QMCHamiltonian>& ham_list)
{
  auto& ham_leader = ham_list.getLeader();
  for (int i_ham_op = 0; i_ham_op < ham_leader.H.size(); ++i_ham_op)
  {
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));
    ham_leader.H[i_ham_op]->acquireResource(collection, HC_list);
  }
}

void QMCHamiltonian::releaseResource(ResourceCollection& collection,
                                     const RefVectorWithLeader<QMCHamiltonian>& ham_list)
{
  auto& ham_leader = ham_list.getLeader();
  for (int i_ham_op = 0; i_ham_op < ham_leader.H.size(); ++i_ham_op)
  {
    const auto HC_list(extract_HC_list(ham_list, i_ham_op));
    ham_leader.H[i_ham_op]->releaseResource(collection, HC_list);
  }
}

std::unique_ptr<QMCHamiltonian> QMCHamiltonian::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  auto myclone = std::make_unique<QMCHamiltonian>(myName);
  for (int i = 0; i < H.size(); ++i)
    H[i]->add2Hamiltonian(qp, psi, *myclone);
  for (int i = 0; i < auxH.size(); ++i)
    auxH[i]->add2Hamiltonian(qp, psi, *myclone);
  //sync indices
  myclone->resetObservables(myIndex, numCollectables);
  //Hamiltonian needs to make sure qp.Collectables are the same as defined by the original Hamiltonian
  if (numCollectables)
  {
    qp.Collectables.clear();
    qp.Collectables.resize(numCollectables);
  }
  //Assume tau is correct for the Kinetic energy operator and assign to the rest of the clones
  //Return_t tau = H[0]->Tau;
  //myclone->setTau(tau);
  return myclone;
}

#ifdef QMC_CUDA
void QMCHamiltonian::evaluate(MCWalkerConfiguration& W, std::vector<RealType>& energyVector)
{
  ScopedTimer local_timer(ham_timer_);
  auto& walkers = W.WalkerList;
  int nw        = walkers.size();
  if (LocalEnergyVector.size() != nw)
  {
    LocalEnergyVector.resize(nw);
    AuxEnergyVector.resize(nw);
  }
  if (energyVector.size() != nw)
    energyVector.resize(nw);
  for (int i = 0; i < LocalEnergyVector.size(); i++)
    LocalEnergyVector[i] = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    H[i]->addEnergy(W, LocalEnergyVector);
    //H[i]->setObservables(Observables);
  }

  for (int iw = 0; iw < walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[WP::LOCALENERGY] = LocalEnergyVector[iw];
    walkers[iw]->getPropertyBase()[WP::LOCALPOTENTIAL] =
        LocalEnergyVector[iw] - walkers[iw]->getPropertyBase()[WP::NUMPROPERTIES];
  }
  energyVector = LocalEnergyVector;
  // P.PropertyList[WP::WP::LOCALENERGY]=LocalEnergy;
  // P.PropertyList[WP::LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for (int i = 0; i < auxH.size(); ++i)
  {
    auxH[i]->addEnergy(W, AuxEnergyVector);
    //auxH[i]->setObservables(Observables);
  }
}


void QMCHamiltonian::evaluate(MCWalkerConfiguration& W,
                              std::vector<RealType>& energyVector,
                              std::vector<std::vector<NonLocalData>>& Txy)
{
  ScopedTimer local_timer(ham_timer_);
  auto& walkers = W.WalkerList;
  int nw        = walkers.size();
  if (LocalEnergyVector.size() != nw)
  {
    LocalEnergyVector.resize(nw);
    AuxEnergyVector.resize(nw);
  }
  if (energyVector.size() != nw)
    energyVector.resize(nw);
  std::fill(LocalEnergyVector.begin(), LocalEnergyVector.end(), 0.0);
  //for (int i=0; i<LocalEnergyVector.size(); i++)
  //  LocalEnergyVector[i] = 0.0;
  for (int i = 0; i < H.size(); ++i)
  {
    ScopedTimer h_timer(my_timers_[i]);
    H[i]->addEnergy(W, LocalEnergyVector, Txy);
  }

  for (int iw = 0; iw < walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[WP::LOCALENERGY] = LocalEnergyVector[iw];
    walkers[iw]->getPropertyBase()[WP::LOCALPOTENTIAL] =
        LocalEnergyVector[iw] - walkers[iw]->getPropertyBase()[WP::NUMPROPERTIES];
  }
  energyVector = LocalEnergyVector;

  if (auxH.size())
  {
    std::fill(AuxEnergyVector.begin(), AuxEnergyVector.end(), 0.0);
    for (int i = 0; i < auxH.size(); ++i)
      auxH[i]->addEnergy(W, AuxEnergyVector);
  }
}
#endif

RefVectorWithLeader<OperatorBase> QMCHamiltonian::extract_HC_list(const RefVectorWithLeader<QMCHamiltonian>& ham_list,
                                                                  int id)
{
  RefVectorWithLeader<OperatorBase> HC_list(*ham_list.getLeader().H[id]);
  HC_list.reserve(ham_list.size());
  for (QMCHamiltonian& H : ham_list)
    HC_list.push_back(*(H.H[id]));
  return HC_list;
}

} // namespace qmcplusplus

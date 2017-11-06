//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file QMCHamiltonianBase.h
 *@brief Declaration of QMCHamiltonianBase
 */
#ifndef QMCPLUSPLUS_HAMILTONIANBASE_H
#define QMCPLUSPLUS_HAMILTONIANBASE_H

#include <Particle/ParticleSet.h>
#include <OhmmsData/RecordProperty.h>
#include <Utilities/RandomGenerator.h>
#include <QMCHamiltonians/observable_helper.h>
#if !defined(REMOVE_TRACEMANAGER)
#include <Estimators/TraceManager.h>
#endif
#include <QMCWaveFunctions/OrbitalSetTraits.h>
#include <bitset>

namespace qmcplusplus
{
class MCWalkerConfiguration;

/**@defgroup hamiltonian Hamiltonian group
 * @brief QMCHamiltonian and its component, QMCHamiltonianBase
 *
 */
class DistanceTableData;
class TrialWaveFunction;
class QMCHamiltonian;

struct NonLocalData: public QMCTraits
{
  IndexType PID;
  RealType Weight;
  PosType Delta;
  inline NonLocalData():PID(-1),Weight(1.0) {}
  inline NonLocalData(IndexType id, RealType w, const PosType& d)
    :PID(id),Weight(w),Delta(d) {}
};

/** @ingroup hamiltonian
 * @brief An abstract class for Local Energy operators
 *
 * Return_t is defined as RealTye.
 * The types should be checked when using complex wave functions.
 */
struct QMCHamiltonianBase: public QMCTraits
{

  /** type of return value of evaluate
   */
  typedef EstimatorRealType Return_t;
  /** typedef for the serialized buffer
   *
   * PooledData<RealType> is used to serialized an anonymous buffer
   */
  typedef ParticleSet::Buffer_t  BufferType;
  ///typedef for the walker
  typedef ParticleSet::Walker_t  Walker_t;
  ///typedef for the ParticleScalar
  typedef ParticleSet::Scalar_t  ParticleScalar_t;

  ///enum to denote energy domain of operators
  enum energy_domains {kinetic=0,potential,no_energy_domain};

  enum quantum_domains {no_quantum_domain=0,classical,quantum,
                        classical_classical,quantum_classical,quantum_quantum};

  ///quantum_domain of the (particle) operator, default = no_quantum_domain
  quantum_domains quantum_domain;
  ///energy domain of the operator (kinetic/potential), default = no_energy_domain
  energy_domains energy_domain;
  ///enum for UpdateMode
  enum {PRIMARY=0, 
    OPTIMIZABLE=1, 
    RATIOUPDATE=2, 
    PHYSICAL=3, 
    COLLECTABLE=4, 
    NONLOCAL=5,
    VIRTUALMOVES=6
  };
  
  ///set the current update mode
  std::bitset<8> UpdateMode;
  ///starting index of this object
  int myIndex;
  ///number of dependents: to be removed
  int Dependants;
  ///current value
  Return_t Value;
  ///a new value for a proposed move
  Return_t NewValue;
  /// This is used to store the value for force on the source
  /// ParticleSet.  It is accumulated if setComputeForces(true).
  ParticleSet::ParticlePos_t IonForce;
  ///reference to the current walker
  Walker_t* tWalker;
  //Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
  ///name of this object
  std::string myName;
  ///name of dependent object: to be removed
  std::string depName;

#if !defined(REMOVE_TRACEMANAGER)
  ///whether traces are being collected
  TraceRequest request;
  bool streaming_scalars;
  bool streaming_particles;
  bool have_required_traces;
  std::vector<RealType> ValueVector;

  ///array to store sample value
  Array<RealType,1>* value_sample;
#endif

  ///constructor
  QMCHamiltonianBase();

  ///virtual destructor
  virtual ~QMCHamiltonianBase() { }

  ///set energy domain
  void set_energy_domain(energy_domains edomain);

  ///return whether the energy domain is valid
  inline bool energy_domain_valid(energy_domains edomain) const
  {
    return edomain!=no_energy_domain;
  }

  ///return whether the energy domain is valid
  inline bool energy_domain_valid() const
  {
    return energy_domain_valid(energy_domain);
  }

  ///set quantum domain
  void set_quantum_domain(quantum_domains qdomain);

  ///set quantum domain for one-body operator
  void one_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P);

  ///set quantum domain for two-body operator
  void two_body_quantum_domain(const ParticleSet& P1,const ParticleSet& P2);

  ///return whether the quantum domain is valid
  bool quantum_domain_valid(quantum_domains qdomain);

  ///return whether the quantum domain is valid
  inline bool quantum_domain_valid()
  {
    return quantum_domain_valid(quantum_domain);
  }

  inline bool is_classical()
  {
    return quantum_domain==classical;
  }
  inline bool is_quantum()
  {
    return quantum_domain==quantum;
  }
  inline bool is_classical_classical()
  {
    return quantum_domain==classical_classical;
  }
  inline bool is_quantum_classical()
  {
    return quantum_domain==quantum_classical;
  }
  inline bool is_quantum_quantum()
  {
    return quantum_domain==quantum_quantum;
  }

  /** return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  inline bool getMode(int i)
  {
    return UpdateMode[i];
  }

  inline bool isNonLocal() const 
  {
    return UpdateMode[NONLOCAL];
  }

  /** named values to  the property list
   * @param plist RecordNameProperty
   *
   * Previously addObservables but it is renamed and a non-virtial function.
   */
  inline void addValue(PropertySetType& plist)
  {
    if(!UpdateMode[COLLECTABLE])
      myIndex=plist.add(myName.c_str());
  }

  /** named values to  the property list
   * @param plist RecordNameProperty
   * @param collectables Observables that are accumulated by evaluate
   *
   * Default implementaton uses addValue(plist)
   */
  virtual void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addValue(plist);
  }

  /*** add to observable descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation is to register a scalar for this->Value
   */
  virtual void registerObservables(std::vector<observable_helper*>& h5desc
                                   , hid_t gid) const ;

  /*** add to collectables descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerCollectables(std::vector<observable_helper*>& h5desc
                                    , hid_t gid) const
  {}

  /** set the values evaluated by this object to plist
   * @param plist RecordNameProperty
   *
   * Default implementation is to assign Value which is updated
   * by evaluate  function using myIndex.
   */
  virtual void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
  }

  virtual void setParticlePropertyList(PropertySetType& plist
                                       , int offset)
  {
    plist[myIndex+offset]=Value;
  }

  //virtual void setHistories(Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker)
  virtual void setHistories(Walker_t& ThisWalker)
  {
    tWalker = &(ThisWalker);
  }

  /** reset the data with the target ParticleSet
   * @param P new target ParticleSet
   */
  virtual void resetTargetParticleSet(ParticleSet& P) = 0;

  /** Evaluate the local energies of an N-particle configuration
   *@param P input configuration containing N particles
   *@return the value of the Hamiltonian
   */
  virtual Return_t evaluate(ParticleSet& P) = 0;
  virtual Return_t rejectedMove(ParticleSet& P)
  {
    return 0;
  }
  virtual Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy) = 0;
  
  /** evaluate value and derivatives wrt the optimizables
   *
   * Default uses evaluate
   */
  virtual Return_t evaluateValueAndDerivatives(ParticleSet& P,
      const opt_variables_type& optvars,
      const std::vector<RealType>& dlogpsi,
      std::vector<RealType>& dhpsioverpsi)
  {
    return evaluate(P);
  }

  /** update data associated with a particleset
   * @param s source particle set
   *
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.
   */
  virtual void update_source(ParticleSet& s) { }
   
  /** return an average value by collective operation
   */
  virtual Return_t getEnsembleAverage()
  {
    return 0.0;
  }

  /** read the input parameter
   * @param cur xml node for a QMCHamiltonianBase object
   */
  virtual bool put(xmlNodePtr cur)=0;

  /** write about the class */
  virtual bool get(std::ostream& os) const =0;

  virtual QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)=0;

  virtual void setRandomGenerator(RandomGenerator_t* rng)
  {
    //empty
  }

  virtual void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi
                               , QMCHamiltonian& targetH);
  //virtual QMCHamiltonianBase* makeDependants(ParticleSet& qp )
  //{
  //  return 0;
  //}

  virtual void setComputeForces(bool compute)
  {
    // empty
  }

#if !defined(REMOVE_TRACEMANAGER)
  ///make trace quantities available
  inline void contribute_trace_quantities()
  {
    contribute_scalar_quantities();
    contribute_particle_quantities();
  }

  ///checkout trace arrays
  inline void checkout_trace_quantities(TraceManager& tm)
  {
    //derived classes must guard individual checkouts using request info
    checkout_scalar_quantities(tm);
    checkout_particle_quantities(tm);
  }

  ///collect scalar trace data
  inline void collect_scalar_traces()
  {
    //app_log()<<"QMCHamiltonianBase::collect_scalar_traces"<< std::endl;
    collect_scalar_quantities();
  }

  ///delete trace arrays
  inline void delete_trace_quantities()
  {
    delete_scalar_quantities();
    delete_particle_quantities();
    streaming_scalars    = false;
    streaming_particles  = false;
    have_required_traces = false;
    request.reset();
  }

  virtual void contribute_scalar_quantities()
  {
    request.contribute_scalar(myName);
  }

  virtual void checkout_scalar_quantities(TraceManager& tm)
  {
    streaming_scalars = request.streaming_scalar(myName);
    if( streaming_scalars)
      value_sample = tm.checkout_real<1>(myName);
  }

  virtual void collect_scalar_quantities()
  {
    if( streaming_scalars)
      (*value_sample)(0) = Value;
  }

  virtual void delete_scalar_quantities()
  {
    if( streaming_scalars)
      delete value_sample;
  }

  virtual void contribute_particle_quantities() {};
  virtual void checkout_particle_quantities(TraceManager& tm) {};
  virtual void delete_particle_quantities() {};
  virtual void get_required_traces(TraceManager& tm) {};
#endif

  virtual void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy);

  virtual void addEnergy(MCWalkerConfiguration &W,
                         std::vector<RealType> &LocalEnergy,
                         std::vector<std::vector<NonLocalData> > &Txy)
  {
    addEnergy (W, LocalEnergy);
  }

};
}
#endif



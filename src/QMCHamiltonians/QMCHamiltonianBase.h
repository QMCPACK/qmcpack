//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file QMCHamiltonianBase.h
 *@brief Declaration of QMCHamiltonianBase
 */
#ifndef QMCPLUSPLUS_HAMILTONIANBASE_H
#define QMCPLUSPLUS_HAMILTONIANBASE_H

#include <Particle/ParticleSet.h>
#include <OhmmsData/RecordProperty.h>
#include <Utilities/RandomGenerator.h>
#include <QMCHamiltonians/observable_helper.h>
#include <Estimators/TraceManager.h>
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
  typedef RealType Return_t;
  /** typedef for the serialized buffer
   *
   * PooledData<RealType> is used to serialized an anonymous buffer
   */
  typedef ParticleSet::Buffer_t  BufferType;
  ///typedef for the walker
  typedef ParticleSet::Walker_t  Walker_t;

  ///enum for UpdateMode
  enum {PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL, COLLECTABLE};
  ///set the current update mode
  bitset<8> UpdateMode;
  ///starting index of this object
  int myIndex;
  ///number of dependents: to be removed
  int Dependants;
  ///current value
  RealType Value;
  ///a new value for a proposed move
  RealType NewValue;
  /// This is used to store the value for force on the source
  /// ParticleSet.  It is accumulated if setComputeForces(true).
  ParticleSet::ParticlePos_t IonForce;
  ///reference to the current walker
  Walker_t* tWalker;
  //Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
  ///name of this object
  string myName;
  ///name of dependent object: to be removed
  string depName;
  ///whether traces are being collected
  TraceRequest trace_request;
  bool tracing;
  bool tracing_scalar_quantities;
  bool tracing_particle_quantities;
  bool have_required_traces;
  ///array to store sample value
  Array<RealType,1>* value_sample;

  ///constructor
  QMCHamiltonianBase()
    :myIndex(-1),Value(0.0),Dependants(0),tWalker(0)
  {
    tracing = false;
    tracing_scalar_quantities = false;
    tracing_particle_quantities = false;
    have_required_traces = false;
    UpdateMode.set(PRIMARY,1);
  }

  ///virtual destructor
  virtual ~QMCHamiltonianBase() { }

  /** return the mode i
   * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
   */
  inline bool getMode(int i)
  {
    return UpdateMode[i];
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
  virtual void registerObservables(vector<observable_helper*>& h5desc
                                   , hid_t gid) const ;

  /*** add to collectables descriptor for hdf5
   * @param h5desc contains a set of hdf5 descriptors for a scalar observable
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerCollectables(vector<observable_helper*>& h5desc
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
  virtual Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) = 0;

  /** update data associated with a particleset
   * @param s source particle set
   *
   * Default implementation does nothing. Only A-A interactions for s needs to implement its own method.
   */
  virtual void update_source(ParticleSet& s) { }
   

  /*@{
   * @brief Functions to handle particle-by-particle move
   *
   * Default implementations use evaluate.
   */
  virtual Return_t registerData(ParticleSet& P, BufferType& buffer)
  {
    return evaluate(P);
  }
  virtual Return_t updateBuffer(ParticleSet& P, BufferType& buf)
  {
    return evaluate(P);
  }
  virtual void copyFromBuffer(ParticleSet& P, BufferType& buf)
  {
    Value=evaluate(P);
  }
  virtual void copyToBuffer(ParticleSet& P, BufferType& buf)
  {
  }
  virtual Return_t evaluatePbyP(ParticleSet& P, int active)
  {
    APP_ABORT(myName + " missing evaluatePbyP");
    return NewValue;
  }
  virtual void acceptMove(int active)
  {
    Value=NewValue;
  }
  virtual void rejectMove(int active)
  {
  }
  /*@}*/

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



  ///request that traces be made available
  inline void request_traces(TraceManager& tm)
  {
    tm.add_trace_request(trace_request);
  }


  ///checkout trace arrays
  inline void initialize_traces(TraceManager& tm)
  {
    TraceRequest& traces_requested = tm.get_trace_request(myName);
    tracing_scalar_quantities   = traces_requested.scalars;
    tracing_particle_quantities = traces_requested.particles;
    tracing = tracing_scalar_quantities || tracing_particle_quantities;
    if(tracing_scalar_quantities)
      checkout_scalar_arrays(tm);
    if(tracing_particle_quantities)
      checkout_particle_arrays(tm);
  }

  ///collect scalar trace data
  inline void collect_scalar_traces()
  {
    //app_log()<<"QMCHamiltonianBase::collect_scalar_traces"<<endl;
    if(tracing_scalar_quantities)
      collect_scalar_samples();
  }

  ///delete trace arrays
  inline void finalize_traces()
  {
    if(tracing_scalar_quantities)
      delete_scalar_arrays();
    if(tracing_particle_quantities)
      delete_particle_arrays();
    tracing_scalar_quantities = false;
    tracing_particle_quantities = false;
    have_required_traces = false;
  }

  virtual void checkout_scalar_arrays(TraceManager& tm)
  {
    value_sample = tm.checkout_real<1>(myName);
  }

  virtual void collect_scalar_samples()
  {
    (*value_sample)(0) = Value;
  }

  virtual void delete_scalar_arrays()
  {
    delete value_sample;
  }

  virtual void checkout_particle_arrays(TraceManager& tm) {};
  virtual void delete_particle_arrays() {};
  virtual void get_required_traces(TraceManager& tm) {};


  ////////////////////////////////////
  // Vectorized evaluation on GPUs  //
  ////////////////////////////////////
#ifdef QMC_CUDA
  vector<RealType> ValueVector;

  virtual void addEnergy(MCWalkerConfiguration &W,
                         vector<RealType> &LocalEnergy)
  {
    app_error() << "Need specialization for " << myName
                << "::addEnergy(MCWalkerConfiguration &W).\n";
  }

  virtual void addEnergy(MCWalkerConfiguration &W,
                         vector<RealType> &LocalEnergy,
                         vector<vector<NonLocalData> > &Txy)
  {
    addEnergy (W, LocalEnergy);
  }
#endif


};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/


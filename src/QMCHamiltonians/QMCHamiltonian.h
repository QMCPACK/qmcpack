//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file QMCHamiltonian.h
 *@brief Declaration of QMCHamiltonian
 */
#ifndef QMCPLUSPLUS_HAMILTONIAN_H
#define QMCPLUSPLUS_HAMILTONIAN_H
#include "Configuration.h"
#include <QMCHamiltonians/OperatorBase.h>
#include "QMCHamiltonians/NonLocalECPotential.h"
#if !defined(REMOVE_TRACEMANAGER)
#include <Estimators/TraceManager.h>
#endif
#include <QMCWaveFunctions/OrbitalSetTraits.h>
namespace qmcplusplus
{
class MCWalkerConfiguration;
class NewTimer;
class HamiltonianFactory;

/**  Collection of Local Energy Operators
 *
 * Note that QMCHamiltonian is not derived from QMCHmailtonianBase.
 */
class QMCHamiltonian
{
  friend class HamiltonianFactory;

public:
  typedef OperatorBase::RealType RealType;
  typedef OperatorBase::ValueType ValueType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  typedef OperatorBase::PropertySetType PropertySetType;
  typedef OperatorBase::BufferType BufferType;
  typedef OperatorBase::Walker_t Walker_t;

  enum
  {
    DIM = OHMMS_DIM
  };

  ///constructor
  QMCHamiltonian();

  ///destructor
  ~QMCHamiltonian();

  ///add an operator
  void addOperator(OperatorBase* h, const std::string& aname, bool physical = true);

  ///record the name-type pair of an operator
  void addOperatorType(const std::string& name, const std::string& type);

  ///return type of named H element or fail
  const std::string& getOperatorType(const std::string& name);

  ///return the number of Hamiltonians
  inline int size() const { return H.size(); }

  ///return the total number of Hamiltonians (physical + aux)
  inline int total_size() const { return H.size() + auxH.size(); }

  /** return OperatorBase with the name aname
   * @param aname name of a OperatorBase
   * @return 0 if aname is not found.
   */
  OperatorBase* getHamiltonian(const std::string& aname);

  /** return i-th OperatorBase
   * @param i index of the OperatorBase
   * @return H[i]
   */
  OperatorBase* getHamiltonian(int i) { return H[i]; }

#if !defined(REMOVE_TRACEMANAGER)
  ///initialize trace data
  void initialize_traces(TraceManager& tm, ParticleSet& P);

  // ///collect scalar trace data
  //void collect_scalar_traces();

  ///collect walker trace data
  void collect_walker_traces(Walker_t& walker, int step);

  ///finalize trace data
  void finalize_traces();
#endif

  /**
   * \defgroup Functions to get/put observables
   */
  /*@{*/
  /** add each term to the PropertyList for averages
   * @param plist a set of properties to which this Hamiltonian add the observables
   */
  //void addObservables(PropertySetType& plist);

  /** add each term to P.PropertyList and P.mcObservables
   * @param P particle set to which observables are added
   * @return the number of observables
   */
  int addObservables(ParticleSet& P);

  /** register obsevables so that their averages can be dumped to hdf5
   * @param h5desc has observable_helper* for each h5 group
   * @param gid h5 group id to which the observable groups are added.
   */
  void registerObservables(std::vector<observable_helper*>& h5desc, hid_t gid) const;
  /** register collectables so that their averages can be dumped to hdf5
   * @param h5desc has observable_helper* for each h5 group
   * @param gid h5 group id to which the observable groups are added.
   *
   * Add observable_helper information for the data stored in ParticleSet::mcObservables.
   */
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const;
  ///retrun the starting index
  inline int startIndex() const { return myIndex; }
  ///return the size of observables
  inline int sizeOfObservables() const { return Observables.size(); }
  ///return the size of collectables
  inline int sizeOfCollectables() const { return numCollectables; }
  ///return the value of the i-th observable
  inline RealType getObservable(int i) const { return Observables.Values[i]; }
  ///return the value of the observable with a set name if it exists
  inline int getObservable(std::string Oname) const
  {
    int rtval(-1);
    for (int io = 0; io < Observables.size(); io++)
    {
      if (Observables.Names[io] == Oname)
        return io;
    }
    return rtval;
  }
  ///return the name of the i-th observable
  inline std::string getObservableName(int i) const { return Observables.Names[i]; }
  ///save the values of Hamiltonian elements to the Properties
  template<class IT>
  inline void saveProperty(IT first)
  {
    first[LOCALPOTENTIAL] = LocalEnergy - KineticEnergy;
    copy(Observables.begin(), Observables.end(), first + myIndex);
  }
  /*@}*/

  template<class IT>
  inline void setProperty(IT first)
  {
    //       LocalEnergy=first[LOCALENERGY];
    //       KineticEnergy=LocalEnergy-first[LOCALPOTENTIAL];
    copy(first + myIndex, first + myIndex + Observables.size(), Observables.begin());
  }

  void update_source(ParticleSet& s);

  ////return the LocalEnergy \f$=\sum_i H^{qmc}_{i}\f$
  inline FullPrecRealType getLocalEnergy() { return LocalEnergy; }
  ////return the LocalPotential \f$=\sum_i H^{qmc}_{i} - KE\f$
  inline FullPrecRealType getLocalPotential() { return LocalEnergy - KineticEnergy; }
  void auxHevaluate(ParticleSet& P);
  void auxHevaluate(ParticleSet& P, Walker_t& ThisWalker);
  void auxHevaluate(ParticleSet& P, Walker_t& ThisWalker, bool do_properties, bool do_collectables);
  void rejectedMove(ParticleSet& P, Walker_t& ThisWalker);
  ///** set Tau for each Hamiltonian
  // */
  //inline void setTau(RealType tau)
  //{
  //  for(int i=0; i< H.size();i++)H[i]->setTau(tau);
  //}

  /** set PRIMARY bit of all the components
   */
  inline void setPrimary(bool primary)
  {
    for (int i = 0; i < H.size(); i++)
      H[i]->UpdateMode.set(OperatorBase::PRIMARY, primary);
  }

  /////Set Tau inside each of the Hamiltonian elements
  //void setTau(Return_t tau)
  //{
  //  for(int i=0; i<H.size(); i++) H[i]->setTau(tau);
  //  for(int i=0; i<auxH.size(); i++) auxH[i]->setTau(tau);
  //}

  ///** return if WaveFunction Ratio needs to be evaluated
  // *
  // * This is added to handle orbital-dependent OperatorBase during
  // * orbital optimizations.
  // */
  //inline bool needRatio() {
  //  bool dependOnOrbital=false;
  //  for(int i=0; i< H.size();i++)
  //    if(H[i]->UpdateMode[OperatorBase::RATIOUPDATE]) dependOnOrbital=true;
  //  return dependOnOrbital;
  //}

  /** evaluate Local Energy
   * @param P ParticleSet
   * @return the local energy
   *
   * P.R, P.G and P.L are used to evaluate the LocalEnergy.
   */
  FullPrecRealType evaluate(ParticleSet& P);

  /** batched version of evaluate for LocalEnergy 
   *
   *  Encapsulation is ignored for H_list hamiltonians method uses its status as QMCHamiltonian to break encapsulation.
   *  ParticleSet is also updated.
   *  Bugs could easily be created by accessing this scope.
   *  This should be set to static and fixed.
   */
  static std::vector<QMCHamiltonian::FullPrecRealType> flex_evaluate(const RefVector<QMCHamiltonian>& H_list,
                                                                     const RefVector<ParticleSet>& P_list);

  /** evaluate Local energy with Toperators updated.
   * @param P ParticleSEt
   * @return Local energy
   */
  FullPrecRealType evaluateWithToperator(ParticleSet& P);

  /** batched version of evaluate Local energy with Toperators updated.
   */
  static std::vector<QMCHamiltonian::FullPrecRealType> flex_evaluateWithToperator(RefVector<QMCHamiltonian>& H_list,
                                                                                  RefVector<ParticleSet>& P_list);


  /** evaluate energy and derivatives wrt to the variables
   * @param P ParticleSet
   * @param optvars current optimiable variables
   * @param dlogpsi \f$\partial \ln \Psi({\bf R})/\partial \alpha \f$
   * @param dhpsioverpsi \f$\partial(\hat{h}\Psi({\bf R})/\Psi({\bf R})) /\partial \alpha \f$
   * @param compute_deriv if true, compute dhpsioverpsi of the non-local potential component
   */
  FullPrecRealType evaluateValueAndDerivatives(ParticleSet& P,
                                               const opt_variables_type& optvars,
                                               std::vector<ValueType>& dlogpsi,
                                               std::vector<ValueType>& dhpsioverpsi,
                                               bool compute_deriv);

  /** evaluate local energy and derivatives w.r.t ionic coordinates.  
  * @param P target particle set (electrons)
  * @param ions source particle set (ions)
  * @param psi Trial wave function
  * @param hf_terms  Re [(dH)Psi]/Psi
  * @param pulay_terms Re [(H-E_L)dPsi]/Psi 
  * @param wf_grad  Re (dPsi/Psi)
  * @return Local Energy.
  */
  FullPrecRealType evaluateIonDerivs(ParticleSet& P,
                                     ParticleSet& ions,
                                     TrialWaveFunction& psi,
                                     ParticleSet::ParticlePos_t& hf_terms,
                                     ParticleSet::ParticlePos_t& pulay_terms,
                                     ParticleSet::ParticlePos_t& wf_grad);
  /** set non local moves options
   * @param cur the xml input
   */
  void setNonLocalMoves(xmlNodePtr cur)
  {
    if (nlpp_ptr != nullptr)
      nlpp_ptr->setNonLocalMoves(cur);
  }

  void setNonLocalMoves(const std::string& non_local_move_option,
                                        const double tau,
                                        const double alpha,
                                        const double gamma)
  {
    if (nlpp_ptr != nullptr)
      nlpp_ptr->setNonLocalMoves(non_local_move_option,
                                        tau,
                                        alpha,
                                        gamma);
  }
  /** make non local moves
   * @param P particle set
   * @return the number of accepted moves
   */
  int makeNonLocalMoves(ParticleSet& P)
  {
    if (nlpp_ptr == nullptr)
      return 0;
    else
      return nlpp_ptr->makeNonLocalMovesPbyP(P);
  }

  static std::vector<int> flex_makeNonLocalMoves(RefVector<QMCHamiltonian>& h_list, RefVector<ParticleSet>& p_list);
  /** evaluate energy 
   * @param P quantum particleset
   * @param free_nlpp if true, non-local PP is a variable
   * @return KE + NonLocal potential
   */
  FullPrecRealType evaluateVariableEnergy(ParticleSet& P, bool free_nlpp);

  /** return an average value of the LocalEnergy
   *
   * Introduced to get a collective value
   */
  FullPrecRealType getEnsembleAverage();

  void resetTargetParticleSet(ParticleSet& P);

  /** By mistake, QMCHamiltonian::getName(int i) is used
   * and this is in conflict with the declaration of OhmmsElementBase.
   * For the moment, QMCHamiltonian is not inherited from OhmmsElementBase.
   */
  void setName(const std::string& aname) { myName = aname; }


  std::string getName() const { return myName; }

  bool get(std::ostream& os) const;

  RealType get_LocalEnergy() const { return LocalEnergy; }

  void setRandomGenerator(RandomGenerator_t* rng);

  static void updateNonKinetic(OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset);
  static void updateKinetic(OperatorBase& op, QMCHamiltonian& ham, ParticleSet& pset);
  
  /** return a clone */
  QMCHamiltonian* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

#ifdef QMC_CUDA
  ////////////////////////////////////////////
  // Vectorized evaluation routines for GPU //
  ////////////////////////////////////////////
  void evaluate(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy);
  void evaluate(MCWalkerConfiguration& W,
                std::vector<RealType>& energyVector,
                std::vector<std::vector<NonLocalData>>& Txy);

private:
  /////////////////////
  // Vectorized data //
  /////////////////////
  std::vector<QMCHamiltonian::FullPrecRealType> LocalEnergyVector, KineticEnergyVector, AuxEnergyVector;
#endif

private:
  ///starting index
  int myIndex;
  ///starting index
  int numCollectables;
  ///Current Local Energy
  FullPrecRealType LocalEnergy;
  ///Current Kinetic Energy
  FullPrecRealType KineticEnergy;
  ///Current Local Energy for the proposed move
  FullPrecRealType NewLocalEnergy;
  ///getName is in the way
  std::string myName;
  ///vector of Hamiltonians
  std::vector<OperatorBase*> H;
  ///pointer to NonLocalECP
  NonLocalECPotential* nlpp_ptr;
  ///vector of Hamiltonians
  std::vector<OperatorBase*> auxH;
  ///timers
  std::vector<NewTimer*> myTimers;
  ///types of component operators
  std::map<std::string, std::string> operator_types;
  ///data
  PropertySetType Observables;
  /** reset Observables and counters
   * @param start starting index within PropertyList
   * @param ncollects number of collectables
   */
  void resetObservables(int start, int ncollects);

  // helper function for extracting a list of Hamiltonian components from a list of QMCHamiltonian::H.
  static RefVector<OperatorBase> extract_HC_list(const RefVector<QMCHamiltonian>& H_list, int id);

#if !defined(REMOVE_TRACEMANAGER)
  ///traces variables
  TraceRequest request;
  bool streaming_position;
  Array<TraceInt, 1>* id_sample;
  Array<TraceInt, 1>* pid_sample;
  Array<TraceInt, 1>* step_sample;
  Array<TraceInt, 1>* gen_sample;
  Array<TraceInt, 1>* age_sample;
  Array<TraceInt, 1>* mult_sample;
  Array<TraceReal, 1>* weight_sample;
  Array<TraceReal, 2>* position_sample;
#endif
};
} // namespace qmcplusplus
#endif

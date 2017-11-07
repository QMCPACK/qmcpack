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
#include <QMCHamiltonians/QMCHamiltonianBase.h>
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

  typedef QMCHamiltonianBase::RealType  RealType;
  typedef QMCHamiltonianBase::ValueType ValueType;
  typedef QMCHamiltonianBase::Return_t  Return_t;
  typedef QMCHamiltonianBase::PropertySetType  PropertySetType;
  typedef QMCHamiltonianBase::BufferType  BufferType;
  typedef QMCHamiltonianBase::Walker_t  Walker_t;

  enum {DIM=OHMMS_DIM};

  ///constructor
  QMCHamiltonian();

  ///destructor
  ~QMCHamiltonian();

  ///add an operator
  void addOperator(QMCHamiltonianBase* h, const std::string& aname, bool physical=true);

  ///record the name-type pair of an operator
  void addOperatorType(const std::string& name,const std::string& type);

  ///return type of named H element or fail
  const std::string& getOperatorType(const std::string& name);

  ///return the number of Hamiltonians
  inline int size() const
  {
    return H.size();
  }

  ///return the total number of Hamiltonians (physical + aux)
  inline int total_size() const
  {
    return H.size()+auxH.size();
  }

  /** return QMCHamiltonianBase with the name aname
   * @param aname name of a QMCHamiltonianBase
   * @return 0 if aname is not found.
   */
  QMCHamiltonianBase* getHamiltonian(const std::string& aname);

  /** return i-th QMCHamiltonianBase
   * @param i index of the QMCHamiltonianBase
   * @return H[i]
   */
  QMCHamiltonianBase* getHamiltonian(int i)
  {
    return H[i];
  }

#if !defined(REMOVE_TRACEMANAGER)
  ///initialize trace data
  void initialize_traces(TraceManager& tm,ParticleSet& P);

  // ///collect scalar trace data
  //void collect_scalar_traces();

  ///collect walker trace data
  void collect_walker_traces(Walker_t& walker,int step);

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
  void registerObservables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;
  /** register collectables so that their averages can be dumped to hdf5
   * @param h5desc has observable_helper* for each h5 group
   * @param gid h5 group id to which the observable groups are added.
   *
   * Add observable_helper information for the data stored in ParticleSet::mcObservables.
   */
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;
  ///retrun the starting index
  inline int startIndex() const
  {
    return myIndex;
  }
  ///return the size of observables
  inline int sizeOfObservables() const
  {
    return Observables.size();
  }
  ///return the size of collectables
  inline int sizeOfCollectables() const
  {
    return numCollectables;
  }
  ///return the value of the i-th observable
  inline RealType getObservable(int i) const
  {
    return Observables.Values[i];
  }
  ///return the value of the observable with a set name if it exists
  inline int getObservable( std::string Oname) const
  {
    int rtval(-1);
    for(int io=0; io<Observables.size(); io++)
    {
      if (Observables.Names[io]==Oname)
        return io;
    }
    return rtval;
  }
  ///return the name of the i-th observable
  inline std::string getObservableName(int i) const
  {
    return Observables.Names[i];
  }
  ///save the values of Hamiltonian elements to the Properties
  template<class IT>
  inline
  void saveProperty(IT first)
  {
    first[LOCALPOTENTIAL]= LocalEnergy-KineticEnergy;
    copy(Observables.begin(),Observables.end(),first+myIndex);
  }
  /*@}*/

  template<class IT>
  inline
  void setProperty(IT first)
  {
//       LocalEnergy=first[LOCALENERGY];
//       KineticEnergy=LocalEnergy-first[LOCALPOTENTIAL];
    copy(first+myIndex,first+myIndex+Observables.size(),Observables.begin());
  }

  void update_source(ParticleSet& s);

  ////return the LocalEnergy \f$=\sum_i H^{qmc}_{i}\f$
  inline Return_t getLocalEnergy()
  {
    return LocalEnergy;
  }
  ////return the LocalPotential \f$=\sum_i H^{qmc}_{i} - KE\f$
  inline Return_t getLocalPotential()
  {
    return LocalEnergy-KineticEnergy;
  }
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
    for(int i=0; i< H.size(); i++)
      H[i]->UpdateMode.set(QMCHamiltonianBase::PRIMARY,primary);
  }

  /////Set Tau inside each of the Hamiltonian elements
  //void setTau(Return_t tau)
  //{
  //  for(int i=0; i<H.size(); i++) H[i]->setTau(tau);
  //  for(int i=0; i<auxH.size(); i++) auxH[i]->setTau(tau);
  //}

  ///** return if WaveFunction Ratio needs to be evaluated
  // *
  // * This is added to handle orbital-dependent QMCHamiltonianBase during
  // * orbital optimizations.
  // */
  //inline bool needRatio() {
  //  bool dependOnOrbital=false;
  //  for(int i=0; i< H.size();i++)
  //    if(H[i]->UpdateMode[QMCHamiltonianBase::RATIOUPDATE]) dependOnOrbital=true;
  //  return dependOnOrbital;
  //}

  /** evaluate Local Energy
   * @param P ParticleSet
   * @return the local energy
   *
   * P.R, P.G and P.L are used to evaluate the LocalEnergy.
   */
  Return_t evaluate(ParticleSet& P);

  /** evaluate Local and NonLocal energies
   * @param P ParticleSEt
   * @param Txy transition matrix of nonlocal Hamiltonians
   * @return Local energy
   */
  Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy);
  
  /** evaluate energy and derivatives wrt to the variables
   * @param P ParticleSet
   * @param optvars current optimiable variables
   * @param dlogpsi \f$\partial \ln \Psi({\bf R})/\partial \alpha \f$
   * @param dhpsioverpsi \f$\partial(\hat{h}\Psi({\bf R})/\Psi({\bf R})) /\partial \alpha \f$
   * @param compute_deriv if true, compute dhpsioverpsi of the non-local potential component
   */
  RealType evaluateValueAndDerivatives(ParticleSet& P,
      const opt_variables_type& optvars,
      std::vector<RealType>& dlogpsi,
      std::vector<RealType>& dhpsioverpsi,
      bool compute_deriv);
  
  /** evaluate energy 
   * @param P quantum particleset
   * @param free_nlpp if true, non-local PP is a variable
   * @return KE + NonLocal potential
   */
  RealType evaluateVariableEnergy(ParticleSet& P, bool free_nlpp);

  /** return an average value of the LocalEnergy
   *
   * Introduced to get a collective value
   */
  Return_t getEnsembleAverage();

  void resetTargetParticleSet(ParticleSet& P);

  /** By mistake, QMCHamiltonian::getName(int i) is used
   * and this is in conflict with the declaration of OhmmsElementBase.
   * For the moment, QMCHamiltonian is not inherited from OhmmsElementBase.
   */
  void setName(const std::string& aname)
  {
    myName=aname;
  }


  std::string getName() const
  {
    return myName;
  }

  bool get(std::ostream& os) const;

  void setRandomGenerator(RandomGenerator_t* rng);

  /** return a clone */
  QMCHamiltonian* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  ////////////////////////////////////////////
  // Vectorized evaluation routines for GPU //
  ////////////////////////////////////////////
  void evaluate (MCWalkerConfiguration &W,  std::vector<RealType> &LocalEnergy);
  void evaluate(MCWalkerConfiguration &W, std::vector<RealType> &energyVector,
                std::vector<std::vector<NonLocalData> > &Txy);

#ifdef QMC_CUDA
private:
  /////////////////////
  // Vectorized data //
  /////////////////////
  std::vector<Return_t> LocalEnergyVector, KineticEnergyVector,
         AuxEnergyVector;
#endif

private:
  ///starting index
  int myIndex;
  ///starting index
  int numCollectables;
  ///enable virtual moves 
  bool EnableVirtualMoves;
  ///Current Local Energy
  Return_t LocalEnergy;
  ///Current Kinetic Energy
  Return_t KineticEnergy;
  ///Current Local Energy for the proposed move
  Return_t NewLocalEnergy;
  ///getName is in the way
  std::string myName;
  ///vector of Hamiltonians
  std::vector<QMCHamiltonianBase*> H;
  ///vector of Hamiltonians
  std::vector<QMCHamiltonianBase*> auxH;
  ///timers
  std::vector<NewTimer*> myTimers;
  ///types of component operators
  std::map<std::string,std::string> operator_types;
  ///data
  PropertySetType Observables;
  /** reset Observables and counters
   * @param start starting index within PropertyList
   * @param ncollects number of collectables
   */
  void resetObservables(int start, int ncollects);


#if !defined(REMOVE_TRACEMANAGER)
  ///traces variables
  TraceRequest request;
  bool streaming_position;
  Array<TraceInt,1>*  id_sample;
  Array<TraceInt,1>*  pid_sample;
  Array<TraceInt,1>*  step_sample;
  Array<TraceInt,1>*  gen_sample;
  Array<TraceInt,1>*  age_sample;
  Array<TraceInt,1>*  mult_sample;
  Array<TraceReal,1>* weight_sample;
  Array<TraceReal,2>* position_sample;
#endif
};
}
#endif



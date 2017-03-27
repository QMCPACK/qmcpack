//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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
/** @file SkAllEstimator.h
 * @brief Declare SkAllEstimator
 */
#ifndef QMCPLUSPLUS_SK_ALL_ESTIMATOR_H
#define QMCPLUSPLUS_SK_ALL_ESTIMATOR_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <vector>
namespace qmcplusplus
{

/** SkAllEstimator evaluate the structure factor of the target particleset
 *
 * <estimator name="sk" type="sk" debug="no"/>
 */
class SkAllEstimator: public QMCHamiltonianBase
{
public:

  SkAllEstimator(ParticleSet& ions, ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void evaluateIonIon();
  
  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

protected:
//  ParticleSet *sourcePtcl;
  ParticleSet *elns;
  ParticleSet *ions;
  /** number of species */
  int NumSpecies;
  int NumeSpecies;
  int NumIonSpecies; 
 /** number of kpoints */
  unsigned int NumK;
  /** number of kshells */
  int MaxKshell;
  /** normalization factor */
  RealType OneOverN;
  /** kshell counters */
  std::vector<int> Kshell;
  /** instantaneous structure factor  */
  std::vector<RealType> Kmag;
  /** 1.0/degenracy for a ksell */
  std::vector<RealType> OneOverDnk;
  /** \f$rho_k = \sum_{\alpha} \rho_k^{\alpha} \f$ for species index \f$\alpha\f$ */
#if defined(USE_REAL_STRUCT_FACTOR)
  Vector<RealType> RhokTot_r,RhokTot_i;
#else
  Vector<ComplexType> RhokTot;
#endif
  Vector<RealType> values;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();

  bool hdf5_out;
};

}
#endif


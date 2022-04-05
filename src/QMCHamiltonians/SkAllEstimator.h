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
#include "QMCHamiltonians/OperatorBase.h"
#include <vector>
namespace qmcplusplus
{
/** SkAllEstimator evaluate the structure factor of the target particleset
 *
 * <estimator name="sk" type="sk" debug="no"/>
 */
class SkAllEstimator : public OperatorBase
{
public:
  SkAllEstimator(ParticleSet& ions, ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void evaluateIonIon();

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

protected:
  //  ParticleSet *sourcePtcl;
  ParticleSet* elns;
  ParticleSet* ions;
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
  /** 1.0/degenracy for a kshell */
  std::vector<RealType> OneOverDnk;
  /** \f$rho_k = \sum_{\alpha} \rho_k^{\alpha} \f$ for species index \f$\alpha\f$ */
  Vector<RealType> RhokTot_r, RhokTot_i;
  Vector<RealType> values;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();

  bool hdf5_out;
};

} // namespace qmcplusplus
#endif

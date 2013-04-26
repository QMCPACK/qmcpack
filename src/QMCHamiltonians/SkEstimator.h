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
/** @file SkEstimator.h
 * @brief Declare SkEstimator
 */
#ifndef QMCPLUSPLUS_SK_ESTIMATOR_H
#define QMCPLUSPLUS_SK_ESTIMATOR_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
namespace qmcplusplus
{

/** SkEstimator evaluate the structure factor of the target particleset
 *
 * <estimator name="sk" type="sk" debug="no"/>
 */
class SkEstimator: public QMCHamiltonianBase
{
public:

  SkEstimator(ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

protected:
  ParticleSet *sourcePtcl;
  /** number of species */
  int NumSpecies;
  /** number of kpoints */
  int NumK;
  /** number of kshells */
  int MaxKshell;
  /** normalization factor */
  RealType OneOverN;
  /** kshell counters */
  vector<int> Kshell;
  /** instantaneous structure factor  */
  vector<RealType> Kmag;
  /** 1.0/degenracy for a ksell */
  vector<RealType> OneOverDnk;
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

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

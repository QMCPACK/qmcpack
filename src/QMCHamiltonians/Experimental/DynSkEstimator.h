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
/** @file DynSkEstimator.h
 * @brief Declare DynSkEstimator
 */
#ifndef QMCPLUSPLUS_DYN_SK_ESTIMATOR_H
#define QMCPLUSPLUS_DYN_SK_ESTIMATOR_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
namespace qmcplusplus
{

/** DynSkEstimator evaluate the dynamic structure factor of the target particleset
 *
 * <estimator name="dsk" type="dynsk" debug="no"/>
 */
class DynSkEstimator: public QMCHamiltonianBase
{
public:

  DynSkEstimator(ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  Return_t calculate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t rejectedMove(ParticleSet& P)
  {
    int lastindex = tWalker->PHindex[pindx]-1;
    if (lastindex<0)
      lastindex += NumT;
    for (int i=0; i<NumK*2; i++)
      tWalker->addPropertyHistoryPoint(pindx+i,  tWalker->PropertyHistory[pindx+i][lastindex]);
    calculate(P);
    return 0.0;
  }

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool putSpecial(xmlNodePtr cur, ParticleSet& P);
  bool put(xmlNodePtr cur)
  {
    return true;
  }
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
  Vector<ComplexType> RhokTot;
  Vector<RealType> values;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();

  bool hdf5_out;

//     storing old rhoK
  /**     starting index in walker of stored rho values, last index(circular queue **/
  int pindx,cindx;
  /**    length of stored rhoK **/
  int NumT,MinT;

};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

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
#ifndef QMCPLUSPLUS_PAIRCOOR_HAMILTONIAN_H
#define QMCPLUSPLUS_PAIRCOOR_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsMatrix.h>

namespace qmcplusplus
{

/** gofr estimator
 *
 * Compute pair correlation function for the target particle set and optionally any source particles
 */
class PairCorrEstimator: public QMCHamiltonianBase
{
public:

  /** constructor
   * @param elns target particle set
   * @param sources list of source particle sets
   *
   * Use the elns.DistTables to evaluate the pair correlation functions.
   */
  PairCorrEstimator(ParticleSet& elns, string& sources);

  void resetTargetParticleSet(ParticleSet& P);

  /* evaluate the pair correlation functions */
  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(vector<observable_helper*>& h5list, hid_t gid) const;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void report();

private:
  ///number of bins
  int NumBins;
  /// maximum distance
  RealType Dmax;
  /// bin size
  RealType Delta;
  /// one of bin size
  RealType DeltaInv;
  /// volume of the cell
  RealType Volume;
  ///save pair indices
  vector<int> pair_ids;
  /// table indexs for other type
  vector<int> other_ids;
  /// offset of the gofr's associated with others_id
  vector<int> other_offsets;
  /////save source indices
  //vector<int> source_ids;
  ///normalization factor
  Matrix<RealType> norm_factor;
  int num_species,N_e;
  std::vector<RealType> n_vec;
  /////data
  //Matrix<RealType> gof_r;
  ///prefix of each gof_r
  vector<string> gof_r_prefix;
  /** resize the internal data
   * @param nbins number of bins for the historgram
   */
  void resize(int nbins);
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  PairCorrEstimator(ParticleSet& elns, std::string& sources);

  void resetTargetParticleSet(ParticleSet& P);

  /* evaluate the pair correlation functions */
  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void registerCollectables(std::vector<observable_helper*>& h5list, hid_t gid) const;
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
  std::vector<int> pair_ids;
  /// table indexs for other type
  std::vector<int> other_ids;
  /// offset of the gofr's associated with others_id
  std::vector<int> other_offsets;
  /////save source indices
  //vector<int> source_ids;
  ///normalization factor
  Matrix<RealType> norm_factor;
  int num_species,N_e;
  std::vector<RealType> n_vec;
  /////data
  //Matrix<RealType> gof_r;
  ///prefix of each gof_r
  std::vector<std::string> gof_r_prefix;
  /** resize the internal data
   * @param nbins number of bins for the historgram
   */
  void resize(int nbins);
};

}
#endif


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
#include "QMCHamiltonians/OperatorBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{
/** gofr estimator
 *
 * Compute pair correlation function for the target particle set and optionally any source particles
 */
class PairCorrEstimator : public OperatorBase
{
public:
  /** constructor
   * @param elns target particle set
   * @param sources list of source particle sets
   *
   * Use the elns.DistTables to evaluate the pair correlation functions.
   */
  PairCorrEstimator(ParticleSet& elns, std::string& sources);

  void resetTargetParticleSet(ParticleSet& P) override;

  /* evaluate the pair correlation functions */
  Return_t evaluate(ParticleSet& P) override;

  /// generate the unique pair id from the group ids of particle i and j and the number of species
  static int gen_pair_id(const int ig, const int jg, const int ns);

  void addObservables(PropertySetType& plist) {}
  void addObservables(PropertySetType& plist, BufferType& collectables) override;
  void registerCollectables(std::vector<ObservableHelper>& h5list, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  void set_norm_factor();
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
  int num_species, N_e;
  std::vector<RealType> n_vec;
  // AA table ID
  const int d_aa_ID_;
  /////data
  //Matrix<RealType> gof_r;
  ///prefix of each gof_r
  std::vector<std::string> gof_r_prefix;
  /** resize the internal data
   * @param nbins number of bins for the historgram
   */
  void resize(int nbins);
};

} // namespace qmcplusplus
#endif

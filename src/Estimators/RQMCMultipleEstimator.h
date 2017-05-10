//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_RQMC_ESTIMATOR_H
#define QMCPLUSPLUS_RQMC_ESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"
#include "QMCDrivers/SpaceWarp.h"
#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus
{

class QMCHamiltonian;
class TrialWaveFunction;

struct RQMCEstimator: public ScalarEstimatorBase
{

  enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, LE_INDEX};

  ///index to keep track how many times accumulate is called
  int CurrentWalker;

  ///total number of walkers
  int NumWalkers;

  ///starting column index of this estimator
  int FirstColumnIndex;

  ///index of the first Hamiltonian element
  int FirstHamiltonian;

  /** the number of duplicated wavefunctions/hamiltonians
   *
   * NumCopies is the row size of Walker_t::Properties
   */
  int NumCopies;

  /** the number of components of a Hamiltonian
   *
   * It also works as the column index for the local energy. For instacne,
   * UmbrellaSamples[walker_index]->operator()(copy_index,NumOperators)
   * is the local energy of the copy_index wavefunction for the walker_index walker
   */
  int NumOperators;

  // containers for observables (copied from RQMCMultiple.cpp)
  std::vector<double>AveEloc,AveWeight;

  ///The Reptile: a chain of beads
  MultiChain* Reptile;

  int Directionless, nSteps;
  double Tau;

  ///the index of weight
  //int WeightIndex;

  /** the column size of EnergyContaner_t for UmbrellaSamples
   *
   * NumCols = NumOperators+2
   * 2 is for the local energy and weight
   int NumCols;
   */

  /** local energy data for each walker.
  */
  //Matrix<RealType> UmbrellaEnergy, UmbrellaWeight, RatioIJ;
  Matrix<RealType> RatioIJ;

  /** accumulated total local energies, weights
   *
   * The dimension of esum is NumCopies by LE_INDEX
   */
  Matrix<RealType> esum;

  /** accumulated local operators, e.g., Kinetic energy
   *
   * The dimension of elocal is NumCopies by NumOperators.
   */
  Matrix<RealType> elocal;

  /** the names of esum
   *
   * The dimension of esum_name is LE_INDEX by NumCopies. The columns
   * of the same field are bunched together.
   */
  Matrix<std::string> esum_name;

  /** the names of elocal
   *
   * The dimension of elocal_name is NumOperators by NumCopies. The columns
   * of the same field are bunched together.
   */
  Matrix<std::string> elocal_name;

  /** the names of energy differences
  */
  std::vector<std::string> ediff_name;

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  RQMCEstimator(QMCHamiltonian& h, int hcopy=1);

  RQMCEstimator(const RQMCEstimator& rest);

  void accumulate(const Walker_t& awalker, RealType wgt);

  /*@{*/
  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    std::accumulate(**first,1.0);
  }
  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();
  /*@}*/


  //RealType std::accumulate(WalkerIterator first, WalkerIterator last) {
  //  int nw=0;
  //  while(first != last) {
  //    std::accumulate(**first,(*first)->Weight);
  //    ++first;++nw;
  //  }
  //  return nw;
  //}

  void initialize(MultiChain* reptileRef, int setDirect, double setTau, int setSteps);

  ///** initialize the multi-configuration data
  // *
  // * @param W MCWalkerConfiguration
  // * @param h Collection of QMCHamiltonians*
  // * @param psi Collection of TrialWaveFunctions*
  // * @param tau time step for the initial weighted drift
  // * @param require_register if true, use buffer for particle-by-particle
  // */
  //void initialize(MCWalkerConfiguration& W,
  //    std::vector<QMCHamiltonian*>& h, std::vector<TrialWaveFunction*>& psi,
  //    RealType tau, std::vector<RealType>& Norm, bool require_register=false);
  //void initialize(MCWalkerConfiguration& W, std::vector<ParticleSet*>& WW, SpaceWarp& Warp,
  //    std::vector<QMCHamiltonian*>& h, std::vector<TrialWaveFunction*>& psi,
  //    RealType tau, std::vector<RealType>& Norm, bool require_register=false);
};

}
#endif

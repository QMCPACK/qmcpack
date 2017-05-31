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
    
    
    
    
    
    
    
    
#ifndef QMCPLUSPLUS_STRUCTUREFACTOR_ESTIMATOR_H
#define QMCPLUSPLUS_STRUCTUREFACTOR_ESTIMATOR_H
#include "Estimators/CompositeEstimators.h"
#include "Estimators/VectorEstimatorImpl.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus
{

class SkEstimator: public CompositeEstimatorBase
{
  /** number of species */
  int NumSpecies;
  /** number of kpoints */
  int NumK;
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
  Vector<ComplexType> RhokTot;
  /** instantaneous structure factor  */
  Vector<RealType> SkInst;
public:

  /** constructor
   * @param source particleset
   */
  SkEstimator(ParticleSet& source);

  /** copy constructor
   */
  SkEstimator(const SkEstimator& a);

  /** virtal destrctor */
  ~SkEstimator();

  //@{
  ///implement virtual functions
  /** create gofr group */
  hid_t createGroup(hid_t gid);
  void resetTargetParticleSet(ParticleSet& p);
  void startAccumulate();
  void accumulate(ParticleSet& p);
  void stopAccumulate(RealType wgtnorm);
  ///TODO: not implemented
  CompositeEstimatorBase* clone();
};
}

#endif

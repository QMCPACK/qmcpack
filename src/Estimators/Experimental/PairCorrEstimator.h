//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_PAIRCORRELATION_ESTIMATORBASE_H
#define QMCPLUSPLUS_PAIRCORRELATION_ESTIMATORBASE_H
#include "Estimators/CompositeEstimators.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
//#include <boost/numeric/ublas/matrix.hpp>

namespace qmcplusplus
{

//temporary
class AppendData;

class PairCorrEstimator: public CompositeEstimatorBase
{
  ///typedef for the storage
  typedef Matrix<RealType> StorageType;
  //typedef boost::numeric::ublas::matrix<RealType> StorageType;
  ///true if source == target
  bool Symmetric;
  ///save the source particleset
  const ParticleSet& sourcePtcl;
  /** distance table */
  const DistanceTableData*  myTable;
  /** number of centers */
  int Centers;
  /** number of distinct pair types */
  int NumPairTypes;
  /** maximum distance */
  RealType Dmax;
  /** bin size */
  RealType Delta;
  /** one of bin size */
  RealType DeltaInv;
  /** local copy of pair index */
  vector<int> PairID;
  /** histogram of the distances per measurement*/
  StorageType gofrInst;
  /** gofr */
  StorageType gofr;
  /** (gofr)^2 */
  StorageType gofr2;
  /** running error of gofr */
  StorageType gofrerr;
  ///debug only
  vector<ofstream*> fout;

  AppendData* v_h;
  AppendData* v2_h;
public:

  /** constructor
   * @param source particleset
   */
  PairCorrEstimator(ParticleSet& source);

  /** constructor
   * @param source particleset
   * @param target particleset
   */
  PairCorrEstimator(const ParticleSet& source, ParticleSet& target);

  /** virtal destrctor */
  ~PairCorrEstimator();

  void resetTargetParticleSet(ParticleSet& p);
  void open(hid_t hroot);
  void close();
  /** prepare data collect */
  void startAccumulate();
  /** accumulate the observables */
  void accumulate(ParticleSet& p, RealType wgt);
  /** reweight of the current cummulative  values */
  void stopAccumulate(RealType wgtinv);
  void startBlock(int steps);
  void stopBlock(RealType wgtnorm, RealType errnorm);
  void setBound(RealType rmax, RealType dr);

private:
  PairCorrEstimator(const PairCorrEstimator& pc): sourcePtcl(pc.sourcePtcl) {}
};

//class PairCorrEstimator: public CompositeEstimatorBase
//{
//  ///true if source == target
//  bool Symmetric;
//  ///save the source particleset
//  const ParticleSet& sourcePtcl;
//  /** distance table */
//  const DistanceTableData*  myTable;
//  /** number of distinc pairs */
//  int NumPairTypes;
//  /** maximum distance */
//  RealType Dmax;
//  /** bin size */
//  RealType Delta;
//  /** one of bin size */
//  RealType DeltaInv;
//  /** histogram of the distances per measurement*/
//  Vector<RealType> dCInst;
//  /** histogram of the distances per block */
//  Vector<RealType> dCBlock;
//  ///debug only
//  ofstream* fout;
//  public:

//  /** constructor
//   * @param source particleset
//   */
//  PairCorrEstimator(ParticleSet& source);

//  /** constructor
//   * @param source particleset
//   * @param target particleset
//   */
//  PairCorrEstimator(const ParticleSet& source, ParticleSet& target);

//  /** virtal destrctor */
//  ~PairCorrEstimator();

//  void resetTargetParticleSet(ParticleSet& p);
//  /** prepare data collect */
//  void startAccumulate();
//  /** accumulate the observables */
//  void accumulate(ParticleSet& p, RealType wgt);
//  /** reweight of the current cummulative  values */
//  void stopAccumulate(RealType wgtinv);
//  void startBlock(int steps);
//  void stopBlock(RealType wgtinv);

//  void setBound(RealType rmax, RealType dr);

//  private:
//  PairCorrEstimator(const PairCorrEstimator& pc): sourcePtcl(pc.sourcePtcl) {}
//};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $
 ***************************************************************************/

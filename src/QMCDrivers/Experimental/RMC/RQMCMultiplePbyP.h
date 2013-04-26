//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_REPMULTIPLE_PbyPH
#define QMCPLUSPLUS_REPMULTIPLE_PbyPH

#include "QMCDrivers/QMCDriver.h"
#include "OhmmsPETE/OhmmsVector.h"


namespace qmcplusplus
{
class Bead_ParticleSet;
class Bead;
class MultiChain;
class PolymerEstimator;

/** @ingroup QMCDrivers MultiplePsi
 * @brief Implements the RMC algorithm for energy differences
 */
class RQMCMultiplePbyP: public QMCDriver
{

public:

  enum { MinusDirection=0, PlusDirection=1, Directionless=2};
  double KEcut;
  /// Constructor.
  RQMCMultiplePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  /// Destructor
  ~RQMCMultiplePbyP();

  bool run();
  bool put(xmlNodePtr q);

protected:

  /** boolean for initialization
   *
   *\if true,
   *use clones for a chain.
   *\else
   *use drift-diffusion to form a chain
   *\endif
   */
  ///The length of polymers
  int ReptileLength;

  PolymerEstimator *multiEstimator;

  ///
  int forward,backward,itail,inext;
  vector<NewTimer*> myTimers;

  ///the number of turns per block
  int NumTurns;

  int nptcl,nObs;
  ///use scaled drift "true" = yes (default) else no
  string scaleBeadDrift;
  ///the number of H/Psi pairs
  int nPsi;
  ParticleSet::ParticlePos_t grand_transProb;//a
  ///The Reptile: a chain of beads
  MultiChain* Reptile;

  ///The new bead
  Bead *NewBead;

  ///move polymers
  void moveReptile();
  void ChooseSlices(int ReptileLength,int &startSlice,int &endSlice);
  int startSlice;
  int endSlice;
  int sliceToThrow;
  ///initialize polymers
  void initReptile();
  void initReptile_new();

  ///for the first run starting with a point, set reference properties (sign)
  void setReptileProperties();

  void checkReptileProperties();
  void checkBeadInfo(int i,bool dontDie=false);

  ///Working arrays
  Vector<RealType> NewGlobalAction,DeltaG;
  Vector<int> NewGlobalSignWgt,WeightSign;
  Vector<int> DeltaSign;
  void resizeArrays(int n);

  ///overwrite recordBlock
  void recordBlock(int block);

private:

  /// Copy Constructor (disabled)
  RQMCMultiplePbyP(const RQMCMultiplePbyP& a): QMCDriver(a) { }

  /// Copy operator (disabled).
  RQMCMultiplePbyP& operator=(const RQMCMultiplePbyP&)
  {
    return *this;
  }

  ParticleSet::ParticlePos_t gRand;
  RealType MSS,Tauoverm,sqrtTauoverm;
  RealType lambda;
  int MaxLevel;
  ParticleSet::ParticleGradient_t dG;
  ParticleSet::ParticleLaplacian_t dL;
  void InitBisection(Buffer_t &w_buffer);
  vector<Bead_ParticleSet*> tempReptile;
  vector<Bead_ParticleSet*> tempReptile2;
  vector<vector<TrialWaveFunction*> > psiReptile;
  vector<vector<QMCHamiltonian*> > hReptile;
  vector<vector<TrialWaveFunction*> > psiReptile2;
  vector<vector<QMCHamiltonian*> > hReptile2;
  void ActionChange_displace(bool oldData);
  void moveReptile_displace();


  void RegisterMe();


  RealType LogSampleProb(vector<Bead_ParticleSet*> &tempReptile,
                         int startSlice, int endSlice,
                         vector<int> &particles, int level);

  RealType BisectionSample(vector<Bead_ParticleSet*> &tempReptile,
                           int particleToMove,
                           int level);

  /*    double BisectionSample(vector<Bead*> &tempReptile,  */
  /* 			   int particleToMove, */
  /* 			   int level); */
  /*     double LogSampleProb(vector<Bead*> &tempReptile, */
  /* 			 int startSlice, int endSlice,  */
  /* 			 vector<int> particles, int level); */



  void  PutInBox (PosType &v);


  void moveReptile_bisection();
  void moveReptile_bisection_end();

  void ActionChange(bool oldData);
  void ActionChange_wf(bool oldData,int sliceToThrow);


};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3019 $   $Date: 2008-08-18 16:49:32 -0500 (Mon, 18 Aug 2008) $
 * $Id: RQMCMultiplePbyP.h 3019 2008-08-18 21:49:32Z jnkim $
 ***************************************************************************/

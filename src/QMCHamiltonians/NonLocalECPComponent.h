//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

/** Contains a set of radial grid potentials around a center.
*/
struct NonLocalECPComponent: public QMCTraits
{

  typedef std::vector<PosType>  SpherGridType;
  typedef OneDimGridBase<RealType> GridType;
  typedef OneDimCubicSpline<RealType> RadialPotentialType;

  ///index of the distance table
  int myTableIndex;
  ///Non Local part: angular momentum, potential and grid
  int lmax;
  ///the number of non-local channels
  int nchannel;
  ///the number of nknot
  int nknot;
  ///Maximum cutoff the non-local pseudopotential
  RealType Rmax;
  ///random number generator
  RandomGenerator_t* myRNG;
  ///Angular momentum map
  aligned_vector<int> angpp_m;
  ///Weight of the angular momentum
  aligned_vector<RealType> wgt_angpp_m;
  /// Lfactor1[l]=(2*l+1)/(l+1)
  RealType Lfactor1[8];
  /// Lfactor1[l]=(l)/(l+1)
  RealType Lfactor2[8];
  ///Non-Local part of the pseudo-potential
  std::vector<RadialPotentialType*> nlpp_m;
  ///fixed Spherical Grid for species
  SpherGridType sgridxyz_m;
  ///randomized spherical grid
  SpherGridType rrotsgrid_m;
  ///weight of the spherical grid
  std::vector<RealType> sgridweight_m;
  ///Working arrays
  std::vector<RealType> psiratio,vrad,dvrad,wvec,Amat,dAmat;
  std::vector<PosType> psigrad, psigrad_source;
  std::vector<RealType> lpol, dlpol;

  // For Pulay correction to the force
  std::vector<RealType> WarpNorm;
  ParticleSet::ParticleGradient_t dG;
  ParticleSet::ParticleLaplacian_t dL;
  /// First index is knot, second is electron
  Matrix<PosType> Gnew;
  ///The gradient of the wave function w.r.t. the ion position
  ParticleSet::ParticleGradient_t Gion;

  ///virtual particle set: delay initialization
  VirtualParticleSet* VP;

  //DistanceTableData* myTable;

#if !defined(REMOVE_TRACEMANAGER)
  ///pointers to trace data of containing NonLocalECPotential object
  Array<TraceReal,1>* Ve_sample;
  Array<TraceReal,1>* Vi_sample;
  bool streaming_particles;
#endif


  NonLocalECPComponent();

  ///destructor
  ~NonLocalECPComponent();

  NonLocalECPComponent* makeClone();

  ///add a new Non Local component
  void add(int l, RadialPotentialType* pp);

  ///add knots to the spherical grid
  void addknot(const PosType& xyz, RealType weight)
  {
    sgridxyz_m.push_back(xyz);
    sgridweight_m.push_back(weight);
  }

  void resize_warrays(int n,int m,int l);

  void randomize_grid(ParticleSet::ParticlePos_t& sphere, bool randomize);
  template<typename T> void randomize_grid(std::vector<T> &sphere);

  RealType evaluateOne(ParticleSet& W, int iat, TrialWaveFunction& Psi, 
      int iel, RealType r, const PosType& dr, bool Tmove, std::vector<NonLocalData>& Txy) const;

  RealType evaluate(ParticleSet& W, int iat, TrialWaveFunction& Psi,
                    PosType &force_iat);

  RealType evaluate(ParticleSet& W, ParticleSet &ions, int iat, TrialWaveFunction& Psi,
                    PosType &force_iat, PosType &pulay_iat);

  RealType
  evaluate(ParticleSet& W, TrialWaveFunction& Psi,int iat, std::vector<NonLocalData>& Txy,
           PosType &force_iat);

  /** compute with virtual moves */
  RealType evaluateVP(const ParticleSet& W, int iat, TrialWaveFunction& Psi);
  RealType evaluateVP(const ParticleSet& W, int iat, TrialWaveFunction& Psi,std::vector<NonLocalData>& Txy);

  RealType
  evaluateValueAndDerivatives(ParticleSet& P,
      int iat, TrialWaveFunction& psi,
      const opt_variables_type& optvars,
      const std::vector<RealType>& dlogpsi,
      std::vector<RealType>& dhpsioverpsi);

  void print(std::ostream& os);

  void initVirtualParticle(ParticleSet* qp);

  void setRandomGenerator(RandomGenerator_t* rng)
  {
    myRNG=rng;
  }

  // For space-warp transformation used in Pulay correction
  inline RealType WarpFunction (RealType r)
  {
    return 1.0/(r*r*r*r);
  }


}; //end of RadialPotentialSet

}
#endif



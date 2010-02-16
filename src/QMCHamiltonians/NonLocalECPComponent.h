//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus {

  /** Contains a set of radial grid potentials around a center.  
  */
  struct NonLocalECPComponent: public QMCTraits { 

    typedef vector<PosType>  SpherGridType;
    typedef OneDimGridBase<RealType> GridType;
    typedef OneDimCubicSpline<RealType> RadialPotentialType;

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
    vector<int> angpp_m;
    ///Weight of the angular momentum
    vector<RealType> wgt_angpp_m;
    /// Lfactor1[l]=(2*l+1)/(l+1)
    vector<RealType> Lfactor1;
    /// Lfactor1[l]=(l)/(l+1)
    vector<RealType> Lfactor2;
    ///Non-Local part of the pseudo-potential
    vector<RadialPotentialType*> nlpp_m;
    ///fixed Spherical Grid for species
    SpherGridType sgridxyz_m;
    ///randomized spherical grid
    SpherGridType rrotsgrid_m;
    ///weight of the spherical grid
    vector<RealType> sgridweight_m;
    ///Working arrays
    vector<RealType> psiratio,vrad,dvrad,wvec,Amat,dAmat;
    vector<PosType> psigrad, psigrad_source;
    vector<RealType> lpol, dlpol;

    // For Pulay correction to the force
    vector<RealType> WarpNorm;
    ParticleSet::ParticleGradient_t dG;
    ParticleSet::ParticleLaplacian_t dL;
    // First index is knot, second is electron
    Matrix<PosType> Gnew;
    // The gradient of the wave function w.r.t. the ion position
    // at each quadrature point.
    ParticleSet::ParticleGradient_t Gion;
    // For space-warp transformation used in Pulay correction
    inline RealType WarpFunction (RealType r) { return 1.0/(r*r*r*r); }


    DistanceTableData* myTable;

    NonLocalECPComponent();

    ///destructor
    ~NonLocalECPComponent();

    NonLocalECPComponent* makeClone();

    ///add a new Non Local component
    void add(int l, RadialPotentialType* pp); 

    ///add knots to the spherical grid
    void addknot(PosType xyz, RealType weight){
      sgridxyz_m.push_back(xyz);
      sgridweight_m.push_back(weight);
    }

    void resize_warrays(int n,int m,int l);

    void randomize_grid(ParticleSet::ParticlePos_t& sphere, bool randomize);
    template<typename T> void randomize_grid(vector<T> &sphere);

    RealType 
      evaluate(ParticleSet& W, int iat, TrialWaveFunction& Psi);

    RealType evaluate(ParticleSet& W, int iat, TrialWaveFunction& Psi,
		      PosType &force_iat);

    RealType evaluate(ParticleSet& W, ParticleSet &ions, int iat, TrialWaveFunction& Psi,
		      PosType &force_iat, PosType &pulay_iat);


    RealType 
    evaluate(ParticleSet& W, TrialWaveFunction& Psi,int iat, vector<NonLocalData>& Txy);

    RealType 
    evaluate(ParticleSet& W, TrialWaveFunction& Psi,int iat, vector<NonLocalData>& Txy,
	     PosType &force_iat);


    void print(std::ostream& os);

    void setRandomGenerator(RandomGenerator_t* rng) 
    {
      myRNG=rng;
    }

  }; //end of RadialPotentialSet

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


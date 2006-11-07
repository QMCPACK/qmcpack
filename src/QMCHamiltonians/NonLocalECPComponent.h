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
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus {

  /** Contains a set of radial grid potentials around a center.  
  */
  struct NonLocalECPComponent: public QMCTraits { 

    typedef vector<PosType>  SpherGridType;
    typedef OneDimGridBase<RealType> GridType;
    typedef OneDimLinearSpline<RealType> RadialPotentialType;

    ///Non Local part: angular momentum, potential and grid
    int lmax;
    ///the number of non-local channels
    int nchannel;
    ///the number of nknot
    int nknot;
    ///Maximum cutoff the non-local pseudopotential
    RealType Rmax;
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
    vector<ValueType> psiratio,vrad,wvec,Amat;
    vector<RealType> lpol;

    DistanceTableData* myTable;

    NonLocalECPComponent();

    ///destructor
    ~NonLocalECPComponent();

    ///add a new Non Local component
    void add(int l, RadialPotentialType* pp); 

    ///add knots to the spherical grid
    void addknot(PosType xyz, RealType weight){
      sgridxyz_m.push_back(xyz);
      sgridweight_m.push_back(weight);
    }

    void resize_warrays(int n,int m,int l);

    ///Randomly rotate sgrid_m
    template <class PA>
      inline void randomize_grid(PA& sphere, bool randomize) {
        if(randomize) {
          const RealType twopi(6.28318530718);
          RealType phi(twopi*Random()),psi(twopi*Random()),cth(Random()-0.5),
                   sph(std::sin(phi)),cph(std::cos(phi)),sth(std::sqrt(1-cth*cth)),sps(std::sin(psi)),
                   cps(std::cos(psi));
          Tensor<double,3> rmat( cph*cth*cps-sph*sps, sph*cth*cps+cph*sps,-sth*cps,
              -cph*cth*sps-sph*cps,-sph*cth*sps+cph*cps, sth*sps,
              cph*sth,             sph*sth,             cth     );
          SpherGridType::iterator it(sgridxyz_m.begin());
          SpherGridType::iterator it_end(sgridxyz_m.end());
          SpherGridType::iterator jt(rrotsgrid_m.begin());
          int ic=0;
          while(it != it_end) {*jt = dot(rmat,*it); ++it; ++jt;}
          //copy the radomized grid to sphere
          std::copy(rrotsgrid_m.begin(), rrotsgrid_m.end(), sphere.begin());
        } else {
          //copy sphere to the radomized grid
          std::copy(sphere.begin(), sphere.end(), rrotsgrid_m.begin());
        }
      }

    ValueType 
      evaluate(ParticleSet& W, int iat, TrialWaveFunction& Psi);

    ValueType 
    evaluate(ParticleSet& W, TrialWaveFunction& Psi,int iat, vector<NonLocalData>& Txy);
  }; //end of RadialPotentialSet

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


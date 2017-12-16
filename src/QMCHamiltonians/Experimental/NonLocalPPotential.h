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
    
    
#ifndef QMCPLUSPLUS_NONLOCALPPOTENTIAL_H
#define QMCPLUSPLUS_NONLOCALPPOTENTIAL_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

class WalkerSetRef;

/** @ingroup hamiltonian
 * \brief Evaluate the semi local potentials
 */
struct NonLocalPPotential: public QMCHamiltonianBase
{

  typedef std::vector<PosType>  SpherGridType;
  typedef OneDimGridBase<ValueType> GridType;
  typedef OneDimCubicSpline<ValueType> LocalPotentialType;

  /** Contains a set of radial grid potentials around a center.
   */
  struct RadialPotentialSet
  {
    ///true, if this species has RadialPotential
    bool HasNonLocalPP;
    ///Non Local part: angular momentum, potential and grid
    int lmax;
    ///the number of non-local channels
    int nchannel;
    ///the number of nknot
    int nknot;
    ///Maximum cutoff the non-local pseudopotential
    RealType Rmax;
    ///Angular momentum map
    std::vector<int> angpp_m;
    ///Weight of the angular momentum
    std::vector<RealType> wgt_angpp_m;
    /// Lfactor1[l]=(2*l+1)/(l+1)
    std::vector<RealType> Lfactor1;
    /// Lfactor1[l]=(l)/(l+1)
    std::vector<RealType> Lfactor2;

    ///Local part of the pseudo-potential
    LocalPotentialType* lpp_m;
    ///Grid of the local potential
    GridType* lgrid_m;
    ///Non-Local part of the pseudo-potential
    std::vector<LocalPotentialType*> nlpp_m;
    ///Grid of the non-local potential
    std::vector<GridType*> nlgrid_m;
    ///fixed Spherical Grid for species
    SpherGridType sgridxyz_m;
    ///randomized spherical grid
    SpherGridType rrotsgrid_m;
    ///weight of the spherical grid
    std::vector<ValueType> sgridweight_m;

    ///Working arrays
    std::vector<ValueType> psiratio,vrad,wvec,Amat,lpol;

    inline RadialPotentialSet():
      HasNonLocalPP(false), lmax(0), nchannel(0), nknot(0),
      Rmax(-1), lpp_m(0), lgrid_m(0) {}

    ///destructor
    ~RadialPotentialSet();

    ///add a new Local component
    void add(GridType* agrid, LocalPotentialType* pp)
    {
      lgrid_m=agrid;
      lpp_m=pp;
    }

    ///add a new Non Local component
    void add(int angmom, GridType* agrid, LocalPotentialType* pp)
    {
      angpp_m.push_back(angmom);
      wgt_angpp_m.push_back(static_cast<RealType>(2*angmom+1));
      nlgrid_m.push_back(agrid);
      nlpp_m.push_back(pp);
    }

    ///add knots to the spherical grid
    void addknot(PosType xyz, ValueType weight)
    {
      sgridxyz_m.push_back(xyz);
      sgridweight_m.push_back(weight);
    }

    void resize_warrays(int n,int m,int l);

    ///Randomly rotate sgrid_m
    template <class PA>
    inline void randomize_grid(PA& sphere, bool randomize)
    {
      if(randomize)
      {
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
        while(it != it_end)
        {
          *jt = dot(rmat,*it);
          ++it;
          ++jt;
        }
        //copy the radomized grid to sphere
        copy(rrotsgrid_m.begin(), rrotsgrid_m.end(), sphere.begin());
      }
      else
      {
        //copy sphere to the radomized grid
        copy(sphere.begin(), sphere.end(), rrotsgrid_m.begin());
      }
    }

    ValueType
    evaluate(ParticleSet& W, DistanceTableData* d_table,
             int iat, TrialWaveFunction& Psi, bool randomize);
  }; //end of RadialPotentialSet
  /// Maximum number of channels (all centers)
  int maxnonloc;
  /// Maximum number of spherical grid points (all centers)
  int maxsgridpts;
  /// Highest angular momentum channel (all centers)
  int maxangmom;
  ///the distance table containing electron-nuclei distances
  DistanceTableData* d_table;
  ///the set of local-potentials (one for each ion)
  std::vector<RadialPotentialSet*> PP;
  ///reference to the center ion
  ParticleSet& IonConfig;
  ///unique index for each ion
  ParticleSet::ParticleIndex_t Centers;

  TrialWaveFunction* Psi;

  NonLocalPPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

  ~NonLocalPPotential();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "NoLocalPPotential: " << d_table->origin().getName();
    return true;
  }


  QMCHamiltonianBase* clone();
};
}
#endif



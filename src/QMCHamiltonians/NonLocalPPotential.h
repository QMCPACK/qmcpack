//////////////////////////////////////////////////////////////////
// (c) Copyright 2004 by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_QMC_NONLOCALPPOTENTIAL_H
#define OHMMS_QMC_NONLOCALPPOTENTIAL_H
#include <fstream>
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/HDFNumericAttrib.h"

namespace ohmmsqmc {


  class WalkerSetRef;

  /**
   * \brief Evaluate the local potentials (either pseudo or full
   core) around each ion.
   */

  struct NonLocalPPotential: public QMCHamiltonianBase {

  /**
   *\brief Contains a set of radial grid potentials around a 
   center.  
   *
   Evaluates the grid potential given the electron-
   ion distances.
   */
  struct RadialPotentialSet: public QMCTraits {

    typedef OneDimGridBase<ValueType> GridType;
    typedef OneDimGridFunctor<ValueType> LocalPotentialType;
    typedef vector<PosType> SpherGridType; 

    ///Local part of the pseudo-potential
    vector<LocalPotentialType*> lpp_m;
    vector<GridType*> lgrid_m;

    ///Non Local part: angular momentum, potential and grid
    int lmax;
    RealType Rmax;
    vector<int> angpp_m;
    vector<LocalPotentialType*> nlpp_m;
    vector<GridType*> nlgrid_m;
     
    ///Spherical Grid
    SpherGridType sgridxyz_m;
    vector<double> sgridweight_m;

    ///destructor
    ~RadialPotentialSet();

    ///add a new Local component
    void add(GridType* agrid, LocalPotentialType* pp) {
      lgrid_m.push_back(agrid);
      lpp_m.push_back(pp);
    }

    ///add a new Non Local component
    void add(int angmom, GridType* agrid, LocalPotentialType* pp) {
      angpp_m.push_back(angmom);
      nlgrid_m.push_back(agrid);
      nlpp_m.push_back(pp);
      if(lmax<angmom)lmax=angmom;
    }

    ///add knots to the spherical grid
    void addknot(PosType xyz, double weight){
      sgridxyz_m.push_back(xyz);
      sgridweight_m.push_back(weight);
    }

    ///Randomly rotate sgrid_m
    void randomize_grid();

    PosType getknot(int i){
      return sgridxyz_m[i];
    }

    /*!
     *\param d_table the distance table containing electron-
     ion distances  
     *\param iat the index of the nucleus
     *\return \f$ \sum_{I=1}^{N_C}\sum_{i=1}^{N_P} V(r_{iI}) \f$
     *\brief Evaluates the set of radial grid potentials around
     a single ion.
     */
    inline ValueType evaluate(DistanceTableData* d_table, int iat) {
      ValueType esum = 0.0;     
      for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++){
	for(int ig=0; ig<lgrid_m.size(); ig++) 
	  lgrid_m[ig]->index(d_table->r(nn));
	for(int ip=0;ip<lpp_m.size(); ip++) 
	  esum += lpp_m[ip]->evaluate(d_table->r(nn),d_table->rinv(nn));
      }
      return esum;
    }
  };

    ///the distance table containing electron-nuclei distances  
    DistanceTableData* d_table;
    ///the set of local-potentials (one for each ion)
    vector<RadialPotentialSet*> PP;
    ///unique index for each ion
    const ParticleSet::ParticleIndex_t& Centers;
    /// Maximum number of channels (all centers)
    int maxnonloc;
    /// Maximum number of spherical grid points (all centers)
    int maxsgridpts;
    /// Highest angular momentum channel (all centers)
    int maxangmom;
    /// Working arrays
    ValueType *psiratio,*vrad,*wvec,*Amat,*lpol;

    TrialWaveFunction& Psi;

    typedef RadialPotentialSet::GridType GridType;
    typedef RadialPotentialSet::LocalPotentialType LocalPotentialType;

    NonLocalPPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

    ~NonLocalPPotential();
     
    ValueType evaluate(ParticleSet& W);

    inline ValueType evaluate(ParticleSet& P, RealType& x) {
      return x=evaluate(P);
    }

#ifdef WALKERSTORAGE_COLUMNMAJOR
    inline void 
      evaluate(WalkerSetRef& P, ValueVectorType& LE) { }
#else
    inline void 
      evaluate(WalkerSetRef& P, ValueVectorType& LE) { }
#endif
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


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
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace ohmmsqmc {


  class WalkerSetRef;

  /**
   * \brief Evaluate the local potentials (either pseudo or full
   core) around each ion.
   */

  struct NonLocalPPotential: public QMCHamiltonianBase {

    typedef OneDimGridBase<ValueType> GridType;
    typedef OneDimGridFunctor<ValueType> LocalPotentialType;


    /**
     *\brief Contains a set of radial grid potentials around a 
     center.  
     *
     Evaluates the grid potential given the electron-
     ion distances.
     */
    struct RadialPotentialSet { 

      typedef vector<PosType> SpherGridType; 

      ///Local part of the pseudo-potential
      vector<LocalPotentialType*> lpp_m;
      vector<GridType*> lgrid_m;

      ///Non Local part: angular momentum, potential and grid
      int lmax;
      RealType Rmax;
      vector<int> angpp_m, wgt_angpp_m;
      vector<LocalPotentialType*> nlpp_m;
      vector<GridType*> nlgrid_m;

      ///Spherical Grid
      SpherGridType sgridxyz_m, rrotsgrid_m;
      vector<ValueType> sgridweight_m;

      ///Working arrays
      vector<ValueType> psiratio,vrad,wvec,Amat,lpol;

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
	wgt_angpp_m.push_back(2*angmom+1);
	nlgrid_m.push_back(agrid);
	nlpp_m.push_back(pp);
      }

      ///add knots to the spherical grid
      void addknot(PosType xyz, ValueType weight){
	sgridxyz_m.push_back(xyz);
	sgridweight_m.push_back(weight);
      }
      
      void resize_warrays(int n,int m,int l){
	psiratio.resize(n);
	vrad.resize(m);
	wvec.resize(m);
	Amat.resize(n*m);
	lpol.resize(l+1);
	rrotsgrid_m.resize(n);
      }

      ///Randomly rotate sgrid_m
      void randomize_grid();

      ///
      ValueType evaluate(ParticleSet& W, DistanceTableData* d_table, 
	  int iat, TrialWaveFunction& Psi);

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
	  //for(int ig=0; ig<lgrid_m.size(); ig++) 
	  //  lgrid_m[ig]->index(d_table->r(nn));
	  for(int ip=0;ip<lpp_m.size(); ip++) 
	    esum += lpp_m[ip]->evaluate(d_table->r(nn),d_table->rinv(nn));
	}
	return esum;
      }
    }; //end of RadialPotentialSet


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

    TrialWaveFunction& Psi;

    NonLocalPPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi);

    ~NonLocalPPotential();

    inline ValueType evaluate(ParticleSet& W) {
      
      RealType esum=0.0;
      //loop over all the ions
      for(int iat=0; iat<Centers.size(); iat++) {
	esum += PP[Centers[iat]]->evaluate(d_table,iat);
	if(PP[Centers[iat]]->nlpp_m.size()==0) continue;
	esum += PP[Centers[iat]]->evaluate(W,d_table,iat,Psi);
      }
      return esum;
    }

    inline ValueType evaluate(ParticleSet& P, RealType& x) {
      return x=evaluate(P);
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_PARTICLEBASIC_FUNCTIONS_H
#define OHMMS_PARTICLEBASIC_FUNCTIONS_H

#include "ParticleBase/ParticleUtility.h"

/**@file ParticleFunctions.h
 *@ingroup ptclutil
 * Utility functions that operate on Particle class
 */

/**@fn template<class PL, class GIM> void PartitionGrid(PL& lattice, GIM& grid)
 *@param lattice Particle Layout class, e.g., CrystalLattice<T,D>
 *@param grid three by D-dimensional index array
 *
 *@brief Make three-level partition of the layout.
 *
 *grid(i,j) for i the level index and j the direction index
 *
 *Requirement for template GIM: access operator [i][j]
 *
 *Recommended template GIM: vector<TinyVector<int,D> > 
 *
 *The levels are:
 *<ul>
 *<il> i=0 for MPI level
 *<il> i=1 for OMP level
 *<il> i=2 for grid level
 *</ul>
 */
template<class PL, class GIM>
void PartitionGrid(PL& lattice, GIM& grid) {
  
  const int ndim = PL::DIM;
  lattice.clear();
    
  //==============================================================
  // grid partitions are created
  //==============================================================
  int mppgrid_id = lattice.makeGrid(&grid[0][0],-1);
  int ompgrid_id = lattice.makeGrid(&grid[1][0],mppgrid_id);

  //==============================================================
  // \warning makeGrid function refines a subgrid based on the parent grid
  //==============================================================
  int g[ndim];
  for(int idim=0; idim<ndim; idim++) {
    if(grid[1][idim] > grid[2][idim]) grid[1][idim] = grid[2][idim];
    g[idim] = grid[0][idim]*grid[1][idim];
    grid[2][idim] /= g[idim];
  }

  lattice.makeGrid(&grid[2][0],ompgrid_id);
  for(int idim=0; idim<ndim; idim++) {
    //lattice.Grid[idim] = grid[2][idim]*grid[1][idim];
    lattice.Grid[idim] = grid[2][idim]*g[idim];
  }

}
  

/**@fn template<class PT, class GIV, class GIM> 
 * void ExpandSuperCell(PT& in_, PT& out_, GIV& uc_grid, GIM& grid)
 *@param in_ the original Particle object to be expanded.
 *@param out_ the output Particle object
 *@param uc_grid D-dimensional index array
 *@param grid three by D-dimensional index array
 *
 *@brief Expand the input Particle object by uc_grid. 
 *
 <ul>
 <li>uc_grid[i] for the number of unit cell in the first direction to 
 be duplicated.

 Requirement for template GIV: access operator [i]

 Recommended template GIV: TinyVector<int,D>  

 <li> grid(i,j) for i the level index and j the direction index

 Requirement for template GIM: access operator [i][j]

 Recommended template GIM: vector<TinyVector<int,D> > 
 
 The levels are:
 <ul>
 <il> i=0 for MPI level
 <il> i=1 for OMP level
 <il> i=2 for grid level
 </ul>
 </ul> 
*/
template<class PT, class GIV, class GIM>
bool ExpandSuperCell(PT& in_, 
		     PT& out_, 
		     GIV& uc_grid,
		     GIM& grid) {
  
  ///Prepare lattice
  out_.Lattice.set(in_.Lattice,&uc_grid[0]);
  typename PT::ParticleLayout_t& lattice = out_.Lattice;

  ///assign the grid information
  //for(int i=0; i<grid.size(); i++)  lattice.Grid[i] = grid[i];

  lattice.makeGrid(grid);

  ///\todogoing to use selectedset
  int natom = in_.getLocalNum();
  int ntot = uc_grid[0]*uc_grid[1]*uc_grid[2]*natom;
  
  ///first convert the input position to unit, if not in unit
  in_.Lattice.convert2Unit(in_.R);

  typedef typename PT::Scalar_t Scalar_t;
  typename PT::SingleParticlePos_t del;
  const typename PT::ParticlePos_t& R = in_.R;
  const typename PT::ParticleIndex_t& gid = in_.GroupID;

  out_.create(ntot);
  out_.R.setUnit(PosUnit::CartesianUnit);
  int nacc = 0;
  for(int ic=0; ic<uc_grid[0]; ic++) {
    del[0] = static_cast<Scalar_t>(ic);
    for(int jc=0; jc<uc_grid[1]; jc++) {
      del[1] = static_cast<Scalar_t>(jc);
      for(int kc=0; kc<uc_grid[2]; kc++) {
	del[2] = static_cast<Scalar_t>(kc);
	for(int nb=0; nb<in_.getLocalNum(); nb++) {
	  out_.R[nacc] = in_.Lattice.toCart(R[nb]+del);
	  out_.ID[nacc] = nacc;
	  out_.GroupID[nacc] = gid[nb];
	  nacc++;
	}
      }
    }
  }

  return true;
}

template<class PT, class GIV>
inline 
bool ExpandSuperCell(PT& inout_, GIV& uc_grid) {

  PT temp(inout_);
  vector<GIV> grid(3);
  grid[0] = GIV(1);
  grid[1] = GIV(1);
  grid[2] = uc_grid;

  inout_.clear();
  return ExpandSuperCell(temp,inout_,uc_grid, grid);
}

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

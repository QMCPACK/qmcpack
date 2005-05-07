//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_UNIFORMGRIDLAYOUT_H
#define OHMMS_UNIFORMGRIDLAYOUT_H

#include "Lattice/CrystalLattice.h"
#include "Lattice/UniformCartesianGrid.h"
#include <iostream>

/**@file UniformGridLayout.h
 *@brief Declaration of UniformGridLayout<T,D>
 */

/** generic UniformGridLayout. Do nothing. 
*/
template<class T, unsigned D> struct UniformGridLayout{ };

/** specialization of UniformGridLayout<T,3> for 3-Dim layout 
 *
 *This class represents a supercell and is one of the layout classes 
 *to define a trait class PT for ParticleBase<PT> (as a matter of fact,
 *the only layout class implemented in ohmms).
 *
 *It is inherited from CrystalLattice<T,3> and
 *UniformCartesianGrid<T,3>: CrystalLattice defines a general
 *supercell and UniformCartesianGrid defines a cubic-cell view of a
 *general supercell. This layout class is intended for typical
 *condensed matter or molecular systems, where a uniform grid partition
 *works fairly well.
 *
 *The additional interfaces of UniformGridLayout class are primarily
 *for NNEngine classes, which use the parition information to evaluate
 *the neighbor lists efficiently. 
 *
 */
template<class T>
class UniformGridLayout<T,3>: public CrystalLattice<T,3>,
			      public UniformCartesianGrid<T,3> { 

public:

  /**enumeration for grid levels*/
  enum {MPI_GRID = 0, /*!< mpi level */
	OMP_GRID,     /*!< open mp level */
	SPATIAL_GRID, /*!< spatial level */
	GridLevel = 3 /*!< counter used to tell the number of levels */
  };

  typedef CrystalLattice<T,3>       Base_t;
  typedef UniformCartesianGrid<T,3> Grid_t;
  typedef UniformGridLayout<T,3>    This_t;

  typedef typename Base_t::SingleParticlePos_t   SingleParticlePos_t;
  typedef typename Base_t::SingleParticleIndex_t SingleParticleIndex_t;

  /**default constructor
   *
   *Create 1x1x1 cubic view for a supercell assuming that no
   *parallelization is used.
   */
  inline UniformGridLayout(){
    for(int i=0; i<GridLevel; i++) Grid[i] = SingleParticleIndex_t(1);
    SuperGrid.reserve(GridLevel);
    for(int i=0; i<GridLevel; i++) SuperGrid.push_back(NULL);
  } 


  inline ~UniformGridLayout(){
    for(int i=0; i<SuperGrid.size(); i++) 
      if(SuperGrid[i]) delete SuperGrid[i];
  } 

  ///return the first cell connected to the ig cell
  inline int first_connected(int ig) const { return c_offset[ig];}
    
  ///return the last cell connected to the ig cell
  inline int last_connected(int ig) const { return c_offset[ig+1];}
    
  /// return the cell index
  inline int id_connected(int j) const { return c_id[j];}
    
  /// return the correction vector according to the boundary conditions
  inline SingleParticlePos_t bc(int j) const { return c_bc[j];}

  /// return the total number of paritions
  inline int ncontexts() const { return c_offset.size()-1;}
    
  /// return the maximum number of connected cells
  int connectGrid(T rmax) {
    ///create the spatial grid
    setGrid(Grid[SPATIAL_GRID]);

    SingleParticlePos_t u0(this->Delta[0],0.0,0.0);
    SingleParticlePos_t u1(0.0,this->Delta[1],0.0);
    SingleParticlePos_t u2(0.0,0.0,this->Delta[2]);
    T RmaxSq = rmax*rmax;

    ///calculate the extend of linked cells
    int nx = static_cast<int>(sqrt(RmaxSq/Dot(u0,u0)))+1;
    int ny = static_cast<int>(sqrt(RmaxSq/Dot(u1,u1)))+1;
    int nz = static_cast<int>(sqrt(RmaxSq/Dot(u2,u2)))+1;


    c_offset.resize(this->NumGrids+1);
    int ntot = this->NumGrids*(2*nx+1)*(2*ny+1)*(2*nz+1);
    if(c_id.capacity() < ntot) c_id.reserve(ntot);
    if(c_bc.capacity() < ntot) c_bc.reserve(ntot);
    if(u_bc.capacity() < ntot) u_bc.reserve(ntot);

    int maxnc = 0, gtot = 0;
    SingleParticlePos_t dx(this->Delta[0],this->Delta[1],this->Delta[2]),org,d;
    c_offset[0] = 0;
    for(int ig=0; ig<this->NP[0]; ig++) {
      org[0] = (static_cast<T>(ig)+0.5)*dx[0];
      for(int jg=0; jg<this->NP[1]; jg++) {
        org[1] = (static_cast<T>(jg)+0.5)*dx[1];
        for(int kg=0; kg<this->NP[2]; kg++) {
          org[2] = (static_cast<T>(kg)+0.5)*dx[2];
          T x,y,z;
          int nconnected = 0;
          for(int ix=-nx; ix<=nx; ix++) {
            d[0] = x = org[0]+T(ix)*dx[0];
            if(this->BoxBConds[0]) {
              x = fmod(d[0],1.0);
              if(x<0.0) x += 1.0;
            }
            if(x<0 || x>=1) continue;
            d[0] -= x;
            for(int jx=-ny; jx<=ny; jx++) {
              d[1] = y = org[1]+T(jx)*dx[1];
              if(this->BoxBConds[1]) {
                y = fmod(d[1],1.0); if(y<0) y += 1.0;
              }
              if(y<0 || y>=1) continue;
              d[1] -= y;
              for(int kx=-nz; kx<=nz; kx++) {
                d[2] = z = org[2]+T(kx)*dx[2];
                if(this->BoxBConds[2]) {
          	z = fmod(d[2],1.0); if(z<0) z += 1.0;
                }
                if(z<0 || z>=1) continue;
                d[2] -= z;
                int iloc = loc(x,y,z);
                if(iloc == gtot && ix == 0 && jx == 0 && kx == 0) continue;
                c_id.push_back(iloc);
                u_bc.push_back(d);
                c_bc.push_back(toCart(d));
                nconnected++;
              }
            }
          }
          c_offset[gtot+1] = c_offset[gtot]+nconnected; gtot++;
          maxnc = max(maxnc,nconnected);
        }
      }
    }  
    return maxnc; // return the maxmimum number of connected cells
  }
    
  template<class GIM>
  inline void makeGrid(const GIM& mgrid) {
    if(mgrid.size() < GridLevel) {
      Grid[SPATIAL_GRID] = mgrid.back();
    } else {
      for(int ig=0; ig<GridLevel; ig++)  Grid[ig] = mgrid[ig];
    }
    ///first build the spatial grid
    setGrid(Grid[SPATIAL_GRID]);
  }

  inline Grid_t* restrict 
  addGrid(int glevel, const SingleParticleIndex_t& agrid) {
    if(SuperGrid[glevel]) {
      SuperGrid[glevel]->setGrid(agrid);
      return SuperGrid[glevel];
    } else {
      Grid_t* g = new Grid_t(agrid);    
      SuperGrid[glevel] = g;
      return g;
    }
  }

  inline Grid_t* restrict getGrid(int glevel) { return SuperGrid[glevel];}
  inline const Grid_t* restrict getGrid(int glevel) const { 
    return SuperGrid[glevel];
  }

  inline int ngrid(int ig) const { return Grid[SPATIAL_GRID][ig];}
  inline int ngrid(int glevel, int ig) const { return Grid[glevel][ig];}
  void print(std::ostream& os) const;

  inline void update() {
    for(int i=0; i<u_bc.size(); i++)  c_bc.push_back(u_bc[i]);
  }

  /** set the lattice vector with a tensor
   *@param lat a tensor representing a supercell
   *
   *In addition to setting the internal parameters, it updates the Cartesian displacement vectors.
   */
  inline void update(const Tensor<T,3>& lat) {
    CrystalLattice<T,3>::set(lat);
    for(int i=0; i<u_bc.size(); i++)  c_bc.push_back(u_bc[i]);
  }

private:

  ///Uniform grid partitions at MPI_GRID, OPENMP_GRID, and SPATIAL_GRID levels.
  SingleParticleIndex_t Grid[GridLevel];

  ///UniformCartesianGrid for multi levels.
  std::vector<Grid_t*> SuperGrid;
  ///offsets to determine cell conection
  std::vector<int> c_offset;
  ///cell index connected to each cell
  std::vector<int> c_id;
  ///displacement vectors due to the boundary conditions
  std::vector<SingleParticlePos_t> c_bc;
  ///displacement vectors due to the boundary conditions in the unit-cell unit
  std::vector<SingleParticlePos_t> u_bc;
};


template<class T>
void UniformGridLayout<T,3>::print(std::ostream& os) const {
  Base_t::print(os);
  Grid_t::printGrid(os);
  for(int ig=0; ig<c_offset.size()-1; ig++) {
    os << ig << " has neighboring cell "  
       << c_offset[ig+1]-c_offset[ig]<< std::endl;
    for(int ii=c_offset[ig]; ii<c_offset[ig+1]; ii++) {
      os << c_id[ii] << " " << c_bc[ii] << std::endl;
    }
  }
}
#endif  
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

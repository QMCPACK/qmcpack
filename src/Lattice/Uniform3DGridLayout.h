//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_UNIFORMGRIDLAYOUT_3D_NONTEMP_H
#define OHMMS_UNIFORMGRIDLAYOUT_3D_NONTEMP_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <vector>
#include <iostream>
#include "Lattice/CrystalLattice.h"
#include "Lattice/UniformCartesianGrid.h"

/**@file Uniform3DGridLayout.h
 *@brief Concrete uniform grid layout for 3D
 */

namespace qmcplusplus
{
/** A concrete, specialized class equivalent to UniformGridLayout<T,3>
 *
 *This class represents a supercell and is one of the layout classes
 *to define a trait class PT for ParticleBase<PT> (as a matter of fact,
 *the only layout class implemented in ohmms).
 *
 *It is inherited from CrystalLattice<OHMMS_PRECISION,3> and
 *UniformCartesianGrid<OHMMS_PRECISION,3>: CrystalLattice defines a general
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
class Uniform3DGridLayout: public CrystalLattice<OHMMS_PRECISION,3,OHMMS_ORTHO>,
  public UniformCartesianGrid<OHMMS_PRECISION,3>
{

public:

  enum {IsOrthogonal = OHMMS_ORTHO};

  /**enumeration for grid levels*/
  enum {MPI_GRID = 0, /*!< mpi level */
        OMP_GRID,     /*!< open mp level */
        SPATIAL_GRID, /*!< spatial level */
        GridLevel = 3 /*!< counter used to tell the number of levels */
       };

  typedef OHMMS_PRECISION                          value_type;
  typedef CrystalLattice<value_type,3,OHMMS_ORTHO> Base_t;
  typedef UniformCartesianGrid<value_type,3>       Grid_t;
  typedef Uniform3DGridLayout                      This_t;

  typedef Base_t::SingleParticlePos_t   SingleParticlePos_t;
  typedef Base_t::SingleParticleIndex_t SingleParticleIndex_t;

  /**default constructor
   *
   *Create 1x1x1 cubic view for a supercell assuming that no
   *parallelization is used.
   */
  inline Uniform3DGridLayout():MaxConnections(0), LR_dim_cutoff(15.0), NCMax(0), Status(0),
    LR_rc(1e6), LR_kc(0.0)
  {
    for(int i=0; i<GridLevel; i++)
      Grid[i] = SingleParticleIndex_t(1);
    SuperGrid.reserve(GridLevel);
    for(int i=0; i<GridLevel; i++)
      SuperGrid.push_back(NULL);
  }

  ///**copy constructor
  // */
  //inline Uniform3DGridLayout(const Uniform3DGridLayout& pl){
  //  copy(pl);
  //}


  /** copy function
   */
  inline void copy(const Uniform3DGridLayout& pl)
  {
    MaxConnections=0;
    LR_dim_cutoff=pl.LR_dim_cutoff;
    LR_kc=pl.LR_kc;
    LR_rc=pl.LR_rc;
    Base_t::set(pl);
    Grid_t::makeCopy(pl);
  }

  inline ~Uniform3DGridLayout()
  {
    for(int i=0; i<SuperGrid.size(); i++)
      if(SuperGrid[i])
        delete SuperGrid[i];
  }

  ///return the first cell connected to the ig cell
  inline int first_connected(int ig) const
  {
    return c_offset[ig];
  }

  ///return the last cell connected by the interaction radius to the ig cell
  inline int last_connected(int ig) const
  {
    return c_max[ig];
  }

  ///return the last cell connected by the connectivity to the ig cell
  inline int max_connected(int ig) const
  {
    return c_offset[ig+1];
  }

  /// return the cell index
  inline int id_connected(int j) const
  {
    return c_id[j];
  }

  /// return the correction vector according to the boundary conditions
  inline SingleParticlePos_t bc(int j) const
  {
    return c_bc[j];
  }

  /// return the total number of paritions
  inline int ncontexts() const
  {
    return c_offset.size()-1;
  }

  /// return the maximum number of connected cells
  int connectGrid(value_type int_rad, value_type con_rad=-1);

  void makeShell(std::vector<SingleParticleIndex_t>& RS,
                 std::vector<bool>& FS,
                 const SingleParticleIndex_t& nc,
                 value_type rr);

  template<class GIM>
  inline void makeGrid(const GIM& mgrid)
  {
    if(mgrid.size() < GridLevel)
    {
      Grid[SPATIAL_GRID] = mgrid.back();
    }
    else
    {
      for(int ig=0; ig<GridLevel; ig++)
        Grid[ig] = mgrid[ig];
    }
    ///first build the spatial grid
    setGrid(Grid[SPATIAL_GRID]);
  }

  inline Grid_t*
  addGrid(int glevel, const SingleParticleIndex_t& agrid)
  {
    if(SuperGrid[glevel])
    {
      SuperGrid[glevel]->setGrid(agrid);
      return SuperGrid[glevel];
    }
    else
    {
      Grid_t* g = new Grid_t(agrid);
      SuperGrid[glevel] = g;
      return g;
    }
  }

  inline Grid_t* getGrid(int glevel)
  {
    return SuperGrid[glevel];
  }
  inline const Grid_t* getGrid(int glevel) const
  {
    return SuperGrid[glevel];
  }
  void initGrid(Grid_t& agrid, const Grid_t& subgrid);

  void checkGrid(value_type int_rad);

  inline int ngrid(int ig) const
  {
    return Grid[SPATIAL_GRID][ig];
  }
  inline int ngrid(int glevel, int ig) const
  {
    return Grid[glevel][ig];
  }
  void print(std::ostream& os, int level=2) const;

  inline void update()
  {
    for(int i=0; i<u_bc.size(); i++)
      c_bc[i]=toCart(u_bc[i]);
  }

  ///Set LR_rc = radius of smallest sphere inside box and kc=dim/rc
  void SetLRCutoffs();

  void printCutoffs();

  /** set the lattice vector with a tensor
   *@param lat a tensor representing a supercell
   *
   *In addition to setting the internal parameters,
   *it updates the Cartesian displacement vectors.
   */
  inline void update(const Tensor<value_type,3>& lat)
  {
    Base_t::set(lat);
    connectGrid(InteractionRadius,ConnectionRadius);
  }

  ///Uniform grid partitions at MPI_GRID, OPENMP_GRID, and SPATIAL_GRID levels.
  SingleParticleIndex_t Grid[GridLevel];

  ///integer to define a status
  unsigned int Status;

  ///The maximum number of connections per cell
  int MaxConnections;

  ///The interaction radius
  value_type InteractionRadius;

  ///The radius to define a connected cell
  value_type ConnectionRadius;

  ///Dimensionless cutoff radius for G/R breakups
  value_type LR_dim_cutoff;
  value_type LR_rc;
  value_type LR_kc;

  SingleParticleIndex_t NCMax;

  ///UniformCartesianGrid for multi levels.
  std::vector<Grid_t*> SuperGrid;

  ///offsets to determine cell conection
  std::vector<int> c_offset, c_max;

  ///cell index connected to each cell
  std::vector<int> c_id;

  ///displacement vectors due to the boundary conditions
  std::vector<SingleParticlePos_t> c_bc;

  ///displacement vectors due to the boundary conditions in the unit-cell unit
  std::vector<SingleParticlePos_t> u_bc;

  //Possible that this is more efficient
  //struct ConnectedCell {
  //  int ID;
  //  int Size[2];
  //  SingleParticlePos_t u_bc;
  //  SingleParticlePos_t c_bc;
  //};
  //typedef std::vector<ConnectedCell> CCContainer_t;
  //std::vector<CCContainer_t*> Connections;
};
}
#endif

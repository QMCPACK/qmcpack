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
    
    



#ifndef OHMMS_MULTIGRID_PARTICLE_LAYOUT_H
#define OHMMS_MULTIGRID_PARTICLE_LAYOUT_H
#include <string>
#include <vector>

#include "Lattice/Region.h"
#include "Lattice/UniformCartesianGrid.h"

//forward declaration
template<class T, unsigned D> class MakeMultiGrid;

/*! \class MultiGridParticleLayout<class T, unsigned D>
 *  \brief A class for grid layout
 */
template<class T, unsigned D>
struct MultiGridParticleLayout
{

  typedef UniformCartesianGrid<T,D> PtclGrid_t;

  std::vector<PtclGrid_t*> dGrid;

  MultiGridParticleLayout() { }
  ~MultiGridParticleLayout()
  {
    for(int i=0; i<dGrid.size(); i++)
      delete dGrid[i];
  }

  inline const PtclGrid_t* getGrid(int level) const
  {
    return dGrid[level]; // return a grid at level
  }


  /*! \fn makeGrid(int *grid, int igrid=-1)
   *  \param grid[D] a grid partition
   *  \param igrid index of an existing grid to make a sub partition
   *  \return the index of a new grid
   *  \note if(igrid == -1) a new grid is created
   *        else  a new subgrid based on igrid-th grid is created
   */
  int makeGrid(int *grid, int igrid=-1)
  {
    int n=dGrid.size();
    dGrid.push_back(new PtclGrid_t);
    if(igrid >= 0  && n>0)
      MakeMultiGrid<T,D>::refineGrid(*dGrid[n],*dGrid[igrid],grid);
    else
      MakeMultiGrid<T,D>::makeGrid(*dGrid[n],grid);
    return n;
  }

  void resetGrid(int igrid, int* grid)
  {
    if(igrid == dGrid.size())
      //need to create a new grid
    {
      dGrid.push_back(new PtclGrid_t);
    }
    if(igrid > 0)
    {
      MakeMultiGrid<T,D>::refineGrid(*dGrid[igrid],*dGrid[igrid-1],grid);
    }
    else
    {
      MakeMultiGrid<T,D>::makeGrid(*dGrid[igrid],grid);
    }
  }

  void printGrid(std::ostream& os)
  {
    for(int i=0; i<dGrid.size(); i++)
      dGrid[i]->print(os);
  }

  template<class P>  void update(P*, int);
};


template<class T, unsigned D>
struct MakeMultiGrid {};

template<class T>
struct MakeMultiGrid<T,3>
{

  typedef UniformCartesianGrid<T,3> ThisGrid_t;

  static void makeGrid(ThisGrid_t& grid, int* ng)
  {
    int ngtot=ng[0]*ng[1]*ng[2];
    grid.set(ng);
    grid.NodeID.resize(ngtot);
    grid.NodeDist.resize(ngtot+1);
    grid.Node.resize(ngtot);
    for(int i=0; i<=ngtot; i++)
      grid.NodeDist[i] = i;
    T ri[3];
    int ig=0;
    for(int ic=0; ic<grid.NP[0]; ic++)
    {
      ri[0] = grid.Delta[0]*static_cast<T>(ic);
      for(int jc=0; jc<grid.NP[1]; jc++)
      {
        ri[1] = grid.Delta[1]*static_cast<T>(jc);
        for(int kc=0; kc<grid.NP[2]; kc++)
        {
          ri[2] = grid.Delta[2]*static_cast<T>(kc);
          grid.NodeID[grid.key(ic,jc,kc)] = ig;
          grid.Node[ig].set(ri,grid.Delta);
          ig++;
        }
      }
    }
  }

  static void
  refineGrid(ThisGrid_t& subgrid, const ThisGrid_t& biggrid, int *ng)
  {
    int ngbig = biggrid.NP[0]*biggrid.NP[1]*biggrid.NP[2];
    int ngsub = ng[0]*ng[1]*ng[2];
    int ngtot = ngbig*ngsub;
    int ngnew[3];
    ngnew[0] = ng[0]*biggrid.NP[0];
    ngnew[1] = ng[1]*biggrid.NP[1];
    ngnew[2] = ng[2]*biggrid.NP[2];
    subgrid.set(ngnew);
    subgrid.NodeID.resize(ngtot);
    subgrid.Node.resize(ngtot);
    T ri[3],orig[3];
    subgrid.NodeDist.resize(ngbig+1);
    subgrid.NodeDist[0] = 0;
    int igrid=0;
    for(int ig=0; ig<biggrid.Node.size(); ig++)
    {
      //origin of a mesh
      orig[0] = biggrid.Node[ig].Ri[0];
      orig[1] = biggrid.Node[ig].Ri[1];
      orig[2] = biggrid.Node[ig].Ri[2];
      //integer indices of a mesh
      subgrid.getcoord(orig,ngnew);
      //ngnew[0]*=ng[0]; ngnew[1]*=ng[1]; ngnew[2]*=ng[2];
      for(int isub=0; isub<ng[0]; isub++)
      {
        ri[0] = subgrid.Delta[0]*static_cast<T>(isub)+orig[0];
        for(int jsub=0; jsub<ng[1]; jsub++)
        {
          ri[1] = subgrid.Delta[1]*static_cast<T>(jsub)+orig[1];
          for(int ksub=0; ksub<ng[2]; ksub++)
          {
            ri[2] = subgrid.Delta[2]*static_cast<T>(ksub)+orig[2];
            subgrid.NodeID[subgrid.key(ngnew[0]+isub,
                                       ngnew[1]+jsub,
                                       ngnew[2]+ksub)] = igrid;
            subgrid.Node[igrid].set(ri,subgrid.Delta);
            igrid++;
          }
        }
      }
      subgrid.NodeDist[ig+1] = subgrid.NodeDist[ig]+ngsub;
    }
  }

//    //Intended for MPI_Cartesian functions but
//  #ifdef USE_MPI
//    static void makeMPIGrid(ThisGrid_t& grid, const int* ng, const int* period) {
//      int ngtot=ng[0]*ng[1]*ng[2];
//      grid.NodeID.resize(ngtot);
//      grid.NodeDist.resize(ngtot+1);
//      grid.Node.resize(ngtot);
//      for(int i=0; i<=ngtot; i++) grid.NodeDist[i] = i;
//      gird.set(ng);
//      Communicate* Comm = CommCreate::get();
//      int numnode = Comm->getNumNodes();
//      int nodeid  = Comm->getNodeID(),newnode;
//      bool reorder = true;
//      MPI_Comm cartcomm;
//      int MyCoord[3];
//      // remapping a MPI grid using MPI_Cart
//      MPI_Cart_create(MPI_COMM_WORLD,3,grid.NP,period,reorder,&cartcomm);
//      MPI_Cart_map(cartcomm,3,grid.NP,period,&newnode);
//      MPI_Cart_get(cartcomm,3,grid.NP,period,MyCoord);
//      Comm->setNodeID(newnode);
//      Comm->setCommID(cartcomm);
//      int MyID = newnode;
//      int coord[3],ii=0;
//      double ri[3];
//      for(int i=0;i<ng[0]; i++) {
//        ri[0] = grid.Delta[0]*static_cast<double>(i);
//        for(int j=0; j<ng[1]; j++) {
//          ri[1] = grid.Delta[1]*static_cast<double>(j);
//          for(int k=0; k<ng[2]; k++) {
//            coord[0] = i; coord[1] = j; coord[2] = k;
//            MPI_Cart_rank(cartcomm,coord,&nodeid);
//            ri[2] = grid.Delta[2]*static_cast<double>(k);
//            grid.Node[nodeid].set(ri,grid.Delta);
//            grid.NodeID[grid.key(i,j,k)] = nodeid;
//  	}
//        }
//      }
//    }
//  #endif
};

template<class T, unsigned D>
template<class P>
void MultiGridParticleLayout<T,D>::update(P* ptcl, int imode)
{
//  std::cout << "Calling MultiGridParticleLayout::update" << std::endl;
//    typedef P::ParticlePos_t ParticlePos_t;
//    typedef P::ParticleIndex_t ParticleIndex_t;
//    ParticlePos_t r;
//    ParticlePos_t id;
//    r.resize(ptcl->getLocalNum());
//    id.resize(ptcl->getLocalNum());
//    r = ptcl->R;
//    id = ptcl->IonID;
}

#endif

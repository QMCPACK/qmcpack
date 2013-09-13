//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SPACEGRID_H
#define QMCPLUSPLUS_SPACEGRID_H

#include <Configuration.h>
#include <OhmmsPETE/Tensor.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <Utilities/PooledData.h>
#include <QMCHamiltonians/observable_helper.h>
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{


class SpaceGrid : public QMCTraits, public PtclOnLatticeTraits
{
public:
  typedef TinyVector<RealType,DIM> Point;
  typedef PooledData<RealType>     BufferType;
  typedef Matrix<RealType>         Matrix_t;

  SpaceGrid(int& nvalues);
  bool put(xmlNodePtr cur, map<string,Point>& points,ParticlePos_t& R,vector<RealType>& Z,int ndp,bool abort_on_fail=true)
  {
    Rptcl = &R;
    Zptcl = &Z;
    ndparticles = ndp;
    return put(cur,points,abort_on_fail);
  }
  bool put(xmlNodePtr cur, map<string,Point>& points,bool abort_on_fail=true);
  bool initialize_rectilinear(xmlNodePtr cur, string& coord, map<string,Point>& points);
  bool initialize_voronoi(map<string,Point>& points);
  void write_description(ostream& os, string& indent);
  int allocate_buffer_space(BufferType& buf);
  void registerCollectables(vector<observable_helper*>& h5desc,
                            hid_t gid, int grid_index) const;
  void evaluate(const ParticlePos_t& R, const Matrix<RealType>& values,
                BufferType& buf,vector<bool>& particles_outside,
                const DistanceTableData& dtab);

  bool check_grid(void);
  inline int nDomains(void)
  {
    return ndomains;
  }

  void sum(const BufferType& buf,RealType* vals);

  int buffer_start;
  int buffer_end;

  //private:

  //same for all spacegrids
  enum {cartesian=0,cylindrical,spherical,voronoi,ncoordinates} coordinate;
  int buffer_offset;
  int ndomains;
  int nvalues_per_domain;
  Matrix<RealType> domain_volumes;
  Matrix<RealType> domain_centers;

  //in use if sorting by particle count
  bool chempot;
  int npmin,npmax;
  int npvalues;
  Matrix<RealType> cellsamples;
  enum {vacuum,neutral,noref} reference;
  vector<int> reference_count;

  //really only used for cartesian-like grids
  Point origin;
  Tensor<RealType,DIM> axes;
  Tensor<RealType,DIM> axinv;
  RealType volume;
  Matrix<RealType> domain_uwidths;
  string axlabel[DIM];
  vector<int> gmap[DIM];
  RealType odu[DIM];
  RealType umin[DIM];
  RealType umax[DIM];
  int dimensions[DIM];
  int dm[DIM];

  //voronoi grids
  ParticlePos_t* Rptcl;
  vector<RealType>* Zptcl;
  struct irpair
  {
    RealType r;
    int i;
  };
  vector<irpair> nearcell;
  int ndparticles;

  //used only in evaluate
  Point u,ub;
};


}

#endif

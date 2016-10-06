//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  bool put(xmlNodePtr cur, std::map<std::string,Point>& points,ParticlePos_t& R,std::vector<RealType>& Z,int ndp,bool abort_on_fail=true)
  {
    Rptcl = &R;
    Zptcl = &Z;
    ndparticles = ndp;
    return put(cur,points,abort_on_fail);
  }
  bool put(xmlNodePtr cur, std::map<std::string,Point>& points,bool abort_on_fail=true);
  bool initialize_rectilinear(xmlNodePtr cur, std::string& coord, std::map<std::string,Point>& points);
  bool initialize_voronoi(std::map<std::string,Point>& points);
  void write_description(std::ostream& os, std::string& indent);
  int allocate_buffer_space(BufferType& buf);
  void registerCollectables(std::vector<observable_helper*>& h5desc,
                            hid_t gid, int grid_index) const;
  void evaluate(const ParticlePos_t& R, const Matrix<RealType>& values,
                BufferType& buf,std::vector<bool>& particles_outside,
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
  std::vector<int> reference_count;

  //really only used for cartesian-like grids
  Point origin;
  Tensor<RealType,DIM> axes;
  Tensor<RealType,DIM> axinv;
  RealType volume;
  Matrix<RealType> domain_uwidths;
  std::string axlabel[DIM];
  std::vector<int> gmap[DIM];
  RealType odu[DIM];
  RealType umin[DIM];
  RealType umax[DIM];
  int dimensions[DIM];
  int dm[DIM];

  //voronoi grids
  ParticlePos_t* Rptcl;
  std::vector<RealType>* Zptcl;
  struct irpair
  {
    RealType r;
    int i;
  };
  std::vector<irpair> nearcell;
  int ndparticles;

  //used only in evaluate
  Point u,ub;
};


}

#endif

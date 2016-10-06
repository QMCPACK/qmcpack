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
    
    
#ifndef QMCPLUSPLUS_REFERENCE_POINTS_H
#define QMCPLUSPLUS_REFERENCE_POINTS_H

#include <Configuration.h>
#include <OhmmsData/OhmmsElementBase.h>
#include <Particle/ParticleSet.h>
#include <QMCHamiltonians/observable_helper.h>
#include <OhmmsPETE/Tensor.h>

namespace qmcplusplus
{


class ReferencePoints: public QMCTraits
{
public:
  typedef TinyVector<RealType,DIM> Point;
  typedef Tensor<RealType,DIM>     Tensor_t;

  std::map<std::string,Point> points;
  Tensor_t axes;

  bool put(xmlNodePtr cur, ParticleSet& P, std::vector<ParticleSet*>& Pref);
  bool put(ParticleSet& P, std::vector<ParticleSet*>& Pref);
  void write_description(std::ostream& os, std::string& indent);
  void save(std::vector<observable_helper*>& h5desc, hid_t gid) const;

private:
  enum Coordinate {cellC=0,cartesianC,ndirections,nodir};
  Coordinate coordinate;
};



}


#endif

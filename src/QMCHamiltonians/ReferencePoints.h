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

  map<string,Point> points;
  Tensor_t axes;

  bool put(xmlNodePtr cur, ParticleSet& P, vector<ParticleSet*>& Pref);
  bool put(ParticleSet& P, vector<ParticleSet*>& Pref);
  void write_description(ostream& os, string& indent);
  void save(vector<observable_helper*>& h5desc, hid_t gid) const;

private:
  enum Coordinate {cellC=0,cartesianC,ndirections,nodir};
  Coordinate coordinate;
};



}


#endif

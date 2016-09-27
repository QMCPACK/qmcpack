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
    
    

#ifndef QMCPLUSPLUS_STATIC_STRUCTURE_FACTOR_H
#define QMCPLUSPLUS_STATIC_STRUCTURE_FACTOR_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>

namespace qmcplusplus
{

class StaticStructureFactor : public QMCHamiltonianBase
{
 public:
 
  typedef std::vector<RealType> k2_t;
  typedef std::vector<RealType> dens_t;
  typedef std::vector<PosType>  pts_t;

  //data members
  int nspecies;
  std::vector<std::string> species_name;
  RealType ecut;
  int nkpoints;
  const ParticleSet& Pinit;

  //constructor/destructor
  StaticStructureFactor(ParticleSet& P);
  ~StaticStructureFactor() { }

  //standard interface
  QMCHamiltonianBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);
  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P); 
  }

  //required for Collectables interface
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P)                      { }
  void setObservables(PropertySetType& plist)                      { }
  void setParticlePropertyList(PropertySetType& plist, int offset) { }

#if !defined(REMOVE_TRACEMANAGER)
  void checkout_scalar_arrays(TraceManager& tm)                    { }
  void collect_scalar_samples()                                    { }
  void delete_scalar_arrays()                                      { }
#endif

  //obsolete?
  bool get(std::ostream& os) const { return false; }

  //local functions
  void reset();
  void report(const std::string& pad);
  void postprocess_density(const std::string& infile,const std::string& species,
                           pts_t& points,dens_t& density,dens_t& density_err);

};

}

#endif

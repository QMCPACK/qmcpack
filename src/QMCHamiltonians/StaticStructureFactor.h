//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_STATIC_STRUCTURE_FACTOR_H
#define QMCPLUSPLUS_STATIC_STRUCTURE_FACTOR_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>

namespace qmcplusplus
{

class StaticStructureFactor : public QMCHamiltonianBase
{
 public:
 
  typedef vector<RealType> k2_t;
  typedef vector<RealType> dens_t;
  typedef vector<PosType>  pts_t;

  //data members
  int nspecies;
  vector<string> species_name;
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
  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P); 
  }

  //required for Collectables interface
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P)                      { }
  void setObservables(PropertySetType& plist)                      { }
  void setParticlePropertyList(PropertySetType& plist, int offset) { }
  void checkout_scalar_arrays(TraceManager& tm)                    { }
  void collect_scalar_samples()                                    { }
  void delete_scalar_arrays()                                      { }

  //obsolete?
  bool get(std::ostream& os) const { return false; }

  //local functions
  void reset();
  void report(const string& pad);
  void postprocess_density(const string& infile,const string& species,
                           pts_t& points,dens_t& density,dens_t& density_err);

};

}

#endif

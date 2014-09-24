//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPIN_DENSITY_H
#define QMCPLUSPLUS_SPIN_DENSITY_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus
{

class SpinDensity : public QMCHamiltonianBase
{
 public:

  typedef map<string,ParticleSet*> PSPool;
  typedef ParticleSet::ParticleLayout_t Lattice_t;
  typedef DistanceTableData::ripair ripair;
  typedef vector<RealType> dens_t;
  typedef vector<PosType>  pts_t;
  typedef ParticleSet::ParticlePos_t ParticlePos_t;

  ParticleSet* Ptmp;
  PSPool& psetpool;
  ParticlePos_t Robs;

  //data members
  int nspecies;
  vector<int>    species_size;
  vector<string> species_name;
  Lattice_t cell;
  PosType   corner;
  TinyVector<int,DIM> grid;
  TinyVector<int,DIM> gdims;
  int npoints;
  bool voronoi_grid;
  int dtable_index;
  vector<ripair> nearest;

  //constructor/destructor
  SpinDensity(ParticleSet& P,PSPool& psp);
  ~SpinDensity() { }

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
  void evaluate_voronoi(ParticleSet& P);
  void evaluate_cartesian(ParticleSet& P);
  void test(int moves,ParticleSet& P);
  Return_t test_evaluate(ParticleSet& P,int& pmin,int& pmax);
  void postprocess_density(const string& infile,const string& species,
                           pts_t& points,dens_t& density,dens_t& density_err);

};

}

#endif

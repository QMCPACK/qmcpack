//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_ENERGY_DENSITY_ESTIMATOR_H
#define QMCPLUSPLUS_ENERGY_DENSITY_ESTIMATOR_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <QMCHamiltonians/ReferencePoints.h>
#include <QMCHamiltonians/SpaceGrid.h>
#include <map>
#include <vector>

namespace qmcplusplus
{

class EnergyDensityEstimator: public QMCHamiltonianBase, public PtclOnLatticeTraits
{
public:

  typedef ReferencePoints::Point Point;
  typedef map<string,ParticleSet*> PSPool;

  EnergyDensityEstimator(PSPool& PSP, const string& defaultKE);
  ~EnergyDensityEstimator();

  void resetTargetParticleSet(ParticleSet& P);
  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }
  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool put(xmlNodePtr cur, ParticleSet& P);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void write_description(ostream& os);

  void write_EDValues(void);
  void write_nonzero_domains(const ParticleSet& P);
  void write_Collectables(string& label,int& cnt,ParticleSet& P);

private:

  //original xml
  xmlNodePtr input_xml;
  //system information
  string defKE;
  PSPool& psetpool;
  ParticleSet* Pdynamic;
  ParticleSet* Pstatic;
  ParticleSet* get_particleset(string& psname);
  int dtable_index;
  int nparticles;
  //collection of points from which to build spacegrid origin and axes
  ReferencePoints ref;
  //EnergyDenstity quantities
  enum {W=0,T,V,nEDValues};
  Matrix<RealType> EDValues;
  //for EnergyDensity of particles falling outside any spacegrid
  int outside_buffer_offset;
  vector<bool> particles_outside;
  //spacegrids are used to find which cell domain
  //  contains the Energy information of particles
  vector<SpaceGrid*> spacegrids;
  //particle positions
  ParticlePos_t R;
  //number of samples accumulated
  int nsamples;

  //needed (temporarily) for chempot
  //ParticleSet should carry Zptcl so it doesn't have
  // to be computed everywhere from species
  vector<RealType> Zptcl;
  ParticlePos_t Rptcl;
  void set_ptcl(void);
  void unset_ptcl(void);

  TraceSample<TraceReal>*         w_trace;
  TraceSample<TraceReal>*         Td_trace;
  CombinedTraceSample<TraceReal>* Vd_trace;
  CombinedTraceSample<TraceReal>* Vs_trace;

  virtual void get_required_traces(TraceManager& tm);

  virtual void contribute_scalar_quantities()               { }
  virtual void checkout_scalar_quantities(TraceManager& tm) { }
  virtual void collect_scalar_quantities()                  { }
  virtual void delete_scalar_quantities()                   { }
};


}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 1 $   $Date: 2011-02-02 16:21:33 -0500 (Wed, 02 Feb 2011) $
 * $Id: EnergyDensityEstimator.h 1 2011-02-02 16:21:33Z jtkrogel $
 ***************************************************************************/

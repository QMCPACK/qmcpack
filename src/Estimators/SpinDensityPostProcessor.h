//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPIN_DENSITY_POST_PROCESSOR_H
#define QMCPLUSPLUS_SPIN_DENSITY_POST_PROCESSOR_H

#include <Estimators/PostProcessorBase.h>

namespace qmcplusplus
{

class QMCHamitonian;

struct SpinDensityPostProcessor : public PostProcessorBase
{
  typedef vector<RealType> dens_t;
  typedef vector<PosType>  pts_t;

  ParticleSet&    Pq;
  ParticleSet&    Pc;
  QMCHamiltonian& H;

  RealType dV;
  int nspecies;
  vector<string> species_name;
  vector<int>    species_size;

  vector<string> sources;
  string format;
  string format_ext;
  Lattice_t cell;
  PosType   corner;
  TinyVector<int,DIM> grid;
  pts_t gridpoints;
  int npoints;
  string normalization;

  SpinDensityPostProcessor(ParticleSet& pq,ParticleSet& pc,QMCHamiltonian& h);

  ~SpinDensityPostProcessor() { }

  void put(xmlNodePtr cur);

  void postprocess();

  template<typename SDO>
  void get_density(const string& infile,const string& species,
                   QMCHamiltonianBase* h,dens_t& density,dens_t& density_err);

  void normalize(int nparticles,dens_t& density,dens_t& density_err);

  void write_density(const string& outfile,dens_t& density);

  void write_density_xsf(const string& outfile,dens_t& density);

  /** get species of particle i from a ParticleSet
   *    wouldn't it be nice if this were in ParticleSet?
   */
  inline const string& pname(ParticleSet& P,int i)
  {
    return P.mySpecies.speciesName[P.GroupID[i]];
  }
};

}

#endif

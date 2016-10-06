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
    
    



#ifndef QMCPLUSPLUS_SPIN_DENSITY_POST_PROCESSOR_H
#define QMCPLUSPLUS_SPIN_DENSITY_POST_PROCESSOR_H

#include <Estimators/PostProcessorBase.h>

namespace qmcplusplus
{

class QMCHamitonian;

struct SpinDensityPostProcessor : public PostProcessorBase
{
  typedef std::vector<RealType> dens_t;
  typedef std::vector<PosType>  pts_t;

  ParticleSet&    Pq;
  ParticleSet&    Pc;
  QMCHamiltonian& H;

  RealType dV;
  int nspecies;
  std::vector<std::string> species_name;
  std::vector<int>    species_size;

  std::vector<std::string> sources;
  std::string format;
  std::string format_ext;
  Lattice_t cell;
  PosType   corner;
  TinyVector<int,DIM> grid;
  pts_t gridpoints;
  int npoints;
  std::string normalization;

  SpinDensityPostProcessor(ParticleSet& pq,ParticleSet& pc,QMCHamiltonian& h);

  ~SpinDensityPostProcessor() { }

  void put(xmlNodePtr cur);

  void postprocess();

  template<typename SDO>
  void get_density(const std::string& infile,const std::string& species,
                   QMCHamiltonianBase* h,dens_t& density,dens_t& density_err);

  void normalize(int nparticles,dens_t& density,dens_t& density_err);

  void write_density(const std::string& outfile,dens_t& density);

  void write_density_xsf(const std::string& outfile,dens_t& density);

  /** get species of particle i from a ParticleSet
   *    wouldn't it be nice if this were in ParticleSet?
   */
  inline const std::string& pname(ParticleSet& P,int i)
  {
    return P.mySpecies.speciesName[P.GroupID[i]];
  }
};

}

#endif

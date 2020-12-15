//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: SpinDensity.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "SpinDensityNew.h"
namespace qmcplusplus
{
SpinDensityNew::SpinDensityNew(SpinDensityInput&& input, const SpeciesSet& species) : input_(input), species_(species)
{
  myName = "SpinDensity";

  // This code is quite suspect.
  // I think it is checking the membersize is either the last attribute or adding it.  If its already there but not
  // not last it fails. Not sure why we care yet.
  int isize = species_.findAttribute("membersize");
  // That there will be a number of particles of a particular species is an invariant
  // but the SpeciesSet fails to say that so we have this
  if (!isize)
    throw std::runtime_error("SpinDensity(P) Species set does not have the required attribute 'membersize'");
  for (int s = 0; s < species_.size(); ++s)
    species_size.push_back(species_(isize, s));

}

OperatorEstBase* SpinDensityNew::makeClone(ParticleSet& P, TrialWaveFunction& Psi) { return new SpinDensityNew(*this); }

QMCTraits::FullPrecRealType SpinDensityNew::evaluate(ParticleSet& P)
{
  QMCT::RealType w = tWalker->Weight;
  int p            = 0;
  int offset       = myIndex;
  for (int s = 0; s < species_.size(); ++s, offset += input_.get_npoints())
    for (int ps = 0; ps < species_.attribName.size(); ++ps, ++p)
    {
      QMCT::PosType u = input_.get_cell().toUnit(P.R[p] - input_.get_corner());
      //bool inside = true;
      //for(int d=0;d<DIM;++d)
      //  inside &= u[d]>0.0 && u[d]<1.0;
      //if(inside)
      //{
      int point = offset;
      for (int d = 0; d < QMCT::DIM; ++d)
        point += input_.get_gdims()[d] * ((int)(input_.get_grid()[d] * (u[d] - std::floor(u[d])))); //periodic only
      P.Collectables[point] += w;
      //}
    }
  return 0.0;
};

void SpinDensityNew::report(const std::string& pad)
{
  app_log() << pad << "SpinDensity report" << std::endl;
  app_log() << pad << "  dim     = " << QMCT::DIM << std::endl;
  app_log() << pad << "  npoints = " << input_.get_npoints() << std::endl;
  app_log() << pad << "  grid    = " << input_.get_grid() << std::endl;
  app_log() << pad << "  gdims   = " << input_.get_gdims() << std::endl;
  app_log() << pad << "  corner  = " << input_.get_corner() << std::endl;
  app_log() << pad << "  center  = " << input_.get_corner() + input_.get_cell().Center << std::endl;
  app_log() << pad << "  cell " << std::endl;
  for (int d = 0; d < QMCT::DIM; ++d)
    app_log() << pad << "    " << d << " " << input_.get_cell().Rv[d] << std::endl;
  app_log() << pad << "  end cell " << std::endl;
  app_log() << pad << "  nspecies = " << species_.size() << std::endl;
  for (int s = 0; s < species_.size(); ++s)
    app_log() << pad << "    species[" << s << "]"
       << " = " << species_.speciesName[s] << " " << species_.attribName.size() << std::endl;
  app_log() << pad << "end SpinDensity report" << std::endl;
}

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "SpinDensity.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{
SpinDensity::SpinDensity(ParticleSet& P)
{
  // get particle information
  SpeciesSet& species = P.getSpeciesSet();
  nspecies            = species.size();
  for (int s = 0; s < nspecies; ++s)
    species_size.push_back(P.groupsize(s));
  for (int s = 0; s < nspecies; ++s)
    species_name.push_back(species.speciesName[s]);
  reset();

  //jtk: spin density only works for periodic bc's for now
  //     abort if using open boundary conditions
  bool open_bcs = (P.getLattice().SuperCellEnum == SUPERCELL_OPEN);
  if (open_bcs)
  {
    APP_ABORT("SpinDensity is not implemented for open boundary conditions at present\n  please contact the developers "
              "if you need this feature");
  }

  Ptmp = &P;
}


void SpinDensity::reset()
{
  name_ = "SpinDensity";
  update_mode_.set(COLLECTABLE, 1);
  corner = 0.0;
}


std::unique_ptr<OperatorBase> SpinDensity::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
{
  return std::make_unique<SpinDensity>(*this);
}


bool SpinDensity::put(xmlNodePtr cur)
{
  using std::ceil;
  using std::sqrt;
  reset();
  std::string write_report = "no";
  OhmmsAttributeSet attrib;
  attrib.add(name_, "name");
  attrib.add(write_report, "report");
  attrib.put(cur);

  bool have_dr     = false;
  bool have_grid   = false;
  bool have_center = false;
  bool have_corner = false;
  bool have_cell   = false;

  PosType dr;
  PosType center;
  Tensor<RealType, DIM> axes;

  int test_moves = 0;

  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename((const char*)element->name);
    if (ename == "parameter")
    {
      const std::string name(getXMLAttributeValue(element, "name"));
      if (name == "dr")
      {
        have_dr = true;
        putContent(dr, element);
      }
      else if (name == "grid")
      {
        have_grid = true;
        putContent(grid, element);
      }
      else if (name == "corner")
      {
        have_corner = true;
        putContent(corner, element);
      }
      else if (name == "center")
      {
        have_center = true;
        putContent(center, element);
      }
      else if (name == "cell")
      {
        have_cell = true;
        putContent(axes, element);
      }
      else if (name == "test_moves")
        putContent(test_moves, element);
    }
    element = element->next;
  }

  if (have_dr && have_grid)
  {
    APP_ABORT("SpinDensity::put  dr and grid are provided, this is ambiguous");
  }
  else if (!have_dr && !have_grid)
    APP_ABORT("SpinDensity::put  must provide dr or grid");

  if (have_corner && have_center)
    APP_ABORT("SpinDensity::put  corner and center are provided, this is ambiguous");
  if (have_cell)
  {
    cell.set(axes);
    if (!have_corner && !have_center)
      APP_ABORT("SpinDensity::put  must provide corner or center");
  }
  else
    cell = Ptmp->getLattice();

  if (have_center)
    corner = center - cell.Center;

  if (have_dr)
    for (int d = 0; d < DIM; ++d)
      grid[d] = (int)ceil(sqrt(dot(cell.Rv[d], cell.Rv[d])) / dr[d]);

  npoints = 1;
  for (int d = 0; d < DIM; ++d)
    npoints *= grid[d];
  gdims[0] = npoints / grid[0];
  for (int d = 1; d < DIM; ++d)
    gdims[d] = gdims[d - 1] / grid[d];

  if (write_report == "yes")
    report("  ");
  if (test_moves > 0)
    test(test_moves, *Ptmp);

  return true;
}


void SpinDensity::report(const std::string& pad)
{
  app_log() << pad << "SpinDensity report" << std::endl;
  app_log() << pad << "  dim     = " << DIM << std::endl;
  app_log() << pad << "  npoints = " << npoints << std::endl;
  app_log() << pad << "  grid    = " << grid << std::endl;
  app_log() << pad << "  gdims   = " << gdims << std::endl;
  app_log() << pad << "  corner  = " << corner << std::endl;
  app_log() << pad << "  center  = " << corner + cell.Center << std::endl;
  app_log() << pad << "  cell " << std::endl;
  for (int d = 0; d < DIM; ++d)
    app_log() << pad << "    " << d << " " << cell.Rv[d] << std::endl;
  app_log() << pad << "  end cell " << std::endl;
  app_log() << pad << "  nspecies = " << nspecies << std::endl;
  for (int s = 0; s < nspecies; ++s)
    app_log() << pad << "    species[" << s << "]"
              << " = " << species_name[s] << " " << species_size[s] << std::endl;
  app_log() << pad << "end SpinDensity report" << std::endl;
}


void SpinDensity::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<RealType> tmp(nspecies * npoints);
  collectables.add(tmp.begin(), tmp.end());
}


void SpinDensity::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  std::vector<int> ng(1);
  ng[0] = npoints;

  hdf_path hdf_name{name_};
  for (int s = 0; s < nspecies; ++s)
  {
    h5desc.emplace_back(hdf_name / species_name[s]);
    auto& oh = h5desc.back();
    oh.set_dimensions(ng, my_index_ + s * npoints);
  }
}


SpinDensity::Return_t SpinDensity::evaluate(ParticleSet& P)
{
  RealType w = t_walker_->Weight;
  int p      = 0;
  int offset = my_index_;
  for (int s = 0; s < nspecies; ++s, offset += npoints)
    for (int ps = 0; ps < species_size[s]; ++ps, ++p)
    {
      PosType u = cell.toUnit(P.R[p] - corner);
      //bool inside = true;
      //for(int d=0;d<DIM;++d)
      //  inside &= u[d]>0.0 && u[d]<1.0;
      //if(inside)
      //{
      int point = offset;
      for (int d = 0; d < DIM; ++d)
        point += gdims[d] * ((int)(grid[d] * (u[d] - std::floor(u[d])))); //periodic only
      P.Collectables[point] += w;
      //}
    }
  return 0.0;
}


void SpinDensity::test(int moves, ParticleSet& P)
{
  app_log() << "  SpinDensity test" << std::endl;
  RandomGenerator rng;
  int particles = P.getTotalNum();
  int pmin      = std::numeric_limits<int>::max();
  int pmax      = std::numeric_limits<int>::min();
  for (int m = 0; m < moves; ++m)
  {
    for (int p = 0; p < particles; ++p)
    {
      PosType u;
      for (int d = 0; d < DIM; ++d)
        u[d] = rng();
      P.R[p] = P.getLattice().toCart(u);
    }
    test_evaluate(P, pmin, pmax);
  }
  app_log() << "  end SpinDensity test" << std::endl;
  APP_ABORT("SpinDensity::test  test complete");
}


SpinDensity::Return_t SpinDensity::test_evaluate(ParticleSet& P, int& pmin, int& pmax)
{
  int p      = 0;
  int offset = 0;
  for (int s = 0; s < nspecies; ++s, offset += npoints)
    for (int ps = 0; ps < species_size[s]; ++ps, ++p)
    {
      PosType u   = cell.toUnit(P.R[p] - corner);
      bool inside = true;
      for (int d = 0; d < DIM; ++d)
        inside &= u[d] > 0.0 && u[d] < 1.0;
      if (inside)
      {
        int point = offset;
        for (int d = 0; d < DIM; ++d)
          point += gdims[d] * ((int)(u[d] * grid[d]));
        pmin = std::min(pmin, point - offset);
        pmax = std::max(pmax, point - offset);
      }
    }
  app_log() << "    pmin = " << pmin << " pmax = " << pmax << " npoints = " << npoints << std::endl;
  return 0.0;
}

} // namespace qmcplusplus

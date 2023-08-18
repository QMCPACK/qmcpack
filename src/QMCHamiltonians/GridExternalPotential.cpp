//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Kevin Ryczko, kryczko@uottawa.ca, University of Ottawa
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "GridExternalPotential.h"
#include "OhmmsData/AttributeSet.h"


namespace qmcplusplus
{
GridExternalPotential::GridExternalPotential(ParticleSet& P) : ps_(P)
{
  setEnergyDomain(POTENTIAL);
  oneBodyQuantumDomain(P);
}

bool GridExternalPotential::put(xmlNodePtr cur)
{
  using std::sqrt;


  OhmmsAttributeSet attrib;

  /*
    Future update must include differing grid sizes in the x, y, and z directions. 

    Everything is the same size for now.
  */

  Ugrid grid;
  grid.start = -2.0;
  grid.end   = 2.0;
  grid.num   = 11;
  BCtype_d BC;
  BC.lCode = NATURAL;
  BC.rCode = NATURAL;
  std::string file_name;
  std::string dataset_name;
  bool pbc;

  attrib.add(grid.start, "start");
  attrib.add(grid.end, "end");
  attrib.add(grid.num, "num");
  attrib.add(file_name, "file_name");
  attrib.add(dataset_name, "dataset_name");
  attrib.add(pbc, "pbc");
  attrib.put(cur);

  double delta = (grid.end - grid.start) / (grid.num - 1);


  if (pbc)
  {
    BC.lCode = PERIODIC;
    BC.rCode = PERIODIC;
    delta    = (grid.end - grid.start) / (grid.num);
  }

  Array<double, 3> data(grid.num, grid.num, grid.num);

  hdf_archive hin;
  bool read_okay = hin.open(file_name, H5F_ACC_RDONLY);
  if (!read_okay)
  {
    app_log() << "Failed to open HDF5 file: " << file_name << "." << std::endl;
  }
  else
  {
    app_log() << "    ==============================\n"
              << "    Information of grid:\n"
              << "    Grid start: " << grid.start << std::endl
              << "    Grid end: " << grid.end << std::endl
              << "    Grid num: " << grid.num << std::endl
              << "    Grid delta: " << delta << std::endl
              << "    Grid file_name: " << file_name << std::endl
              << "    Grid dataset_name: " << dataset_name << std::endl
              << "    Periodic: " << pbc << std::endl
              << "    ==============================\n";
  }

  hin.read(data, dataset_name);

  spline_data_.reset(create_UBspline_3d_d(grid, grid, grid, BC, BC, BC, data.data()), destroy_Bspline);

  return true;
}


bool GridExternalPotential::get(std::ostream& os) const
{
  os << "External grid potential" << std::endl;
  return true;
}


std::unique_ptr<OperatorBase> GridExternalPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<GridExternalPotential>(*this);
}


GridExternalPotential::Return_t GridExternalPotential::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles_)
    value_ = evaluate_sp(P);
  else
  {
#endif
    value_ = 0.0;
    for (int i = 0; i < P.getTotalNum(); ++i)
    {
      PosType r = P.R[i];
      P.getLattice().applyMinimumImage(r);
      double val = 0.0;
      eval_UBspline_3d_d(spline_data_.get(), r[0], r[1], r[2], &val);

      value_ += val;
    }
#if !defined(REMOVE_TRACEMANAGER)
  }
#endif
  return value_;
}

GridExternalPotential::Return_t GridExternalPotential::evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
{
  return evaluate(P);
}

#if !defined(REMOVE_TRACEMANAGER)
void GridExternalPotential::contributeParticleQuantities() { request_.contribute_array(name_); }

void GridExternalPotential::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_);
  if (streaming_particles_)
    v_sample_ = tm.checkout_real<1>(name_, ps_);
}

void GridExternalPotential::deleteParticleQuantities()
{
  if (streaming_particles_)
    delete v_sample_;
}

GridExternalPotential::Return_t GridExternalPotential::evaluate_sp(ParticleSet& P)
{
  Array<TraceReal, 1>& V_samp = *v_sample_;
  value_                      = 0.0;
  for (int i = 0; i < P.getTotalNum(); ++i)
  {
    PosType r   = P.R[i];
    RealType v1 = dot(r, r);
    V_samp(i)   = v1;
    value_ += v1;
  }
  return value_;
}
#endif

} // namespace qmcplusplus

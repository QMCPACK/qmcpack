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


#include <QMCHamiltonians/GridExternalPotential.h>
#include <OhmmsData/AttributeSet.h>


namespace qmcplusplus
{
bool GridExternalPotential::put(xmlNodePtr cur)
{
  using std::sqrt;


  OhmmsAttributeSet attrib;
  attrib.put(cur);

  Ugrid grid;
  grid.start = -2.0;
  grid.end = 2.0;
  grid.num = 11;
  BCtype_d BC;
  BC.lCode = NATURAL;
  BC.rCode = NATURAL;

  double delta = (grid.end - grid.start)/(grid.num-1);

  Array<double, 3> data(grid.num, grid.num, grid.num);

  hdf_archive hin;
  bool read_okay = hin.open("sho.h5",H5F_ACC_RDONLY);
  if (!read_okay) {
    app_log() << "Failed to open sho.h5" << std::endl;
  }

  hin.read(data, "pot_data");
  

  spline_data = create_UBspline_3d_d(grid, grid, grid, 
                                     BC, BC, BC, data.data());

  return true;
}


bool GridExternalPotential::get(std::ostream& os) const
{
  os << "External grid potential" << std::endl;
  return true;
}


QMCHamiltonianBase* GridExternalPotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return new GridExternalPotential(*this);
}


GridExternalPotential::Return_t GridExternalPotential::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles)
    Value = evaluate_sp(P);
  else
  {
#endif
    Value              = 0.0;
    for (int i = 0; i < P.getTotalNum(); ++i)
    {
      PosType r = P.R[i];
      double val = 0.0;
      eval_UBspline_3d_d(spline_data, r[0], r[1], r[2], &val);
      Value += val;
    }
#if !defined(REMOVE_TRACEMANAGER)
  }
#endif
  return Value;
}


#if !defined(REMOVE_TRACEMANAGER)
GridExternalPotential::Return_t GridExternalPotential::evaluate_sp(ParticleSet& P)
{
  Array<TraceReal, 1>& V_samp = *V_sample;
  Value                       = 0.0;
  for (int i = 0; i < P.getTotalNum(); ++i)
  {
    PosType r   = P.R[i];
    RealType v1 = dot(r, r);
    V_samp(i)   = v1;
    Value += v1;
  }
  return Value;
}
#endif

} // namespace qmcplusplus

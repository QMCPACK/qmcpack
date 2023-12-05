//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCHamiltonians/ECPComponentBuilder.h"
namespace qmcplusplus
{
void ECPComponentBuilder::buildL2(xmlNodePtr cur)
{
  // report to log and check for global grid (pseudo/grid)
  app_log() << "    ECPComponentBuilder::buildL2 " << std::endl;
  if (grid_global == 0)
  {
    app_error() << "    Global grid needs to be defined." << std::endl;
    myComm->barrier_and_abort("ECPComponentBuilder::buildL2");
  }

  // read <L2> attributes
  std::string eunits = "hartree";
  std::string format = "r*V";
  RealType rcut      = -1.0;
  OhmmsAttributeSet attrib;
  attrib.add(eunits, "units");
  attrib.add(format, "format");
  attrib.add(rcut, "cutoff");
  attrib.put(cur);

  // check validity of <L2> attributes
  RealType Vprefactor = 1.0;
  if (eunits.find("ydberg") < eunits.size())
  {
    app_log() << "    Input energy unit = Rydberg " << std::endl;
    Vprefactor = 0.5;
  }
  else
  {
    app_log() << "    Assuming Hartree unit" << std::endl;
  }
  bool is_r_times_V(true);
  if (format == "r*V")
    is_r_times_V = true;
  else if (format == "V")
    is_r_times_V = false;
  else
  {
    app_error() << "  Unrecognized format \"" << format << "\" in PP file." << std::endl;
    myComm->barrier_and_abort("ECPComponentBuilder::buildL2");
  }
  if (rcut < 0.0)
  {
    app_error() << "  L2 attribute cutoff is missing or negative.  Cutoff must be a positive real number." << std::endl;
    myComm->barrier_and_abort("ECPComponentBuilder::buildL2");
  }
  int npts = grid_global->size();

  // read in the data
  std::vector<RealType> vL2in(npts);
  xmlNodePtr c = cur->children;
  while (c != NULL)
  {
    if (xmlStrEqual(c->name, (const xmlChar*)"radfunc"))
    {
      xmlNodePtr c1 = c->children;
      while (c1 != NULL)
      {
        if (xmlStrEqual(c1->name, (const xmlChar*)"data"))
          putContent(vL2in, c1);
        c1 = c1->next;
      }
    }
    c = c->next;
  }

  // convert to r*V if needed
  if (!is_r_times_V)
  {
    app_log() << "  Input pseudopotential is converted into r*V" << std::endl;
    for (int i = 0; i < npts; i++)
      vL2in[i] *= grid_global->r(i);
  }

  // create a new grid with maximum extent set to the cutoff radius
  RealType rmax        = rcut;
  const int max_points = 100000;
  app_log() << "   Creating a Linear Grid Rmax=" << rmax << std::endl;
  auto grid  = std::make_unique<LinearGrid<RealType>>();
  RealType d = 1e-4;
  int ng;
  if (grid_global->getGridTag() == LINEAR_1DGRID)
  {
    ng = (int)std::ceil(rmax * grid_global->DeltaInv) + 1;
    if (ng <= max_points)
    {
      app_log() << "  Using global grid with delta = " << grid_global->Delta << std::endl;
      rmax = grid_global->Delta * (ng - 1);
      grid->set(0.0, rmax, ng);
    }
    else
      grid->set(0.0, rmax, max_points);
  }
  else
  {
    ng = std::min(max_points, static_cast<int>(rmax / d) + 1);
    grid->set(0, rmax, ng);
  }
  d = grid->Delta;

  // follow doBreakUp and reinterpolate the data
  //   likely necessary for data inputted on other than a linear grid
  //   but seems extraneous for data on linear grid (most common case)
  int ngIn = npts - 2;
  std::vector<RealType> new_vL2(ng);
  std::vector<RealType> new_vL2in(ngIn);
  for (int i = 0; i < ngIn; i++)
    new_vL2in[i] = Vprefactor * vL2in[i];
  OneDimCubicSpline<mRealType> infunc(grid_global->makeClone(), new_vL2in);
  infunc.spline(0, 0.0, ngIn - 1, 0.0);
  for (int i = 1; i < ng - 1; i++)
  {
    RealType r = d * i;
    new_vL2[i] = infunc.splint(r) / r;
  }
  new_vL2[0]      = new_vL2[1];
  new_vL2[ng - 1] = 0.0;
  auto vL2        = std::make_unique<RadialPotentialType>(std::move(grid), new_vL2);
  vL2->spline();

  // save the splined L2 potential
  pp_L2       = std::make_unique<L2RadialPotential>();
  pp_L2->vL2  = std::move(vL2);
  pp_L2->rcut = rcut;

  //app_log()<<std::endl;
  //for(int i=0;i<ng;i++)
  //  app_log()<<d*i<<" "<<pp_L2->vL2->splint(d*i)*d*i<<std::endl;
  //app_log()<<std::endl;
  //app_log()<<pp_L2->rcut;
  //APP_ABORT("here");
}

} // namespace qmcplusplus

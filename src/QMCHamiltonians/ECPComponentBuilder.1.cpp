//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"

namespace qmcplusplus
{
void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur)
{
  std::unique_ptr<mGridType> grid_semilocal;
  RealType rmax = pp_nonloc->Rmax;
  cur           = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "grid")
    {
      grid_semilocal = createGrid(cur);
      rmax           = grid_semilocal->rmax();
    }
    else if (cname == "vps")
    {
      //should be able to overwrite rmax
      int l           = angMon[getXMLAttributeValue(cur, "l")];
      Lmax            = std::max(l, Lmax);
      xmlNodePtr cur1 = cur->children;
      while (cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if (cname1 == "basisGroup")
        {
          pp_nonloc->add(l, createVrWithBasisGroup(cur1, grid_semilocal.get()));
        }
        cur1 = cur1->next;
      }
      NumNonLocal++;
    }
    cur = cur->next;
  }
  pp_nonloc->lmax = Lmax;
  pp_nonloc->Rmax = rmax;
}

ECPComponentBuilder::RadialPotentialType* ECPComponentBuilder::createVrWithBasisGroup(xmlNodePtr cur, mGridType* agrid)
{
  //todo rcut should be reset if necessary
  using InFuncType = GaussianTimesRN<RealType>;
  InFuncType a;
  a.putBasisGroup(cur);
  bool ignore        = true;
  const RealType eps = 1e-4;
  mRealType rout     = agrid->rmax() * 2;
  while (ignore && rout > agrid->rmax())
  {
    ignore = (std::abs(a.f(rout)) < eps);
    rout -= 0.01;
  }
  rout += 0.01;
  app_log() << "  cutoff for non-local pseudopotential = " << agrid->rmax() << std::endl;
  app_log() << "  calculated cutoff for " << eps << " = " << rout << std::endl;
  int ng = agrid->size();
  //all ways reset grid. Ye Luo
  std::unique_ptr<GridType> agrid_new;
  //if(agrid->GridTag != LINEAR_1DGRID)// || rout > agrid->rmax())
  {
    RealType ri = 0.0;
    RealType rf = std::max(rout, agrid->rmax());
    agrid_new   = std::make_unique<LinearGrid<RealType>>();
    agrid_new->set(ri, rf, ng);
    app_log() << "  Reset the grid for SemiLocal component to a LinearGrid. " << std::endl;
  }
  std::vector<RealType> v(ng);
  for (int ig = 0; ig < ng; ig++)
  {
    v[ig] = a.f((*agrid_new)[ig]);
  }
  v[ng - 1]                = 0.0;
  RadialPotentialType* app = new RadialPotentialType(std::move(agrid_new), v);
  app->spline();
  return app;
}

/** build a Local Pseudopotential
 *
 * For numerical stability, the radial function of the local pseudopotential is stored
 * as \f$rV(r)/Z_{eff}/Z_{e}\f$ for the distance r and \f$Z_e = -1\f$.
 * The coulomb factor $Z_{e}Z_{eff}/r$ is applied by LocalECPotential
 * (see LocalECPotential::evaluate).
 */
void ECPComponentBuilder::buildLocal(xmlNodePtr cur)
{
  if (pp_loc)
    return; //something is wrong

  std::string vFormat("V");
  const std::string v_str(getXMLAttributeValue(cur, "format"));
  if (!v_str.empty())
    vFormat = v_str;

  int vPowerCorrection = 1;
  if (vFormat == "r*V")
  {
    app_log() << "  Local pseudopotential format = r*V" << std::endl;
    vPowerCorrection = 0;
  }
  else
  {
    app_log() << "  Local pseudopotential format = V" << std::endl;
  }
  using InFuncType = GaussianTimesRN<RealType>;
  std::unique_ptr<GridType> grid_local;
  std::unique_ptr<mGridType> grid_local_inp;
  InFuncType vr;
  bool bareCoulomb = true;
  cur              = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "grid")
    {
      grid_local_inp = createGrid(cur, true);
    }
    else if (cname == "basisGroup")
    {
      vr.putBasisGroup(cur, vPowerCorrection);
      bareCoulomb = false;
    }
    cur = cur->next;
  }
  if (grid_local_inp == nullptr)
  {
    if (grid_global == nullptr)
      myComm->barrier_and_abort("ECPComponentBuilder::buildLocal Missing grid information. ");
    grid_local = std::make_unique<LinearGrid<RealType>>();
    grid_local->set(grid_global->rmin(), grid_global->rmax(), grid_global->size());
  }
  else
  {
    grid_local = std::make_unique<LinearGrid<RealType>>();
    grid_local->set(grid_local_inp->rmin(), grid_local_inp->rmax(), grid_local_inp->size());
  }
  if (grid_local->GridTag == CUSTOM_1DGRID)
    myComm->barrier_and_abort("ECPComponentBuilder::buildLocal Custom grid is used. Need to recast to the linear grid");
  else
  {
    std::vector<RealType> v;
    if (bareCoulomb)
    {
      app_log() << "   Bare Coulomb potential is used." << std::endl;
      grid_local->set(0.0, 1., 3);
      v.resize(3);
      for (int ig = 0; ig < 3; ig++)
        v[ig] = 1.0;
      pp_loc = std::make_unique<RadialPotentialType>(std::move(grid_local), v);
      pp_loc->spline(0, 0.0, 2, 0.0);
    }
    else
    {
      app_log() << "   Guassian basisGroup is used: base power " << vr.basePower << std::endl;
      const RealType eps = 1e-12; //numeric_limits<RealType>::epsilon();//1e-12;
      RealType zinv      = 1.0 / Zeff;
      RealType r         = 10.;
      bool ignore        = true;
      int last           = grid_local->size() - 1;
      while (ignore && last)
      {
        r      = (*grid_local)[last];
        ignore = (std::abs(zinv * vr.f(r)) < eps);
        --last;
      }
      if (last == 0)
        myComm->barrier_and_abort("ECPComponentBuilder::buildLocal. Illegal Local Pseudopotential");
      //Add the reset values here
      int ng = static_cast<int>(r / 1e-3) + 1;
      app_log() << "     Use a Linear Grid: [0," << r << "] Number of points = " << ng << std::endl;
      grid_local->set(0.0, r, ng);
      v.resize(ng);
      for (int ig = 1; ig < ng - 1; ig++)
      {
        double r = (*grid_local)[ig];
        v[ig]    = 1.0 - zinv * vr.f(r);
      }
      v[0]      = 2.0 * v[1] - v[2];
      v[ng - 1] = 1.0;
      pp_loc    = std::make_unique<RadialPotentialType>(std::move(grid_local), v);
      pp_loc->spline(); //use the fixed conditions
    }
  }
}

std::unique_ptr<ECPComponentBuilder::mGridType> ECPComponentBuilder::createGrid(xmlNodePtr cur, bool useLinear)
{
  mRealType ri     = 1e-6;
  mRealType rf     = 100.0;
  mRealType ascale = -1.0e0;
  mRealType astep  = -1.0;
  //mRealType astep = 1.25e-2;
  int npts = 1001;
  std::string gridType("log");
  std::string gridID("global");
  OhmmsAttributeSet radAttrib;
  radAttrib.add(gridType, "type");
  radAttrib.add(gridID, "grid_id");
  radAttrib.add(gridID, "grid_def");
  radAttrib.add(gridID, "name");
  radAttrib.add(gridID, "id");
  radAttrib.add(npts, "npts");
  radAttrib.add(ri, "ri");
  radAttrib.add(rf, "rf");
  radAttrib.add(ascale, "ascale");
  radAttrib.add(astep, "astep");
  radAttrib.add(ascale, "scale");
  radAttrib.add(astep, "step");
  radAttrib.put(cur);
  auto git = grid_inp.find(gridID);
  if (git != grid_inp.end())
  {
    return (*git).second->makeClone(); //use the same grid
  }
  //overwrite the grid type to linear starting at 0.0
  if (useLinear)
  {
    gridType = "linear";
    ri       = 0.0;
  }
  std::unique_ptr<mGridType> agrid;
  if (gridType == "log")
  {
    if (ascale > 0.0)
    {
      app_log() << "   Log grid scale=" << ascale << " step=" << astep << " npts=" << npts << std::endl;
      agrid = std::make_unique<LogGridZero<mRealType>>();
      agrid->set(astep, ascale, npts);
    }
    else
    {
      if (ri < std::numeric_limits<mRealType>::epsilon())
      {
        ri = std::numeric_limits<mRealType>::epsilon();
      }
      agrid = std::make_unique<LogGrid<mRealType>>();
      agrid->set(ri, rf, npts);
    }
  }
  else if (gridType == "linear")
  {
    agrid = std::make_unique<LinearGrid<mRealType>>();
    if (astep > 0.0)
    {
      npts = static_cast<int>((rf - ri) / astep) + 1;
    }
    agrid->set(ri, rf, npts);
    app_log() << "   Linear grid  ri=" << ri << " rf=" << rf << " npts = " << npts << std::endl;
  }
  else
  //accept numerical data
  {
    xmlNodePtr cur1 = cur->children;
    while (cur1 != NULL)
    {
      std::string cname((const char*)cur1->name);
      if (cname == "data")
      {
        std::vector<double> gIn(npts);
        putContent(gIn, cur1);
        agrid = std::make_unique<NumericalGrid<mRealType>>(gIn);
        app_log() << "   Numerical grid npts = " << gIn.size() << std::endl;
      }
      cur1 = cur1->next;
    }
  }
  return agrid;
}

} // namespace qmcplusplus

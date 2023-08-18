//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
// Modifications Copyright (C) 2021 Advanced Micro Devices, Inc. All rights reserved.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jakub Kurzak, jakurzak@amd.com, Advanced Micro Devices, Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"
#include "Numerics/OneDimCubicSpline.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/SimpleParser.h"

namespace qmcplusplus
{
//this is abinit/siesta format
void ECPComponentBuilder::buildSemiLocalAndLocal(std::vector<xmlNodePtr>& semiPtr)
{
  app_log() << "    ECPComponentBuilder::buildSemiLocalAndLocal " << std::endl;
  if (grid_global == 0)
    myComm->barrier_and_abort("ECPComponentBuilder::buildSemiLocalAndLocal. Global grid needs to be defined.\n");
  // There should only be one semilocal tag
  if (semiPtr.size() > 1)
  {
    std::stringstream err_msg;
    err_msg << "ECPComponentBuilder::buildSemiLocalAndLocal. "
            << "We have more than one semilocal sections in the PP xml file";
    myComm->barrier_and_abort(err_msg.str());
  }
  RealType rmax = -1;
  //attributes: initailize by defaults
  std::string eunits("hartree");
  std::string format("r*V");
  std::string lloc;
  int ndown = 1;
  int nup   = 0;
  int nso   = 0;
  OhmmsAttributeSet aAttrib;
  int quad_rule     = -1;
  int local_channel = -1;
  aAttrib.add(eunits, "units");
  aAttrib.add(format, "format");
  aAttrib.add(ndown, "npots-down");
  aAttrib.add(nup, "npots-up");
  aAttrib.add(local_channel, "l-local");
  aAttrib.add(quad_rule, "nrule");
  aAttrib.add(Srule, "srule");
  aAttrib.add(nso, "npots-so");

  xmlNodePtr cur_semilocal = semiPtr[0];
  aAttrib.put(cur_semilocal);

  // settle Nrule. Priority: current value (from input file) > PP XML file > lmax derived
  if (quad_rule > -1 && Nrule > -1)
  {
    app_warning() << " Nrule setting found in both qmcpack input (Nrule = " << Nrule
                  << ") and pseudopotential file (Nrule = " << quad_rule << ")."
                  << " Using nrule setting in qmcpack input file." << std::endl;
  }
  else if (quad_rule > -1 && Nrule == -1)
  {
    app_log() << " Nrule setting found in pseudopotential file and used." << std::endl;
    Nrule = quad_rule;
  }
  else if (quad_rule == -1 && Nrule > -1)
    app_log() << " Nrule setting found in qmcpack input file and used." << std::endl;
  else
  {
    //Sperical quadrature rules set by exact integration up to lmax of
    //nonlocal channels.
    //From J. Chem. Phys. 95 (3467) (1991)
    //Keeping Nrule = 4 as default for lmax <= 5.
    switch (pp_nonloc->lmax)
    {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      Nrule = 4;
      break;
    case 6:
    case 7:
      Nrule = 6;
      break;
    case 8:
    case 9:
    case 10:
    case 11:
      Nrule = 7;
      break;
    default:
      myComm->barrier_and_abort("Default value for pseudopotential nrule not determined.");
      break;
    }
    app_warning() << "Nrule was not determined from qmcpack input or pseudopotential file. Setting sensible default."
                  << std::endl;
  }


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
    std::stringstream err_msg;
    err_msg << "ECPComponentBuilder::buildSemiLocalAndLocal."
            << "Unrecognized format \"" << format << "\" in PP file.\n";
    myComm->barrier_and_abort(err_msg.str());
  }
  // We cannot construct the potentials as we construct them since
  // we may not know which one is local yet.

  std::vector<int> angList;
  std::vector<int> angListSO; //For spin-orbit, if it exists
  std::vector<xmlNodePtr> vpsPtr;
  std::vector<xmlNodePtr> vpsoPtr; //For spin-orbit, if it exists.
  Lmax   = -1;
  LmaxSO = -1;
  // Now read vps sections
  xmlNodePtr cur_vps = cur_semilocal->children;
  while (cur_vps != NULL)
  {
    std::string vname((const char*)cur_vps->name);
    if (vname == "vps")
    {
      OhmmsAttributeSet aAttrib;
      std::string lstr("s");
      RealType rc = -1.0;
      aAttrib.add(lstr, "l");
      aAttrib.add(rc, "cutoff");
      aAttrib.put(cur_vps);
      rmax = std::max(rmax, rc);
      if (angMon.find(lstr) == angMon.end())
      {
        std::stringstream err_msg;
        err_msg << "ECPComponentBuilder::buildSemiLocalAndLocal. "
                << "Requested angular momentum " << lstr << " not available.\n";
        myComm->barrier_and_abort(err_msg.str());
      }
      int l = angMon[lstr];
      angList.push_back(l);
      vpsPtr.push_back(cur_vps);
      Lmax = std::max(Lmax, l); //count the maximum L
    }
    else if (vname == "vps_so") //This accumulates the spin-orbit corrections, if defined.
    {
      OhmmsAttributeSet aAttrib;
      std::string lstr("s");
      RealType rc = -1.0;
      aAttrib.add(lstr, "l");
      aAttrib.add(rc, "cutoff");
      aAttrib.put(cur_vps);
      rmax = std::max(rmax, rc);
      if (angMon.find(lstr) == angMon.end())
      {
        std::stringstream err_msg;
        err_msg << "ECPComponentBuilder::buildSemiLocalAndLocal. "
                << "Requested angular momentum " << lstr << " not available for SO.\n";
        myComm->barrier_and_abort(err_msg.str());
      }
      int l = angMon[lstr];
      angListSO.push_back(l);
      vpsoPtr.push_back(cur_vps);
      LmaxSO = std::max(LmaxSO, l); //count the maximum L
    }
    cur_vps = cur_vps->next;
  }

  if (rmax < 0)
    rmax = 1.8;

  // settle Llocal. Priority: current value (from input file) > PP XML file
  if (local_channel > -1 && Llocal > -1)
  {
    app_warning() << " l-local setting found in both qmcpack input (l-local = " << Llocal
                  << ") and pseudopotential file (l-local = " << local_channel << ")."
                  << " Using l-local setting in qmcpack input file." << std::endl;
  }
  else if (local_channel > -1 && Llocal == -1)
  {
    app_log() << " l-local setting found in pseudopotential file and used." << std::endl;
    Llocal = local_channel;
  }
  else if (local_channel == -1 && Llocal > -1)
    app_log() << " l-local setting found in qmcpack input file and used." << std::endl;
  else if (angList.size() == 1)
  {
    Llocal = Lmax;
    app_log() << "    Only one vps is found. Set the local component=" << Lmax << std::endl;
  }
  else
  {
    app_error() << "The local channel is specified in neither the pseudopotential file nor the input file.\n"
                << "Please add \'l-local=\"n\"\' attribute to either file.\n";
    myComm->barrier_and_abort("ECPComponentBuilder::doBreakUp");
  }

  if (angListSO.size() != nso)
  {
    std::stringstream ssout;
    ssout << "Error. npots-so=" << angListSO.size() << " while declared number of SO channels is " << nso << std::endl;
    std::string outstring("");
    outstring = ssout.str();

    myComm->barrier_and_abort(outstring.c_str());
  }
  int npts = grid_global->size();
  Matrix<mRealType> vnn(angList.size(), npts);
  for (int l = 0; l < angList.size(); l++)
  {
    std::vector<mRealType> vt(npts);
    xmlNodePtr c = vpsPtr[l]->children;
    while (c != NULL)
    {
      if (xmlStrEqual(c->name, (const xmlChar*)"radfunc"))
      {
        xmlNodePtr c1 = c->children;
        while (c1 != NULL)
        {
          if (xmlStrEqual(c1->name, (const xmlChar*)"data"))
            putContent(vt, c1);
          c1 = c1->next;
        }
      }
      c = c->next;
    }
    //copy the numerical data with the correct map
    copy(vt.begin(), vt.end(), vnn[angList[l]]);
  }

  //Grabbing the spin-orbit functions from XML.
  Matrix<mRealType> vnnso(angListSO.size(), npts);
  for (int l = 0; l < angListSO.size(); l++)
  {
    std::vector<mRealType> vtso(npts);
    xmlNodePtr c = vpsoPtr[l]->children;
    while (c != NULL)
    {
      if (xmlStrEqual(c->name, (const xmlChar*)"radfunc"))
      {
        xmlNodePtr c1 = c->children;
        while (c1 != NULL)
        {
          if (xmlStrEqual(c1->name, (const xmlChar*)"data"))
            putContent(vtso, c1);
          c1 = c1->next;
        }
      }
      c = c->next;
    }
    //copy the numerical data with the correct map
    //So this is weird, but I feel like l should be the proper index for vnnso,
    //with angListSO[l] being the actual angular momentum channel referred to by l.
    //This differs from the parsing of the nonlocal pseudopotential piece, but whatever.
    copy(vtso.begin(), vtso.end(), vnnso[l]);
  }


  ////rather stupid to do this but necessary
  //vector<RealType> temp(npts);
  //for(int i=0; i<npts; i++) temp[i]=grid_global->r(i);
  if (!is_r_times_V)
  {
    app_log() << "  Input pseudopotential is converted into r*V" << std::endl;
    for (int i = 0; i < vnn.rows(); i++)
      for (int j = 0; j < npts; j++)
        vnn[i][j] *= grid_global->r(j);
    for (int i = 0; i < vnnso.rows(); i++)
      for (int j = 0; j < npts; j++)
        vnnso[i][j] *= grid_global->r(j);
  }
  app_log() << "   Number of angular momentum channels " << angList.size() << std::endl;
  app_log() << "   Maximum angular momentum channel (Lmax) " << Lmax << std::endl;
  doBreakUp(angList, vnn, rmax, Vprefactor);

  //If any spinorbit terms are found...
  if (nso > 0)
    buildSO(angListSO, vnnso, rmax, 1.0);
  else
  {
    //No SO channels found. Delete pp_so
    pp_so.reset();
  }
}

//Most of this is copied directly from doBreakUp, but is separated to prevent from cluttering doBreakUp.
//This function takes the input grid and input SO potential table, interpolates it onto a linear grid if necessary
//via a cubic spline, and then adds the final splined RadialPotentialType object to the SOECPComponent object.
void ECPComponentBuilder::buildSO(const std::vector<int>& angList,
                                  const Matrix<mRealType>& vnnso,
                                  RealType rmax,
                                  mRealType Vprefactor)
{
  const int max_points = 100000;
  app_log() << "   Creating a Linear Grid Rmax=" << rmax << std::endl;
  //this is a new grid
  mRealType d = 1e-4;
  auto agrid  = std::make_unique<LinearGrid<RealType>>();
  // If the global grid is already linear, do not interpolate the data
  int ng;
  if (grid_global->getGridTag() == LINEAR_1DGRID)
  {
    ng = (int)std::ceil(rmax * grid_global->DeltaInv) + 1;
    if (ng <= max_points)
    {
      app_log() << "  Using global grid with delta = " << grid_global->Delta << std::endl;
      rmax = grid_global->Delta * (ng - 1);
      agrid->set(0.0, rmax, ng);
    }
    else
      agrid->set(0.0, rmax, max_points);
  }
  else
  {
    ng = std::min(max_points, static_cast<int>(rmax / d) + 1);
    agrid->set(0, rmax, ng);
  }
  // This is critical!!!
  // If d is not reset, we generate an error in the interpolated PP!
  d        = agrid->Delta;
  int ngIn = vnnso.cols() - 2;
  std::vector<RealType> newP(ng);
  std::vector<mRealType> newPin(ngIn);
  for (int l = 0; l < angList.size(); l++)
  {
    const mRealType* restrict vp = vnnso[l];
    for (int i = 0; i < ngIn; i++)
      newPin[i] = Vprefactor * vp[i];

    OneDimCubicSpline<mRealType> infunc(grid_global->makeClone(), newPin);
    infunc.spline(0, 0.0, ngIn - 1, 0.0);
    for (int i = 1; i < ng - 1; i++)
    {
      mRealType r = d * i;
      newP[i]     = infunc.splint(r) / r;
    }
    newP[0]                  = newP[1];
    newP[ng - 1]             = 0.0;
    RadialPotentialType* app = new RadialPotentialType(agrid->makeClone(), newP);
    app->spline();
    pp_so->add(angList[l], app);
  }
  NumSO = angList.size();
  pp_so->setRmax(rmax);
}

bool ECPComponentBuilder::parseCasino(const std::string& fname, xmlNodePtr cur)
{
  app_log() << "   Start ECPComponentBuilder::parseCasino" << std::endl;
  RealType rmax = 2.0;
  Llocal        = -1;
  Lmax          = -1;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(rmax, "cutoff");
  aAttrib.add(Llocal, "l-local");
  aAttrib.add(Lmax, "lmax");
  aAttrib.add(Nrule, "nrule");
  aAttrib.put(cur);

  std::ifstream fin(fname.c_str(), std::ios_base::in);
  if (!fin)
  {
    app_error() << "Could not open file " << fname << std::endl;
    myComm->barrier_and_abort("ECPComponentBuilder::parseCasino");
  }
  if (!pp_nonloc)
    pp_nonloc = std::make_unique<NonLocalECPComponent>();
  OhmmsAsciiParser aParser;
  int npts = 0, idummy;
  std::string eunits("rydberg");
  app_log() << "    ECPComponentBuilder::parseCasino" << std::endl;
  aParser.skiplines(fin, 1); //Header
  aParser.skiplines(fin, 1); //Atomic number and pseudo-charge
  aParser.getValue(fin, AtomicNumber, Zeff);
  app_log() << "      Atomic number = " << AtomicNumber << "  Zeff = " << Zeff << std::endl;
  aParser.skiplines(fin, 1); //Energy units (rydberg/hartree/ev):
  aParser.getValue(fin, eunits);
  app_log() << "      Unit of the potentials = " << eunits << std::endl;
  mRealType Vprefactor = (eunits == "rydberg") ? 0.5 : 1.0;
  aParser.skiplines(fin, 1); //Angular momentum of local component (0=s,1=p,2=d..)
  aParser.getValue(fin, idummy);
  if (Lmax < 0)
    Lmax = idummy;
  aParser.skiplines(fin, 1); //NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)
  aParser.skiplines(fin, 1); //0 0, not sure what to do yet
  aParser.skiplines(fin, 1); //Number of grid points
  aParser.getValue(fin, npts);
  app_log() << "      Input Grid size = " << npts << std::endl;
  std::vector<mRealType> temp(npts);
  aParser.skiplines(fin, 1); //R(i) in atomic units
  aParser.getValues(fin, temp.begin(), temp.end());
  //create a global grid of numerical type
  grid_global = std::make_unique<NumericalGrid<mRealType>>(temp);
  Matrix<mRealType> vnn(Lmax + 1, npts);
  for (int l = 0; l <= Lmax; l++)
  {
    aParser.skiplines(fin, 1);
    aParser.getValues(fin, vnn[l], vnn[l] + npts);
  }
  std::vector<int> angList(Lmax + 1);
  for (int l = 0; l <= Lmax; l++)
    angList[l] = l;
  // Now, check to see what maximum cutoff should be
  if (vnn.size() > 1)
  {
    const double tolerance = 1.0e-5;
    double rc_check        = grid_global->r(npts - 1);
    for (int j = npts - 1; j > 0; j--)
    {
      bool closeEnough = true;
      for (int i = 0; i < vnn.rows(); i++)
        for (int k = i + 1; k < vnn.rows(); k++)
          if (std::abs(vnn[i][j] - vnn[k][j]) > tolerance)
            closeEnough = false;
      if (!closeEnough)
      {
        rc_check = grid_global->r(j);
        break;
      }
    }
    app_log() << "  Maxium cutoff for non-local pseudopotentials " << rc_check << std::endl;
  }
  doBreakUp(angList, vnn, rmax, Vprefactor);
  SetQuadratureRule(Nrule);
  app_log() << "    Non-local pseudopotential parameters" << std::endl;
  pp_nonloc->print(app_log());
  return true;
}


/** Separate local from non-local potentials
 * @param angList angular momentum list
 * @param vnn semilocal tables of size (angList.size(),global_grid->size())
 * @param rmax cutoff radius
 * @param Vprefactor conversion factor to Hartree
 *
 * Note that local pseudopotential is r*V !!!
 */
void ECPComponentBuilder::doBreakUp(const std::vector<int>& angList,
                                    const Matrix<mRealType>& vnn,
                                    RealType rmax,
                                    mRealType Vprefactor)
{
  //ALERT magic number
  const int max_points = 100000;
  app_log() << "   Creating a Linear Grid Rmax=" << rmax << std::endl;
  //this is a new grid
  mRealType d = 1e-4;
  auto agrid  = std::make_unique<LinearGrid<RealType>>();
  // If the global grid is already linear, do not interpolate the data
  int ng;
  if (grid_global->getGridTag() == LINEAR_1DGRID)
  {
    ng = (int)std::ceil(rmax * grid_global->DeltaInv) + 1;
    if (ng <= max_points)
    {
      app_log() << "  Using global grid with delta = " << grid_global->Delta << std::endl;
      rmax = grid_global->Delta * (ng - 1);
      agrid->set(0.0, rmax, ng);
    }
    else
      agrid->set(0.0, rmax, max_points);
  }
  else
  {
    ng = std::min(max_points, static_cast<int>(rmax / d) + 1);
    agrid->set(0, rmax, ng);
  }
  // This is critical!!!
  // If d is not reset, we generate an error in the interpolated PP!
  d        = agrid->Delta;
  int ngIn = vnn.cols() - 2;

  assert(angList.size() > 0 && "ECPComponentBuilder::doBreakUp angList cannot be empty!");
  //find the index of local
  int iLlocal = -1;
  for (int l = 0; l < angList.size(); l++)
    if (angList[l] == Llocal)
      iLlocal = l;
  std::vector<RealType> newP(ng);
  std::vector<mRealType> newPin(ngIn);
  for (int l = 0; l < angList.size(); l++)
  {
    if (angList[l] == Llocal)
      continue;
    const mRealType* restrict vp    = vnn[angList[l]];
    const mRealType* restrict vpLoc = vnn[iLlocal];
    int ll                          = angList[l];
    for (int i = 0; i < ngIn; i++)
      newPin[i] = Vprefactor * (vp[i] - vpLoc[i]);
    OneDimCubicSpline<mRealType> infunc(grid_global->makeClone(), newPin);
    infunc.spline(0, 0.0, ngIn - 1, 0.0);
    for (int i = 1; i < ng - 1; i++)
    {
      mRealType r = d * i;
      newP[i]     = infunc.splint(r) / r;
    }
    newP[0]      = newP[1];
    newP[ng - 1] = 0.0;
    auto app     = new RadialPotentialType(agrid->makeClone(), newP);
    app->spline();
    pp_nonloc->add(angList[l], app);
  }
  NumNonLocal = Lmax;
  if (Llocal == Lmax)
    Lmax--;
  if (NumNonLocal)
  {
    pp_nonloc->lmax = Lmax;
    pp_nonloc->Rmax = rmax;
  }
  else
  {
    //only one component, remove non-local potentials
    pp_nonloc.reset();
  }
  {
    // Spline local potential on original grid
    newPin[0]                       = 0.0;
    mRealType vfac                  = -Vprefactor / Zeff;
    const mRealType* restrict vpLoc = vnn[iLlocal];
    for (int i = 0; i < ngIn; i++)
      newPin[i] = vfac * vpLoc[i];
    double dy0 = (newPin[1] - newPin[0]) / ((*grid_global)[1] - (*grid_global)[0]);
    OneDimCubicSpline<mRealType> infunc(grid_global->makeClone(), newPin);
    infunc.spline(0, dy0, ngIn - 1, 0.0);
    int m          = grid_global->size();
    double loc_max = grid_global->r(m - 1);
    int nloc       = (int)std::floor(loc_max / d);
    loc_max        = (nloc - 1) * d;
    auto grid_loc  = std::make_unique<LinearGrid<RealType>>();
    grid_loc->set(0.0, loc_max, nloc);
    app_log() << "   Making L=" << Llocal << " a local potential with a radial cutoff of " << loc_max << std::endl;
    std::vector<RealType> newPloc(nloc);
    for (int i = 1; i < nloc - 1; i++)
    {
      mRealType r = d * i;
      newPloc[i]  = infunc.splint(r);
    }
    newPloc[0]        = 0.0;
    newPloc[nloc - 1] = 1.0;
    pp_loc            = std::make_unique<RadialPotentialType>(std::move(grid_loc), newPloc);
    pp_loc->spline(0, dy0, nloc - 1, 0.0);
    // for (double r=0.0; r<3.50001; r+=0.001)
    //   fprintf (stderr, "%10.5f %10.5f\n", r, pp_loc->splint(r));
  }
}

} // namespace qmcplusplus

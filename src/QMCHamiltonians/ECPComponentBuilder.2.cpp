//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/OneDimCubicSpline.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/SimpleParser.h"
//#include "QMCHamiltonians/FSAtomPseudoPot.h"
//#include "Utilities/IteratorUtility.h"
#ifdef QMC_CUDA
#include <cuda_runtime_api.h>
#endif

namespace qmcplusplus
{
//this is abinit/siesta format
void ECPComponentBuilder::buildSemiLocalAndLocal(std::vector<xmlNodePtr>& semiPtr)
{
  app_log() << "    ECPComponentBuilder::buildSemiLocalAndLocal " << std::endl;
  if(grid_global == 0)
  {
    app_error() << "    Global grid needs to be defined." << std::endl;
    APP_ABORT("ECPComponentBuilder::buildSemiLocalAndLocal");
  }
  // There should only be one semilocal tag
  if (semiPtr.size()> 1)
  {
    app_error() << "    We have more than one semilocal sections in the PP xml file." << std::endl;
    APP_ABORT("ECPComponentBuilder::buildSemiLocalAndLocal");
  }
  RealType rmax = -1;
  //attributes: initailize by defaults
  std::string eunits("hartree");
  std::string format("r*V");
  std::string lloc;
  int ndown=1;
  int nup=0;
  Llocal = -1;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(eunits,"units");
  aAttrib.add(format,"format");
  aAttrib.add(ndown,"npots-down");
  aAttrib.add(nup,"npots-up");
  aAttrib.add(Llocal,"l-local");
  aAttrib.add(Nrule,"nrule");
  xmlNodePtr cur_semilocal = semiPtr[0];
  aAttrib.put(cur_semilocal);
  RealType Vprefactor=1.0;
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
  else
    if (format == "V")
      is_r_times_V = false;
    else
    {
      app_error() << "Unrecognized format \"" << format << "\" in PP file." << std::endl;
      APP_ABORT("ECPComponentBuilder::buildSemiLocalAndLocal");
    }
  // We cannot construct the potentials as we construct them since
  // we may not know which one is local yet.
  std::vector<int> angList;
  std::vector<xmlNodePtr> vpsPtr;
  int iLocal=-1;
  Lmax=-1;
  // Now read vps sections
  xmlNodePtr cur_vps = cur_semilocal->children;
  while (cur_vps != NULL)
  {
    std::string vname ((const char*)cur_vps->name);
    if (vname == "vps")
    {
      OhmmsAttributeSet aAttrib;
      std::string lstr("s");
      RealType rc=-1.0;
      aAttrib.add(lstr,"l");
      aAttrib.add(rc,"cutoff");
      aAttrib.put(cur_vps);
      rmax=std::max(rmax,rc);
      int l=angMon[lstr];
      angList.push_back(l);
      vpsPtr.push_back(cur_vps);
      Lmax=std::max(Lmax,l); //count the maximum L
    }
    cur_vps = cur_vps->next;
  }
  if(rmax<0)
    rmax=1.8;
  if(angList.size()==1)
  {
    Llocal=Lmax;
    app_log() << "    Only one vps is found. Set the local component=" << Lmax << std::endl;
  }
  int npts=grid_global->size();
  Matrix<mRealType> vnn(angList.size(),npts);
  for(int l=0; l<angList.size(); l++)
  {
    std::vector<mRealType>  vt(npts);
    xmlNodePtr c=vpsPtr[l]->children;
    while(c != NULL)
    {
      if(xmlStrEqual(c->name,(const xmlChar*)"radfunc"))
      {
        xmlNodePtr c1=c->children;
        while(c1 != NULL)
        {
          if(xmlStrEqual(c1->name,(const xmlChar*)"data"))
            putContent(vt,c1);
          c1=c1->next;
        }
      }
      c=c->next;
    }
    //copy the numerical data with the correct map
    copy(vt.begin(),vt.end(),vnn[angList[l]]);
  }
  ////rather stupid to do this but necessary
  //vector<RealType> temp(npts);
  //for(int i=0; i<npts; i++) temp[i]=grid_global->r(i);
  if(!is_r_times_V)
  {
    app_log() << "  Input pseudopotential is converted into r*V" << std::endl;
    for(int i=0; i<vnn.rows(); i++)
      for(int j=0; j <npts; j++)
        vnn[i][j] *= grid_global->r(j);
  }
  app_log() << "   Number of angular momentum channels " << angList.size() << std::endl;
  app_log() << "   Maximum angular momentum channel " << Lmax << std::endl;
  doBreakUp(angList,vnn,rmax,Vprefactor);
}

bool
ECPComponentBuilder::parseCasino(const std::string& fname, xmlNodePtr cur)
{
  app_log() << "   Start ECPComponentBuilder::parseCasino" << std::endl;
  RealType rmax=2.0;
  Llocal=-1;
  Lmax=-1;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(rmax,"cutoff");
  aAttrib.add(Llocal,"l-local");
  aAttrib.add(Lmax,"lmax");
  aAttrib.add(Nrule,"nrule");
  aAttrib.put(cur);
  //const xmlChar* rptr=xmlGetProp(cur,(const xmlChar*)"cutoff");
  //if(rptr != NULL) rmax = atof((const char*)rptr);
  //app_log() << "   Creating a Linear Grid Rmax=" << rmax << std::endl;
  //const RealType d=5e-4;
  //LinearGrid<RealType>* agrid = new LinearGrid<RealType>;
  //int ng=static_cast<int>(rmax/d)+1;
  //agrid->set(0,rmax,ng);
  std::ifstream fin(fname.c_str(),std::ios_base::in);
  if(!fin)
  {
    app_error() << "Could not open file " << fname << std::endl;
    APP_ABORT("ECPComponentBuilder::parseCasino");
  }
  if(pp_nonloc==0)
    pp_nonloc=new NonLocalECPComponent;
  OhmmsAsciiParser aParser;
  int atomNumber=0;
  int npts=0, idummy;
  std::string eunits("rydberg");
  app_log() << "    ECPComponentBuilder::parseCasino" << std::endl;
  aParser.skiplines(fin,1);//Header
  aParser.skiplines(fin,1);//Atomic number and pseudo-charge
  aParser.getValue(fin,atomNumber,Zeff);
  app_log() << "      Atomic number = " << atomNumber << "  Zeff = " << Zeff << std::endl;
  aParser.skiplines(fin,1);//Energy units (rydberg/hartree/ev):
  aParser.getValue(fin,eunits);
  app_log() << "      Unit of the potentials = " << eunits << std::endl;
  mRealType Vprefactor = (eunits == "rydberg")?0.5:1.0;
  aParser.skiplines(fin,1);//Angular momentum of local component (0=s,1=p,2=d..)
  aParser.getValue(fin,idummy);
  if(Lmax<0)
    Lmax=idummy;
  aParser.skiplines(fin,1);//NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)
  aParser.skiplines(fin,1);//0 0, not sure what to do yet
  aParser.skiplines(fin,1);//Number of grid points
  aParser.getValue(fin,npts);
  app_log() << "      Input Grid size = " << npts << std::endl;
  std::vector<mRealType> temp(npts);
  aParser.skiplines(fin,1);//R(i) in atomic units
  aParser.getValues(fin,temp.begin(),temp.end());
  //create a global grid of numerical type
  grid_global= new NumericalGrid<mRealType>(temp);
  Matrix<mRealType> vnn(Lmax+1,npts);
  for(int l=0; l<=Lmax; l++)
  {
    aParser.skiplines(fin,1);
    aParser.getValues(fin,vnn[l],vnn[l]+npts);
  }
  std::vector<int> angList(Lmax+1);
  for(int l=0; l<=Lmax; l++)
    angList[l]=l;
  // Now, check to see what maximum cutoff should be
  if(vnn.size()>1)
  {
    const double tolerance=1.0e-5;
    double rc_check = grid_global->r(npts-1);
    for (int j=npts-1; j>0; j++)
    {
      bool closeEnough = true;
      for (int i=0; i<vnn.rows(); i++)
        for (int k=i+1; k<vnn.rows(); k++)
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
  doBreakUp(angList,vnn,rmax,Vprefactor);
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
void
ECPComponentBuilder::doBreakUp(const std::vector<int>& angList,
                               const Matrix<mRealType>& vnn,
                               RealType rmax, mRealType Vprefactor)
{
#ifdef QMC_CUDA
  int device;
  cudaGetDevice(&device);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device);
  const int max_points = deviceProp.maxTexture1D-1;
#else
  const int max_points = 100000;
#endif
  app_log() << "   Creating a Linear Grid Rmax=" << rmax << std::endl;
  //this is a new grid
  mRealType d=1e-4;
  LinearGrid<RealType>* agrid = new LinearGrid<RealType>;
  // If the global grid is already linear, do not interpolate the data
  int ng;
  if (grid_global->getGridTag() == LINEAR_1DGRID)
  {
    ng = (int)std::ceil(rmax*grid_global->DeltaInv) + 1;
    if (ng <= max_points)
    {
      app_log() << "  Using global grid with delta = "
                << grid_global->Delta << std::endl;
      rmax = grid_global->Delta * (ng-1);
      agrid->set(0.0,rmax,ng);
    }
    else
      agrid->set(0.0,rmax,max_points);
  }
  else
  {
    ng = std::min(max_points, static_cast<int>(rmax/d)+1);
    agrid->set(0,rmax,ng);
  }
  // This is critical!!!
  // If d is not reset, we generate an error in the interpolated PP!
  d = agrid->Delta;
  int ngIn=vnn.cols()-2;
  if (Llocal == -1 && Lmax > 0)
  {
    app_error() << "The local channel is not specified in the pseudopotential file.\n"
                << "Please add \'l-local=\"n\"\' attribute the semilocal section of the fsatom XML file.\n";
    APP_ABORT("ECPComponentBuilder::doBreakUp");
    // Llocal = Lmax;
  }
  //find the index of local
  int iLlocal=-1;
  for(int l=0; l<angList.size(); l++)
    if(angList[l] == Llocal)
      iLlocal=l;
  std::vector<RealType> newP(ng);
  std::vector<mRealType> newPin(ngIn);
  for(int l=0; l<angList.size(); l++)
  {
    if(angList[l] == Llocal)
      continue;
    const mRealType* restrict vp=vnn[angList[l]];
    const mRealType* restrict vpLoc=vnn[iLlocal];
    int ll=angList[l];
    for(int i=0; i<ngIn; i++)
      newPin[i]=Vprefactor*(vp[i]-vpLoc[i]);
    OneDimCubicSpline<mRealType> infunc(grid_global,newPin);
    infunc.spline(0,0.0,ngIn-1,0.0);
    for(int i=1; i<ng-1; i++)
    {
      mRealType r=d*i;
      newP[i]=infunc.splint(r)/r;
    }
    newP[0]=newP[1];
    newP[ng-1]=0.0;
    RadialPotentialType *app = new RadialPotentialType(agrid,newP);
    app->spline();
    pp_nonloc->add(angList[l],app);
  }
  NumNonLocal=Lmax;
  if (Llocal == Lmax)
    Lmax--;
  if(NumNonLocal)
  {
    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;
  }
  else
  {
    //only one component, remove non-local potentials
    delete pp_nonloc;
    pp_nonloc=0;
  }
  {
    // Spline local potential on original grid
    newPin[0]=0.0;
    mRealType vfac=-Vprefactor/Zeff;
    const mRealType* restrict vpLoc=vnn[iLlocal];
    for(int i=0; i<ngIn; i++)
      newPin[i]=vfac*vpLoc[i];
    double dy0 = (newPin[1] - newPin[0])/((*grid_global)[1]-(*grid_global)[0]);
    OneDimCubicSpline<mRealType> infunc(grid_global,newPin);
    infunc.spline(0,dy0,ngIn-1,0.0);
    int m = grid_global->size();
    double loc_max = grid_global->r(m-1);
    int nloc = (int)std::floor(loc_max / d);
    loc_max = (nloc-1) * d;
    LinearGrid<RealType>* grid_loc = new LinearGrid<RealType>;
    grid_loc->set(0.0, loc_max, nloc);
    app_log() << "   Making L=" << Llocal
              << " a local potential with a radial cutoff of "
              << loc_max << std::endl;
    std::vector<RealType> newPloc(nloc);
    for(int i=1; i<nloc-1; i++)
    {
      mRealType r=d*i;
      newPloc[i]=infunc.splint(r);
    }
    newPloc[0]      = 0.0;
    newPloc[nloc-1] = 1.0;
    pp_loc = new RadialPotentialType(grid_loc,newPloc);
    pp_loc->spline(0, dy0, nloc-1, 0.0);
    // for (double r=0.0; r<3.50001; r+=0.001)
    //   fprintf (stderr, "%10.5f %10.5f\n", r, pp_loc->splint(r));
  }
}

} // namespace qmcPlusPlus

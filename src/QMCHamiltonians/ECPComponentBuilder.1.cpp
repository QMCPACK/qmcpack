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

namespace qmcplusplus
{

void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur)
{
  mGridType* grid_semilocal=0;
  RealType rmax= pp_nonloc->Rmax;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "grid")
    {
      grid_semilocal=createGrid(cur);
      rmax=grid_semilocal->rmax();
    }
    else if(cname == "vps")
    {
      //should be able to overwrite rmax
      int l=angMon[(const char*)xmlGetProp(cur,(const xmlChar*)"l")];
      Lmax = std::max(l,Lmax);
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)cur1->name);
        if(cname1 == "basisGroup")
        {
          pp_nonloc->add(l,createVrWithBasisGroup(cur1,grid_semilocal));
        }
        //else if(cname1 == "data")
        //{
        //  pp_nonloc->add(l,createVrWithData(cur1,grid_semilocal));
        //}
        cur1=cur1->next;
      }
      NumNonLocal++;
    }
    cur=cur->next;
  }
  pp_nonloc->lmax=Lmax;
  pp_nonloc->Rmax=rmax;
}

ECPComponentBuilder::RadialPotentialType*
ECPComponentBuilder::createVrWithBasisGroup(xmlNodePtr cur, mGridType* agrid)
{
  //todo rcut should be reset if necessary
  typedef GaussianTimesRN<RealType> InFuncType;
  InFuncType a;
  a.putBasisGroup(cur);
  bool ignore=true;
  const RealType eps=1e-4;
  mRealType rout=agrid->rmax()*2;
  while(ignore&&rout>agrid->rmax())
  {
    ignore=(std::abs(a.f(rout))<eps);
    rout-=0.01;
  }
  rout += 0.01;
  app_log() << "  cutoff for non-local pseudopotential = " << agrid->rmax() << std::endl;
  app_log() << "  calculated cutoff for " << eps << " = " << rout << std::endl;
  int ng = agrid->size();
  //all ways reset grid. Ye Luo
  GridType* agrid_new;
  //if(agrid->GridTag != LINEAR_1DGRID)// || rout > agrid->rmax())
  {
    RealType ri=0.0;
    RealType rf=std::max(rout,agrid->rmax());
    agrid_new = new LinearGrid<RealType>;
    agrid_new->set(ri,rf,ng);
    app_log() << "  Reset the grid for SemiLocal component to a LinearGrid. " << std::endl;
  }
  std::vector<RealType> v(ng);
  for(int ig=0; ig<ng; ig++)
  {
    v[ig]=a.f((*agrid_new)[ig]);
  }
  v[ng-1]=0.0;
  RadialPotentialType *app=new RadialPotentialType(agrid_new,v);
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
  if(pp_loc)
    return; //something is wrong
  std::string vFormat("V");
  const xmlChar* vptr=xmlGetProp(cur,(const xmlChar*)"format");
  if(vptr != NULL)
  {
    vFormat=(const char*)vptr;
  }
  int vPowerCorrection=1;
  if(vFormat == "r*V")
  {
    app_log() << "  Local pseudopotential format = r*V" << std::endl;
    vPowerCorrection=0;
  }
  else
  {
    app_log() << "  Local pseudopotential format = V" << std::endl;
  }
  typedef GaussianTimesRN<RealType> InFuncType;
  GridType* grid_local=0;
  mGridType* grid_local_inp=0;
  InFuncType vr;
  bool bareCoulomb=true;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "grid")
    {
      grid_local_inp=createGrid(cur,true);
    }
    else
      if(cname == "basisGroup")
      {
        vr.putBasisGroup(cur,vPowerCorrection);
        bareCoulomb=false;
      }
      else
        if(cname == "data")
        {
          pp_loc=createVrWithData(cur,grid_local_inp,vPowerCorrection);
          app_log() << "  Local pseduopotential in a <data/>" << std::endl;
          return;
        }
    cur=cur->next;
  }
  if(grid_local_inp == 0)
  {
    if(grid_global == 0)
    {
      APP_ABORT("ECPComponentBuilder::buildLocal Missing grid information. ");
    }
    grid_local = new LinearGrid<RealType>;
    grid_local->set(grid_global->rmin(),grid_global->rmax(),grid_global->size());
  }
  else
  {
    grid_local = new LinearGrid<RealType>;
    grid_local->set(grid_local_inp->rmin(),grid_local_inp->rmax(),grid_local_inp->size());
  }
  if(grid_local->GridTag == CUSTOM_1DGRID)
  {
    APP_ABORT("ECPComponentBuilder::buildLocal Custom grid is used. Need to recast to the linear grid");
  }
  else
  {
    std::vector<RealType> v;
    if(bareCoulomb)
    {
      app_log() << "   Bare Coulomb potential is used." << std::endl;
      grid_local->set(0.0,1.,3);
      v.resize(3);
      for(int ig=0; ig<3; ig++)
        v[ig]=1.0;
      pp_loc=new RadialPotentialType(grid_local,v);
      pp_loc->spline(0,0.0,2,0.0);
    }
    else
    {
      app_log() << "   Guassian basisGroup is used: base power " << vr.basePower << std::endl;
      const RealType eps=1e-12;//numeric_limits<RealType>::epsilon();//1e-12;
      RealType zinv=1.0/Zeff;
      RealType r=10.;
      bool ignore=true;
      int last=grid_local->size()-1;
      while(ignore&&last)
      {
        r=(*grid_local)[last];
        ignore=(std::abs(zinv*vr.f(r))<eps);
        --last;
      }
      if(last ==0)
      {
        app_error() << "  Illegal Local Pseudopotential " << std::endl;
      }
      //Add the reset values here
      int ng=static_cast<int>(r/1e-3)+1;
      app_log() << "     Use a Linear Grid: [0,"<< r << "] Number of points = " << ng << std::endl;
      grid_local->set(0.0,r,ng);
      v.resize(ng);
      for(int ig=1; ig<ng-1; ig++)
      {
        double r=(*grid_local)[ig];
        v[ig]=1.0-zinv*vr.f(r);
      }
      v[0]=2.0*v[1]-v[2];
      v[ng-1]=1.0;
      pp_loc=new RadialPotentialType(grid_local,v);
      pp_loc->spline(); //use the fixed conditions
    }
  }
}

ECPComponentBuilder::mGridType* ECPComponentBuilder::createGrid(xmlNodePtr cur, bool useLinear)
{
  mRealType ri = 1e-6;
  mRealType rf = 100.0;
  mRealType ascale = -1.0e0;
  mRealType astep = -1.0;
  //mRealType astep = 1.25e-2;
  int npts = 1001;
  std::string gridType("log");
  std::string gridID("global");
  OhmmsAttributeSet radAttrib;
  radAttrib.add(gridType,"type");
  radAttrib.add(gridID,"grid_id");
  radAttrib.add(gridID,"grid_def");
  radAttrib.add(gridID,"name");
  radAttrib.add(gridID,"id");
  radAttrib.add(npts,"npts");
  radAttrib.add(ri,"ri");
  radAttrib.add(rf,"rf");
  radAttrib.add(ascale,"ascale");
  radAttrib.add(astep,"astep");
  radAttrib.add(ascale,"scale");
  radAttrib.add(astep,"step");
  radAttrib.put(cur);
  std::map<std::string,mGridType*>::iterator git(grid_inp.find(gridID));
  if(git != grid_inp.end())
  {
    return (*git).second; //use the same grid
  }
  //overwrite the grid type to linear starting at 0.0
  if(useLinear)
  {
    gridType="linear";
    ri=0.0;
  }
  mGridType *agrid=0;
  if(gridType == "log")
  {
    if(ascale>0.0)
    {
      app_log() << "   Log grid scale=" << ascale << " step=" << astep << " npts=" << npts << std::endl;
      agrid = new LogGridZero<mRealType>;
      agrid->set(astep,ascale,npts);
    }
    else
    {
      if(ri<std::numeric_limits<mRealType>::epsilon())
      {
        ri= std::numeric_limits<mRealType>::epsilon();
      }
      agrid = new LogGrid<mRealType>;
      agrid->set(ri,rf,npts);
    }
  }
  else
    if(gridType == "linear")
    {
      agrid = new LinearGrid<mRealType>;
      if(astep>0.0)
      {
        npts = static_cast<int>((rf-ri)/astep)+1;
      }
      agrid->set(ri,rf,npts);
      app_log() << "   Linear grid  ri=" << ri << " rf=" << rf << " npts = " << npts << std::endl;
    }
    else
      //accept numerical data
    {
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname((const char*)cur1->name);
        if(cname == "data")
        {
          std::vector<double> gIn(npts);
          putContent(gIn,cur1);
          agrid = new NumericalGrid<mRealType>(gIn);
          app_log() << "   Numerical grid npts = " <<  gIn.size() << std::endl;
        }
        cur1 = cur1->next;
      }
    }
  return agrid;
}

/** Disable pseudo/semilocal/vps/data */
ECPComponentBuilder::RadialPotentialType*
ECPComponentBuilder::createVrWithData(xmlNodePtr cur, mGridType* agrid, int rCorrection)
{
  return 0;
  //  RealType rcIn = agrid->rmax();
  //  //use the maximum value of the grid
  //  if(RcutMax<0) RcutMax=rcIn;
  //  //create a new linear grid if the input grid is not good enough
  //  GridType *newgrid=0;
  //  if(agrid->GridTag != LINEAR_1DGRID || RcutMax < rcIn)
  //  {
  //    const RealType delta=1000.; // use 1/000
  //    newgrid = new LinearGrid<RealType>;
  //    newgrid->set(0.0,RcutMax,static_cast<int>(RcutMax*delta)+1);
  //  }
  //  //read the numerical data
  //  std::vector<RealType> pdata;
  //  putContent(pdata,cur);
  //  if(pdata.size() != agrid->size())
  //  {
  //    app_error() << "  ECPComponentBuilder::createVrWithData vsp/data size does not match." << std::endl;
  //    abort(); //FIXABORT
  //  }
  //  if(rCorrection == 1)
  //  {
  //    for(int i=0; i<agrid->size(); i++) pdata[i] *= (*agrid)[i];
  //  }
  //  if(newgrid)
  //  {
  //    OneDimCubicSpline<RealType> inFunc(grid_global,pdata);
  //    inFunc.spline();
  //    int ng=newgrid->size();
  //    pdata.resize(ng);
  //    for(int i=0; i<ng; i++)
  //    {
  //      RealType r((*agrid)[i]);
  //      pdata[i]=inFunc.splint(r);
  //    }
  //    if(agrid->rmin()>0.0) pdata[0]=pdata[1];
  //    RadialPotentialType *app = new RadialPotentialType(newgrid,pdata);
  //    app->spline();
  //    return app;
  //  }
  //  else
  //  {//use Radial potential with the input grid
  //    RadialPotentialType *app = new RadialPotentialType(agrid,pdata);
  //    app->spline();
  //    return app;
  //  }
}


} // namespace qmcPlusPlus

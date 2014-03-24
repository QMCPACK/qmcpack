//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"

namespace qmcplusplus
{

void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur)
{
  GridType* grid_semilocal=0;
  RealType rmax= pp_nonloc->Rmax;
  cur=cur->children;
  while(cur != NULL)
  {
    string cname((const char*)cur->name);
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
        string cname1((const char*)cur1->name);
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
ECPComponentBuilder::createVrWithBasisGroup(xmlNodePtr cur, GridType* agrid)
{
  //todo rcut should be reset if necessary
  typedef GaussianTimesRN<RealType> InFuncType;
  InFuncType a;
  a.putBasisGroup(cur);
  bool ignore=true;
  const RealType eps=1e-4;
  RealType rout=agrid->rmax()*2;
  while(ignore&&rout>agrid->rmax())
  {
    ignore=(abs(a.f(rout))<eps);
    rout-=0.01;
  }
  rout += 0.01;
  app_log() << "  cutoff for non-local pseudopotential = " << agrid->rmax() << endl;
  app_log() << "  calculated cutoff for " << eps << " = " << rout << endl;
  int ng = agrid->size();
  if(agrid->GridTag != LINEAR_1DGRID)// || rout > agrid->rmax())
  {
    RealType ri=0.0;
    RealType rf=std::max(rout,agrid->rmax());
    delete agrid;
    agrid = new LinearGrid<RealType>;
    agrid->set(ri,rf,ng);
    app_log() << "  Reset the grid for SemiLocal component to a LinearGrid. " << endl;
  }
  std::vector<RealType> v(ng);
  for(int ig=0; ig<ng; ig++)
  {
    v[ig]=a.f((*agrid)[ig]);
  }
  v[ng-1]=0.0;
  RadialPotentialType *app=new RadialPotentialType(agrid,v);
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
  string vFormat("V");
  const xmlChar* vptr=xmlGetProp(cur,(const xmlChar*)"format");
  if(vptr != NULL)
  {
    vFormat=(const char*)vptr;
  }
  int vPowerCorrection=1;
  if(vFormat == "r*V")
  {
    app_log() << "  Local pseudopotential format = r*V" << endl;
    vPowerCorrection=0;
  }
  else
  {
    app_log() << "  Local pseudopotential format = V" << endl;
  }
  typedef GaussianTimesRN<RealType> InFuncType;
  GridType* grid_local=0;
  InFuncType vr;
  bool bareCoulomb=true;
  cur=cur->children;
  while(cur != NULL)
  {
    string cname((const char*)cur->name);
    if(cname == "grid")
    {
      grid_local=createGrid(cur,true);
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
          pp_loc=createVrWithData(cur,grid_local,vPowerCorrection);
          app_log() << "  Local pseduopotential in a <data/>" << endl;
          return;
        }
    cur=cur->next;
  }
  if(grid_local == 0)
  {
    if(grid_global == 0)
    {
      APP_ABORT("ECPComponentBuilder::buildLocal Missing grid information. ");
    }
    grid_local=grid_global;
  }
  if(grid_local->GridTag == CUSTOM_1DGRID)
  {
    APP_ABORT("ECPComponentBuilder::buildLocal Custom grid is used. Need to recast to the linear grid");
  }
  else
  {
    vector<RealType> v;
    if(bareCoulomb)
    {
      app_log() << "   Bare Coulomb potential is used." << endl;
      grid_local->set(0.0,1.,3);
      v.resize(3);
      for(int ig=0; ig<3; ig++)
        v[ig]=1.0;
      pp_loc=new RadialPotentialType(grid_local,v);
      pp_loc->spline(0,0.0,2,0.0);
    }
    else
    {
      app_log() << "   Guassian basisGroup is used: base power " << vr.basePower << endl;
      const RealType eps=1e-12;//numeric_limits<RealType>::epsilon();//1e-12;
      RealType zinv=1.0/Zeff;
      RealType r=10.;
      bool ignore=true;
      int last=grid_local->size()-1;
      while(ignore&&last)
      {
        r=(*grid_local)[last];
        ignore=(abs(zinv*vr.f(r))<eps);
        --last;
      }
      if(last ==0)
      {
        app_error() << "  Illegal Local Pseudopotential " <<endl;
      }
      //Add the reset values here
      int ng=static_cast<int>(r/1e-3)+1;
      app_log() << "     Use a Linear Grid: [0,"<< r << "] Number of points = " << ng << endl;
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

ECPComponentBuilder::GridType* ECPComponentBuilder::createGrid(xmlNodePtr cur, bool useLinear)
{
  RealType ri = 1e-6;
  RealType rf = 100.0;
  RealType ascale = -1.0e0;
  RealType astep = -1.0;
  //RealType astep = 1.25e-2;
  int npts = 1001;
  string gridType("log");
  string gridID("global");
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
  map<string,GridType*>::iterator git(grid_inp.find(gridID));
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
  GridType *agrid=0;
  if(gridType == "log")
  {
    if(ascale>0.0)
    {
      app_log() << "   Log grid scale=" << ascale << " step=" << astep << " npts=" << npts << endl;
      agrid = new LogGridZero<RealType>;
      agrid->set(astep,ascale,npts);
    }
    else
    {
      if(ri<numeric_limits<RealType>::epsilon())
      {
        ri=numeric_limits<RealType>::epsilon();
      }
      agrid = new LogGrid<RealType>;
      agrid->set(ri,rf,npts);
    }
  }
  else
    if(gridType == "linear")
    {
      agrid = new LinearGrid<RealType>;
      if(astep>0.0)
      {
        npts = static_cast<int>((rf-ri)/astep)+1;
      }
      agrid->set(ri,rf,npts);
      app_log() << "   Linear grid  ri=" << ri << " rf=" << rf << " npts = " << npts << endl;
    }
    else
      //accept numerical data
    {
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        string cname((const char*)cur1->name);
        if(cname == "data")
        {
          std::vector<double> gIn(npts);
          putContent(gIn,cur1);
          agrid = new NumericalGrid<RealType>(gIn);
          app_log() << "   Numerical grid npts = " <<  gIn.size() << endl;
        }
        cur1 = cur1->next;
      }
    }
  return agrid;
}

/** Disable pseudo/semilocal/vps/data */
ECPComponentBuilder::RadialPotentialType*
ECPComponentBuilder::createVrWithData(xmlNodePtr cur, GridType* agrid, int rCorrection)
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
  //  vector<RealType> pdata;
  //  putContent(pdata,cur);
  //  if(pdata.size() != agrid->size())
  //  {
  //    app_error() << "  ECPComponentBuilder::createVrWithData vsp/data size does not match." << endl;
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1551 $   $Date: 2006-12-02 09:32:17 -0600 (Sat, 02 Dec 2006) $
 * $Id: ECPComponentBuilder.cpp 1551 2006-12-02 15:32:17Z jnkim $
 ***************************************************************************/

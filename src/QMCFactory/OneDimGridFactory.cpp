//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCFactory/OneDimGridFactory.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{

//initialize the static data
OneDimGridFactory::GridObjectMapType OneDimGridFactory::GridObjects;

OneDimGridFactory::GridType* OneDimGridFactory::createGrid(xmlNodePtr cur)
{
  GridType *agrid=0;
  RealType ri = 1e-5;
  RealType rf = 100.0;
  RealType ascale = -1.0e0;
  RealType astep = 1.25e-2;
  IndexType npts = 1001;
  string gridType("log");
  string gridID("invalid");
  OhmmsAttributeSet radAttrib;
  radAttrib.add(gridType,"type");
  radAttrib.add(npts,"npts");
  radAttrib.add(ri,"ri");
  radAttrib.add(rf,"rf");
  radAttrib.add(ascale,"ascale");
  radAttrib.add(astep,"astep");
  radAttrib.add(ascale,"scale");
  radAttrib.add(astep,"step");
  radAttrib.add(gridID,"id");
  radAttrib.add(gridID,"name");
  radAttrib.add(gridID,"ref");
  if(cur != NULL)
    radAttrib.put(cur);
  //return for the same grid
  bool hasName= (gridID != "invalid");
  if(hasName)
  {
    std::map<string,GridType*>::iterator it(GridObjects.find(gridID));
    if(it != GridObjects.end())
    {
      app_log() << "  Reuse " << gridID << " grid" << endl;
      return (*it).second;
    }
  }
  if(gridType == "log")
  {
    if(ascale>0.0)
    {
      LOGMSG("Using log grid with default values: scale = " << ascale << " step = " << astep << " npts = " << npts)
      agrid = new LogGridZero<RealType>;
      agrid->set(astep,ascale,npts);
    }
    else
    {
      LOGMSG("Using log grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
      if(ri<numeric_limits<RealType>::epsilon())
      {
        ri=numeric_limits<RealType>::epsilon();
        app_error() << "   LogGrid cannot accept r=0 for the initial point. Using ri=" << ri << endl;
      }
      agrid = new LogGrid<RealType>;
      agrid->set(ri,rf,npts);
    }
  }
  else
    if(gridType == "linear")
    {
      LOGMSG("Using linear grid with default values: ri = " << ri << " rf = " << rf << " npts = " << npts)
      agrid = new LinearGrid<RealType>;
      agrid->set(ri,rf,npts);
    }
  if(hasName)
  {
    GridObjects[gridID]=agrid;
  }
  else
  {
    char gname[16];
    int s=GridObjects.size();
    sprintf(gname,"g1_%d",s);
    GridObjects[gname]=agrid;
  }
  return agrid;
}

OneDimGridFactory::RealType
OneDimGridFactory::setSmoothCutoff(GridType* agrid, xmlNodePtr cur)
{
  //first create one if none
  if(agrid == 0)
    agrid = createGrid(cur);
  RealType rmax=agrid->rmax();
  RealType rmin=agrid->rmin();
  //This should check the targetPtcl::Lattice
  RealType rcut=rmax+1.0; //rcut is set to larget than
  if(cur != NULL)
  {
    const xmlChar* rcPtr = xmlGetProp(cur,(const xmlChar*)"rc");
    if(rcPtr != NULL)
    {
      rcut=atof((const char*)rcPtr);
    }
  }
  if(rcut>rmax)
    //set it to 0.99
  {
    rcut = (rmax-rmin)*0.99+rmin;
  }
  return rcut;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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

#include <QMCHamiltonians/ReferencePoints.h>
#include <Utilities/string_utils.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

bool ReferencePoints::put(xmlNodePtr cur, ParticleSet& P, vector<ParticleSet*>& Pref)
{
  app_log()<<"  Entering ReferencePoints::put"<<endl;
  bool succeeded=true;
  put(P,Pref);
  OhmmsAttributeSet ra;
  string coord="";
  ra.add(coord,"coord");
  ra.put(cur);
  for(int i=0; i<DIM; i++)
    for(int d=0; d<DIM; d++)
      axes(d,i)=P.Lattice.a(i)[d];
  Tensor_t crd;
  if(coord=="cell")
  {
    coordinate=cellC;
    crd = axes;
  }
  else
    if(coord=="cartesian")
    {
      coordinate=cartesianC;
      for(int i=0; i<DIM; i++)
        for(int d=0; d<DIM; d++)
          if(d==i)
            crd(i,i)=1.0;
          else
            crd(d,i)=0.0;
    }
    else
    {
      app_log() << endl;
      app_log() << "    Valid coordinates must be provided for element reference_points." << endl;
      app_log() << "      You provided: "<<coord<<endl;
      app_log() << "      Options are cell or cartesian." << endl;
      app_log() << endl;
      succeeded=false;
    }
  //read in the point contents
  app_log()<<"    reading reference_points contents"<<endl;
  string contents = (const char*)(xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1));
  vector<string> lines = split(strip(contents),"\n");
  for(int i=0; i<lines.size(); i++)
  {
    vector<string> tokens = split(strip(lines[i])," ");
    if(tokens.size()!=DIM+1)
    {
      app_log() << "  reference point has 4 entries, given "<< tokens.size() <<": "<<lines[i]<< endl;
      succeeded=false;
    }
    else
    {
      Point rp;
      for(int d=0; d<DIM; d++)
      {
        rp[d]=string2real(tokens[d+1]);
      }
      rp = dot(crd,rp);
      points[tokens[0]]=rp;
    }
  }
  return succeeded;
}

bool ReferencePoints::put(ParticleSet& P, vector<ParticleSet*>& Psets)
{
  //get axes and origin information from the ParticleSet
  points["zero"]  = 0*P.Lattice.a(0);
  points["a1"]    =   P.Lattice.a(0);
  points["a2"]    =   P.Lattice.a(1);
  points["a3"]    =   P.Lattice.a(2);
  //points["center"]= .5*(P.Lattice.a(0)+P.Lattice.a(1)+P.Lattice.a(2))
  //set points on face centers
  points["f1p"] = points["zero"]+.5*points["a1"];
  points["f1m"] = points["zero"]-.5*points["a1"];
  points["f2p"] = points["zero"]+.5*points["a2"];
  points["f2m"] = points["zero"]-.5*points["a2"];
  points["f3p"] = points["zero"]+.5*points["a3"];
  points["f3m"] = points["zero"]-.5*points["a3"];
  //set points on cell corners
  points["cmmm"]=points["zero"]+.5*( -1*points["a1"] -points["a2"] -points["a3"]);
  points["cpmm"]=points["zero"]+.5*(    points["a1"] -points["a2"] -points["a3"]);
  points["cmpm"]=points["zero"]+.5*( -1*points["a1"] +points["a2"] -points["a3"]);
  points["cmmp"]=points["zero"]+.5*( -1*points["a1"] -points["a2"] +points["a3"]);
  points["cmpp"]=points["zero"]+.5*( -1*points["a1"] +points["a2"] +points["a3"]);
  points["cpmp"]=points["zero"]+.5*(    points["a1"] -points["a2"] +points["a3"]);
  points["cppm"]=points["zero"]+.5*(    points["a1"] +points["a2"] -points["a3"]);
  points["cppp"]=points["zero"]+.5*(    points["a1"] +points["a2"] +points["a3"]);
  //get points from requested particle sets
  int cshift=1;
  for(int i=0; i<Psets.size(); i++)
  {
    ParticleSet& PS = *Psets[i];
    for(int p=0; p<PS.getTotalNum(); p++)
    {
      stringstream ss;
      ss<<p+cshift;
      points[PS.getName()+ss.str()]=PS.R[p];
    }
  }
  return true;
}


void ReferencePoints::write_description(ostream& os, string& indent)
{
  os<< indent+"reference_points"  <<endl;
  map<string,Point>::const_iterator it,end=points.end();
  for(it=points.begin(); it!=end; ++it)
  {
    os<< indent+"  " << it->first << ": " << it->second <<endl;
  }
  os<< indent+"end reference_points"  <<endl;
  return;
}

void ReferencePoints::save(vector<observable_helper*>& h5desc, hid_t gid) const
{
  observable_helper* oh = new observable_helper("reference_points");
  oh->open(gid);
  map<string,Point>::const_iterator it;
  for (it=points.begin(); it != points.end(); ++it)
  {
    oh->addProperty(const_cast<Point&>(it->second),it->first);
  }
  h5desc.push_back(oh);
  return;
}


}

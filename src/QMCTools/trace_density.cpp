//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file hyperslab.cpp
 * @brief Test code to read a chunk of data
 */
#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/FastParticleOperators.h>
#include <ParticleIO/ParticleLayoutIO.h>
//#include <OhmmsPETE/OhmmsVector.h>
//#include <OhmmsPETE/OhmmsMatrix.h>
//#include <algorithm>
#include <io/hdf_archive.h>
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"

/** xml input file for density
<?xml version="1.0"?>
<qmcsystem>
  <simulationcell>
    <parameter name="scale">30 </parameter>
    <parameter name="lattice">
      1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0
    </parameter>
    <parameter name="bconds">p p p</parameter>
  </simulationcell>
  <density file_root="XYZ" first="0" last="1">
  </density>
</qmcsystem>
*/

using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  Communicate* comm=OHMMS::Controller;
  OhmmsInfo Welcome("density",comm->rank());
  Random.init(0,1,11);

  ParticleSet::ParticleLayout_t simcell;
  string trace_root;
  int node_start=0, node_end=0;
  {
    Libxml2Document adoc;
    bool success = adoc.parse(argv[1]);
    xmlNodePtr cur=adoc.getRoot()->children;
    while(cur != NULL)
    {
      string cname((const char*)cur->name);
      if(cname.find("cell")<cname.size())
      {
        LatticeParser a(simcell);
        a.put(cur);
      }
      else if(cname =="density")
      {
        OhmmsAttributeSet a;
        a.add(trace_root,"file_root");
        a.add(node_start,"start");
        a.add(node_end,"end");
        a.put(cur);
      }
      cur=cur->next;
    }
  }

  //use proper root and extensition convention
  string h5fname=trace_root+".h5";

  int numptcls=0;
  vector<int> id;
  vector<double> pos;
  vector<double> weight;

  { //will be made into a class with id,weight and the data to analyze, e.g. pos
    hdf_archive hout(comm);
    hout.open(h5fname,H5F_ACC_RDONLY);

    ///read id
    hsize_t dims[2];
    bool g=hout.push("/int_data",false);
    if(g)
    {
      hid_t gid=hout.top();
      bool gotit=get_space(gid,"traces",2,dims);
    }
    {
      id.resize(dims[0]*dims[1]);
      TinyVector<int,2> gcounts(dims[0],dims[1]);
      hyperslab_proxy<vector<int>, 2> id_slab(id,gcounts);
      hout.read(id_slab,"traces");
    }
    hout.pop();

    ///read weight
    int nblocks=dims[0];
    int blockoffset=0;
    weight.resize(nblocks);
    {
      TinyVector<int,2> gcounts(dims[0],dims[1]);
      TinyVector<int,2> counts(dims[0],1);
      TinyVector<int,2> offsets(blockoffset,0);
      hyperslab_proxy<vector<double>, 2> w_slab(weight,gcounts,counts,offsets);
      hout.read(w_slab,"traces");
    }

    g=hout.push("/real_data",false);
    if(g)
    {
      hid_t gid=hout.top();
      bool gotit=get_space(gid,"traces",2,dims);
    }
    int pos_start=3;
    int pos_end=7;
    hout.read(pos_start,"layout/e/position/row_start");
    hout.read(pos_end,"layout/e/position/row_end");
    int npos=pos_end-pos_start;
    numptcls=npos/3;
    Timer now;
    now.restart();
    pos.resize(nblocks*npos);
    {
      TinyVector<int,2> gcounts(dims[0],dims[1]);
      TinyVector<int,2> counts(nblocks,npos);
      TinyVector<int,2> offsets(blockoffset,pos_start);
      hyperslab_proxy<vector<double>, 2> pos_slab(pos,gcounts,counts,offsets);
      hout.read(pos_slab,"traces");
    }
    int mbytes=(pos_start*pos_end*sizeof(double))>>20;

    app_log() << "Time to read positions " << now.elapsed() << endl;
    app_log() << "pos_start = " << pos_start << " pos_end = " << pos_end << " mbytes " << mbytes<< endl;

    hout.pop();
  }

  typedef ParticleSet::ParticlePos_t ParticlePos_t;
  typedef ParticleSet::Tensor_t Tensor_t;
  ParticlePos_t R(numptcls); //cartesian coordiates
  ParticlePos_t Ru(numptcls);//unit coordinates [0,1)
  std::copy(pos.data(),pos.data()+3*numptcls, &(R[0][0]));
  //convert Cartesian to unit in a box
  ApplyBConds<ParticlePos_t,Tensor_t,3,false>::Cart2Unit(R,simcell.G,Ru,0,numptcls);

  for(int i=0; i<8; ++i)
    app_log() << R[i] << " " << Ru[i] << endl;

  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

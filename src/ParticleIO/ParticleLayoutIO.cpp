//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include "OhmmsData/FileUtility.h"
using namespace std;
#include "Utilities/OhmmsInfo.h"
#include "ParticleIO/ParticleLayoutIO.h"

namespace ohmmsqmc {

  bool LatticeParser::put(xmlNodePtr cur){

    double a0 = 1.0;

    cur = cur->xmlChildrenNode;
    typedef ParticleLayout_t::SingleParticleIndex_t SingleParticleIndex_t;
    vector<SingleParticleIndex_t> grid(3,SingleParticleIndex_t(1));

    TinyVector<string,OHMMS_DIM> bconds("p");
    //int mpigrid = ParticleLayout_t::MPI_GRID;
    //int ompgrid = ParticleLayout_t::OMP_GRID;
    //int finegrid =  ParticleLayout_t::SPATIAL_GRID;
    int mpigrid = 0;
    int ompgrid = 1;
    int finegrid = 2;

    while (cur != NULL) {
      if(!xmlStrcmp(cur->name, (const xmlChar *)"parameter")) {
        string aname((const char*)(xmlGetProp(cur, (const xmlChar *) "name")));
	if(aname == "scale") {
	  putContent(a0,cur);
	} else if(aname == "lattice") {
	  putContent(ref_.R,cur);
	} else if(aname == "grid") {
	  putContent(grid[2],cur);
	} else if(aname == "omgrid") {
	  putContent(grid[1],cur);
	} else if(aname == "mpigrid") {
	  putContent(grid[0],cur);
	} else if(aname == "bconds") {
	  putContent(bconds,cur);
	}
      }
      cur = cur->next;
    }

    for(int idir=0;idir<OHMMS_DIM; idir++) {
      char b = bconds[idir][0];
      if(b == 'n' || b == 'N') {
	ref_.BoxBConds[idir] = false;
      }
    }

    ref_.R *= a0;
    ref_.reset();
    //@todo need to add UniformGridLayout features
    //ref_.makeGrid(grid);
    return true;
  }


  bool LatticeXMLWriter::get(ostream& os) const {
    os<< "<unitcell>" << endl;
    os<< "<parameter name=\"lattice\" datatype=\"tensor\">" << endl;
    os << ref_.R << "</parameter>" << endl;
    os << "<parameter name=\"bconds\">";
    for(int idir=0; idir<OHMMS_DIM; idir++) {
      if(ref_.BoxBConds[idir])
	os << "p ";
      else
	os << "n ";
    }
    os << "</parameter>" << endl;
    ///only write the spatial grid but may choose to write mpi and openmp
    os << "<parameter name=\"grid\">";
    //@todo write grid-partition
    for(int idir=0; idir<OHMMS_DIM; idir++) {
      os << "1 "; //os <<ref_.size(idir) << " "; 
    }
    os << "</parameter>" << endl;
    os << "</unitcell>" << endl;
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

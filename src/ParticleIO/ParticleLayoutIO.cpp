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

    typedef ParticleLayout_t::SingleParticleIndex_t SingleParticleIndex_t;
    vector<SingleParticleIndex_t> grid(3,SingleParticleIndex_t(1));

    TinyVector<string,OHMMS_DIM> bconds("p");
    std::size_t finegrid =  ParticleLayout_t::SPATIAL_GRID;
    std::size_t ompgrid = ParticleLayout_t::OMP_GRID;
    std::size_t mpigrid = ParticleLayout_t::MPI_GRID;

    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "parameter") {
        string aname((const char*)(xmlGetProp(cur, (const xmlChar *) "name")));
	if(aname == "scale") {
	  putContent(a0,cur);
	} else if(aname == "lattice") {
	  putContent(ref_.R,cur);
	} else if(aname == "grid") {
	  putContent(grid[finegrid],cur);
	} else if(aname == "ompgrid") {
	  putContent(ref_.Grid[ompgrid],cur);
	  //putContent(grid[ompgrid],cur);
          //ref_.Grid[ompgrid]=grid[ompgrid];
	} else if(aname == "mpigrid") {
	  putContent(ref_.Grid[mpigrid],cur);
	  //putContent(grid[mpigrid],cur);
          //ref_.Grid[mpigrid]=grid[mpigrid];
	} else if(aname == "bconds") {
	  putContent(bconds,cur);
	} else if(aname == "LR_dim_cutoff") {
	  putContent(ref_.LR_dim_cutoff,cur);
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
    ref_.makeGrid(grid);
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
    std::size_t finegrid =  ParticleLayout_t::SPATIAL_GRID;
    os << "<parameter name=\"grid\">";
    for(int idir=0; idir<OHMMS_DIM; idir++) {
      os <<ref_.getGrid(finegrid)->size(idir) << " "; 
    }
    os << "</parameter>" << endl;
    //os << "<parameter name=\"omega\">" << ref_.ABC << "</parameter>" << endl;
    os << "</unitcell>" << endl;
    return true;
  }

  xmlNodePtr LatticeXMLWriter::createNode() {
    xmlNodePtr cur = xmlNewNode(NULL,(const xmlChar*)"unitcell");
    std::ostringstream l;
    l.setf(ios_base::scientific);
    l.precision(12);
    l  << ref_.R;
    xmlNodePtr p=xmlNewTextChild(cur,NULL,
        (const xmlChar*)"parameter", (const xmlChar*)l.str().c_str());
    xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)"lattice");
    return cur;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

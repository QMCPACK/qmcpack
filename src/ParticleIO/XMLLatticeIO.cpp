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
using namespace std;
#include "Configuration.h"
#include "OhmmsData/FileUtility.h"
#include "OhmmsData/libxmldefs.h"
#include "Utilities/SpeciesCollection.h"
#include "ParticleIO/XMLLatticeIO.h"

using namespace xmlpp;
using namespace ohmmsqmc;

///reading from a file
bool XMLLatticeIO::parse(const char* fname) {

  DomParser parser;
  parser.parse_file(fname);

  const Node* root
    = parser.get_document()->get_root_node(); //deleted by DomParser.

  NodeSet pset = root->find("//UnitCell");

  for(int i=0; i<pset.size(); i++) put(pset[i]);

  return true;
}

bool XMLLatticeIO::put(Node* pnode) {

  RealType a0 = 1.0;

  Node::NodeList param 
    = (dynamic_cast<Element*>(pnode))->get_children("PARAMETER");

  typedef TinyVector<IndexType,DIM> SingleParticleIndex_t;

  vector<SingleParticleIndex_t> grid(3,SingleParticleIndex_t(1));
  TinyVector<string,DIM> bconds("p");

  Node::NodeList::iterator it = param.begin();
  while(it != param.end()) {
    Element* a = dynamic_cast<Element*>(*it);
    const string& aname = a->get_attribute("name")->get_value();
    _xmlNode* cur = a->cobj();
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
    it++;  
  }

  for(int idir=0;idir<DIM; idir++) {
    char b = bconds[idir][0];
    if(b == 'n' || b == 'N') {
      //ref_.BConds[idir] = ParticleNoBCond;
      //ref_.BConds.wrapper(idir) = ParticleNoBCond;
      ref_.BoxBConds[idir] = false;
    }
  }

  ref_.R *= a0;
  ref_.reset();

  //PartitionGrid(ref_,grid);
  return true;
}

bool XMLLatticeIO::get(ostream& os) {
  os<< "<UnitCell>" << endl;
  os<< "<PARAMETER name=\"lattice\">" << endl;
  os << ref_.R << "</PARAMETER>" << endl;
  os << "<PARAMETER name=\"bconds\">";
  for(int idir=0; idir<DIM; idir++) {
    if(ref_.BoxBConds[idir])
      os << "p ";
    else
      os << "n ";
  }
  os << "</PARAMETER>" << endl;
  os << "</UnitCell>" << endl;
  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

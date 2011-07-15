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
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"

namespace qmcplusplus {

  bool LatticeParser::put(xmlNodePtr cur){

    const int DIM=ParticleLayout_t::SingleParticlePos_t::Size;
    double a0 = 1.0;

    double rs=-1.0;

    typedef ParticleLayout_t::SingleParticleIndex_t SingleParticleIndex_t;
    vector<SingleParticleIndex_t> grid(DIM,SingleParticleIndex_t(1));

    TinyVector<string,DIM> bconds("p");
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
	} else if(aname == "mpigrid") {
	  putContent(ref_.Grid[mpigrid],cur);
	} else if(aname == "bconds") {
	  putContent(bconds,cur);
          int boxsum=0;
          for(int idir=0;idir<DIM; idir++) {
            char b = bconds[idir][0];
            if(b == 'n' || b == 'N') {
              ref_.BoxBConds[idir] = false;
            } else {
              ref_.BoxBConds[idir] = true;
              boxsum++;
            }
          }
          if(boxsum>0 && boxsum<DIM)
          {
            APP_ABORT(" LatticeParser::put \n   Mixed boundary is not supported. Set \n   <parameter name=\"bconds\">p p p </parameter>\n");
          }
	} else if(aname == "LR_dim_cutoff") {
	  putContent(ref_.LR_dim_cutoff,cur);
	} else if(aname == "rs") {
          int nptcl=0;
          int nsh=0;
          int pol=0;
          OhmmsAttributeSet rAttrib;
          rAttrib.add(nptcl,"condition");
          rAttrib.add(pol,"polarized");
          rAttrib.add(nsh,"shell");
          rAttrib.put(cur);
	  putContent(rs,cur);
          HEGGrid<double,OHMMS_DIM> heg(ref_);
          if(pol==0)
          {
            if(nsh>0)
              nptcl=2*heg.getNumberOfKpoints(nsh);
            else
              nsh=heg.getShellIndex(nptcl/2);
          }
          else
          {
//             spin polarized
            if(nsh>0)
              nptcl=heg.getNumberOfKpoints(nsh);
            else
              nsh=heg.getShellIndex(nptcl);
          }
          double acubic=heg.getCellLength(nptcl,rs);
          //double acubic=pow(4.0*M_PI*nptcl/3.0,1.0/3.0)*rs;
          app_log() << "  " << OHMMS_DIM << "D HEG system"
            << "\n     rs  = " << rs;
            
            if(pol==0)
            {
              app_log() << "\n     number of up particles = " << nptcl/2
              << "\n     number of dn particles = " << nptcl/2 ;
            }
            else
            {
              app_log() << "\n     number of up particles = " << nptcl;
            }
            app_log()<< "\n     filled kshells      = " << nsh 
            << "\n     lattice constant    = " << acubic << " bohr"<< endl;
          ref_.R=0.0;
          for(int idim=0; idim<DIM; idim++)
            for(int jdim=0; jdim<DIM; jdim++)
              if (idim==jdim) ref_.R(idim,jdim)=acubic;
              else ref_.R(idim,jdim)=0.0;
          a0=1.0;
        }
      } 
      cur = cur->next;
    }

    ref_.R *= a0;
    ref_.reset();
    ref_.makeGrid(grid);
    app_log() << std::fixed;
    app_log() << "  Simulation cell radius = " << ref_.SimulationCellRadius << endl;
    app_log() << "  Wigner-Seitz    radius = " << ref_.WignerSeitzRadius    << endl;

    return true;
  }


  bool LatticeXMLWriter::get(ostream& os) const {
    os<< "<unitcell>" << endl;
    os<< "<parameter name=\"lattice\" datatype=\"tensor\">" << endl;
    os << ref_.R << "</parameter>" << endl;
    os << "<parameter name=\"bconds\">";
    const int DIM=ParticleLayout_t::SingleParticlePos_t::Size;
    for(int idir=0; idir<DIM; idir++) {
      if(ref_.BoxBConds[idir])
	os << "p ";
      else
	os << "n ";
    }
    os << "</parameter>" << endl;
    ///only write the spatial grid but may choose to write mpi and openmp
    std::size_t finegrid =  ParticleLayout_t::SPATIAL_GRID;
    os << "<parameter name=\"grid\">";
    for(int idir=0; idir<DIM; idir++) {
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

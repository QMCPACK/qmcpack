//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include "OhmmsData/FileUtility.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"

namespace qmcplusplus
{

bool LatticeParser::put(xmlNodePtr cur)
{
  const int DIM=ParticleLayout_t::SingleParticlePos_t::Size;
  double a0 =  1.0;
  double rs = -1.0;
  int nptcl =  0;  
  int nsh   =  0; //for backwards compatibility w/ odd heg initialization style
  int pol   =  0;
  typedef ParticleLayout_t::SingleParticleIndex_t SingleParticleIndex_t;
  std::vector<SingleParticleIndex_t> grid(DIM,SingleParticleIndex_t(1));
  TinyVector<std::string,DIM> bconds("p");
  std::size_t finegrid =  ParticleLayout_t::SPATIAL_GRID;
  std::size_t ompgrid = ParticleLayout_t::OMP_GRID;
  std::size_t mpigrid = ParticleLayout_t::MPI_GRID;

  Tensor<OHMMS_PRECISION_FULL,DIM> lattice_in;
  bool lattice_defined=false;
  bool bconds_defined=false;
  int boxsum=0;

  app_log() << " Lattice" << std::endl;
  app_log() << " -------" << std::endl;
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "parameter")
    {
      std::string aname((const char*)(xmlGetProp(cur, (const xmlChar *) "name")));
      if(aname == "scale")
      {
        putContent(a0,cur);
      }
      else if(aname == "lattice")
      {
        const char *units_prop = (const char *)(xmlGetProp(cur, (const xmlChar *) "units"));
        if (units_prop && std::string(units_prop) != "bohr") {
          APP_ABORT("LatticeParser::put. Only atomic units (bohr) supported for lattice units. Input file uses: " << std::string(units_prop));
        }

        putContent(lattice_in,cur);
        lattice_defined=true;
        //putContent(ref_.R,cur);
      }
      else if(aname == "grid")
      {
        putContent(grid[finegrid],cur);
      }
      else if(aname == "ompgrid")
      {
        putContent(ref_.Grid[ompgrid],cur);
      }
      else if(aname == "mpigrid")
      {
        putContent(ref_.Grid[mpigrid],cur);
      }
      else if(aname == "bconds")
      {
        putContent(bconds,cur);
        bconds_defined=true;
        for(int idir=0; idir<DIM; idir++)
        {
          char b = bconds[idir][0];
          if(b == 'n' || b == 'N')
          {
            ref_.BoxBConds[idir] = false;
          }
          else
          {
            ref_.BoxBConds[idir] = true;
            boxsum++;
          }
        }
      }
      else if(aname == "vacuum")
      {
        putContent(ref_.VacuumScale,cur);
      }
      else if(aname == "LR_dim_cutoff")
      {
        putContent(ref_.LR_dim_cutoff,cur);
      }
      else if(aname == "rs")
      {
        lattice_defined=true;
        OhmmsAttributeSet rAttrib;
        rAttrib.add(nptcl,"condition");
        rAttrib.add(pol,"polarized");
        rAttrib.add(nsh,"shell");
        rAttrib.put(cur);
        putContent(rs,cur);
      }
      else if(aname == "nparticles")
      {
        putContent(nptcl,cur);
      }
    }
    cur = cur->next;
  }
  // checking boundary conditions
  if(lattice_defined)
  {
    if(!bconds_defined)
    {
      app_log() << "  Lattice is specified but boundary conditions are not. Assuming PBC." << std::endl;
      ref_.BoxBConds = true;
    }
  }
  else
  {
    if(boxsum==0)
    {
      app_log() << "  Lattice is not specified for the Open BC. Add a huge box." << std::endl;
      lattice_in=0;
      for(int idir=0; idir<DIM; idir++)
        lattice_in(idir,idir)=1e5;
    }
    else
    {
      APP_ABORT(" LatticeParser::put \n   Mixed boundary is supported only when a lattice is specified!");
    }
  }
  //special heg processing
  if(rs>0.0)
  {
    HEGGrid<ParticleLayout_t::Scalar_t,OHMMS_DIM> heg(ref_);
    if(pol==0)
    {
      if(nsh>0)
        nptcl=2*heg.getNumberOfKpoints(nsh);
      else
        nsh=heg.getShellIndex(nptcl/2);
    }
    else
    { //             spin polarized
      if(nsh>0)
        nptcl=heg.getNumberOfKpoints(nsh);
      else
        nsh=heg.getShellIndex(nptcl);
    }
    ParticleLayout_t::Scalar_t acubic=heg.getCellLength(nptcl,rs);
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
             << "\n     lattice constant    = " << acubic << " bohr"<< std::endl;
    lattice_in=0.0;
    for(int idim=0; idim<DIM; idim++)
      lattice_in(idim,idim)=acubic;
    a0=1.0;
  }

  lattice_in *= a0;
  ref_.set(lattice_in);
  ref_.makeGrid(grid);
  if(ref_.SuperCellEnum == SUPERCELL_OPEN)
    ref_.WignerSeitzRadius=ref_.SimulationCellRadius;
  std::string unit_name = "bohr";
  app_log() << std::fixed;
  app_log() << "  Simulation cell radius   = " << ref_.SimulationCellRadius << " " << unit_name << std::endl;
  app_log() << "  Wigner-Seitz cell radius = " << ref_.WignerSeitzRadius    << " " << unit_name << std::endl;
  app_log() << std::endl;

  //initialize the global cell
  //qmc_common.theSuperCell=lattice_in;
  return lattice_defined;
}


bool LatticeXMLWriter::get(std::ostream& os) const
{
  os<< "<unitcell>" << std::endl;
  os<< "<parameter name=\"lattice\" datatype=\"tensor\">" << std::endl;
  os << ref_.R << "</parameter>" << std::endl;
  os << "<parameter name=\"bconds\">";
  const int DIM=ParticleLayout_t::SingleParticlePos_t::Size;
  for(int idir=0; idir<DIM; idir++)
  {
    if(ref_.BoxBConds[idir])
      os << "p ";
    else
      os << "n ";
  }
  os << "</parameter>" << std::endl;
  ///only write the spatial grid but may choose to write mpi and openmp
  std::size_t finegrid =  ParticleLayout_t::SPATIAL_GRID;
  os << "<parameter name=\"grid\">";
  for(int idir=0; idir<DIM; idir++)
  {
    os <<ref_.getGrid(finegrid)->size(idir) << " ";
  }
  os << "</parameter>" << std::endl;
  //os << "<parameter name=\"omega\">" << ref_.ABC << "</parameter>" << std::endl;
  os << "</unitcell>" << std::endl;
  return true;
}

xmlNodePtr LatticeXMLWriter::createNode()
{
  xmlNodePtr cur = xmlNewNode(NULL,(const xmlChar*)"unitcell");
  std::ostringstream l;
  l.setf(std::ios_base::scientific);
  l.precision(12);
  l  << ref_.R;
  xmlNodePtr p=xmlNewTextChild(cur,NULL,
                               (const xmlChar*)"parameter", (const xmlChar*)l.str().c_str());
  xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)"lattice");
  return cur;
}
}

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
#include "LatticeIO.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "LongRange/LRCoulombSingleton.h"
#include "ModernStringUtils.hpp"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
bool LatticeParser::put(xmlNodePtr cur)
{
  const int DIM             = ParticleLayout::SingleParticlePos::Size;
  double a0                 = 1.0;
  double rs                 = -1.0;
  int nptcl                 = 0;
  int nsh                   = 0; //for backwards compatibility w/ odd heg initialization style
  int pol                   = 0;
  using SingleParticleIndex = ParticleLayout::SingleParticleIndex;
  TinyVector<std::string, DIM> bconds("p");

  Tensor<OHMMS_PRECISION_FULL, DIM> lattice_in;
  bool lattice_defined = false;
  bool bconds_defined  = false;
  int boxsum           = 0;

  app_summary() << std::endl;
  app_summary() << " Lattice" << std::endl;
  app_summary() << " -------" << std::endl;
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "parameter")
    {
      const std::string aname(getXMLAttributeValue(cur, "name"));
      if (aname == "scale")
      {
        putContent(a0, cur);
      }
      else if (aname == "lattice")
      {
        const std::string units_prop(getXMLAttributeValue(cur, "units"));
        if (!units_prop.empty() && units_prop != "bohr")
        {
          std::ostringstream err_msg;
          err_msg << "LatticeParser::put. Only atomic units (bohr) supported for lattice units. Input file uses: "
                  << units_prop;
          throw UniformCommunicateError(err_msg.str());
        }

        putContent(lattice_in, cur);
        lattice_defined = true;
        //putContent(ref_.R,cur);
      }
      else if (aname == "bconds")
      {
        putContent(bconds, cur);
        bconds_defined = true;
        for (int idir = 0; idir < DIM; idir++)
        {
          char b = bconds[idir][0];
          if (b == 'n' || b == 'N')
          {
            ref_.BoxBConds[idir] = false;
          }
          else if (b == 'p' || b == 'P')
          {
            ref_.BoxBConds[idir] = true;
            boxsum++;
          }
          else
          {
            std::ostringstream err_msg;
            err_msg << "LatticeParser::put. Unknown label '" + bconds[idir] +
                    "' used for periodicity. Only 'p', 'P', 'n' and 'N' are valid!";
            throw UniformCommunicateError(err_msg.str());
          }

          // Protect BCs which are not implemented.
          if (idir > 0 && !ref_.BoxBConds[idir - 1] && ref_.BoxBConds[idir])
          {
            std::ostringstream err_msg;
            err_msg
                << "LatticeParser::put. In \"bconds\", non periodic directions must be placed after the periodic ones.";
            throw UniformCommunicateError(err_msg.str());
          }
        }
      }
      else if (aname == "vacuum")
      {
        putContent(ref_.VacuumScale, cur);
      }
      else if (aname == "LR_dim_cutoff")
      {
        putContent(ref_.LR_dim_cutoff, cur);
      }
      else if (aname == "LR_handler")
      {
        std::string handler_type("opt_breakup");
        //This chops whitespace so the simple str == comparisons work
        putContent(handler_type, cur);
        handler_type = lowerCase(handler_type);
        if (handler_type == "ewald")
          LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
        else if (handler_type == "opt_breakup")
          LRCoulombSingleton::this_lr_type = LRCoulombSingleton::ESLER;
        else if (handler_type == "opt_breakup_original")
          LRCoulombSingleton::this_lr_type = LRCoulombSingleton::NATOLI;
        else if (handler_type == "ewald_strict2d")
        {
          LRCoulombSingleton::this_lr_type = LRCoulombSingleton::STRICT2D;
          ref_.ndim                        = 2;
        }
        else if (handler_type == "ewald_quasi2d")
          LRCoulombSingleton::this_lr_type = LRCoulombSingleton::QUASI2D;
        else
          throw UniformCommunicateError("LatticeParser::put. Long range breakup handler not recognized.");
      }
      else if (aname == "LR_tol")
      {
        putContent(ref_.LR_tol, cur);
      }
      else if (aname == "rs")
      {
        lattice_defined = true;
        OhmmsAttributeSet rAttrib;
        rAttrib.add(nptcl, "condition");
        rAttrib.add(pol, "polarized");
        rAttrib.add(nsh, "shell");
        rAttrib.put(cur);
        putContent(rs, cur);
      }
      else if (aname == "nparticles")
      {
        putContent(nptcl, cur);
      }
    }
    cur = cur->next;
  }

  // checking boundary conditions
  if (lattice_defined)
  {
    if (!bconds_defined)
    {
      app_log() << "  Lattice is specified but boundary conditions are not. Assuming PBC." << std::endl;
      ref_.BoxBConds = true;
    }
  }
  else if (boxsum == 0)
    app_log() << "  Lattice is not specified for the Open BC. Add a huge box." << std::endl;
  else
    throw UniformCommunicateError("LatticeParser::put. Mixed boundary is supported only when a lattice is specified!");

  //special heg processing
  if (rs > 0.0)
  {
    HEGGrid<ParticleLayout::Scalar_t> heg(ref_);
    if (pol == 0)
    {
      if (nsh > 0)
        nptcl = 2 * heg.getNumberOfKpoints(nsh);
      else
        nsh = heg.getShellIndex(nptcl / 2);
    }
    else
    { //             spin polarized
      if (nsh > 0)
        nptcl = heg.getNumberOfKpoints(nsh);
      else
        nsh = heg.getShellIndex(nptcl);
    }
    ParticleLayout::Scalar_t acubic = heg.getCellLength(nptcl, rs);
    app_log() << "  " << OHMMS_DIM << "D HEG system"
              << "\n     rs  = " << rs;
    if (pol == 0)
    {
      app_log() << "\n     number of up particles = " << nptcl / 2 << "\n     number of dn particles = " << nptcl / 2;
    }
    else
    {
      app_log() << "\n     number of up particles = " << nptcl;
    }
    app_log() << "\n     filled kshells      = " << nsh << "\n     lattice constant    = " << acubic << " bohr"
              << std::endl;
    lattice_in = 0.0;
    for (int idim = 0; idim < DIM; idim++)
      lattice_in(idim, idim) = acubic;
    a0 = 1.0;
  }

  if (lattice_defined)
  {
    lattice_in *= a0;
    ref_.set(lattice_in);
  }

  if (ref_.SuperCellEnum != SUPERCELL_SLAB && LRCoulombSingleton::isQuasi2D())
    throw UniformCommunicateError("LatticeParser::put. Quasi 2D Ewald only works with boundary condition 'p p n'!");

  if (ref_.SuperCellEnum == SUPERCELL_OPEN)
    ref_.WignerSeitzRadius = ref_.SimulationCellRadius;

  std::string unit_name = "bohr";
  app_log() << std::fixed;
  app_log() << "  Simulation cell radius   = " << ref_.SimulationCellRadius << " " << unit_name << std::endl;
  app_log() << "  Wigner-Seitz cell radius = " << ref_.WignerSeitzRadius << " " << unit_name << std::endl;
  app_log() << std::endl;

  return lattice_defined;
}


bool LatticeXMLWriter::get(std::ostream& os) const
{
  os << "<unitcell>" << std::endl;
  os << R"(<parameter name="lattice" datatype="tensor">)" << std::endl;
  os << ref_.R << "</parameter>" << std::endl;
  os << "<parameter name=\"bconds\">";
  const int DIM = ParticleLayout::SingleParticlePos::Size;
  for (int idir = 0; idir < DIM; idir++)
  {
    if (ref_.BoxBConds[idir])
      os << "p ";
    else
      os << "n ";
  }
  os << "</parameter>" << std::endl;
  os << "</unitcell>" << std::endl;
  return true;
}

xmlNodePtr LatticeXMLWriter::createNode()
{
  xmlNodePtr cur = xmlNewNode(NULL, (const xmlChar*)"unitcell");
  std::ostringstream l;
  l.setf(std::ios_base::scientific);
  l.precision(12);
  l << ref_.R;
  xmlNodePtr p = xmlNewTextChild(cur, NULL, (const xmlChar*)"parameter", (const xmlChar*)l.str().c_str());
  xmlNewProp(p, (const xmlChar*)"name", (const xmlChar*)"lattice");
  return cur;
}
} // namespace qmcplusplus

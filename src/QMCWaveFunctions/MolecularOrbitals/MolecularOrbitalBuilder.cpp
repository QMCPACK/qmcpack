//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"

namespace ohmmsqmc {

  bool MolecularOrbitalBuilder::put(xmlNodePtr cur) {
    
    bool usegrid = true;
    if(xmlHasProp(cur, (const xmlChar*)"usegrid")) {
      string a((const char*)(xmlGetProp(cur, (const xmlChar *)"usegrid")));
      if(a == "false" || a == "0") {
	usegrid = false;
      }
    }

    if(usegrid) {
      LOGMSG("Using radial grids for molecular orbitals")
      GridMolecularOrbitals a(wfs_ref,Ions,Els);
      a.put(cur);
    } else {
      LOGMSG("Using STO for molecular orbitals")
      STOMolecularOrbitals a(wfs_ref,Ions,Els);
      a.put(cur);
    }

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

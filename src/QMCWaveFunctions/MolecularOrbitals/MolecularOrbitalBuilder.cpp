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
    if(xmlHasProp(cur, (const xmlChar*)"src")){
      string fname_in = (const char*)(xmlGetProp(cur, (const xmlChar *)"src"));
      LOGMSG("Opening external file " << fname_in)
      putOpen(fname_in);
    } else {
      putSpecial(cur);
    }
    return true;
  }

  bool MolecularOrbitalBuilder::putSpecial(xmlNodePtr cur) {
    
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

  bool MolecularOrbitalBuilder::putOpen(const string& fname_in){

      xmlDocPtr doc=NULL;
      xmlNsPtr ns;
      xmlNodePtr cur;
      
      // build an XML tree from a the file;
      doc = xmlParseFile(fname_in.c_str());
      if (doc == NULL) {
	ERRORMSG(fname_in << " does not exist")
	  return false;
      } 

      ///using XPath instead of recursive search
      xmlXPathContextPtr context;
      // xmlXPathObjectPtr result;
      context = xmlXPathNewContext(doc);
      const char* name = OrbitalBuilderBase::detset_tag.c_str();
      string sdname("//");
      sdname.append(OrbitalBuilderBase::detset_tag);
      xmlXPathObjectPtr result
	= xmlXPathEvalExpression((const xmlChar*)(sdname.c_str()),context);
      //      result = xmlXPathEvalExpression((const xmlChar*)"//determinantset",context);
  
      if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
	
	ERRORMSG(fname_in << " does not contain any Data!")
	  
	  } else {
	    xmlNodePtr cur = result->nodesetval->nodeTab[0];
	    putSpecial(cur);
	  }

      //free local objects
      xmlXPathFreeObject(result);
      xmlFreeDoc(doc);
      
      return true;
    }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

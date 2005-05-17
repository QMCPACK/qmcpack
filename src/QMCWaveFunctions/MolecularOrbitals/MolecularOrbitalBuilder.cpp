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
#include "OhmmsData/AttributeSet.h"

namespace ohmmsqmc {

  bool MolecularOrbitalBuilder::put(xmlNodePtr cur) {
    if(xmlHasProp(cur, (const xmlChar*)"href")){
      string fname_in = (const char*)(xmlGetProp(cur, (const xmlChar *)"href"));
      LOGMSG("Opening external file " << fname_in)
      putOpen(fname_in);
    } else {
      putSpecial(cur);
    }
    return true;
  }

  /** process cur xml node
   *
   * Check attribute list
   * - transform if yes, using numerical grid
   * - source    the particle set providing the nuclei default is "i"
   *\xmlonly
   *  <determinantset type="MolecularOrbital" source="i" transform="yes"/>
   *\endxmlonly
   */
  bool MolecularOrbitalBuilder::putSpecial(xmlNodePtr cur) {
    
    bool usegrid = true;
    string transform("yes"), source("i");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(transform,"transform");
    aAttrib.add(source,"source");
    aAttrib.put(cur);

    if(transform == "no") 
      usegrid=false;
    else 
      usegrid=true;

    //ionic system for the molecualr orbitals
    ParticleSet* ions=0;

    //initialize with the source tag
    PtclPoolType::iterator pit(ptclPool.find(source));
    if(pit != ptclPool.end()) {
      LOGMSG("Molecular orbital uses an ionic system " << source)
      ions=(*pit).second; 
    }

    //check if distance table is used
    xmlNodePtr tcur=cur->children;
    while(ions == 0 && tcur != NULL) {
      if(xmlStrEqual(tcur->name,(const xmlChar*)"distancetable")) {
        const xmlChar*  a=xmlGetProp(tcur,(const xmlChar*)"source");
        if(a) { source = (const char*)a;}
        pit = ptclPool.find(source);
        if(pit != ptclPool.end()) ions = (*pit).second;
      }
      tcur=tcur->next;
    }

    if(ions == 0){
       ERRORMSG("Any molecular orbital needs an ionic system. Missing " << source)
       return false;
    }

    if(usegrid) {
      LOGMSG("Using radial grids for molecular orbitals for " << targetPtcl.getName() << " with " << ions->getName())
      GridMolecularOrbitals a(targetPtcl,targetPsi,*ions);
      a.put(cur);
    } else {
      LOGMSG("Using analytic STO molecular orbitals for " << targetPtcl.getName() << " with " << ions->getName())
      STOMolecularOrbitals a(targetPtcl,targetPsi,*ions);
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
      putSpecial(result->nodesetval->nodeTab[0]);
    }

    //free local objects
    xmlXPathFreeObject(result);
    xmlXPathFreeContext(context);
    xmlFreeDoc(doc);
    
    return true;
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

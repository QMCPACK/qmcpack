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
#include "OhmmsData/FileUtility.h"
#include "Utilities/OhmmsInfo.h"
//#include "Utilities/SpeciesCollection.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/HDFParticleIO.h"
#include "ParticleBase/ParticleFunctions.h"
using namespace ohmmsqmc;


XMLParticleParser::XMLParticleParser(Particle_t& aptcl, bool donotresize):
  AssignmentOnly(donotresize),
  ref_(aptcl)
{ 
}

/**reading particleset node from a file
 *@param fname_in a file name to open
 *@param pformat_in the format of the file, not used
 *@return true, if successful
 *
 *Check the type of the external source to work on.
 *The external source itself can be an xml file.
 */
bool XMLParticleParser::put(const string& fname_in, 
			    const string& fext_in) {

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
  xmlXPathObjectPtr result;
  context = xmlXPathNewContext(doc);
  result = xmlXPathEvalExpression((const xmlChar*)"//particleset",context);
  
  if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {

    ERRORMSG(fname_in << " does not contain any ParticleSet")

  } else {

    cur = result->nodesetval->nodeTab[0];
    string fname_ext(null_tag), pformat(null_tag);

    bool use_ext = false;
    ///process attributes: type or format
    xmlAttrPtr att = cur->properties;
    while(!use_ext && att != NULL) {
      string aname((const char*)(att->name));
      if(aname == "src" || aname=="href") {
	fname_ext = (const char*)(att->children->content);
	pformat = getExtension(fname_ext);
	use_ext = true;
      } else if(aname == "srctype") {
	pformat =(const char*)(att->children->content);
	use_ext = true;
      }
      att = att->next;
    }

    if(use_ext) {
      if(pformat == "h5") {
	HDFParticleParser ahandle(ref_);
	ahandle.put(cur);
      } else {
	ERRORMSG("Unknown file extension " << pformat << " of " << fname_ext)
      }
    } else {
      putSpecial(cur);
    }
  }

  //free local objects
  xmlXPathFreeObject(result);
  xmlXPathFreeContext(context);
  xmlFreeDoc(doc);

  return true;
}

/** process xmlnode &lt;particleset/&gt;
 *@param cur the xmlnode to work on
 *
 *If the node has src or href attribute, use an external file.
 */
bool XMLParticleParser::put(xmlNodePtr cur) {

  ///process attributes: type or format
  string fname("none"), pformat("xml");
  xmlAttrPtr att = cur->properties;
  while(att != NULL) {
    string aname((const char*)(att->name));
    if(aname == "src" || aname=="href") {
      fname = (const char*)(att->children->content);
      pformat = getExtension(fname);
    } else if(aname == "srctype") {
      pformat =(const char*)(att->children->content);
    }
    att = att->next;
  }

  if(fname == "none") 
    return putSpecial(cur);
  else
    return put(fname,pformat);
}


/** process xmlnode &lt;particleset/&gt; which contains everything about the particle set to initialize
 *@param cur the xmlnode to work on
 *
 */
bool XMLParticleParser::putSpecial(xmlNodePtr cur) {

  string pname("none");
  xmlDocPtr doc = cur->doc;

  //the number of particles that are initialized by <attrib/>
  int nat = 0;

  //process attributes of particleset
  xmlAttrPtr att = cur->properties;
  while(att != NULL) {
    string aname((const char*)(att->name));
    if(aname == "name") {
      pname = (const char*)(att->children->content);
    } else if(aname == "num" || aname == "size") {
      nat = atoi((const char*)(att->children->content));
    } 
    att = att->next;
  }
  
  ///count the number of atom added one by one
  xmlNodePtr cur0 = cur->xmlChildrenNode;

  //total count of the particles to be created
  int ntot = 0;
  int ng = 0;
  vector<int> nat_group;

  vector<xmlNodePtr> atom_ptr;

  //pre-process the nodes to count the number of particles to be added
  while(cur0 != NULL) {
    string cname((const char*)cur0->name);
    if(cname == "atom") {
      ntot++;	
      atom_ptr.push_back(cur0);
    } else if(cname == "group") {
      nat_group.push_back(0);
      if(xmlHasProp(cur0, (const xmlChar *) "size")) {
	nat_group[ng] = atoi((const char*)(xmlGetProp(cur0, (const xmlChar *) "size")));  
	ntot += nat_group[ng];
      }
      ng++;
    } else if(cname == attrib_tag) {
      int size_att = 0;
      if(xmlHasProp(cur0, (const xmlChar *) "size")) {
	size_att = atoi((const char*)(xmlGetProp(cur0, (const xmlChar *) "size")));  
      }
      if(size_att) {
	if(size_att != nat) {
	  WARNMSG("\tOverwriting the size of the particle by //particleset/attrib/@size=" << size_att)
	  nat = size_att;
	}
      }
    }
    cur0 = cur0->next;
  }

  ntot += nat;  
  ref_.setName(pname.c_str());
  int nloc = ref_.getLocalNum();

  //treat assignment only differently
  if(AssignmentOnly) {
    ntot = 0;
    nloc = 0;
    for(int iat=0;iat<ref_.getTotalNum(); iat++) ref_.ID[iat]=iat;
  }
   
  if(ntot) {
    LOGMSG("\tCreating " << ntot << " particles for " << pname << " set")
    ref_.create(ntot);
    //assign default ID
    int nloci=nloc;
    for(int iat=0;iat<ntot; iat++,nloci++) ref_.ID[iat]=nloci;
  }
  
  bool expand = false;
  TinyVector<int,OHMMS_DIM> uc_grid(1);      
  
  SpeciesSet* tspecies = &(ref_.Species); //SpeciesCollection::getSpecies();
  cur = cur->xmlChildrenNode;

  //reset the group counter
  ng = 0;  
  while (cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "UnitCell" || cname == "unitcell") {
      LatticeParser lat(ref_.Lattice);
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL) {
	if (!xmlStrcmp(tcur->name, (const xmlChar *)"parameter") &&
	    !xmlStrcmp(xmlGetProp(tcur, (const xmlChar *) "name"),
		       (const xmlChar*)"uc_grid")) {
	  expand = true;
	  putContent(uc_grid,tcur);
	}
	tcur = tcur->next;
      }
      lat.put(cur);
    } else if (cname == attrib_tag) {
      getPtclAttrib(cur,nat,nloc);
    } else  if (cname == "group") {
      //adding a group consiting of a number of atoms of the same species
      int sid
	= tspecies->addSpecies((const char*)
			       (xmlGetProp(cur, (const xmlChar *)"name")));

      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL) {
	string tcname((const char*)tcur->name);
	if(tcname == "parameter") {
	  int iproperty 
	    = tspecies->addAttribute((const char*)(xmlGetProp(tcur,(const xmlChar*)"name")));
	  double ap;
	  putContent(ap,tcur);
	  tspecies->operator()(iproperty,sid) = ap;
	} else if(nat_group[ng] && tcname == attrib_tag) { 
	  cout << "Adding groups " << nat_group[ng] << endl;
	  getPtclAttrib(tcur,nat_group[ng],nloc);
	}
	tcur = tcur->next;
      }
      for(int iat=0; iat<nat_group[ng]; iat++, nloc++)  ref_.GroupID[nloc] = sid;
      ng++;
    }
    cur = cur->next;
  }

  //have read from <attrib/>'s and <group/>'s. Time to add <atom/>'s
  nloc += nat;
  for(int ia=0; ia<atom_ptr.size(); ia++,nloc++) {
    cur = atom_ptr[ia];
    int sid
      = tspecies->addSpecies((const char*)
			      (xmlGetProp(cur, (const xmlChar *)"name")));
    att = cur->properties;
    int inunit = 0; //ref_.R.getUnit();
    while(att!= NULL) {
      string aname((const char*)(att->name));
      if(aname == condition_tag) {
	inunit =atoi((const char*)(att->children->content));
      }
      att = att->next;
    }
    
    Particle_t::SingleParticlePos_t pos;
    istringstream stream((const char*)
			 (xmlNodeListGetString(doc, cur->xmlChildrenNode, 1)));
    if(inunit == ref_.R.getUnit()) {
      stream >> ref_.R[nloc]; 
    } else {
      stream >> pos;
      if(inunit) 
	ref_.R[nloc] = ref_.Lattice.toCart(pos);
      else 
	ref_.R[nloc] = ref_.Lattice.toUnit(pos);
    }
    ref_.ID[nloc] = nloc; 
    ref_.GroupID[nloc] = sid;
  }

  if(expand) {
    ExpandSuperCell(ref_,uc_grid);
    ref_.Lattice.print(cout);
  }
  
  return true;
}

void XMLParticleParser::getPtclAttrib(xmlNodePtr cur, int nat, int nloc) {

  xmlDocPtr doc = cur->doc;
  string oname, otype;
  int utype = 0;
  xmlAttrPtr att = cur->properties;
  while(att != NULL) {
    string aname((const char*)(att->name));
    const char* vname = (const char*)(att->children->content);
    if(aname == "name") { oname = vname;} 
    else if(aname == datatype_tag) {otype = vname;}
    else if(aname == condition_tag) {utype = atoi(vname);}
    att = att->next;
  }
  if(oname.empty() || otype.empty()) {
    ERRORMSG("No default value is set for name and type attributes")
    return;
  }
      
  istringstream  stream((const char*)
			(xmlNodeListGetString(doc, cur->xmlChildrenNode, 1)));
  int t_id = ref_.getAttribType(otype);
  if(oname == ionid_tag) { 
    if(otype == stringtype_tag) {
      int nloci = nloc;
      string a;   
      for(int iat=0; iat<nat; iat++,nloci++) {
	stream >> a;
	ref_.GroupID[nloci] = ref_.Species.addSpecies(a);
      }
    } else {
      int nloci = nloc;
      for(int iat=0; iat<nat; iat++,nloci++) stream >> ref_.GroupID[nloci];
    }
  } else {
    if(t_id == PA_IndexType) {
      ParticleIndex_t& itmp = *(ref_.getIndexAttrib(oname));
      itmp.setUnit(utype);
      int nloci = nloc;
      for(int iat=0; iat<nat; iat++,nloci++) stream >> itmp[nloci];
    } else if(t_id == PA_ScalarType) {
      ParticleScalar_t& stmp = *(ref_.getScalarAttrib(oname));
      stmp.setUnit(utype);
      int nloci = nloc;
      for(int iat=0; iat<nat; iat++,nloci++) stream >> stmp[nloci];
    } else if(t_id == PA_PositionType) {
      ParticlePos_t& ptmp = *(ref_.getVectorAttrib(oname));
      ptmp.setUnit(utype);
      int nloci = nloc;
      for(int iat=0; iat<nat; iat++,nloci++) stream >> ptmp[nloci];
    } else if(t_id == PA_TensorType) {
      ParticleTensor_t& itmp = *(ref_.getTensorAttrib(oname));
      itmp.setUnit(utype);
      int nloci = nloc;
      for(int iat=0; iat<nat; iat++,nloci++) stream >> itmp[nloci];
    }
  }
}


XMLSaveParticle::XMLSaveParticle(Particle_t& pin): ref_(pin){ }

XMLSaveParticle::~XMLSaveParticle() { 
}

void XMLSaveParticle::reset(const char* fileroot){

  FileRoot = fileroot;
  FileName = fileroot;
  FileName.append(".xml");
  SpeciesName = ref_.Species.speciesName;
}

void XMLSaveParticle::report(int iter) {

  // writing a meta file
  ofstream fxml(FileName.c_str()); // always overwrite
  fxml << "<?xml version=\"1.0\"?>" << endl;
  get(fxml,1);
  fxml.close();
}

void XMLSaveParticle::get(ostream& fxml, int olevel) const  {

  ref_.begin_node(fxml);

  fxml.setf(ios::scientific);
  fxml.precision(15);
  LatticeXMLWriter latticeout(ref_.Lattice);
  latticeout.get(fxml);

  for(int i=0; i<SpeciesName.size(); i++) {
    fxml << "<group name=\"" << SpeciesName[i] << "\"/>" << endl;
  }

  //only write the local particles
  int nloc = ref_.getLocalNum();

  if(olevel) {
    Particle_t::PAListIterator it = ref_.first_attrib();
    while(it != ref_.last_attrib()) {
      OhmmsObject* ooref= (*it).second;
//       if(ooref->objName() == ionid_tag) {
// 	IonName.begin_node(fxml); 
//  	for(int iat=0; iat<nloc; iat++) {
//  	  fxml << IonName[iat] << " ";
//  	  if(iat%20 == 19) fxml << endl;
//  	}
// 	IonName.end_node(fxml); 
//       } else {
	int t_id = ref_.getAttribType(ooref->typeName());
	int o_id = ooref->id(); 
        ooref->begin_node(fxml);
	if(t_id == PA_IndexType) {
	  const ParticleIndex_t& itmp = *(ref_.getIndexAttrib(o_id));
	  for(int iat=0; iat<nloc; iat++) {
	    fxml << itmp[iat] << " ";
	    if(iat%20 == 19) fxml << endl;
	  }
	} else if(t_id == PA_ScalarType) {
	  fxml.precision(6);
	  const ParticleScalar_t& stmp =*(ref_.getScalarAttrib(o_id));
	  for(int iat=0; iat<nloc; iat++) {
	    fxml << stmp[iat] << " ";
	    if(iat%5 == 4) fxml << endl;
	  }
	  if(nloc%5 != 0) fxml<< endl;
	} else if (t_id == PA_PositionType) {
	  fxml.precision(15);
	  const ParticlePos_t& rtmp =*(ref_.getVectorAttrib(o_id));

	  for(int iat=0; iat<nloc; iat++) {
	    fxml << rtmp[iat] << endl;
	  }
	} else if (t_id == PA_TensorType) {
	  fxml.precision(15);
	  const ParticleTensor_t& ttmp =*(ref_.getTensorAttrib(o_id));
	  for(int iat=0; iat<nloc; iat++) {
	    fxml << ttmp[iat];
	  }
	}
        ooref->end_node(fxml);
	//      }
      it++;
    }
  } else {

    ref_.R.begin_node(fxml);
    for(int iat=0; iat<nloc; iat++) fxml << ref_.R[iat] << endl;
    ref_.R.end_node(fxml);

    ref_.GroupID.begin_node(fxml);
    for(int iat=0; iat<nloc; iat++) {
      fxml << ref_.GroupID[iat] << " "; if(iat%20 == 19) fxml << endl;
    }
    if(nloc%20 != 19) fxml << endl;
    ref_.GroupID.end_node(fxml);
  }
  ref_.end_node(fxml);
}

bool XMLSaveParticle::put(xmlNodePtr cur) {
  
  return true;
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

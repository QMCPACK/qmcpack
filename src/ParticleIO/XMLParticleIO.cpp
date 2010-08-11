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
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/OhmmsInfo.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/ParticleIOUtility.h"
//#include "ParticleIO/HDFParticleIO.h"
#include "ParticleBase/ParticleFunctions.h"
#include "ParticleBase/RandomSeqGenerator.h"
using namespace qmcplusplus;


XMLParticleParser::XMLParticleParser(Particle_t& aptcl, Tensor<int,OHMMS_DIM>& tmat
    , bool donotresize):
  AssignmentOnly(donotresize),ref_(aptcl),TileMatrix(tmat)
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
    app_error() << fname_in << " does not contain any ParticleSet" << endl;
  } else {
    xmlNodePtr cur=result->nodesetval->nodeTab[0];
    string fname, pformat;
    OhmmsAttributeSet pAttrib;
    pAttrib.add(fname,"src"); pAttrib.add(fname,"href");
    pAttrib.add(pformat,"srctype");
    pAttrib.put(cur);

    if(fname.size()) pformat=getExtension(fname);

    if(pformat.empty())
      putSpecial(cur);
    //else
    //{
    //  if(pformat == "h5") {
    //    HDFParticleParser ahandle(ref_);
    //    ahandle.put(cur);
    //  } else {
    //    app_error() << "  Unknown file extension " << pformat << endl;
    //  }
    //} 
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
  string fname, pformat("xml");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(fname,"src");
  pAttrib.add(fname,"href");
  pAttrib.add(pformat,"srctype");
  pAttrib.put(cur);

  if(fname.empty()) 
    return putSpecial(cur);
  else
  {//overwrite the format
    pformat = getExtension(fname);
    return put(fname,pformat);
  }
}


/** process xmlnode &lt;particleset/&gt; which contains everything about the particle set to initialize
 *@param cur the xmlnode to work on
 *
 */
bool XMLParticleParser::putSpecial(xmlNodePtr cur) {

  string pname("none");
  //xmlDocPtr doc = cur->doc;

  //the number of particles that are initialized by <attrib/>
  int nat = 0;
  string randomizeR("no");
  string randomsrc("");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(randomizeR,"random");
  pAttrib.add(randomsrc, "randomsrc");
  pAttrib.add(nat,"size");
  pAttrib.add(pname,"name");
  pAttrib.put(cur);

  ///count the number of atom added one by one
  xmlNodePtr cur0 = cur->xmlChildrenNode;

  //total count of the particles to be created
  int ntot = 0;
  int ng = 0, ng_in=0;
  vector<int> nat_group;

  vector<xmlNodePtr> atom_ptr;

  //pre-process the nodes to count the number of particles to be added
  while(cur0 != NULL) {
    string cname((const char*)cur0->name);
    if(cname == "atom") {
      ntot++;	
      atom_ptr.push_back(cur0);
    } else if(cname == "group") {
      int nat_per_group=0;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(nat_per_group,"size");
      gAttrib.put(cur0);
      nat_group.push_back(nat_per_group);
      ng_in += nat_per_group;
      ntot += nat_per_group;
      ng++;
    } else if(cname == attrib_tag) {
      int size_att = 0;
      OhmmsAttributeSet aAttrib;
      aAttrib.add(size_att,"size");
      aAttrib.put(cur0);
      if(size_att) {
	if(size_att != nat) {
	  app_warning() << "\tOverwriting the size of the particle by //particleset/attrib/@size=" 
            << size_att << endl;
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
    if(ng_in) {
      ref_.create(nat_group);
    } else {
      ref_.create(ntot);
    }
    //assign default ID
    int nloci=nloc;
    for(int iat=0;iat<ntot; iat++,nloci++) ref_.ID[iat]=nloci;
  }
  
  //TinyVector<int,OHMMS_DIM> uc_grid(1);      
  
  SpeciesSet& tspecies(ref_.getSpeciesSet()); //SpeciesCollection::getSpecies();

  cur = cur->xmlChildrenNode;
  //reset the group counter
  ng = 0;  
  while (cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname.find("ell")< cname.size())//accept UnitCell, unitcell, supercell
    { //if(cname == "UnitCell" || cname == "unitcell") {
      LatticeParser lat(ref_.Lattice);
      lat.put(cur);
      //ParameterSet params;
      //params.add(uc_grid,"uc_grid","int");
      //params.put(cur);
    } else if (cname == attrib_tag) {
      getPtclAttrib(cur,nat,nloc);
    } else  if (cname == "group") { //found group
      string sname;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(sname,"name");
      gAttrib.put(cur);
      if(sname.size()) //only if name is found
      {
        int sid=tspecies.addSpecies(sname);
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
          string tcname((const char*)tcur->name);
          if(tcname == "parameter") {
            string pname;
            OhmmsAttributeSet pAttrib;
            pAttrib.add(pname,"name");
            pAttrib.put(tcur);
            if(pname.size())
            {
              int iproperty=tspecies.addAttribute(pname);
              double ap;
              putContent(ap,tcur);
              tspecies(iproperty,sid) = ap;
            }
          } else if(nat_group[ng] && tcname == attrib_tag) { 
            getPtclAttrib(tcur,nat_group[ng],nloc);
          }
          tcur = tcur->next;
        }
        for(int iat=0; iat<nat_group[ng]; iat++, nloc++)  ref_.GroupID[nloc] = sid;
        ng++;
      }
    }
    cur = cur->next;
  }

    expandSuperCell(ref_,TileMatrix);

  //Disable atom 
  ////have read from <attrib/>'s and <group/>'s. Time to add <atom/>'s
  //nloc += nat;
  //for(int ia=0; ia<atom_ptr.size(); ia++,nloc++) {
  //  cur = atom_ptr[ia];
  //  int inunit=0;
  //  string sname;
  //  OhmmsAttributeSet aAttrib;
  //  aAttrib.add(inunit,condition_tag);
  //  aAttrib.add(sname,"name");
  //  aAttrib.put(cur);
  //  if(sname.empty())
  //  {
  //    app_error() << "  Missing atom/@name. Fatal error." << endl;
  //    return false;
  //  }
  //  int sid= tspecies.addSpecies(sname);
  //  ParticleSet::SingleParticlePos_t pos;
  //  putContent(pos,cur);
  //  if(inunit == ref_.R.getUnit())
  //    ref_.R[nloc]=pos;
  //  else
  //  {
  //    if(inunit) 
  //      ref_.R[nloc]=ref_.Lattice.toCart(pos);
  //    else
  //      ref_.R[nloc]=ref_.Lattice.toUnit(pos);
  //  }
  //  ref_.ID[nloc] = nloc; 
  //  ref_.GroupID[nloc] = sid;
  //}

  //disable uc_grid
  //int ngtot=uc_grid[0];
  //for(int idim=1; idim<OHMMS_DIM; idim++) ngtot*=uc_grid[idim];
  //if(ngtot>1) {
  //  ExpandSuperCell(ref_,uc_grid);
  //  //ref_.Lattice.print(cout);
  //}
  
  ref_.RandomSource=randomsrc;

  if(randomizeR == "yes") {
    if(ref_.Lattice.SuperCellEnum)
    { 
      makeUniformRandom(ref_.R);
      ref_.R.setUnit(PosUnit::LatticeUnit);
      ref_.convert2Cart(ref_.R);
    }  
    else if (randomsrc == "") {
      ostringstream o;
      o << "Not know how to randomize R of an open system.\n"
        << "Use randomsrc=\"ion\" instead.\n";
      APP_ABORT(o.str());
    }
  }

  vector<int> numPerGroup(tspecies.getTotalNum(),0);
  for(int iat=0; iat<ref_.GroupID.size(); iat++) 
    numPerGroup[ref_.GroupID[iat]]++;


  int membersize= tspecies.addAttribute("membersize");
  for(int ig=0; ig<tspecies.getTotalNum(); ++ig) {
    tspecies(membersize,ig)=numPerGroup[ig];
  }

  int beforemass=tspecies.numAttributes();
  int massind= tspecies.addAttribute("mass");
  if(beforemass == massind)
  {
    app_log() << "  XMLParticleParser setting mass of  " << ref_.getName() << " to 1.0" << endl;
    for(int ig=0; ig<tspecies.getTotalNum(); ++ig) 
      tspecies(massind,ig)=1.0; 
  }
  //Check the unit of ParticleSet::R and PBC
  ref_.createSK();

  return true;
}

template<typename PAT>
struct ParticleAttribXmlNode
{
  PAT& ref_;

  inline ParticleAttribXmlNode(PAT& a, int utype):ref_(a)
  {
    ref_.InUnit=utype;
  }

  inline bool put(xmlNodePtr cur, int n_in, int start)
  {
    typedef typename PAT::Type_t data_type;
    std::vector<data_type> data_in(n_in);
    putContent(data_in,cur);
    std::copy(data_in.begin(),data_in.end(),ref_.begin()+start);
    return true;
  }
};

void XMLParticleParser::getPtclAttrib(xmlNodePtr cur, int nat, int nloc) {

  string oname, otype;
  int utype=0;
  int size_in=0;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(otype,datatype_tag);//datatype
  pAttrib.add(oname,"name");//name
  pAttrib.add(utype,condition_tag);//condition
  pAttrib.add(size_in,"size");//size
  pAttrib.put(cur);
  
  if(oname.empty() || otype.empty()) {
    app_error() << "   Missing attrib/@name or attrib/@datatype " << endl;
    app_error() << "     <attrib name=\"aname\"  datatype=\"atype\"/>" << endl;
    return;
  }

  int t_id = ref_.getAttribType(otype);
  if(oname == ionid_tag) { 
    if(otype == stringtype_tag) {
      int nloci = nloc;
      vector<string> d_in(nat);
      putContent(d_in,cur);
      for(int iat=0; iat<d_in.size(); iat++,nloci++) {
        ref_.GroupID[nloci] = ref_.getSpeciesSet().addSpecies(d_in[iat]);
      }
    } else {
      ParticleAttribXmlNode<ParticleIndex_t> a(ref_.GroupID,utype);
      a.put(cur,nat,nloc);
    }
  } else {
    if(t_id == PA_IndexType) {
      ParticleAttribXmlNode<ParticleIndex_t> a(*(ref_.getIndexAttrib(oname)),utype);
      a.put(cur,nat,nloc);
    } else if(t_id == PA_ScalarType) {
      ParticleAttribXmlNode<ParticleScalar_t> a(*(ref_.getScalarAttrib(oname)),utype);
      a.put(cur,nat,nloc);
    } else if(t_id == PA_PositionType) {
      ParticleAttribXmlNode<ParticlePos_t> a(*(ref_.getVectorAttrib(oname)),utype);
      a.put(cur,nat,nloc);
    } else if(t_id == PA_TensorType) {
      ParticleAttribXmlNode<ParticleTensor_t> a(*(ref_.getTensorAttrib(oname)),utype);
      a.put(cur,nat,nloc);
    }
  }
}


XMLSaveParticle::XMLSaveParticle(Particle_t& pin): ref_(pin){ }

XMLSaveParticle::~XMLSaveParticle() { 
}

void XMLSaveParticle::reset(const char* fileroot, bool append){
  //append is ignored

  FileRoot = fileroot;
  FileName = fileroot;
  FileName.append(".xml");
  SpeciesName = ref_.getSpeciesSet().speciesName;
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

/** create particleset node
 * @param addlattice if true, add unitcell
 */
xmlNodePtr XMLSaveParticle::createNode(bool addlattice) {

  if(SpeciesName.size()!=ref_.getSpeciesSet().getTotalNum()) {
    SpeciesName = ref_.getSpeciesSet().speciesName;
  }

  //if(addlattice) {
  //  ref_.Lattice.print(cout);
  //}

  xmlNodePtr cur = xmlNewNode(NULL,(const xmlChar*)"particleset");
  xmlNewProp(cur,(const xmlChar*)"name",(const xmlChar*)ref_.getName().c_str());

  if(ref_.groups()>1) {
    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for(int is=0; is<nspecies; is++) {
      std::ostringstream ng;
      ng << ref_.last(is)-ref_.first(is);
      xmlNodePtr g=xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(g,(const xmlChar*)"name",(const xmlChar*)SpeciesName[is].c_str());
      xmlNewProp(g,(const xmlChar*)"size",(const xmlChar*)ng.str().c_str());

      for(int item=0; item<nitem; item++) {
        std::ostringstream prop;
        prop<<mySpecies(item,is);
        xmlNodePtr p=xmlNewTextChild(g,NULL,
            (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)mySpecies.attribName[item].c_str());
      }
      std::ostringstream pos;
      pos.setf(ios_base::scientific);
      pos <<"\n";
      for(int iat=ref_.first(is); iat<ref_.last(is); iat++) {
        pos <<  ref_.R[iat] << endl;
      }
      xmlNodePtr posPtr=xmlNewTextChild(g,NULL,
              (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
      xmlNewProp(posPtr,(const xmlChar*)"name",(const xmlChar*)"position");
      xmlNewProp(posPtr,(const xmlChar*)"datatype",(const xmlChar*)"posArray");

      xmlAddChild(cur,g);
    }
  } else {
    std::ostringstream nat;
    nat<<ref_.getTotalNum();
    xmlNewProp(cur,(const xmlChar*)"size",(const xmlChar*)nat.str().c_str());

    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for(int is=0; is<nspecies; is++) {
      xmlNodePtr g=xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(g,(const xmlChar*)"name",(const xmlChar*)SpeciesName[is].c_str());
      for(int item=0; item<nitem; item++) {
        std::ostringstream prop;
        prop<<mySpecies(item,is);
        xmlNodePtr p=xmlNewTextChild(g,NULL,
            (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)mySpecies.attribName[item].c_str());
      }
      xmlAddChild(cur,g);
    }

    std::ostringstream pos,gid;
    pos.setf(ios_base::scientific);
    pos <<"\n";
    for(int iat=0; iat<ref_.getTotalNum(); iat++) {
      pos << ref_.R[iat] << endl;
    }
    xmlNodePtr posPtr=xmlNewTextChild(cur,NULL,
            (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
    xmlNewProp(posPtr,(const xmlChar*)"name",(const xmlChar*)"position");
    xmlNewProp(posPtr,(const xmlChar*)"datatype",(const xmlChar*)"posArray");

    gid <<"\n ";
    for(int iat=0; iat<ref_.getTotalNum(); iat++) {
      gid << SpeciesName[ref_.GroupID[iat]] << " ";
    }
    gid << endl;
    xmlNodePtr gPtr=xmlNewTextChild(cur,NULL,
            (const xmlChar*)"attrib", (const xmlChar*)gid.str().c_str());
    xmlNewProp(gPtr,(const xmlChar*)"name",(const xmlChar*)"ionid");
    xmlNewProp(gPtr,(const xmlChar*)"datatype",(const xmlChar*)"stringArray");
  }

  return cur;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

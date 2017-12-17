//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
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
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/ParticleIOUtility.h"
//#include "ParticleIO/HDFParticleIO.h"
#include "ParticleBase/ParticleFunctions.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Utilities/ProgressReportEngine.h"
namespace qmcplusplus
{

/** set the property of a SpeciesSet
 * @param tspecies SpeciesSet
 * @param sid id of the species whose properties to be set
 *
 * Need unit handlers but for now, everything is in AU:
 * m_e=1.0, hartree and bohr and unit="au"
 * Example to define C(arbon)
 * <group name="C">
 *   <parameter name="mass" unit="amu">12</parameter>
 *   <parameter name="charge">-6</parameter>
 * </group>
 * Note that unit="amu" is given for the mass.
 * When mass is not given, they are set to the electron mass.
 */
void setSpeciesProperty(SpeciesSet& tspecies, int sid, xmlNodePtr cur)
{
  const double proton_mass=1822.888530063;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "parameter")
    {
      std::string pname;
      std::string unit_name("au");//hartree, bohr & me=1
      OhmmsAttributeSet pAttrib;
      pAttrib.add(pname,"name");
      pAttrib.add(unit_name,"unit");
      pAttrib.put(cur);
      if(pname.size())
      {
        int iproperty=tspecies.addAttribute(pname);
        double ap=0.0;
        double unit_conv=1.0;
        putContent(ap,cur);
        if(pname == "mass")
        {
          if(unit_name=="amu")
            unit_conv=proton_mass;
        }
        tspecies(iproperty,sid) = ap*unit_conv;
      }
    }
    cur=cur->next;
  }
}


XMLParticleParser::XMLParticleParser(Particle_t& aptcl, Tensor<int,OHMMS_DIM>& tmat
                                     , bool donotresize):
  AssignmentOnly(donotresize),ref_(aptcl),TileMatrix(tmat)
{
  //add ref particle attributes
  ref_.createAttributeList(ref_AttribList);
}

/**reading particleset node from a file
 *@param fname_in a file name to open
 *@param pformat_in the format of the file, not used
 *@return true, if successful
 *
 *Check the type of the external source to work on.
 *The external source itself can be an xml file.
 */
bool XMLParticleParser::put(const std::string& fname_in,
                            const std::string& fext_in)
{
  xmlDocPtr doc=NULL;
  xmlNsPtr ns;
  // build an XML tree from a the file;
  doc = xmlParseFile(fname_in.c_str());
  if (doc == NULL)
  {
    ERRORMSG(fname_in << " does not exist")
    return false;
  }
  ///using XPath instead of recursive search
  xmlXPathContextPtr context;
  xmlXPathObjectPtr result;
  context = xmlXPathNewContext(doc);
  result = xmlXPathEvalExpression((const xmlChar*)"//particleset",context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    app_error() << fname_in << " does not contain any ParticleSet" << std::endl;
  }
  else
  {
    xmlNodePtr cur=result->nodesetval->nodeTab[0];
    std::string fname, pformat;
    OhmmsAttributeSet pAttrib;
    pAttrib.add(fname,"src");
    pAttrib.add(fname,"href");
    pAttrib.add(pformat,"srctype");
    pAttrib.put(cur);
    if(fname.size())
      pformat=getExtension(fname);
    if(pformat.empty())
      putSpecial(cur);
    //else
    //{
    //  if(pformat == "h5") {
    //    HDFParticleParser ahandle(ref_);
    //    ahandle.put(cur);
    //  } else {
    //    app_error() << "  Unknown file extension " << pformat << std::endl;
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
bool XMLParticleParser::put(xmlNodePtr cur)
{
  ///process attributes: type or format
  std::string fname, pformat("xml");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(fname,"src");
  pAttrib.add(fname,"href");
  pAttrib.add(pformat,"srctype");
  pAttrib.put(cur);
  if(fname.empty())
    return putSpecial(cur);
  else
  {
    //overwrite the format
    pformat = getExtension(fname);
    return put(fname,pformat);
  }
}


/** process xmlnode &lt;particleset/&gt; which contains everything about the particle set to initialize
 *@param cur the xmlnode to work on
 *
 */
bool XMLParticleParser::putSpecial(xmlNodePtr cur)
{
  ReportEngine PRE("XMLParticleParser","putSpecial");
  std::string pname("none");
  //xmlDocPtr doc = cur->doc;
  //the number of particles that are initialized by <attrib/>
  int nat = 0;
  std::string randomizeR("no");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(randomizeR,"random");
  pAttrib.add(nat,"size");
  pAttrib.add(pname,"name");
  pAttrib.put(cur);
  ///count the number of atom added one by one
  xmlNodePtr cur0 = cur->xmlChildrenNode;
  //total count of the particles to be created
  int ntot = 0;
  int ng = 0, ng_in=0;
  std::vector<int> nat_group;
  std::vector<xmlNodePtr> atom_ptr;
  //pre-process the nodes to count the number of particles to be added
  while(cur0 != NULL)
  {
    std::string cname((const char*)cur0->name);
    if(cname == "atom")
    {
      ntot++;
      atom_ptr.push_back(cur0);
    }
    else if(cname == "group")
    {
      int nat_per_group=0;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(nat_per_group,"size");
      gAttrib.put(cur0);
      nat_group.push_back(nat_per_group);
      ng_in += nat_per_group;
      ntot += nat_per_group;
      ng++;
    }
    else if(cname == attrib_tag)
    {
      int size_att = 0;
      OhmmsAttributeSet aAttrib;
      aAttrib.add(size_att,"size");
      aAttrib.put(cur0);
      if(size_att)
      {
        if(size_att != nat)
        {
          app_warning() << "\tOverwriting the size of the particle by //particleset/attrib/@size="
            << size_att << std::endl;
          nat = size_att;
        }
      }
    }
    cur0 = cur0->next;
  }
  ntot += nat;
  ref_.setName(pname.c_str());
  int nloc = ref_.getTotalNum();
  //treat assignment only differently
  if(AssignmentOnly)
  {
    ntot = 0;
    nloc = 0;
    for(int iat=0; iat<ref_.getTotalNum(); iat++)
      ref_.ID[iat]=iat;
  }
  if(ntot)
  {
    if(ng_in)
    {
      ref_.create(nat_group);
    }
    else
    {
      ref_.create(ntot);
    }
    //assign default ID
    int nloci=nloc;
    for(int iat=0; iat<ntot; iat++,nloci++)
      ref_.ID[iat]=nloci;
  }
  //TinyVector<int,OHMMS_DIM> uc_grid(1);
  SpeciesSet& tspecies(ref_.getSpeciesSet()); //SpeciesCollection::getSpecies();
  cur = cur->xmlChildrenNode;
  //reset the group counter
  ng = 0;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname.find("ell")< cname.size())//accept UnitCell, unitcell, supercell
    {
      //if(cname == "UnitCell" || cname == "unitcell") {
      LatticeParser lat(ref_.Lattice);
      lat.put(cur);
      //ParameterSet params;
      //params.add(uc_grid,"uc_grid","int");
      //params.put(cur);
    }
    else if (cname == attrib_tag)
    {
      getPtclAttrib(cur,nat,nloc);
    }
    else if (cname == "group")
      //found group
    {
      std::string sname;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(sname,"name");
      gAttrib.put(cur);
      if(sname.size()) //only if name is found
      {
        int sid=tspecies.addSpecies(sname);
        setSpeciesProperty(tspecies,sid,cur);
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL)
        {
          std::string tcname((const char*)tcur->name);
          if(nat_group[ng] && tcname == attrib_tag)
          {
            getPtclAttrib(tcur,nat_group[ng],nloc);
          }
          tcur = tcur->next;
        }
        for(int iat=0; iat<nat_group[ng]; iat++, nloc++)
          ref_.GroupID[nloc] = sid;
        ng++;
      }
    }
    cur = cur->next;
  }

    //copy ID -> PCID
    ref_.PCID=ref_.ID;

  expandSuperCell(ref_,TileMatrix);
  if(ref_.Lattice.SuperCellEnum)
  {
    if(randomizeR == "yes")
    {
      makeUniformRandom(ref_.R);
      ref_.R.setUnit(PosUnit::LatticeUnit);
      ref_.convert2Cart(ref_.R);
    }
    else  // put them [0,1) in the cell 
      ref_.applyBC(ref_.R);
  } 

  //this sets Mass, Z
  ref_.resetGroups();
  ref_.createSK();

  return true;
}

/** process xmlnode to reset the properties of a particle set
 * @param cur current node
 * @return true, if successful
 *
 * This resets or adds new attributes to a particle set.
 * It cannot modify the size of the particle set.
 */
bool XMLParticleParser::reset(xmlNodePtr cur)
{
  ReportEngine PRE("XMLParticleParser","reset");
  SpeciesSet& tspecies(ref_.getSpeciesSet());
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "group")
    {
      std::string sname;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(sname,"name");
      gAttrib.put(cur);
      if(sname.size())
      {
        int sid=tspecies.addSpecies(sname);
        setSpeciesProperty(tspecies,sid,cur);
      }
    }
    cur = cur->next;
  }
//  //@todo Will add a member function to ParticleSet to handle these
//  int massind=tspecies.addAttribute("mass");
//  for(int iat=0; iat<ref_.getTotalNum(); iat++)
//    ref_.Mass[iat]=tspecies(massind,ref_.GroupID[iat]);
//
//  int qind=tspecies.addAttribute("charge");
//  for(int iat=0; iat<ref_.getTotalNum(); iat++)
//    ref_.Z[iat]=tspecies(qind,ref_.GroupID[iat]);
//
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
    copy(data_in.begin(),data_in.end(),ref_.begin()+start);
    return true;
  }
};

void XMLParticleParser::getPtclAttrib(xmlNodePtr cur, int nat, int nloc)
{
  std::string oname, otype;
  int utype=0;
  int size_in=0;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(otype,datatype_tag);//datatype
  pAttrib.add(oname,"name");//name
  pAttrib.add(utype,condition_tag);//condition
  pAttrib.add(size_in,"size");//size
  pAttrib.put(cur);
  if(oname.empty() || otype.empty())
  {
    app_error() << "   Missing attrib/@name or attrib/@datatype " << std::endl;
    app_error() << "     <attrib name=\"aname\"  datatype=\"atype\"/>" << std::endl;
    return;
  }
  int t_id = ref_AttribList.getAttribType(otype);

  if(oname == ionid_tag)
  {
    if(otype == stringtype_tag)
    {
      int nloci = nloc;
      std::vector<std::string> d_in(nat);
      putContent(d_in,cur);
      for(int iat=0; iat<d_in.size(); iat++,nloci++)
      {
        ref_.GroupID[nloci] = ref_.getSpeciesSet().addSpecies(d_in[iat]);
      }
    }
    else
    {
      ParticleAttribXmlNode<ParticleIndex_t> a(ref_.GroupID,utype);
      a.put(cur,nat,nloc);
    }
  }
  else
  {
    //very permissive in that a unregistered attribute will be created and stored by ParticleSet
    //cloning is not going to work
    if(t_id == PA_IndexType)
    {
      ParticleIndex_t* obj=nullptr;
      obj=ref_AttribList.getAttribute(otype,oname,obj);
      ParticleAttribXmlNode<ParticleIndex_t> a(*obj,utype);
      a.put(cur,nat,nloc);
    }
    else if(t_id == PA_ScalarType)
    {
      ParticleScalar_t* obj=nullptr;
      obj=ref_AttribList.getAttribute(otype,oname,obj);
      ParticleAttribXmlNode<ParticleScalar_t> a(*obj,utype);
      a.put(cur,nat,nloc);
    }
    else if(t_id == PA_PositionType)
    {
      ParticlePos_t* obj=nullptr;
      obj=ref_AttribList.getAttribute(otype,oname,obj);
      ParticleAttribXmlNode<ParticlePos_t> a(*obj,utype);
      a.put(cur,nat,nloc);
    }
    else if(t_id == PA_TensorType)
    {
      ParticleTensor_t* obj=nullptr;
      obj=ref_AttribList.getAttribute(otype,oname,obj);
      ParticleAttribXmlNode<ParticleTensor_t> a(*obj,utype);
      a.put(cur,nat,nloc);
    }
  }
}


XMLSaveParticle::XMLSaveParticle(Particle_t& pin): ref_(pin) { }

XMLSaveParticle::~XMLSaveParticle()
{
}

void XMLSaveParticle::reset(const char* fileroot, bool append)
{
  //append is ignored
  FileRoot = fileroot;
  FileName = fileroot;
  FileName.append(".xml");
  SpeciesName = ref_.getSpeciesSet().speciesName;
}

void XMLSaveParticle::report(int iter)
{
  // writing a meta file
  std::ofstream fxml(FileName.c_str()); // always overwrite
  fxml << "<?xml version=\"1.0\"?>" << std::endl;
  get(fxml,1);
  fxml.close();
}

void XMLSaveParticle::get(std::ostream& fxml, int olevel) const
{
  ref_.begin_node(fxml);
  fxml.setf(std::ios::scientific);
  fxml.precision(15);
  LatticeXMLWriter latticeout(ref_.Lattice);
  latticeout.get(fxml);
  for(int i=0; i<SpeciesName.size(); i++)
  {
    fxml << "<group name=\"" << SpeciesName[i] << "\"/>" << std::endl;
  }
  //only write the local particles
  int nloc = ref_.getTotalNum();
  if(olevel)
  {
    /*
    Particle_t::PAListIterator it = ref_.first_attrib();
    while(it != ref_.last_attrib())
    {
      OhmmsObject* ooref= (*it).second;
//       if(ooref->objName() == ionid_tag) {
// 	IonName.begin_node(fxml);
//  	for(int iat=0; iat<nloc; iat++) {
//  	  fxml << IonName[iat] << " ";
//  	  if(iat%20 == 19) fxml << std::endl;
//  	}
// 	IonName.end_node(fxml);
//       } else {
      int t_id = ref_AttribList.getAttribType(otype);
      int o_id = ooref->id();
      ooref->begin_node(fxml);
      if(t_id == PA_IndexType)
      {
        const ParticleIndex_t* itmp=dynamic_cast<ParticleIndex_t*>(ooref);
        for(int iat=0; iat<nloc; iat++)
        {
          fxml << (*itmp)[iat] << " ";
          if(iat%20 == 19)
            fxml << std::endl;
        }
      }
      else if(t_id == PA_ScalarType)
      {
        fxml.precision(6);
        const ParticleScalar_t* stmp=dynamic_cast<ParticleScalar_t*>(ooref);
        for(int iat=0; iat<nloc; iat++)
        {
          fxml << (*stmp)[iat] << " ";
          if(iat%5 == 4)
            fxml << std::endl;
        }
        if(nloc%5 != 0)
          fxml<< std::endl;
      }
      else if (t_id == PA_PositionType)
      {
        fxml.precision(15);
        const ParticlePos_t* rtmp=dynamic_cast<ParticlePos_t*>(ooref);
        for(int iat=0; iat<nloc; iat++)
        {
          fxml << (*rtmp)[iat] << std::endl;
        }
      }
      else if (t_id == PA_TensorType)
      {
        fxml.precision(15);
        const ParticleTensor_t* ttmp=dynamic_cast<ParticleTensor_t*>(ooref);
        for(int iat=0; iat<nloc; iat++)
        {
          fxml << (*ttmp)[iat];
        }
      }
      ooref->end_node(fxml);
      //      }
      it++;
    }
    */
  }
  else
  {
    ref_.R.begin_node(fxml);
    for(int iat=0; iat<nloc; iat++)
      fxml << ref_.R[iat] << std::endl;
    ref_.R.end_node(fxml);
    ref_.GroupID.begin_node(fxml);
    for(int iat=0; iat<nloc; iat++)
    {
      fxml << ref_.GroupID[iat] << " ";
      if(iat%20 == 19)
        fxml << std::endl;
    }
    if(nloc%20 != 19)
      fxml << std::endl;
    ref_.GroupID.end_node(fxml);
  }
  ref_.end_node(fxml);
}

bool XMLSaveParticle::put(xmlNodePtr cur)
{
  return true;
}

/** create particleset node
 * @param addlattice if true, add unitcell
 */
xmlNodePtr XMLSaveParticle::createNode(bool addlattice)
{
  if(SpeciesName.size()!=ref_.getSpeciesSet().getTotalNum())
  {
    SpeciesName = ref_.getSpeciesSet().speciesName;
  }
  //if(addlattice) {
  //  ref_.Lattice.print(std::cout);
  //}
  xmlNodePtr cur = xmlNewNode(NULL,(const xmlChar*)"particleset");
  xmlNewProp(cur,(const xmlChar*)"name",(const xmlChar*)ref_.getName().c_str());
  if(ref_.groups()>1)
  {
    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for(int is=0; is<nspecies; is++)
    {
      std::ostringstream ng;
      ng << ref_.last(is)-ref_.first(is);
      xmlNodePtr g=xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(g,(const xmlChar*)"name",(const xmlChar*)SpeciesName[is].c_str());
      xmlNewProp(g,(const xmlChar*)"size",(const xmlChar*)ng.str().c_str());
      for(int item=0; item<nitem; item++)
      {
        std::ostringstream prop;
        prop<<mySpecies(item,is);
        xmlNodePtr p=xmlNewTextChild(g,NULL,
                                     (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)mySpecies.attribName[item].c_str());
      }
      std::ostringstream pos;
      pos.setf(std::ios_base::scientific);
      pos <<"\n";
      for(int iat=ref_.first(is); iat<ref_.last(is); iat++)
      {
        pos <<  ref_.R[iat] << std::endl;
      }
      xmlNodePtr posPtr=xmlNewTextChild(g,NULL,
                                        (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
      xmlNewProp(posPtr,(const xmlChar*)"name",(const xmlChar*)"position");
      xmlNewProp(posPtr,(const xmlChar*)"datatype",(const xmlChar*)"posArray");
      xmlAddChild(cur,g);
    }
  }
  else
  {
    std::ostringstream nat;
    nat<<ref_.getTotalNum();
    xmlNewProp(cur,(const xmlChar*)"size",(const xmlChar*)nat.str().c_str());
    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for(int is=0; is<nspecies; is++)
    {
      xmlNodePtr g=xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(g,(const xmlChar*)"name",(const xmlChar*)SpeciesName[is].c_str());
      for(int item=0; item<nitem; item++)
      {
        std::ostringstream prop;
        prop<<mySpecies(item,is);
        xmlNodePtr p=xmlNewTextChild(g,NULL,
                                     (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p,(const xmlChar*)"name",(const xmlChar*)mySpecies.attribName[item].c_str());
      }
      xmlAddChild(cur,g);
    }
    std::ostringstream pos,gid;
    pos.setf(std::ios_base::scientific);
    pos <<"\n";
    for(int iat=0; iat<ref_.getTotalNum(); iat++)
    {
      pos << ref_.R[iat] << std::endl;
    }
    xmlNodePtr posPtr=xmlNewTextChild(cur,NULL,
                                      (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
    xmlNewProp(posPtr,(const xmlChar*)"name",(const xmlChar*)"position");
    xmlNewProp(posPtr,(const xmlChar*)"datatype",(const xmlChar*)"posArray");
    gid <<"\n ";
    for(int iat=0; iat<ref_.getTotalNum(); iat++)
    {
      gid << SpeciesName[ref_.GroupID[iat]] << " ";
    }
    gid << std::endl;
    xmlNodePtr gPtr=xmlNewTextChild(cur,NULL,
                                    (const xmlChar*)"attrib", (const xmlChar*)gid.str().c_str());
    xmlNewProp(gPtr,(const xmlChar*)"name",(const xmlChar*)"ionid");
    xmlNewProp(gPtr,(const xmlChar*)"datatype",(const xmlChar*)"stringArray");
  }
  return cur;
}
}


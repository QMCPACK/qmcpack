//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
#include "ParticleIO/LatticeIO.h"
#include "XMLParticleIO.h"
#include "ParticleBase/RandomSeqGeneratorGlobal.h"
#include "Utilities/ProgressReportEngine.h"
#include <Message/UniformCommunicateError.h>

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
  const double proton_mass = 1822.888530063;
  cur                      = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "parameter")
    {
      std::string pname;
      std::string unit_name("au"); //hartree, bohr & me=1
      OhmmsAttributeSet pAttrib;
      pAttrib.add(pname, "name");
      pAttrib.add(unit_name, "unit");
      pAttrib.put(cur);
      if (pname.size())
      {
        int iproperty    = tspecies.addAttribute(pname);
        double ap        = 0.0;
        double unit_conv = 1.0;
        putContent(ap, cur);
        if (pname == "mass")
        {
          if (unit_name == "amu")
            unit_conv = proton_mass;
        }
        tspecies(iproperty, sid) = ap * unit_conv;
      }
    }
    cur = cur->next;
  }
}


XMLParticleParser::XMLParticleParser(Particle_t& aptcl) : ref_(aptcl)
{
  //add ref particle attributes
  ref_.createAttributeList(ref_AttribList);
}

/** process xmlnode &lt;particleset/&gt; which contains everything about the particle set to initialize
 *@param cur the xmlnode to work on
 *
 */
bool XMLParticleParser::readXML(xmlNodePtr cur)
{
  ReportEngine PRE("XMLParticleParser", "readXML");

  if (ref_.getTotalNum())
    throw UniformCommunicateError("The ParticleSet object to load XML input was not empty. Report a bug!");

  SpeciesSet& tspecies(ref_.getSpeciesSet());
  if (tspecies.size() != 0)
    throw UniformCommunicateError("The SpeciesSet object to load XML input was not empty. Report a bug!");

  // the total number of particles, once it is set non-zero, always check against it.
  int nat = 0;
  // the number of particles by group, once it is constructed, always check against it.
  std::vector<int> nat_group;

  std::string pname("none");
  std::string randomizeR("no");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(randomizeR, "random");
  pAttrib.add(nat, "size");
  pAttrib.add(pname, "name");
  pAttrib.put(cur);

  ref_.setName(pname.c_str());

  if (nat != 0)
  {
    app_debug() << "Set the total size " << nat
                << " by the 'size' attribute found in 'particleset' XML element node named '" << pname << "'."
                << std::endl;
  }

  bool ionid_found = false;
  { // parse all the 'group's to obtain or verify the total number of particles
    //total count of the particles to be created
    int ntot               = 0;
    int num_non_zero_group = 0;
    bool group_found       = false;

    processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
      if (cname == "atom")
        throw UniformCommunicateError("XML element node atom is no more supported");
      else if (cname.find("ell") < cname.size()) //accept UnitCell, unitcell, supercell
        throw UniformCommunicateError("Constructing cell inside particleset is illegal!");
      else if (cname == "group")
      {
        group_found       = true;
        std::string sname = getXMLAttributeValue(element, "name");
        if (sname.empty())
          throw UniformCommunicateError("'group' element node must include a name attribute!");
        else
        {
          const int sid = tspecies.addSpecies(sname);
          setSpeciesProperty(tspecies, sid, element);
        }

        int nat_per_group = 0;
        OhmmsAttributeSet gAttrib;
        gAttrib.add(nat_per_group, "size");
        gAttrib.put(element);

        nat_group.push_back(nat_per_group);
        ntot += nat_per_group;
        if (nat_per_group > 0)
          num_non_zero_group++;
      }
      else if (cname == attrib_tag && getXMLAttributeValue(element, "name") == ionid_tag)
        ionid_found = true;
    });

    if (!group_found)
      throw UniformCommunicateError("No 'group' XML element node was found. Check XML input!");

    if (nat != 0 && ntot != 0 && nat != ntot)
    {
      std::ostringstream msg;
      msg << "The total number of particles deterimined previously was " << nat
          << "but the sum of the sizes from all the 'group' XML element nodes is " << ntot
          << ". Please check the 'particleset' XML element node!" << std::endl;
      throw UniformCommunicateError(msg.str());
    }

    if (nat == 0 && ntot != 0)
    {
      nat = ntot;
      app_debug() << "Set the total size " << nat << " by the sum of the 'size's on all the 'group' XML element nodes."
                  << std::endl;
    }

    if (ntot > 0 && num_non_zero_group != nat_group.size())
      throw UniformCommunicateError("Some 'group' XML element node doesn't contain a 'size' attribute!");
  }

  { // parse all the 'attrib's to obtain or verify the total number of particles
    processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
      if (cname == attrib_tag)
      {
        std::string sname = getXMLAttributeValue(element, "name");
        if (sname.empty())
          throw UniformCommunicateError("'" + ParticleTags::attrib_tag +
                                        "' XML element node must include a name attribute!");

        int size_att = 0;
        OhmmsAttributeSet aAttrib;
        aAttrib.add(size_att, "size");
        aAttrib.put(element);

        if (nat != 0 && size_att != 0 && nat != size_att)
        {
          std::ostringstream msg;
          msg << "The total number of particles deterimined previously was " << nat
              << " but the 'size' atttribute found on the '" << ParticleTags::attrib_tag
              << "' XML element nodes named '" << sname << "' is " << size_att
              << ". Please check the 'particleset' XML element node!" << std::endl;
          throw UniformCommunicateError(msg.str());
        }

        if (nat == 0 && size_att != 0)
        {
          nat = size_att;
          app_debug() << "Set the total size " << nat << " by the 'size' on the '" << ParticleTags::attrib_tag
                      << "' XML element node named '" << sname << "'." << std::endl;
        }
      }
    });
  }

  if (nat == 0)
    throw UniformCommunicateError("Failed in figuring out the total number of particles. Check XML input!");

  if (ionid_found)
  { // parse ionid and construct input order to stored order
    std::vector<int> map_storage_to_input(nat);
    processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
      if (cname == attrib_tag && getXMLAttributeValue(element, "name") == ionid_tag)
      {
        std::string datatype = getXMLAttributeValue(element, datatype_tag);
        if (datatype != stringtype_tag)
          throw UniformCommunicateError("'ionid' only supports datatype=\"" + stringtype_tag + "\"");
        std::vector<std::string> d_in(nat);
        putContent(d_in, element);
        int storage_index = 0;
        for (int ig = 0; ig < nat_group.size(); ig++)
        {
          const auto& group_species_name = tspecies.getSpeciesName(ig);
          int count_group_size           = 0;
          for (int iat = 0; iat < nat; iat++)
          {
            const int element_index = tspecies.findSpecies(d_in[iat]);
            if (element_index == tspecies.size())
              throw UniformCommunicateError("Element " + d_in[iat] +
                                            " doesn't match any species from 'group' XML element nodes.");
            if (element_index == ig)
            {
              count_group_size++;
              map_storage_to_input[storage_index++] = iat;
            }
          }

          if (count_group_size == 0)
            throw UniformCommunicateError("Element '" + group_species_name + "' not found in 'ionid'.");

          if (nat_group[ig] == 0)
            nat_group[ig] = count_group_size;
          else if (nat_group[ig] != count_group_size)
          {
            std::ostringstream msg;
            msg << "The number of particles of element '" << group_species_name << "' from 'group' XML elment node was "
                << nat_group[ig] << " but 'ionid' contains " << count_group_size << " entries." << std::endl;
            throw UniformCommunicateError(msg.str());
          }
        }
      }
    });

    checkGrouping(nat, nat_group);
    ref_.create(nat_group);
    // save map_storage_to_input
    ref_.setMapStorageToInput(map_storage_to_input);

    for (int iat = 0; iat < nat; iat++)
    {
      processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
        if (cname == attrib_tag && getXMLAttributeValue(element, "name") != ionid_tag)
          getPtclAttrib(element, map_storage_to_input[iat], 1, iat);
      });
    }
  }
  else
  {
    // fix old input with positions outside 'group'
    if (nat_group.size() == 1 && nat_group[0] == 0)
      nat_group[0] = nat;

    checkGrouping(nat, nat_group);
    ref_.create(nat_group);

    // obtain 'attrib' inside 'group'
    size_t start = 0;
    size_t ig    = 0;
    processChildren(cur, [&](const std::string& cname, const xmlNodePtr child) {
      if (cname == "group")
      {
        processChildren(child, [&](const std::string& cname, const xmlNodePtr element) {
          if (cname == attrib_tag)
            getPtclAttrib(element, 0, nat_group[ig], start);
        });
        start += nat_group[ig];
        ig++;
      }
      else if (cname == attrib_tag)
      {
        if (nat_group.size() > 1)
          throw UniformCommunicateError("An 'attrib' XML element node was found outside 'group'"
                                        " without XML element node named 'ionid'."
                                        " Cannot map particles to more than one species. Check XML input!");
        getPtclAttrib(child, 0, nat, 0);
      }
    });
  }

  if (ref_.getLattice().SuperCellEnum)
  {
    if (randomizeR == "yes")
    {
      makeUniformRandom(ref_.R);
      ref_.R.setUnit(PosUnit::Lattice);
      ref_.convert2Cart(ref_.R);
#if !defined(QMC_CUDA)
      makeUniformRandom(ref_.spins);
      ref_.spins *= 2 * M_PI;
#endif
    }
    else // put them [0,1) in the cell
      ref_.applyBC(ref_.R);
  }

  //this sets Mass, Z
  ref_.resetGroups();
  ref_.createSK();

  return true;
}

void XMLParticleParser::checkGrouping(int nat, const std::vector<int>& nat_group) const
{
  app_debug() << "There are " << nat << " particles in " << nat_group.size() << " species containing:" << std::endl;
  for (int ig = 0; ig < nat_group.size(); ig++)
  {
    const auto& group_species_name = ref_.getSpeciesSet().getSpeciesName(ig);
    if (nat_group[ig] == 0)
      throw UniformCommunicateError("Element '" + group_species_name + "' was provided but never referenced.");
    app_debug() << "    " << nat_group[ig] << " '" << group_species_name << "'" << std::endl;
  }

  if (std::accumulate(nat_group.begin(), nat_group.end(), 0) != nat)
    throw UniformCommunicateError(
        "The total number of particles doesn't match the sum of the particle counts of all the species.");
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
  ReportEngine PRE("XMLParticleParser", "reset");
  SpeciesSet& tspecies(ref_.getSpeciesSet());
  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "group")
    {
      std::string sname;
      OhmmsAttributeSet gAttrib;
      gAttrib.add(sname, "name");
      gAttrib.put(cur);
      if (sname.size())
      {
        int sid = tspecies.addSpecies(sname);
        setSpeciesProperty(tspecies, sid, cur);
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

  inline ParticleAttribXmlNode(PAT& a, PosUnit utype) : ref_(a) { ref_.InUnit = utype; }

  inline bool put(xmlNodePtr cur, int in_offset, int copy_size, int out_offset)
  {
    using data_type = typename PAT::Type_t;
    std::vector<data_type> data_in;
    putContent(data_in, cur);
    if (data_in.size() < in_offset + copy_size)
    {
      std::ostringstream msg;
      msg << "Insufficient data to copy from XML input which holds " << data_in.size() << " entries."
          << " Need to copy from [" << in_offset << ", " << in_offset + copy_size << ")." << std::endl;
      throw UniformCommunicateError(msg.str());
    }
    std::copy_n(data_in.begin() + in_offset, copy_size, ref_.begin() + out_offset);
    return true;
  }
};

void XMLParticleParser::getPtclAttrib(xmlNodePtr cur, int in_offset, int copy_size, int out_offset)
{
  std::string oname, otype;
  int utype   = 0;
  int size_in = 0;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(otype, datatype_tag);  //datatype
  pAttrib.add(oname, "name");        //name
  pAttrib.add(utype, condition_tag); //condition
  pAttrib.add(size_in, "size");      //size
  pAttrib.put(cur);
  if (oname.empty() || otype.empty())
  {
    app_error() << "   Missing attrib/@name or attrib/@datatype " << std::endl;
    app_error() << "     <attrib name=\"aname\"  datatype=\"atype\"/>" << std::endl;
    return;
  }
  int t_id = ref_AttribList.getAttribType(otype);

  if (oname == ionid_tag)
    throw UniformCommunicateError("'ionid' should not be parsed by getPtclAttrib.");
  else
  {
    //very permissive in that a unregistered attribute will be created and stored by ParticleSet
    //cloning is not going to work
    if (t_id == PA_IndexType)
    {
      ParticleIndex* obj = nullptr;
      obj                = ref_AttribList.getAttribute(otype, oname, obj);
      ParticleAttribXmlNode<ParticleIndex> a(*obj, static_cast<PosUnit>(utype));
      a.put(cur, in_offset, copy_size, out_offset);
    }
    else if (t_id == PA_ScalarType)
    {
      ParticleScalar* obj = nullptr;
      obj                 = ref_AttribList.getAttribute(otype, oname, obj);
      ParticleAttribXmlNode<ParticleScalar> a(*obj, static_cast<PosUnit>(utype));
      a.put(cur, in_offset, copy_size, out_offset);
    }
    else if (t_id == PA_PositionType)
    {
      ParticlePos* obj = nullptr;
      obj              = ref_AttribList.getAttribute(otype, oname, obj);
      ParticleAttribXmlNode<ParticlePos> a(*obj, static_cast<PosUnit>(utype));
      a.put(cur, in_offset, copy_size, out_offset);
    }
    else if (t_id == PA_TensorType)
    {
      ParticleTensor* obj = nullptr;
      obj                 = ref_AttribList.getAttribute(otype, oname, obj);
      ParticleAttribXmlNode<ParticleTensor> a(*obj, static_cast<PosUnit>(utype));
      a.put(cur, in_offset, copy_size, out_offset);
    }
  }
}


XMLSaveParticle::XMLSaveParticle(Particle_t& pin) : ref_(pin) {}

XMLSaveParticle::~XMLSaveParticle() {}

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
  get(fxml, 1);
  fxml.close();
}

void XMLSaveParticle::get(std::ostream& fxml, int olevel) const
{
  ref_.begin_node(fxml);
  fxml.setf(std::ios::scientific);
  fxml.precision(15);
  LatticeXMLWriter latticeout(ref_.getLattice());
  latticeout.get(fxml);
  for (int i = 0; i < SpeciesName.size(); i++)
  {
    fxml << "<group name=\"" << SpeciesName[i] << "\"/>" << std::endl;
  }
  //only write the local particles
  int nloc = ref_.getTotalNum();
  if (olevel)
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
        const ParticleIndex* itmp=dynamic_cast<ParticleIndex*>(ooref);
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
        const ParticleScalar* stmp=dynamic_cast<ParticleScalar*>(ooref);
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
        const ParticlePos* rtmp=dynamic_cast<ParticlePos*>(ooref);
        for(int iat=0; iat<nloc; iat++)
        {
          fxml << (*rtmp)[iat] << std::endl;
        }
      }
      else if (t_id == PA_TensorType)
      {
        fxml.precision(15);
        const ParticleTensor* ttmp=dynamic_cast<ParticleTensor*>(ooref);
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
    for (int iat = 0; iat < nloc; iat++)
      fxml << ref_.R[iat] << std::endl;
    ref_.R.end_node(fxml);
    ref_.GroupID.begin_node(fxml);
    for (int iat = 0; iat < nloc; iat++)
    {
      fxml << ref_.GroupID[iat] << " ";
      if (iat % 20 == 19)
        fxml << std::endl;
    }
    if (nloc % 20 != 19)
      fxml << std::endl;
    ref_.GroupID.end_node(fxml);
  }
  ref_.end_node(fxml);
}

bool XMLSaveParticle::put(xmlNodePtr cur) { return true; }

/** create particleset node
 * @param addlattice if true, add unitcell
 */
xmlNodePtr XMLSaveParticle::createNode(bool addlattice)
{
  if (SpeciesName.size() != ref_.getSpeciesSet().getTotalNum())
  {
    SpeciesName = ref_.getSpeciesSet().speciesName;
  }
  //if(addlattice) {
  //  ref_.getLattice().print(std::cout);
  //}
  xmlNodePtr cur = xmlNewNode(NULL, (const xmlChar*)"particleset");
  xmlNewProp(cur, (const xmlChar*)"name", (const xmlChar*)ref_.getName().c_str());
  if (ref_.groups() > 1)
  {
    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for (int is = 0; is < nspecies; is++)
    {
      std::ostringstream ng;
      ng << ref_.last(is) - ref_.first(is);
      xmlNodePtr g = xmlNewNode(NULL, (const xmlChar*)"group");
      xmlNewProp(g, (const xmlChar*)"name", (const xmlChar*)SpeciesName[is].c_str());
      xmlNewProp(g, (const xmlChar*)"size", (const xmlChar*)ng.str().c_str());
      for (int item = 0; item < nitem; item++)
      {
        std::ostringstream prop;
        prop << mySpecies(item, is);
        xmlNodePtr p = xmlNewTextChild(g, NULL, (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p, (const xmlChar*)"name", (const xmlChar*)mySpecies.attribName[item].c_str());
      }
      std::ostringstream pos;
      pos.setf(std::ios_base::scientific);
      pos << "\n";
      for (int iat = ref_.first(is); iat < ref_.last(is); iat++)
      {
        pos << ref_.R[iat] << std::endl;
      }
      xmlNodePtr posPtr = xmlNewTextChild(g, NULL, (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
      xmlNewProp(posPtr, (const xmlChar*)"name", (const xmlChar*)"position");
      xmlNewProp(posPtr, (const xmlChar*)"datatype", (const xmlChar*)"posArray");
      xmlAddChild(cur, g);
    }
  }
  else
  {
    std::ostringstream nat;
    nat << ref_.getTotalNum();
    xmlNewProp(cur, (const xmlChar*)"size", (const xmlChar*)nat.str().c_str());
    const SpeciesSet& mySpecies(ref_.getSpeciesSet());
    int nitem(mySpecies.numAttributes());
    int nspecies(mySpecies.getTotalNum());
    for (int is = 0; is < nspecies; is++)
    {
      xmlNodePtr g = xmlNewNode(NULL, (const xmlChar*)"group");
      xmlNewProp(g, (const xmlChar*)"name", (const xmlChar*)SpeciesName[is].c_str());
      for (int item = 0; item < nitem; item++)
      {
        std::ostringstream prop;
        prop << mySpecies(item, is);
        xmlNodePtr p = xmlNewTextChild(g, NULL, (const xmlChar*)"parameter", (const xmlChar*)prop.str().c_str());
        xmlNewProp(p, (const xmlChar*)"name", (const xmlChar*)mySpecies.attribName[item].c_str());
      }
      xmlAddChild(cur, g);
    }
    std::ostringstream pos, gid;
    pos.setf(std::ios_base::scientific);
    pos << "\n";
    for (int iat = 0; iat < ref_.getTotalNum(); iat++)
    {
      pos << ref_.R[iat] << std::endl;
    }
    xmlNodePtr posPtr = xmlNewTextChild(cur, NULL, (const xmlChar*)"attrib", (const xmlChar*)pos.str().c_str());
    xmlNewProp(posPtr, (const xmlChar*)"name", (const xmlChar*)"position");
    xmlNewProp(posPtr, (const xmlChar*)"datatype", (const xmlChar*)"posArray");
    gid << "\n ";
    for (int iat = 0; iat < ref_.getTotalNum(); iat++)
    {
      gid << SpeciesName[ref_.GroupID[iat]] << " ";
    }
    gid << std::endl;
    xmlNodePtr gPtr = xmlNewTextChild(cur, NULL, (const xmlChar*)"attrib", (const xmlChar*)gid.str().c_str());
    xmlNewProp(gPtr, (const xmlChar*)"name", (const xmlChar*)"ionid");
    xmlNewProp(gPtr, (const xmlChar*)"datatype", (const xmlChar*)"stringArray");
  }
  return cur;
}
} // namespace qmcplusplus

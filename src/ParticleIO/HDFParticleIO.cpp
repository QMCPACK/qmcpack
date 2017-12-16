//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include "ParticleIO/ParticleLayoutIO.h"
#include "ParticleIO/HDFParticleIO.h"
#include "ParticleBase/ParticleUtility.h"
//using namespace qmcplusplus;

#if defined(HAVE_LIBHDF5)
#include "ParticleIO/HDFParticleAttrib.h"
namespace qmcplusplus
{
bool HDFParticleParser::put(const char* fname)
{
  xmlDocPtr doc;
  xmlNsPtr ns;
  xmlNodePtr cur;
  // build an XML tree from a the file;
  doc = xmlParseFile(fname);
  if (doc == NULL)
  {
    ERRORMSG(fname << " does not exist")
    return false;
  }
  ///using XPath instead of recursive search
  xmlXPathContextPtr context;
  xmlXPathObjectPtr result;
  context = xmlXPathNewContext(doc);
  result = xmlXPathEvalExpression((const xmlChar*)"//particleset",context);
  if(!xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    cur = result->nodesetval->nodeTab[0];
    put(cur);
  }
  xmlXPathFreeObject(result);
  return true;
}

bool
HDFParticleParser::put(xmlNodePtr cur)
{
  typedef Particle_t::ParticleIndex_t ParticleIndex_t;
  typedef Particle_t::ParticleScalar_t ParticleScalar_t;
  typedef Particle_t::ParticlePos_t ParticlePos_t;
  std::string hfile;
  if(xmlHasProp(cur, (const xmlChar *) "href"))
  {
    hfile = (const char*)(xmlGetProp(cur, (const xmlChar *) "href"));
  }
  if(xmlHasProp(cur, (const xmlChar *) "src"))
  {
    hfile = (const char*)(xmlGetProp(cur, (const xmlChar *) "src"));
  }
  hid_t dataset, dataspace, ret;
  hid_t hdfFile = H5Fopen(hfile.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  int natom = 0;
  dataset = H5Dopen(hdfFile, "atnum");
  ret = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &natom);
  H5Dclose(dataset);
  ref_.create(natom);
  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)(cur->name));
    /**@warning transition from UnitCell -> unitcell */
    if (cname == "UnitCell" || cname == "unitcell")
    {
      LatticeParser lat(ref_.Lattice);
      lat.put(cur);
    }
    else
      if (cname == "attrib")
      {
        std::string oname(null_tag), otype(null_tag);
        int utype = 0;
        ///process attributes: type or format
        xmlAttrPtr att = cur->properties;
        while(att != NULL)
        {
          std::string aname((const char*)(att->name));
          if(aname == "name")
          {
            oname = (const char*)(att->children->content);
          }
          else
            if(aname == datatype_tag)
            {
              otype = (const char*)(att->children->content);
            }
            else
              if(aname == condition_tag)
              {
                utype = atoi( (const char*)(att->children->content));
              }
          att = att->next;
        }
        if(oname == null_tag || otype == null_tag)
          break;
        if(otype == scalartype_tag)
        {
          ParticleScalar_t& stmp = *(ref_.getScalarAttrib(oname.c_str()));
          stmp.setUnit(utype);
          HDFAttribIO<ParticleScalar_t> hdf_s(stmp);
          hdf_s.read(hdfFile,oname.c_str());
        }
        else
          if(otype == postype_tag)
          {
            ParticlePos_t& ptmp = *(ref_.getVectorAttrib(oname.c_str()));
            ptmp.setUnit(utype);
            HDFAttribIO<ParticlePos_t> hdf_p(ptmp);
            hdf_p.read(hdfFile,oname.c_str());
          }
          else
            if(otype == indextype_tag)
            {
              ParticleIndex_t& itmp = *(ref_.getIndexAttrib(oname.c_str()));
              itmp.setUnit(utype);
              HDFAttribIO<ParticleIndex_t> hdf_i(itmp);
              hdf_i.read(hdfFile,oname.c_str());
            }
      }
    cur = cur->next;
  }
  H5Fclose(hdfFile);
  ref_.update(-1);
  return true;
}

HDFSaveParticle::~HDFSaveParticle()
{
}

void HDFSaveParticle::reset(const char* fileroot, bool append)
{
  FileRoot = fileroot;
}


void HDFSaveParticle::report(int iter)
{
//  typedef Particle_t::ParticleIndex_t ParticleIndex_t;
//  typedef Particle_t::ParticleScalar_t ParticleScalar_t;
//  typedef Particle_t::ParticlePos_t ParticlePos_t;
//
//  std::string metafile = FileRoot;
//  metafile.append(".xml");
//
//  std::string hfile = FileRoot;
//  hfile.append(".h5");
//
//  // writing a meta file
//  std::ofstream ftbmd(metafile.c_str()); // always overwrite
//  ftbmd << "<?xml version=\"1.0\"?>" << std::endl;
//  ftbmd << "<particleset src=\"" << hfile << "\" size=\"" << ref_.getTotalNum()
//	<<"\" format=\"hdf5\">" << std::endl;
//  ftbmd.setf(std::ios_base::scientific);
//  ftbmd.precision(15);
//  LatticeXMLWriter latticeout(ref_.Lattice);
//  latticeout.get(ftbmd);
//  Particle_t::PAListIterator it = ref_.first_attrib();
//  while(it != ref_.last_attrib()) {
//    OhmmsObject& ooref= *((*it).second);
//    ooref.begin_node(ftbmd);
//    ooref.end_node(ftbmd);
//    it++;
//  }
//  ftbmd << "</particleset>" << std::endl;
//  ftbmd.close();
//
//  // creat HDF file handler
//  hid_t hdfFile =
//    H5Fcreate(hfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//
//  hsize_t dim = 1;
//  int natom = ref_.getLocalNum();
//  hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
//  hid_t dataset =
//  H5Dcreate(hdfFile, "atnum", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
//  hid_t ret =
//      H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&natom);
//  H5Dclose(dataset);
//  H5Sclose(dataspace);
//
//  // write all the extra ParticleAttrib's added to the ParticleList
//  it = ref_.first_attrib();
//  while(it != ref_.last_attrib()) {
//    OhmmsObject* ooref= (*it).second;
//    int t_id = ref_.getAttribType(ooref->typeName());
//    int o_id = ooref->id();
//    if(t_id == PA_IndexType) {
//      HDFAttribIO<ParticleIndex_t> hdf_i(*(ref_.getIndexAttrib(o_id)));
//      hdf_i.write(hdfFile,ooref->objName().c_str());
//    } else if(t_id == PA_ScalarType) {
//      HDFAttribIO<ParticleScalar_t> hdf_s(*(ref_.getScalarAttrib(o_id)));
//      hdf_s.write(hdfFile,ooref->objName().c_str());
//    } else if(t_id == PA_PositionType) {
//      HDFAttribIO<ParticlePos_t> hdf_p(*(ref_.getVectorAttrib(o_id)));
//      hdf_p.write(hdfFile,ooref->objName().c_str());
//    } else if(t_id == PA_TensorType) {
//      WARNMSG("Not working yet")
//    }
//    it++;
//  }
//  H5Fclose(hdfFile);
}

bool HDFSaveParticle::put(xmlNodePtr cur)
{
  return true;
}
}
#else
bool HDFParticleParser::put(const char* fname)
{
  return false;
}
bool
HDFParticleParser::put(xmlNodePtr cur)
{
  return false;
}

HDFSaveParticle::~HDFSaveParticle()
{
}
void HDFSaveParticle::reset(const char* fileroot, bool append)
{
}
void HDFSaveParticle::report(int iter)
{
}
bool HDFSaveParticle::put(xmlNodePtr cur)
{
  return true;
}
}
#endif


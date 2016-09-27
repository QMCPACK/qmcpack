//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
template<class T>
inline bool
putContent2(std::vector<T>& a, xmlNodePtr cur)
{
  std::istringstream
  stream((const char*)
         (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
  T temp;
  a.clear();
  while(!stream.eof())
  {
    stream >> temp;
    if ( stream.fail() || stream.bad())
      break;
    else
      a.push_back(temp);
  }
  return a.size() > 0;
}


bool
kSpaceJastrowBuilder::put(xmlNodePtr cur)
{
  xmlNodePtr kids = cur->xmlChildrenNode;
  kSpaceJastrow::SymmetryType oneBodySymm, twoBodySymm;
  RealType kc1, kc2;
  std::string symm1_opt, symm2_opt, id1_opt, id2_opt, spin1_opt("no"), spin2_opt("no");
  std::vector<RealType> oneBodyCoefs, twoBodyCoefs;
  // Initialize options
  kc1 = kc2 = 0.0;
  oneBodySymm = twoBodySymm = kSpaceJastrow::CRYSTAL;
  id1_opt = "cG1";
  id2_opt = "cG2";
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    std::vector<RealType> *coefs(NULL);
    std::string* id_opt(NULL);
    if (kidsname == "correlation")
    {
      std::string type_opt;
      OhmmsAttributeSet attrib;
      attrib.add (type_opt, "type");
      attrib.put(kids);
      if (type_opt == "One-Body")
      {
        attrib.add (symm1_opt, "symmetry");
        attrib.add (kc1, "kc");
        attrib.add (spin1_opt, "spinDependent");
        coefs = &oneBodyCoefs;
        id_opt = &id1_opt;
      }
      else
        if (type_opt == "Two-Body")
        {
          attrib.add (symm2_opt, "symmetry");
          attrib.add (kc2, "kc");
          attrib.add (spin2_opt, "spinDependent");
          coefs = &twoBodyCoefs;
          id_opt = &id2_opt;
        }
        else
          app_warning() << "  Unrecognized kSpace type \"" << type_opt
                        << "\" in kSpaceJastrowBuilder::put(xmlNotPtr cur).\n";
      attrib.put (kids);
      // Read the coefficients
      if (coefs)
      {
        xmlNodePtr xmlCoefs = kids->xmlChildrenNode;
        while (xmlCoefs != NULL)
        {
          std::string cname((const char*)xmlCoefs->name);
          if (cname == "coefficients")
          {
            std::string type("0"), id("0");
            OhmmsAttributeSet cAttrib;
            cAttrib.add(*id_opt, "id");
            cAttrib.add(type, "type");
            cAttrib.put(xmlCoefs);
            if (type != "Array")
            {
              app_error() << "Unknown coefficients type """ << type
                          << """ in kSpaceJastrowBuilder.\n"<< "Resetting to ""Array"".\n";
              xmlNewProp (xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
            }
            //vector<T> can be read by this
            putContent2(*coefs,xmlCoefs);
            app_log() << "  Read " << coefs->size() << " coefficients for type "
                      << type_opt << std::endl;
          }
          xmlCoefs = xmlCoefs->next;
        }
      }
    }
    else
      if (kidsname != "text")
      {
        app_warning() << "Unrecognized section \"" << kidsname
                      << "\" in kSpaceJastrowBuilder.\n";
      }
    kids = kids->next;
  }
  // Now build the kSpaceJastrow
  std::map<std::string,kSpaceJastrow::SymmetryType>::iterator symm1 =
    SymmMap.find(symm1_opt);
  if (symm1 != SymmMap.end())
    oneBodySymm = symm1->second;
  std::map<std::string,kSpaceJastrow::SymmetryType>::iterator symm2 =
    SymmMap.find(symm2_opt);
  if (symm2 != SymmMap.end())
    twoBodySymm = symm2->second;
  kSpaceJastrow *jastrow =
    new kSpaceJastrow(sourcePtcl, targetPtcl,
                      oneBodySymm, kc1, id1_opt, spin1_opt=="yes",
                      twoBodySymm, kc2, id2_opt, spin2_opt=="yes");
  jastrow->setCoefficients (oneBodyCoefs, twoBodyCoefs);
  //jastrow->addOptimizables(targetPsi.VarList);
  targetPsi.addOrbital(jastrow,"kSpace");
  return true;
}
}

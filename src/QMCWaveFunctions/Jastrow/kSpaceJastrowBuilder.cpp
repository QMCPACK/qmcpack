//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: kesler@ciw.edu
//   Tel:    
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  bool
  kSpaceJastrowBuilder::put(xmlNodePtr cur)
  {
    cerr << "In  kSpaceJastrowBuilder::put(xmlNodePtr cur).\n";
    xmlNodePtr kids = cur->xmlChildrenNode;
    kSpaceJastrow::SymmetryType oneBodySymm, twoBodySymm;
    RealType kc1, kc2;
    string symm1_opt, symm2_opt, kc1_opt, kc2_opt;
    // Initialize options
    kc1 = kc2 = 0.0;
    oneBodySymm = twoBodySymm = kSpaceJastrow::CRYSTAL;
    
    while (kids != NULL) {
      std::string kidsname = (char*)kids->name;
      if (kidsname == "correlation") {
	string type_opt;
	OhmmsAttributeSet attrib;
	attrib.add (type_opt, "type");
	attrib.put(kids);
	if (type_opt == "One-Body") {
	  attrib.add (symm1_opt, "symmetry");
	  attrib.add (kc1, "kc");
	}
	else if (type_opt == "Two-Body") {
	  attrib.add (symm2_opt, "symmetry");
	  attrib.add (kc2, "kc");
	}
	else 
	  app_warning() << "  Unrecognized kSpace type \"" << type_opt 
			<< "\" in kSpaceJastrowBuilder::put(xmlNotPtr cur).\n";
	attrib.put (kids);
      }
      else if (kidsname != "text") {
	app_warning() << "Unrecognized section \"" << kidsname 
		      << "\" in kSpaceJastrowBuilder.\n";
      }
      kids = kids->next;
    }
    // Now build the kSpaceJastrow
    std::map<string,kSpaceJastrow::SymmetryType>::iterator symm1 = 
      SymmMap.find(symm1_opt);
    if (symm1 != SymmMap.end())
      oneBodySymm = symm1->second;
    
    std::map<string,kSpaceJastrow::SymmetryType>::iterator symm2 = 
      SymmMap.find(symm2_opt);
    if (symm2 != SymmMap.end())
      twoBodySymm = symm2->second;
    
    kSpaceJastrow *jastrow = new kSpaceJastrow(sourcePtcl, targetPtcl,
					       oneBodySymm, kc1,
					       twoBodySymm, kc2);
    return true;
  }
}

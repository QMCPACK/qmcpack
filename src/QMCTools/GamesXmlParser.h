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
    
    



#ifndef QMCPLUSPLUS_TOOLS_GAMESS_XML_PARSER_H
#define QMCPLUSPLUS_TOOLS_GAMESS_XML_PARSER_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class GamesXmlParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{


  void getGeometry(std::vector<xmlNodePtr>&);
  void getGaussianCenters(std::vector<xmlNodePtr>&);
  void getEigVectors(std::vector<xmlNodePtr>&);

  void getControlParameters(xmlNodePtr);

public:

  GamesXmlParser();

  GamesXmlParser(int argc, char** argv);

  void parse(const std::string& fname);
};
#endif

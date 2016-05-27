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

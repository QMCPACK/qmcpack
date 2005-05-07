#ifndef OHMMS_TOOLS_GAMESS_XML_PARSER_H
#define OHMMS_TOOLS_GAMESS_XML_PARSER_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class GamesXmlParser: public QMCGaussianParserBase, 
                    public OhmmsAsciiParser {


  void getGeometry(vector<xmlNodePtr>&);
  void getGaussianCenters(vector<xmlNodePtr>&);
  void getEigVectors(vector<xmlNodePtr>&);

public:

  GamesXmlParser();

  GamesXmlParser(int argc, char** argv);

  void parse(const std::string& fname);
};
#endif

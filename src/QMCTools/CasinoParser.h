#ifndef OHMMS_TOOLS_CASINOPARSER_H
#define OHMMS_TOOLS_CASINOPARSER_H
#include "Tools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class CasinoParser: public QMCGaussianParserBase, 
                    public OhmmsAsciiParser {

public:

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  //Specialized functions
  void getNumberOfAtoms(std::istream& is);

  void getAtomicPositions(std::istream& is);

  void getAtomicNumbers(std::istream& is);

  void getValenceCharges(std::istream& is);

};
#endif

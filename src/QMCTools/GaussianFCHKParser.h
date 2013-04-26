#ifndef QMCPLUSPLUS_TOOLS_GAUSSIAN_FCHK_H
#define QMCPLUSPLUS_TOOLS_GAUSSIAN_FCHK_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class GaussianFCHKParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  GaussianFCHKParser();

  GaussianFCHKParser(int argc, char** argv);

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);
};
#endif

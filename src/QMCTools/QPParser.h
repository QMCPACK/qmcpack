//#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
//#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class QPParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  QPParser();

  QPParser(int argc, char** argv);

  std::streampos pivot_begin;
  std::vector<std::string> tags;
  bool usingECP;
  std::string MOtype;
  int readtype, numAO;
  int NFZC, NEXT, NTOT, NAC;

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  void getMO(std::istream& is);

  void getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal_alpha);

  void getQPCI(std::istream& is);


};
//#endif

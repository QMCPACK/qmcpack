#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class GamesAsciiParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  GamesAsciiParser();

  GamesAsciiParser(int argc, char** argv);

  streampos pivot_begin;
  vector<std::string> tags;
  bool usingECP;
  std::string MOtype;
  //int nCartMO;
  int readtype;
  int NFZC, NEXT, NTOT, NAC;

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  void getMO(std::istream& is);

  void getCI(std::istream& is);

  void getORMAS(std::istream& is);

  void getCSF(std::istream& is);

  double getCSFSign(vector<int>&);

};
#endif

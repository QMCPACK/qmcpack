//#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
//#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class VSVBParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  VSVBParser();

  VSVBParser(int argc, char** argv);

  streampos pivot_begin;
  vector<std::string> tags;
  bool usingECP;
  std::string MOtype;
  //int nCartMO;
  int readtype;
  int NFZC, NEXT, NTOT, NAC;
  int NumSpinCoupledOrbitals,NumOpenShellOrbitals,NumDoubleOccupiedOrbitals;

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  void getMO(std::istream& is);

  void getMDVSVB(std::istream& is,int NbVbStructures);
   
  std::string getOccup(int val);

};
//#endif

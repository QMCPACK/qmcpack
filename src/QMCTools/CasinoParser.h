//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_TOOLS_CASINOPARSER_H
#define QMCPLUSPLUS_TOOLS_CASINOPARSER_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class CasinoParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

  std::vector<double> BasisCorrection;

public:

  CasinoParser();

  CasinoParser(int argc, char** argv);

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  //Specialized functions
  void getNumberOfAtoms(std::istream& is);

  void getAtomicPositions(std::istream& is);

  void getAtomicNumbers(std::istream& is);

  void getValenceCharges(std::istream& is);

  double contractionCorrection(int shell_id, double alpha);

  void makeCorrections();
};
#endif

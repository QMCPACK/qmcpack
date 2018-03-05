//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//  
//
// File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



//#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
//#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class PyscfParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  PyscfParser();

  PyscfParser(int argc, char** argv);

  std::streampos pivot_begin;
  std::vector<std::string> tags;
  std::string MOtype;
  int readtype, numAO;
  int NFZC, NEXT, NTOT, NAC;

  void parse(const std::string& fname);

  void getGeometry(const std::string& fname);


};
//#endif

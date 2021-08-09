//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//
// File created by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"

class RMGParser : public QMCGaussianParserBase, public OhmmsAsciiParser
{
public:
  RMGParser();

  RMGParser(int argc, char** argv);

  std::vector<std::string> ECP_names;
  //std::vector<std::string> tags;
  int NumberOfSpins;
  void dumpPBC(const std::string& psi_tag, const std::string& ion_tag);


  // void dumpStdInput(const std::string& psi_tag, const std::string& ion_tag);

  //void dumpStdInputProd(const std::string& psi_tag, const std::string& ion_tag);
  void parse(const std::string& fname);
  void getCell(const std::string& fname);
};
//#endif

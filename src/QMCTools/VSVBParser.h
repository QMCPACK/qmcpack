//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
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

class VSVBParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  VSVBParser();

  VSVBParser(int argc, char** argv);

  std::streampos pivot_begin;
  std::vector<std::string> tags;
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

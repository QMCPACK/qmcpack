//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#define QMCPLUSPLUS_TOOLS_GAMESS_OUT_H
#include "QMCTools/QMCGaussianParserBase.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"

class GamesAsciiParser: public QMCGaussianParserBase,
  public OhmmsAsciiParser
{

public:

  GamesAsciiParser();

  GamesAsciiParser(int argc, char** argv);

  std::streampos pivot_begin;
  std::vector<std::string> tags;
  bool usingECP;
  std::string MOtype;
  //int nCartMO;
  int readtype;
  int NFZC, NEXT, NTOT, NAC;

  void parse(const std::string& fname);

  void getGeometry(std::istream& is);

  void getGaussianCenters(std::istream& is);

  void getMO(std::istream& is);

  void getMO_single_set(std::istream& is, Matrix<double> &CartMat, std::vector<value_type>& EigVal_alpha);

  void getCI(std::istream& is);

  void getORMAS(std::istream& is);

  void getCSF(std::istream& is);

  double getCSFSign(std::vector<int>&);

};
#endif

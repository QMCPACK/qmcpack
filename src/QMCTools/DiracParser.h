//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TOOLS_DIRAC_OUT_H
#define QMCPLUSPLUS_TOOLS_DIRAC_OUT_H

#include "QMCTools/QMCGaussianParserBase.h"

typedef std::pair<double, double> primBasis;

struct basisGroup
{
  std::string rid;
  std::string type = "Gaussian";
  int l;
  int n;
  std::vector<primBasis> radfuncs;
};

struct atBasisSet
{
  std::string name       = "Gaussian";
  std::string normalized = "no";
  std::string grid_type  = "log";
  std::string expandYlm  = "Gamess";
  std::string elementType;
  std::string angular = "cartesian";
  double grid_ri      = 1e-6;
  double grid_rf      = 100;
  int grid_npts       = 1001;
  std::vector<basisGroup> basisGroups;
};


class DiracParser : public QMCGaussianParserBase, public OhmmsAsciiParser
{
public:
  DiracParser(int argc, char** argv);
  void parse(const std::string& fname);

private:
  void dumpHDF5(const std::string& fname);
  void getGeometry(std::istream& is);
  void getGaussianCenters(std::istream& is);
  void getSpinors(std::istream& is);
  std::streampos pivot_begin;
  int NumberOfSpecies;
  int version;
  std::string aline;
  std::vector<atBasisSet> basisset;
  std::vector<std::vector<std::complex<double>>> upcoeff;
  std::vector<std::vector<std::complex<double>>> dncoeff;
  int numAO, numMO;
};


#endif

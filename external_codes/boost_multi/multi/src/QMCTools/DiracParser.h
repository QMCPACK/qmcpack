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
#include <string>
#include <map>
#include <cmath>

using primBasis = std::pair<double, double>;

struct basisGroup
{
  int l;
  int n;
  std::vector<primBasis> radfuncs;
};

struct atBasisSet
{
  std::string elementType;
  std::vector<basisGroup> basisGroups;
};

enum class OrbType
{
  CORE,
  ACTIVE,
  VIRTUAL
};

struct fermIrrep
{
  std::string label;
  std::vector<double> energies;
  std::vector<OrbType> orbtypes;
  std::vector<std::vector<std::vector<double>>> spinor_mo_coeffs; //mo,ao,up_r,up_i,dn_r,dn_i
  fermIrrep(std::string label_in, int nSpinors, int numAO);
  fermIrrep generate_kramers_pair();
  int get_num_spinors() { return energies.size(); }
  int get_num_ao() { return spinor_mo_coeffs[0].size(); }
  std::string get_label() { return label; }
};

struct ciState
{
  double energy;
  std::vector<std::string> occstrings;
  std::vector<double> coeffs;
};

class cosciRep
{
public:
  cosciRep(std::string in_label, int nstates);
  std::string label;
  std::vector<ciState> states;
  void printInfo(std::ostream& os, int& tot_state_count);
};

class DiracParser : public QMCGaussianParserBase, public OhmmsAsciiParser
{
  using normMapType = std::map<std::string, double>;

public:
  DiracParser(int argc, char** argv);
  void parse(const std::string& fname) override;

private:
  void dumpHDF5(const std::string& fname);
  void getGeometry(std::istream& is);
  void getGaussianCenters(std::istream& is);
  void getSpinors(std::istream& is);
  void getWF(std::istream& is);
  void getCOSCI(std::istream& is);
  void getSingleDet(std::istream& is);
  std::streampos pivot_begin;
  int NumberOfSpecies;
  int version;
  std::string aline;
  std::vector<atBasisSet> basisset;
  normMapType normMap;

  std::vector<fermIrrep> irreps;
  std::vector<fermIrrep> kp_irreps;
  std::vector<cosciRep> cosciReps;

  void parseCOSCIOrbInfo(std::istream& is, const int irrep_idx, OrbType type);
  int sortAndStoreCOSCIOrbs(OrbType type, const int spinor_component);
};


#endif

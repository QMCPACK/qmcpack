//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TOOLS_EXTERNAL_GAUSSIANPARSERBASE_H
#define QMCPLUSPLUS_TOOLS_EXTERNAL_GAUSSIANPARSERBASE_H
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include "OhmmsData/OhmmsElementBase.h"
#include "Utilities/SimpleParser.h"
#include "Particle/ParticleSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "hdf/hdf_archive.h"

using namespace qmcplusplus;

struct QMCGaussianParserBase
{
  using value_type        = double;
  using SingleParticlePos = ParticleSet::SingleParticlePos;

  bool multideterminant;
  bool multidetH5;
  bool BohrUnit;
  bool SpinRestricted;
  bool Periodicity;
  bool UseHDF5;
  bool PBC;
  bool production;
  bool zeroCI;
  bool orderByExcitation;
  bool addJastrow;
  bool addJastrow3Body;
  bool ECP;
  bool debug;
  bool Structure;
  bool DoCusp;
  // if true, adjust valence electron count output based on gCoreTable
  bool FixValence;
  bool singledetH5;
  bool optDetCoeffs;
  bool usingCSF;
  bool isSpinor;
  int IonChargeIndex;
  int ValenceChargeIndex;
  int AtomicNumberIndex;
  int NumberOfAtoms;
  int NumberOfEls;
  // targeted state number
  int target_state;
  int SpinMultiplicity;
  int NumberOfAlpha, NumberOfBeta;
  int SizeOfBasisSet;
  // mmorales: number of Molecular orbitals, not always equal to SizeOfBasisSet
  int numMO, readNO, readGuess, numMO2print;
  int ci_size, ci_nca, ci_ncb, ci_nea, ci_neb, ci_nstates;
  int NbKpts;
  int nbexcitedstates;
  double ci_threshold;


  std::vector<double> STwist_Coord; //Super Twist Coordinates


  std::string Title;
  std::string basisType;
  std::string basisName;
  std::string Normalized;
  std::string CurrentCenter;
  std::string outputFile;
  std::string angular_type;
  std::string expandYlm;
  std::string h5file;
  std::string multih5file;
  std::string WFS_name;
  std::string CodeName;

  const SimulationCell simulation_cell;
  ParticleSet IonSystem;


  std::vector<std::string> GroupName;

  std::vector<int> gShell, gNumber, gBound;
  std::vector<int> Occ_alpha, Occ_beta;
  std::vector<value_type> Qv;
  std::vector<value_type> gExp, gC0, gC1;
  std::vector<value_type> EigVal_alpha, EigVal_beta;
  std::vector<value_type> EigVec;
  //std::vector<GaussianCombo<value_type> > gExp, gC0, gC1;
  //std::string EigVecU, EigVecD;
  std::unique_ptr<xmlNode, void (*)(xmlNodePtr)> gridPtr =
      std::unique_ptr<xmlNode, void (*)(xmlNodePtr)>(nullptr, nullptr);
  std::vector<std::string> CIalpha, CIbeta;
  std::vector<std::string> CSFocc;
  std::vector<std::vector<std::string>> CSFalpha, CSFbeta;
  std::vector<std::vector<double>> CSFexpansion;
  std::vector<double> CIcoeff;
  std::vector<double> X, Y, Z; //LAttice vectors for PBC
  std::vector<int> Image;

  std::vector<int> CIexcitLVL;

  std::vector<std::pair<int, double>> coeff2csf;

  QMCGaussianParserBase();
  QMCGaussianParserBase(int argc, char** argv);

  virtual ~QMCGaussianParserBase() = default;

  void setOccupationNumbers();

  void createGridNode(int argc, char** argv);

  void createSPOSets(xmlNodePtr, xmlNodePtr);
  void createSPOSetsH5(xmlNodePtr, xmlNodePtr);
  void PrepareSPOSetsFromH5(xmlNodePtr, xmlNodePtr);
  xmlNodePtr createElectronSet(const std::string& ion_tag);
  xmlNodePtr createIonSet();
  xmlNodePtr createCell();
  xmlNodePtr createHamiltonian(const std::string& ion_tag, const std::string& psi_tag);
  xmlNodePtr createBasisSet();
  xmlNodePtr createBasisSetWithHDF5();
  xmlNodePtr createCenter(int iat, int _off);
  void createCenterH5(int iat, int _off, int numelem);
  void createShell(int n, int ig, int off_, xmlNodePtr abasis);
  void createShellH5(int n, int ig, int off_, int numelem);

  xmlNodePtr createDeterminantSet();
  xmlNodePtr createMultiDeterminantSet();
  xmlNodePtr createDeterminantSetWithHDF5();
  xmlNodePtr createMultiDeterminantSetFromH5();
  xmlNodePtr createMultiDeterminantSetCIHDF5();
  xmlNodePtr PrepareDeterminantSetFromHDF5();
  xmlNodePtr createJ3();
  xmlNodePtr createJ2();
  xmlNodePtr createJ1();

  xmlNodePtr parameter(xmlNodePtr Parent, std::string Mypara, std::string a);

  int numberOfExcitationsCSF(std::string&);

  virtual void parse(const std::string& fname) = 0;

  virtual void dumpPBC(const std::string& psi_tag, const std::string& ion_tag);

  virtual void dump(const std::string& psi_tag, const std::string& ion_tag);

  void dumpStdInput(const std::string& psi_tag, const std::string& ion_tag);

  void dumpStdInputProd(const std::string& psi_tag, const std::string& ion_tag);


  //static std::vector<std::string> IonName;
  static std::map<int, std::string> IonName;

  static std::vector<std::string> gShellType;
  static std::vector<int> gShellID;

  static const std::vector<double> gCoreTable;

  static void init();
};
#endif

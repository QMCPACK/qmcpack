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
#include "OhmmsData/HDFStringAttrib.h"
#include "io/hdf_archive.h"

using namespace qmcplusplus;

struct QMCGaussianParserBase
{

  typedef double value_type;
  typedef ParticleSet::SingleParticlePos_t SingleParticlePos_t;

  bool multideterminant;
  bool AllH5;
  bool BohrUnit;
  bool SpinRestricted;
  bool Periodicity;
  bool UseHDF5;
  bool production;
  bool zeroCI;
  bool orderByExcitation;
  bool addJastrow;
  bool addJastrow3Body;
  bool ECP;
  bool debug;
  bool Structure;
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
// benali: Point Charge from FMO ESP 
  int * ESPIonChargeIndex;
  int * ESPValenceChargeIndex;
  int * ESPAtomicNumberIndex;
  int TotNumMonomer;
  ParticleSet *ESPSystem;
  std::vector <std::vector<double> > ESP;
  std::vector<std::vector<std::string> > ESPGroupName;
  xmlNodePtr createESPSet(int iesp);
  static std::map<int,std::string> ESPName;
  int FMOIndexI,FMOIndexJ,FMOIndexK;
  bool FMO, FMO1,FMO2,FMO3,DoCusp,FixValence,QP;
  

  std::string Title;
  std::string basisType;
  std::string basisName;
  std::string Normalized;
  std::string CurrentCenter;
  std::string outputFile;
  std::string angular_type;
  std::string h5file;
  std::string WFS_name;
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
  xmlNodePtr gridPtr;
  std::vector<std::string> CIalpha,CIbeta;
  std::vector<std::string> CSFocc;
  std::vector<std::vector<std::string> > CSFalpha,CSFbeta;
  std::vector<std::vector<double> > CSFexpansion;
  std::vector<double> CIcoeff;
  std::vector<int> CIexcitLVL;
  int ci_size,ci_nca,ci_ncb,ci_nea,ci_neb,ci_nstates;
  double ci_threshold;
  bool usingCSF;
  bool VSVB;

  std::vector<std::pair<int,double> > coeff2csf;

  QMCGaussianParserBase();
  QMCGaussianParserBase(int argc, char** argv);

  void setOccupationNumbers();

  void createGridNode(int argc, char** argv);

  void createSPOSets(xmlNodePtr,xmlNodePtr);
  void createSPOSetsH5(xmlNodePtr,xmlNodePtr);
  xmlNodePtr createElectronSet(const std::string& ion_tag);
  xmlNodePtr createIonSet();
  xmlNodePtr createHamiltonian(const std::string& ion_tag, const std::string& psi_tag);
  xmlNodePtr createBasisSet();
  xmlNodePtr createBasisSetWithHDF5();
  xmlNodePtr createCenter(int iat, int _off);
  void createCenterH5(int iat, int _off,int numelem);
  void createShell(int n, int ig, int off_, xmlNodePtr abasis);
  void createShellH5(int n, int ig, int off_,int numelem);
  xmlNodePtr createDeterminantSet();
  xmlNodePtr createMultiDeterminantSet();
  xmlNodePtr createMultiDeterminantSetVSVB();
  xmlNodePtr createMultiDeterminantSetQP();
  xmlNodePtr createMultiDeterminantSetQPHDF5();
  xmlNodePtr createDeterminantSetWithHDF5();
  xmlNodePtr PrepareDeterminantSetFromHDF5();
  xmlNodePtr createJ3();
  xmlNodePtr createJ2();
  xmlNodePtr createJ1();

  xmlNodePtr parameter(xmlNodePtr Parent, std::string Mypara ,std::string a);

  int numberOfExcitationsCSF( std::string&);

  void map2GridFunctors(xmlNodePtr cur);
  virtual void parse(const std::string& fname) = 0;

  virtual void dump(const std::string& psi_tag,
                    const std::string& ion_tag);

  void dumpStdInput(const std::string& psi_tag,
                    const std::string& ion_tag);

  void dumpStdInputProd(const std::string& psi_tag,
                    const std::string& ion_tag);

  virtual void Fmodump(const std::string& psi_tag,
                                 const std::string& ion_tag,
                                 std::string Mytag);

  //static std::vector<std::string> IonName;
  static std::map<int,std::string> IonName;

  static std::vector<std::string> gShellType;
  static std::vector<int>         gShellID;

  static void init();
};
#endif

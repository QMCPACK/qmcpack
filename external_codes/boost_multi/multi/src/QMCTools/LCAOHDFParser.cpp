//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//
// File created by: Anouar Benali, benali@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "LCAOHDFParser.h"
#include <algorithm>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <sstream>


LCAOHDFParser::LCAOHDFParser()
{
  basisName    = "Gaussian";
  Normalized   = "no";
  BohrUnit     = true;
  MOtype       = "Canonical";
  angular_type = "cartesian";
  readtype     = 0;
  NFZC         = 0;
}

LCAOHDFParser::LCAOHDFParser(int argc, char** argv) : QMCGaussianParserBase(argc, argv)
{
  basisName      = "Gaussian";
  Normalized     = "no";
  BohrUnit       = true;
  MOtype         = "Canonical";
  angular_type   = "cartesian";
  SpinRestricted = true;
  readtype       = 0;
  NFZC           = 0;
}

void LCAOHDFParser::parse(const std::string& fname)
{
  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }

  hin.push("application");
  hin.read(CodeName, "code");
  hin.pop();
  std::cout << "Converting Wavefunction from the " << CodeName << " Code" << std::endl;

  hin.push("PBC");
  hin.read(PBC, "PBC");
  hin.pop();

  std::cout << "Periodic Boundary Conditions: " << (PBC ? ("yes") : ("no")) << std::endl;
  std::cout.flush();


  hin.push("parameters");
  hin.read(ECP, "ECP");
  std::cout << "usingECP: " << (ECP ? ("yes") : ("no")) << std::endl;
  std::cout.flush();

  std::cout << "Multideterminants: " << (multideterminant ? ("yes") : ("no")) << std::endl;
  std::cout.flush();

  hin.read(SpinRestricted, "SpinRestricted");
  hin.read(numAO, "numAO");
  hin.read(numMO, "numMO");


  std::cout << "NUMBER OF AOs: " << numAO << std::endl;
  SizeOfBasisSet = numAO;
  std::cout << "Size of Basis Set: " << SizeOfBasisSet << std::endl;
  std::cout << "NUMBER OF MOs: " << numMO << std::endl;

  bool bohr = true;
  hin.read(bohr, "Unit");
  std::cout << "Unit in Bohr=" << bohr << std::endl;
  BohrUnit = bohr;

  hin.read(NumberOfAlpha, "NbAlpha");
  hin.read(NumberOfBeta, "NbBeta");
  hin.read(NumberOfEls, "NbTotElec");
  int ds;
  hin.read(ds, "spin");
  if (CodeName == "PySCF")
    SpinMultiplicity = ds + 1;
  else
    SpinMultiplicity = ds;

  std::cout << "Number of alpha electrons: " << NumberOfAlpha << std::endl;
  std::cout << "Number of beta electrons: " << NumberOfBeta << std::endl;
  std::cout << "Number of electrons: " << NumberOfEls << std::endl;
  std::cout << "SPIN MULTIPLICITY: " << SpinMultiplicity << std::endl;

  hin.pop();
  hin.push("atoms");
  hin.read(NumberOfAtoms, "number_of_atoms");
  std::cout << "NUMBER OF ATOMS: " << NumberOfAtoms << std::endl;
  hin.pop();

  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  Matrix<double> myvec(1, numMO);

  hin.push("Super_Twist");
  if (hin.readEntry(myvec, "eigenval_0"))
  {
    // eigenval_0 exists on file
    for (int i = 0; i < numMO; i++)
      EigVal_alpha[i] = myvec[0][i];

    //Reading Eigenvals for Spin unRestricted calculation. This section is needed to set the occupation numbers
    if (!SpinRestricted)
    {
      hin.read(myvec, "eigenval_1");
      for (int i = 0; i < numMO; i++)
        EigVal_beta[i] = myvec[0][i];
    }
  }
  else
  {
    app_warning() << "eigenval_0 doesn't exist in h5 file. Treat all values zero." << std::endl;
  }
  hin.close();


  IonSystem.create({NumberOfAtoms});
  GroupName.resize(NumberOfAtoms);
  if (PBC)
  {
    getCell(fname);
    getSuperTwist(fname);
    if (debug)
    {
      getGaussianCenters(fname);
      getMO(fname);
    }
  }
  getGeometry(fname);

  if (multideterminant)
  {
    hdf_archive hin;

    if (!hin.open(outputFile.c_str(), H5F_ACC_RDONLY))
    {
      std::cerr << "Could not open H5 file" << std::endl;
      abort();
    }

    hin.push("MultiDet");

    hin.read(ci_size, "NbDet");
    hin.read(ci_nstates, "nstate");
    hin.read(nbexcitedstates, "nexcitedstate");
    CIcoeff.clear();
    CIalpha.clear();
    CIbeta.clear();
    CIcoeff.resize(ci_size);
    CIalpha.resize(ci_size);
    CIbeta.resize(ci_size);

    int ds  = SpinMultiplicity - 1;
    int neb = (NumberOfEls - ds) / 2;
    int nea = NumberOfEls - NumberOfBeta;
    ci_nea  = NumberOfAlpha;
    ci_neb  = NumberOfBeta;
    ci_nca  = nea - ci_nea;
    ci_ncb  = neb - ci_neb;
    std::cout << " Done reading CIs!!" << std::endl;
    hin.close();
  }
}

void LCAOHDFParser::getCell(const std::string& fname)
{
  X.resize(3);
  Y.resize(3);
  Z.resize(3);

  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }
  hin.push("Cell");
  Matrix<double> LatticeVec(3, 3);
  hin.read(LatticeVec, "LatticeVectors");
  X[0] = LatticeVec[0][0];
  X[1] = LatticeVec[0][1];
  X[2] = LatticeVec[0][2];
  Y[0] = LatticeVec[1][0];
  Y[1] = LatticeVec[1][1];
  Y[2] = LatticeVec[1][2];
  Z[0] = LatticeVec[2][0];
  Z[1] = LatticeVec[2][1];
  Z[2] = LatticeVec[2][2];
  hin.close();
  std::cout << "Lattice parameters in Bohr:" << std::endl;
  std::cout << X[0] << "  " << X[1] << "  " << X[2] << std::endl;
  std::cout << Y[0] << "  " << Y[1] << "  " << Y[2] << std::endl;
  std::cout << Z[0] << "  " << Z[1] << "  " << Z[2] << std::endl;
}
void LCAOHDFParser::getGeometry(const std::string& fname)
{
  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }
  hin.push("atoms");

  //atomic numbers
  std::vector<int> atomic_number;
  std::vector<double> q, pos;
  tags.clear();
  //read atomic info
  //MAP: (Atom Number, Species Index)
  std::vector<int> idx(NumberOfAtoms);
  std::map<int, int> AtomIndexmap;
  hin.read(idx, "species_ids");
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    AtomIndexmap.insert(std::pair<int, int>(i, idx[i]));
  }
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    std::string speciesName("species_");
    speciesName = speciesName + std::to_string(AtomIndexmap[i]);
    hin.push(speciesName);
    int zint;
    double z;
    std::string Name;
    hin.read(zint, "atomic_number");
    atomic_number.push_back(zint);
    hin.read(z, "charge");
    q.push_back(z);
    hin.read(Name, "name");
    tags.push_back(Name);
    hin.pop();
  }

  Matrix<double> IonPos(NumberOfAtoms, 3);
  hin.read(IonPos, "positions");

  SpeciesSet& species(IonSystem.getSpeciesSet());
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    IonSystem.R[i][0]                     = IonPos[i][0];
    IonSystem.R[i][1]                     = IonPos[i][1];
    IonSystem.R[i][2]                     = IonPos[i][2];
    GroupName[i]                          = IonName[atomic_number[i]];
    int speciesID                         = species.addSpecies(GroupName[i]);
    IonSystem.GroupID[i]                  = speciesID;
    species(AtomicNumberIndex, speciesID) = atomic_number[i];
    species(IonChargeIndex, speciesID)    = q[i];
  }
  hin.close();
}

void LCAOHDFParser::getSuperTwist(const std::string& fname)
{
  Matrix<double> MyVec(1, 3);
  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }

  hin.push("Super_Twist");

  STwist_Coord.resize(3);

  hin.read(MyVec, "Coord");

  hin.pop();
  STwist_Coord[0] = MyVec[0][0];
  STwist_Coord[1] = MyVec[0][1];
  STwist_Coord[2] = MyVec[0][2];

  hin.close();
}

void LCAOHDFParser::getMO(const std::string& fname)
{
  EigVal_alpha.resize(numMO);
  EigVal_beta.resize(numMO);
  EigVec.resize(2 * SizeOfBasisSet * numMO);

  Matrix<double> CartMat(numMO, SizeOfBasisSet);

  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }

  std::string setname = "/Super_Twist/eigenset_0";
  hin.read(CartMat, setname);
  setname = "/Super_Twist/eigenval_0";
  hin.read(EigVal_alpha, setname);
  std::copy(CartMat.begin(), CartMat.end(), EigVec.begin());

  if (!SpinRestricted)
  {
    setname = "/Super_Twist/eigenset_1";
    hin.read(CartMat, setname);
    setname = "/Super_Twist/eigenval_1";
    hin.read(EigVal_beta, setname);
  }
  std::copy(CartMat.begin(), CartMat.end(), EigVec.begin() + SizeOfBasisSet * numMO);

  hin.close();
  int btot = numMO * SizeOfBasisSet;
  int n = btot / 4, b = 0;
  int dn = btot - n * 4;
  std::ostringstream eig;

  eig.setf(std::ios::scientific, std::ios::floatfield);
  eig.setf(std::ios::right, std::ios::adjustfield);
  eig.precision(14);
  eig << "\n";
  for (int k = 0; k < n; k++)
  {
    eig << std::setw(22) << EigVec[b] << std::setw(22) << EigVec[b + 1] << std::setw(22) << EigVec[b + 2]
        << std::setw(22) << EigVec[b + 3] << "\n";
    b += 4;
  }
  for (int k = 0; k < dn; k++)
  {
    eig << std::setw(22) << EigVec[b++];
  }
  std::cout << eig.str() << std::endl;
  std::cout << "Finished reading MO." << std::endl;
  hin.close();
}

void LCAOHDFParser::getGaussianCenters(const std::string fname) {}

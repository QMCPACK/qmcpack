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

#include "QMCTools/CrystalAsciiParser.h"
#include "io/hdf_archive.h"

CrystalAsciiParser::CrystalAsciiParser()
{
  PBC            = true;
  ECP            = false;
  BohrUnit       = true;
  FixValence     = true;
  SpinRestricted = true;
  IsComplex      = false;
}

CrystalAsciiParser::CrystalAsciiParser(int argc, char** argv) : QMCGaussianParserBase(argc, argv)
{
  PBC            = true;
  ECP            = false;
  BohrUnit       = true;
  FixValence     = true;
  SpinRestricted = true;
  IsComplex      = false;
}

void CrystalAsciiParser::parse(const std::string& fname)
{
  std::string crystal_output_name = fname + ".out";
  std::ifstream fin(crystal_output_name.c_str());
  if (fin.fail())
  {
    std::cerr << "Error when opening file: " << crystal_output_name << std::endl;
    abort();
  }
  pivot_begin = fin.tellg();

  std::string aline;
  search(fin, "*******************************************************************************", aline);
  getwords(currentWords, fin);
  getwords(currentWords, fin);
  version = currentWords[1];

  search(fin, "CALCULATION", aline);
  parsewords(aline.c_str(), currentWords);
  if (currentWords[0] == "CRYSTAL")
  {
    PBC = true;
  }
  else if (currentWords[0] == "MOLECULAR")
  {
    PBC = false;
  }
  else
  {
    std::cerr << "Currently only converts CRYSTAL or MOLECULAR calculations" << std::endl;
    abort();
  }

  search(fin, "N. OF ATOMS PER CELL", aline);
  parsewords(aline.c_str(), currentWords);
  NumberOfAtoms = atoi(currentWords[5].c_str());
  search(fin, "NUMBER OF AO", aline);
  parsewords(aline.c_str(), currentWords);
  numAO = atoi(currentWords[3].c_str());
  search(fin, "N. OF ELECTRONS PER CELL", aline);
  parsewords(aline.c_str(), currentWords);
  NumberOfEls = atoi(currentWords[5].c_str());

  search(fin, "TYPE OF CALCULATION", aline);
  parsewords(aline.c_str(), currentWords);
  if (currentWords[4] == "RESTRICTED")
  {
    SpinRestricted = true;
    spin           = 0;
  }
  else
  {
    SpinRestricted = false;
    fin.seekg(pivot_begin);
    if (lookFor(fin, "ALPHA-BETA ELECTRONS LOCKED TO", aline))
    {
      parsewords(aline.c_str(), currentWords);
      spin = atoi(currentWords[4].c_str()); //2S
    }
    else
    {
      spin = 0;
    }
  }
  NumberOfAlpha = NumberOfEls / 2 + spin;
  NumberOfBeta  = NumberOfEls - NumberOfAlpha;

  SizeOfBasisSet = numAO;
  numMO          = numAO;


  BohrUnit = true; //will convert everything to Bohr

  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);

  if (PBC)
  {
    getCell(fin);
    getKpts(fin);
  }
  else
  {
    NbKpts = 1;
  }

  getGeometry(fin);

  getGaussianCenters(fin);

  getMO(fin);

  dumpHDF5(fname);
}

void CrystalAsciiParser::getCell(std::istream& is)
{
  X.resize(3);
  Y.resize(3);
  Z.resize(3);

  is.clear();
  is.seekg(pivot_begin);
  std::string aline;
  search(is, "DIRECT LATTICE VECTORS COMPON. (A.U.)", aline);
  getwords(currentWords, is);
  getwords(currentWords, is);
  X[0] = atof(currentWords[0].c_str());
  X[1] = atof(currentWords[1].c_str());
  X[2] = atof(currentWords[2].c_str());
  getwords(currentWords, is);
  Y[0] = atof(currentWords[0].c_str());
  Y[1] = atof(currentWords[1].c_str());
  Y[2] = atof(currentWords[2].c_str());
  getwords(currentWords, is);
  Z[0] = atof(currentWords[0].c_str());
  Z[1] = atof(currentWords[1].c_str());
  Z[2] = atof(currentWords[2].c_str());
  std::cout << "Lattice parameters in Bohr:" << std::endl;
  std::cout << X[0] << "  " << X[1] << "  " << X[2] << std::endl;
  std::cout << Y[0] << "  " << Y[1] << "  " << Y[2] << std::endl;
  std::cout << Z[0] << "  " << Z[1] << "  " << Z[2] << std::endl;
}

void CrystalAsciiParser::getKpts(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);
  std::string aline;
  search(is, "NUMBER OF K POINTS IN THE IBZ", aline);
  parsewords(aline.c_str(), currentWords);
  NbKpts = atoi(currentWords[12].c_str());
  search(is, "K POINTS COORDINATES (OBLIQUE COORDINATES IN UNITS OF IS", aline);
  std::string ignore = ")";
  parsewords(aline.c_str(), currentWords, ignore);
  int is_unit = atoi(currentWords[10].c_str());

  Kpoints_Coord.resize(NbKpts);
  for (int i = 0; i < NbKpts; i++)
  {
    getwords(currentWords, is, ignore);
    int nkpts;
    //Has at max 4 kpoints per line
    nkpts = currentWords.size() / 4;
    for (int j = 0; j < nkpts; j++)
    {
      Kpoints_Coord[i].resize(3);
      int a               = atoi(currentWords[4 * j + 1].c_str());
      int b               = atoi(currentWords[4 * j + 2].c_str());
      int c               = atoi(currentWords[4 * j + 3].c_str());
      Kpoints_Coord[i][0] = a / (double)is_unit;
      Kpoints_Coord[i][1] = b / (double)is_unit;
      Kpoints_Coord[i][2] = c / (double)is_unit;
      i++;
    }
    i--;
  }
}

void CrystalAsciiParser::getGeometry(std::istream& is)
{
  const double ang_to_bohr = 1.0 / 0.529177e0;
  tags.clear();
  is.clear();
  is.seekg(pivot_begin);
  std::string aline;
  search(is, "NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL", aline);
  parsewords(aline.c_str(), currentWords);
  int num_irr_ats = atoi(currentWords[8].c_str());
  search(is, "ATOM AT. N.", aline);

  for (int i = 0; i < num_irr_ats; i++)
  {
    getwords(currentWords, is);
    int at = atoi(currentWords[1].c_str());
    atzeff x(at, at);
    convAtNum.insert(std::pair<int, atzeff>(at, x));
  }

  //Find which of the unique atoms have ECPs to determine the Zeff and atomic number
  if (lookFor(is, "PSEUDOPOTENTIAL INFORMATION", aline))
  {
    ECP = true;
    for (int i = 0; i < num_irr_ats; i++)
    {
      if (lookFor(is, "ATOMIC NUMBER", aline))
      {
        parsewords(aline.c_str(), currentWords);
        int atnum                           = atoi(currentWords[2].c_str());
        int zeff                            = atoi(currentWords[5].c_str());
        convAtNum.at(atnum + 200).atomicNum = atnum;
        convAtNum.at(atnum + 200).zeff      = zeff;
      }
    }
  }
  else
  {
    ECP = false;
  }

  is.clear();
  is.seekg(pivot_begin);
  std::vector<int> atomic_number;
  std::vector<int> idx(NumberOfAtoms);
  Matrix<double> IonPos(NumberOfAtoms, 3);
  if (PBC)
  {
    search(is, "CARTESIAN COORDINATES - PRIMITIVE CELL", aline);
    getwords(currentWords, is);
    getwords(currentWords, is);
    getwords(currentWords, is);
    for (int i = 0; i < NumberOfAtoms; i++)
    {
      getwords(currentWords, is);
      idx[i]           = atoi(currentWords[1].c_str());
      std::string name = IonName[convAtNum.at(idx[i]).atomicNum];
      std::transform(name.begin(), name.end(), name.begin(), ::toupper);
      tags.push_back(name);
      for (int d = 0; d < 3; d++)
      {
        IonPos[i][d] = atof(currentWords[3 + d].c_str()) * ang_to_bohr;
      }
    }
  }
  else
  {
    search(is, "ATOMS IN THE ASYMMETRIC UNIT", aline);
    getwords(currentWords, is);
    getwords(currentWords, is);
    for (int i = 0; i < NumberOfAtoms; i++)
    {
      getwords(currentWords, is);
      idx[i]           = atoi(currentWords[2].c_str());
      std::string name = IonName[convAtNum.at(idx[i]).atomicNum];
      std::transform(name.begin(), name.end(), name.begin(), ::toupper);
      tags.push_back(name);
      for (int d = 0; d < 3; d++)
      {
        IonPos[i][d] = atof(currentWords[4 + d].c_str()) * ang_to_bohr;
      }
    }
  }

  for (int i = 0; i < NumberOfAtoms; i++)
  {
    AtomIndexmap.insert(std::pair<int, int>(i, idx[i]));
  }

  SpeciesSet& species(IonSystem.getSpeciesSet());
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      IonSystem.R[i][d] = IonPos[i][d];
    }
    GroupName[i]                          = IonName[convAtNum.at(AtomIndexmap.at(i)).atomicNum];
    int speciesID                         = species.addSpecies(GroupName[i]);
    IonSystem.GroupID[i]                  = speciesID;
    species(AtomicNumberIndex, speciesID) = convAtNum.at(AtomIndexmap.at(i)).atomicNum;
    species(IonChargeIndex, speciesID)    = convAtNum.at(AtomIndexmap.at(i)).zeff;
  }
}

void CrystalAsciiParser::scaleBasis()
{
  for (int i = 0; i < shID.size(); i++)
  {
    int count = 0;
    for (int j = 0; j < shID[i].size() - 1; j++)
    {
      int l       = str_to_l.at(shID[i][j]);
      int nfunc   = ncoeffpershell[i][j];
      double norm = 0.0;
      for (int a = 0; a < nfunc; a++)
      {
        for (int b = 0; b < nfunc; b++)
        {
          double tmp = 2.0 * std::sqrt(expo[i][count + a] * expo[i][count + b]);
          tmp /= expo[i][count + a] + expo[i][count + b];
          tmp = std::pow(tmp, l + 1.5);
          norm += coef[i][count + a] * coef[i][count + b] * tmp;
        }
      }
      norm = 1.0 / std::sqrt(norm);
      for (int a = 0; a < nfunc; a++)
      {
        coef[i][count + a] *= norm;
      }
      count += nfunc;
    }
  }
}


void CrystalAsciiParser::getGaussianCenters(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);
  gBound.resize(NumberOfAtoms + 1);
  int ng, nx;
  std::string aline;
  int nUniqAt = 0;
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    std::map<std::string, int>::iterator it(basisDataMap.find(tags[i]));
    if (it == basisDataMap.end())
    {
      basisDataMap[tags[i]] = nUniqAt++;
    }
  }
  NumberOfSpecies = basisDataMap.size();


  expo.resize(nUniqAt);
  coef.resize(nUniqAt);
  nshll.resize(nUniqAt, 0);
  ncoeffpershell.resize(nUniqAt);
  shID.resize(nUniqAt);

  gsMap[std::string("S")]  = 1;
  gsMap[std::string("SP")] = 2;
  gsMap[std::string("P")]  = 3;
  gsMap[std::string("D")]  = 4;
  gsMap[std::string("F")]  = 5;
  gsMap[std::string("G")]  = 6;

  str_to_l.insert(std::pair<std::string, int>("S", 0));
  str_to_l.insert(std::pair<std::string, int>("P", 1));
  str_to_l.insert(std::pair<std::string, int>("D", 2));
  str_to_l.insert(std::pair<std::string, int>("F", 3));
  str_to_l.insert(std::pair<std::string, int>("G", 4));


  search(is, "LOCAL ATOMIC FUNCTIONS BASIS SET", aline);
  getwords(currentWords, is);
  getwords(currentWords, is);
  getwords(currentWords, is);

  int currPos = -1;
  while (true)
  {
    getwords(currentWords, is);
    if ((currentWords[0] == "INFORMATION") ||
        (currentWords[0] == "*******************************************************************************"))
    {
      break;
    }
    if (currentWords.size() == 5) //found species
    {
      std::map<std::string, int>::iterator it(basisDataMap.find(currentWords[1]));
      if (it == basisDataMap.end())
      {
        std::cerr << "Error in parser.\n";
        abort();
      }
      currPos       = it->second;
      bool newgroup = (nshll[currPos] == 0);
      if (newgroup)
      {
        ncoeffpershell[currPos].clear();
        ncoeffpershell[currPos].push_back(0);
        shID[currPos].clear();
        shID[currPos].push_back("None");
      }
      std::string type;
      while (true)
      {
        std::streampos pivot = is.tellg();
        getwords(currentWords, is);
        if ((currentWords.size() == 5) || (currentWords[0] == "INFORMATION") ||
            (currentWords[0] == "*******************************************************************************") ||
            (currentWords[0] == "WARNING"))
        {
          is.seekg(pivot);
          break;
        }
        if (gsMap.find(currentWords[1]) != gsMap.end() || gsMap.find(currentWords[2]) != gsMap.end())
        {
          if (newgroup)
          {
            if (gsMap.find(currentWords[1]) != gsMap.end())
            {
              type = currentWords[1];
            }
            else
              type = currentWords[2];

            if (gsMap[type] == 2)
            {
              std::cerr << "Can't handle SP basis states yet. Fix latex. \n";
              abort();
            }
            if (gsMap.find(type) == gsMap.end())
            {
              std::cerr << "Unhandled primitive type " << type << std::endl;
              abort();
            }

            nshll[currPos]++;
            ncoeffpershell[currPos].push_back(0);
            shID[currPos].push_back("None");
          }
        }
        else
        {
          if (newgroup)
          {
            double expt, coeff1, coeff2, coeff3;
            if (currentWords.size() == 3)
            { //The primitive has a - sign between coefficients
              int j;
              for (int i = 0; i < 3; i++)
              {
                if (currentWords[i].length() == 19)
                {
                  j = i;
                }
              }
              std::string str1, str2;
              for (int i = 0; i < 19; i++)
              {
                if (i < 9)
                {
                  str1 += currentWords[j][i];
                }
                else
                {
                  str2 += currentWords[j][i];
                }
              }

              if (j == 0)
              {
                expt   = atof(str1.c_str());
                coeff1 = atof(str2.c_str());
                coeff2 = atof(currentWords[1].c_str());
                coeff3 = atof(currentWords[2].c_str());
              }
              else if (j == 1)
              {
                expt   = atof(currentWords[0].c_str());
                coeff1 = atof(str1.c_str());
                coeff2 = atof(str2.c_str());
                coeff3 = atof(currentWords[2].c_str());
              }
              else
              {
                expt   = atof(currentWords[0].c_str());
                coeff1 = atof(currentWords[1].c_str());
                coeff2 = atof(str1.c_str());
                coeff3 = atof(str2.c_str());
              }
            }
            else
            {
              expt   = atof(currentWords[0].c_str());
              coeff1 = atof(currentWords[1].c_str());
              coeff2 = atof(currentWords[2].c_str());
              coeff3 = atof(currentWords[3].c_str());
            }
            expo[currPos].push_back(expt);
            if (gsMap[type] == 1)
            {
              coef[currPos].push_back(coeff1);
            }
            else if (gsMap[type] == 3)
            {
              coef[currPos].push_back(coeff2);
            }
            else
            {
              coef[currPos].push_back(coeff3);
            }
            ncoeffpershell[currPos][nshll[currPos] - 1]++;
            shID[currPos][nshll[currPos] - 1] = type;
          }
        }
      }
    }
  }

  gShell.clear();
  gNumber.clear();
  gExp.clear();
  gC0.clear();
  gC1.clear();
  int gtot = 0;
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    std::map<std::string, int>::iterator it(basisDataMap.find(tags[i]));
    if (it == basisDataMap.end())
    {
      std::cerr << "Error in parser" << std::endl;
      abort();
    }
    gBound[i] = gtot;
    int indx  = it->second;
    gtot += nshll[indx];
    for (int k = 0; k < nshll[indx]; k++)
    {
      gShell.push_back(gsMap[shID[indx][k]]);
      gNumber.push_back(ncoeffpershell[indx][k]);
    }
    for (int k = 0; k < expo[indx].size(); k++)
    {
      gExp.push_back(expo[indx][k]);
      gC0.push_back(coef[indx][k]);
    }
  }
  gBound[NumberOfAtoms] = gtot;

  scaleBasis();
}

void CrystalAsciiParser::getKMO(std::istream& is, std::vector<std::vector<double>>& Mat)
{
  int currTot = Mat.size();
  getwords(currentWords, is);
  getwords(currentWords, is);

  std::vector<double> tmp;

  while (true)
  {
    getwords(currentWords, is);
    if (currentWords.size() == 0)
    {
      continue;
    }
    else if (Mat.size() == numMO)
    {
      break;
    }
    else
    {
      int nmo = currentWords.size();
      getwords(currentWords, is); //empty line
      for (int i = 0; i < nmo; i++)
      {
        Mat.push_back(tmp);
      }
      for (int i = 0; i < SizeOfBasisSet; i++)
      {
        getwords(currentWords, is);
        for (int j = 0; j < nmo; j++)
        {
          Mat[currTot + j].push_back(atof(currentWords[j + 1].c_str()));
        }
      }
      currTot += nmo;
    }
  }
}

void CrystalAsciiParser::getKMO(std::istream& is, std::vector<std::vector<std::complex<double>>>& CMat)
{
  int currTot = CMat.size();
  getwords(currentWords, is);
  getwords(currentWords, is);

  std::vector<std::complex<double>> tmp;

  while (true)
  {
    getwords(currentWords, is);
    if (currentWords.size() == 0)
    {
      continue;
    }
    else if (CMat.size() == numMO)
    {
      break;
    }
    else if (currentWords[0] == currentWords[1])
    { //actual complex kpoint
      int nmo = currentWords.size() / 2;
      getwords(currentWords, is);
      for (int i = 0; i < nmo; i++)
      {
        CMat.push_back(tmp);
      }
      for (int i = 0; i < SizeOfBasisSet; i++)
      {
        getwords(currentWords, is);
        for (int j = 0; j < nmo; j++)
        {
          std::complex<double> c(atof(currentWords[2 * j + 1].c_str()), atof(currentWords[2 * j + 2].c_str()));
          CMat[currTot + j].push_back(c);
        }
      }
      currTot += nmo;
    }
    else
    { //real kpoint, but storing in complex MO matrix
      int nmo = currentWords.size();
      getwords(currentWords, is);
      for (int i = 0; i < nmo; i++)
      {
        CMat.push_back(tmp);
      }
      for (int i = 0; i < SizeOfBasisSet; i++)
      {
        getwords(currentWords, is);
        for (int j = 0; j < nmo; j++)
        {
          std::complex<double> c(atof(currentWords[j + 1].c_str()), 0.0);
          CMat[currTot + j].push_back(c);
        }
      }
      currTot += nmo;
    }
  }
}

void CrystalAsciiParser::getMO(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);

  std::string aline;
  if (!lookFor(is, "HAMILTONIAN EIGENVALUES", aline))
  {
    std::cerr << "Calculation was not run with eigenvalue printing\n";
    std::cerr << "Rerun with the following input in the computational control section: \n";
    std::cerr << "PRINTOUT\n";
    std::cerr << "EIGENALL\n";
    std::cerr << "EIGENVEC\n";
    std::cerr << "-999\n";
    std::cerr << "END\n";
    abort();
  }
  is.clear();
  is.seekg(pivot_begin);
  std::streampos pivot;
  search(is, "CYC   0 ETOT(AU)", aline);
  while (!is.eof())
  {
    getwords(currentWords, is);
    if (currentWords[0] == "HAMILTONIAN" && currentWords[1] == "EIGENVALUES")
    {
      pivot = is.tellg();
    }
    else if (currentWords[0] == "SCF" && currentWords[1] == "ENDED")
    {
      break;
    }
  }
  if (is.eof())
  {
    std::cerr << "SCF Convergence not found\n";
    abort();
  }
  is.seekg(pivot);

  bool finished = false;
  while (!finished)
  {
    getwords(currentWords, is);
    if (currentWords.size() == 4 && currentWords[0].back() == '(')
    {
      std::vector<double> tmp;
      eigvals.push_back(tmp);
      while (true)
      {
        getwords(currentWords, is);
        if (currentWords.size() == 0)
        {
          break;
        }
        else if (currentWords[0] == "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
        {
          finished = true;
          break;
        }
        else
        {
          for (int i = 0; i < currentWords.size(); i++)
          {
            eigvals.back().push_back(atof(currentWords[i].c_str()));
          }
        }
      }
    }
  }

  if (eigvals.size() != NbKpts && SpinRestricted)
  {
    std::cerr << "Did not read correct number of eigenvalues" << std::endl;
    std::cerr << "#Eigvals: " << eigvals.size() << "  NbKpts: " << NbKpts << std::endl;
    abort();
  }
  else if (eigvals.size() / 2 != NbKpts && !SpinRestricted)
  {
    std::cerr << "Did not read correct number of eigenvalues" << std::endl;
    std::cerr << "#Eigvals Alpha+Beta: " << eigvals.size() << "  NbKpts: " << NbKpts << std::endl;
    abort();
  }

  is.clear();
  is.seekg(pivot_begin);
  std::vector<std::streampos> pivots;
  IsComplex = false;
  if (!lookFor(is, "FINAL EIGENVECTORS", aline))
  {
    std::cerr << "Calculation was not run with eigenvector printing\n";
    std::cerr << "Rerun with the following input in the computational control section: \n";
    std::cerr << "PRINTOUT\n";
    std::cerr << "EIGENALL\n";
    std::cerr << "EIGENVEC\n";
    std::cerr << "-999\n";
    std::cerr << "END\n";
    abort();
  }
  getwords(currentWords, is);
  while (true)
  {
    getwords(currentWords, is);
    if (currentWords[0] == "EIGENVECTORS" && currentWords[1] == "IN" && currentWords[2] == "FORTRAN" &&
        currentWords[3] == "UNIT")
    {
      break;
    }
    else if (currentWords.size() == 5 && currentWords[1] == "(")
    {
      pivots.push_back(is.tellg());
      getwords(currentWords, is);
      getwords(currentWords, is);

      getwords(currentWords, is);
      if (currentWords[0] == currentWords[1])
      {
        IsComplex = true;
      }
    }
  }

  if (pivots.size() != NbKpts && SpinRestricted)
  {
    std::cerr << "Wrong number of Kpoints for SpinRestricted calculations\n";
    abort();
  }
  else if (pivots.size() / 2 != NbKpts && !SpinRestricted)
  {
    std::cerr << "Wrong number of Kpoints for SpinUnrestricted calculation\n";
    abort();
  }

  for (int k = 0; k < pivots.size(); k++)
  {
    is.seekg(pivots[k]);
    std::vector<std::vector<double>> Mat;
    std::vector<std::vector<std::complex<double>>> CMat;
    if (IsComplex)
    {
      getKMO(is, CMat);
      complex_kmos.push_back(CMat);
    }
    else
    {
      getKMO(is, Mat);
      real_kmos.push_back(Mat);
    }
  }
}

void CrystalAsciiParser::dumpHDF5(const std::string& fname)
{
  for (int k = 0; k < NbKpts; k++)
  {
    STwist_Coord.resize(3);
    STwist_Coord = Kpoints_Coord[k];
    std::cout << "Creating ESHDF file for SuperTwist " << k << ": " << Kpoints_Coord[k][0] << " " << Kpoints_Coord[k][1] << " "
              << Kpoints_Coord[k][2] << std::endl;
    h5file = fname + "_" + std::to_string(k) + ".h5";
    hdf_archive hout;
    hout.create(h5file.c_str(), H5F_ACC_TRUNC);

    std::string str;

    hout.push("application", true);
    str = "crystal";
    hout.write(str, "code");
    hout.write(version, "version");
    hout.pop();


    hout.push("PBC", true);
    hout.write(PBC, "PBC");
    hout.pop();

    hout.push("atoms", true);
    hout.write(NumberOfAtoms, "number_of_atoms");
    hout.write(NumberOfSpecies, "number_of_species");

    Matrix<double> Pos(NumberOfAtoms, 3);
    for (int i = 0; i < NumberOfAtoms; i++)
    {
      for (int d = 0; d < 3; d++)
      {
        Pos[i][d] = IonSystem.R[i][d];
      }
    }
    hout.write(Pos, "positions");
    for (int i = 0; i < NumberOfSpecies; i++)
    {
      int j = 0;
      while (j < NumberOfAtoms)
      {
        if (i != IonSystem.GroupID[j])
        {
          j++;
        }
        else
        {
          str = "species_" + std::to_string(i);
          hout.push(str, true);
          int at, zeff, core;
          at   = convAtNum.at(AtomIndexmap.at(j)).atomicNum;
          zeff = convAtNum.at(AtomIndexmap.at(j)).zeff;
          core = at - zeff;
          hout.write(at, "atomic_number");
          hout.write(zeff, "charge");
          hout.write(core, "core");
          hout.write(GroupName[j], "name");
          hout.pop();
          break;
        }
      }
    }
    std::vector<int> ids(NumberOfAtoms);
    for (int i = 0; i < NumberOfAtoms; i++)
    {
      ids[i] = IonSystem.GroupID[i];
    }
    hout.write(ids, "species_ids");
    hout.pop();

    if (PBC)
    {
      hout.push("Cell", true);
      Matrix<double> cell(3, 3);
      for (int d = 0; d < 3; d++)
      {
        cell[0][d] = X[d];
        cell[1][d] = Y[d];
        cell[2][d] = Z[d];
      }
      hout.write(cell, "LatticeVectors");
      hout.pop();
    }

    hout.push("parameters", true);
    hout.write(ECP, "ECP");
    hout.write(IsComplex, "IsComplex");
    hout.write(NumberOfAlpha, "NbAlpha");
    hout.write(NumberOfBeta, "NbBeta");
    hout.write(NumberOfEls, "NbTotElec");
    hout.write(spin, "spin");
    hout.write(SpinRestricted, "SpinRestricted");
    hout.write(BohrUnit, "Unit");
    if (SpinRestricted)
    {
      hout.write(numMO, "numMO");
      hout.write(numAO, "numAO");
    }
    else
    {
      hout.write(numMO, "numMO_up");
      hout.write(numMO, "numMO_dn");
      hout.write(numAO, "numAO_up");
      hout.write(numAO, "numAO_dn");
    }
    hout.pop();

    hout.push("basisset", true);
    int NbElements = shID.size();
    hout.write(NbElements, "NbElements");
    str = "LCAOBSet";
    hout.write(str, "name");

    for (int i = 0; i < shID.size(); i++)
    {
      str = "atomicBasisSet" + std::to_string(i);
      hout.push(str, true);
      int NbBasisGroups = shID[i].size() - 1;
      hout.write(NbBasisGroups, "NbBasisGroups");
      str = "spherical";
      hout.write(str, "angular");
      str = "gaussian";
      hout.write(str, "expandYlm");
      int grid_npts  = 1001;
      double grid_rf = 100.0;
      double grid_ri = 1e-06;
      hout.write(grid_npts, "grid_npts");
      hout.write(grid_rf, "grid_rf");
      hout.write(grid_ri, "grid_ri");
      str = "log";
      hout.write(str, "grid_type");
      str = "name";
      hout.write(str, "name");
      str = "no";
      hout.write(str, "normalized");

      for (std::map<std::string, int>::iterator it = basisDataMap.begin(); it != basisDataMap.end(); it++)
      {
        if (it->second == i)
        {
          str = it->first;
        }
      }
      std::string element = str;
      if (element.size() == 2)
      {
        element[1] = std::tolower(element[1]);
      }
      hout.write(element, "elementType");

      int count = 0;
      for (int j = 0; j < NbBasisGroups; j++)
      {
        str = "basisGroup" + std::to_string(j);
        hout.push(str, true);
        int NbRadFunc = ncoeffpershell[i][j];
        hout.write(NbRadFunc, "NbRadFunc");
        std::vector<int> shell(3, 0);
        hout.write(shell, "Shell_coord");
        int l = str_to_l.at(shID[i][j]);
        hout.write(l, "l");
        hout.write(j, "n");
        hout.push("radfunctions", true);
        for (int b = 0; b < NbRadFunc; b++)
        {
          str = "DataRad" + std::to_string(b);
          hout.push(str, true);
          hout.write(coef[i][count], "contraction");
          hout.write(expo[i][count], "exponent");
          hout.pop(); //end DataRad
          count++;
        }
        hout.pop(); //end radfunctions
        str = element + std::to_string(j) + std::to_string(l);
        hout.write(str, "rid");
        str = "Gaussian";
        hout.write(str, "type");
        hout.pop(); //end  basisGroup
      }
      hout.pop(); //end atomicBasisSet
    }
    hout.pop(); //end basisset

    hout.pop();

    str = "Super_Twist";
    hout.push(str, true);
    if (PBC)
    {
      Matrix<double> Tw(1,3);
      Tw[0][0] = Kpoints_Coord[k][0];
      Tw[0][1] = Kpoints_Coord[k][1];
      Tw[0][2] = Kpoints_Coord[k][2];
      hout.write(Tw, "Coord");
    }
    if (IsComplex)
    {
      Matrix<double> real(SizeOfBasisSet, SizeOfBasisSet);
      Matrix<double> imaj(SizeOfBasisSet, SizeOfBasisSet);
      for (int i = 0; i < SizeOfBasisSet; i++)
      {
        for (int j = 0; j < SizeOfBasisSet; j++)
        {
          real[i][j] = std::real(complex_kmos[k][i][j]);
          imaj[i][j] = std::imag(complex_kmos[k][i][j]);
        }
      }
      hout.write(real, "eigenset_0");
      hout.write(imaj, "eigenset_0_imag");
      if (!SpinRestricted)
      {
        for (int i = 0; i < SizeOfBasisSet; i++)
        {
          for (int j = 0; j < SizeOfBasisSet; j++)
          {
            real[i][j] = std::real(complex_kmos[k + NbKpts][i][j]);
            imaj[i][j] = std::imag(complex_kmos[k + NbKpts][i][j]);
          }
        }
        hout.write(real, "eigenset_1");
        hout.write(imaj, "eigenset_1_imag");
      }
    }
    else
    {
      Matrix<double> real(SizeOfBasisSet, SizeOfBasisSet);
      for (int i = 0; i < SizeOfBasisSet; i++)
      {
        for (int j = 0; j < SizeOfBasisSet; j++)
        {
          real[i][j] = real_kmos[k][i][j];
        }
      }
      hout.write(real, "eigenset_0");
      if (!SpinRestricted)
      {
        for (int i = 0; i < SizeOfBasisSet; i++)
        {
          for (int j = 0; j < SizeOfBasisSet; j++)
          {
            real[i][j] = real_kmos[k + NbKpts][i][j];
          }
        }
        hout.write(real, "eigenset_1");
      }
    }
    if (PBC)
    {
      hout.write(eigvals[k], "eigenval_0");
      if (!SpinRestricted)
      {
        hout.write(eigvals[k + NbKpts], "eigenval_1");
      }
      hout.pop();
    }


    hout.close();
  }
}

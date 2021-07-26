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

#include "QMCTools/DiracParser.h"
#include "io/hdf/hdf_archive.h"
#include <cstdio>

typedef std::pair<double, double> primExpCoeff;
typedef std::vector<primExpCoeff> basisFunc;

DiracParser::DiracParser(int argc, char** argv) : QMCGaussianParserBase(argc, argv)
{
  ECP      = false;
  BohrUnit = true;
  PBC      = false;
  normMap  = {
      {"s", 1.0},
      {"px", 1.0},
      {"py", 1.0},
      {"pz", 1.0},
      {"dxx", std::sqrt(3)},
      {"dyy", std::sqrt(3)},
      {"dzz", std::sqrt(3)},
      {"dxy", 1.0},
      {"dxz", 1.0},
      {"dyz", 1.0},
      {"fxxx", std::sqrt(15)},
      {"fyyy", std::sqrt(15)},
      {"fzzz", std::sqrt(15)},
      {"fxxy", std::sqrt(3)},
      {"fxxz", std::sqrt(3)},
      {"fxyy", std::sqrt(3)},
      {"fyyz", std::sqrt(3)},
      {"fxzz", std::sqrt(3)},
      {"fyzz", std::sqrt(3)},
      {"fxyz", 1},
      {"g400", std::sqrt(105)},
      {"g040", std::sqrt(105)},
      {"g004", std::sqrt(105)},
      {"g310", std::sqrt(15)},
      {"g301", std::sqrt(15)},
      {"g130", std::sqrt(15)},
      {"g031", std::sqrt(15)},
      {"g103", std::sqrt(15)},
      {"g013", std::sqrt(15)},
      {"g220", std::sqrt(9)},
      {"g202", std::sqrt(9)},
      {"g022", std::sqrt(9)},
      {"g211", std::sqrt(3)},
      {"g121", std::sqrt(3)},
      {"g112", std::sqrt(3)},
      {"h500", std::sqrt(945)},
      {"h050", std::sqrt(945)},
      {"h005", std::sqrt(945)},
      {"h410", std::sqrt(105)},
      {"h401", std::sqrt(105)},
      {"h140", std::sqrt(105)},
      {"h041", std::sqrt(105)},
      {"h104", std::sqrt(105)},
      {"h014", std::sqrt(105)},
      {"h320", std::sqrt(45)},
      {"h302", std::sqrt(45)},
      {"h230", std::sqrt(45)},
      {"h032", std::sqrt(45)},
      {"h203", std::sqrt(45)},
      {"h023", std::sqrt(45)},
      {"h311", std::sqrt(15)},
      {"h131", std::sqrt(15)},
      {"h113", std::sqrt(15)},
      {"h221", std::sqrt(9)},
      {"h212", std::sqrt(9)},
      {"h122", std::sqrt(9)},
      {"i600", std::sqrt(10395)},
      {"i060", std::sqrt(10395)},
      {"i006", std::sqrt(10395)},
      {"i510", std::sqrt(945)},
      {"i501", std::sqrt(945)},
      {"i150", std::sqrt(945)},
      {"i051", std::sqrt(945)},
      {"i105", std::sqrt(945)},
      {"i015", std::sqrt(945)},
      {"i420", std::sqrt(315)},
      {"i402", std::sqrt(315)},
      {"i240", std::sqrt(315)},
      {"i042", std::sqrt(315)},
      {"i204", std::sqrt(315)},
      {"i024", std::sqrt(315)},
      {"i411", std::sqrt(105)},
      {"i141", std::sqrt(105)},
      {"i114", std::sqrt(105)},
      {"i330", std::sqrt(225)},
      {"i303", std::sqrt(225)},
      {"i033", std::sqrt(225)},
      {"i321", std::sqrt(45)},
      {"i312", std::sqrt(45)},
      {"i231", std::sqrt(45)},
      {"i132", std::sqrt(45)},
      {"i213", std::sqrt(45)},
      {"i123", std::sqrt(45)},
      {"i222", std::sqrt(27)},
  };
}

void DiracParser::parse(const std::string& fname)
{
  std::string dirac_out = fname + ".out";
  std::ifstream fin(dirac_out.c_str());
  if (fin.fail())
  {
    std::cerr << "Error when opening file: " << dirac_out << std::endl;
    abort();
  }
  pivot_begin = fin.tellg();

  search(fin, "Release DIRAC", aline);
  parsewords(aline.c_str(), currentWords);
  version = std::stoi(currentWords[2].erase(0, 5));

  search(fin, "Number of atom types", aline);
  parsewords(aline.c_str(), currentWords);
  NumberOfSpecies = std::atoi(currentWords.back().c_str());

  search(fin, "Total number of atoms", aline);
  parsewords(aline.c_str(), currentWords);
  NumberOfAtoms = std::atoi(currentWords.back().c_str());

  std::cout << "Found " << NumberOfSpecies << " unique species" << std::endl;
  std::cout << "Found " << NumberOfAtoms << " total number of atoms" << std::endl;

  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  getGeometry(fin);
  getGaussianCenters(fin);
  getSpinors(fin);
  dumpHDF5(fname);
}

void DiracParser::getGeometry(std::istream& is)
{
  const double ang_to_bohr = 1.0 / 0.529177e0;
  is.clear();
  is.seekg(pivot_begin);

  search(is, "label    atoms   charge   prim    cont     basis", aline);
  skiplines(is, 1);
  std::map<std::string, int> zMap;
  std::map<std::string, int> zeffMap;
  for (int uat = 0; uat < NumberOfSpecies; uat++)
  {
    getwords(currentWords, is);
    std::string name = currentWords[0];
    int zeff         = std::stoi(currentWords[2]);
    for (int i = 0; i < IonName.size(); i++)
    {
      if (IonName[i] == name)
      {
        std::pair<std::string, int> x(name, i);
        std::pair<std::string, int> y(name, zeff);
        zMap.insert(x);
        zeffMap.insert(y);
        break;
      }
    }
    skiplines(is, 1);
  }
  search(is, "total: ", aline);
  parsewords(aline.c_str(), currentWords);
  numAO = std::stoi(currentWords[4]);

  search(is, "Cartesian coordinates in XYZ format", aline);
  skiplines(is, 4);
  SpeciesSet& species(IonSystem.getSpeciesSet());
  for (int iat = 0; iat < NumberOfAtoms; iat++)
  {
    getwords(currentWords, is);
    for (int d = 0; d < 3; d++)
    {
      IonSystem.R[iat][d] = std::stod(currentWords[1 + d]) * ang_to_bohr;
    }
    GroupName[iat]                        = currentWords[0];
    int speciesID                         = species.addSpecies(GroupName[iat]);
    IonSystem.GroupID[iat]                = speciesID;
    species(AtomicNumberIndex, speciesID) = zMap.find(GroupName[iat])->second;
    species(IonChargeIndex, speciesID)    = zeffMap.find(GroupName[iat])->second;
  }
}

void DiracParser::getGaussianCenters(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);

  search(is, "Contents of the molecule file", aline);
  skiplines(is, 5);
  getwords(currentWords, is);
  int nat = std::stoi(currentWords[1]);
  assert(nat == NumberOfSpecies);
  basisset.resize(NumberOfSpecies);
  SpeciesSet& species(IonSystem.getSpeciesSet());
  for (int iat = 0; iat < NumberOfSpecies; iat++)
  {
    basisset[iat].elementType = species.speciesName[iat];
    std::vector<int> numPerL;
    search(is, "LARGE", aline);
    parsewords(aline.c_str(), currentWords);
    if (currentWords[1] != "EXPLICIT")
    {
      std::cerr << "Cannot extract basis. Rerun with EXPLICIT basis" << std::endl;
      abort();
    }
    for (int n = 3; n < currentWords.size(); n++)
    {
      numPerL.push_back(std::stoi(currentWords[n]));
    }
    //loop through each angular momentum
    for (int l = 0; l < numPerL.size(); l++)
    {
      for (int il = 0; il < numPerL[l]; il++)
      {
        search(is, "f ", aline);
        parsewords(aline.c_str(), currentWords);
        int nprim = std::stoi(currentWords[1]);
        int ncont = std::stoi(currentWords[2]);
        //each basis in this current block is uncontracted
        if (ncont == 0)
        {
          for (int n = 0; n < nprim; n++)
          {
            getwords(currentWords, is);
            primBasis p;
            p.first  = std::stod(currentWords[0]);
            p.second = 1.0;
            basisGroup bas;
            bas.n   = basisset[iat].basisGroups.size();
            bas.l   = l;
            bas.rid = basisset[iat].elementType + std::to_string(bas.n) + std::to_string(bas.l);
            bas.radfuncs.push_back(p);
            basisset[iat].basisGroups.push_back(bas);
          }
        }
        else
        {
          int bgstart = basisset[iat].basisGroups.size();
          for (int n = 0; n < ncont; n++)
          {
            basisGroup bg;
            bg.n   = basisset[iat].basisGroups.size();
            bg.l   = l;
            bg.rid = basisset[iat].elementType + std::to_string(bg.n) + std::to_string(bg.l);
            basisset[iat].basisGroups.push_back(bg);
          }
          for (int n = 0; n < nprim; n++)
          {
            getwords(currentWords, is);
            assert(currentWords.size() == ncont + 1);
            for (int c = 0; c < ncont; c++)
            {
              primBasis p;
              p.first  = std::stod(currentWords[0]);
              p.second = std::stod(currentWords[c + 1]);
              basisset[iat].basisGroups[bgstart + c].radfuncs.push_back(p);
            }
          }
        }
      }
    }
  }
}

void DiracParser::getSpinors(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);
  search(is, "Coefficients from DFCOEF", aline);
  std::streampos startspinors = is.tellg();
  search(is, "************************", aline);
  std::streampos endspinors = is.tellg();

  is.clear();
  is.seekg(startspinors);


  while (lookFor(is, "Electronic eigenvalue no.", aline))
  {
    skiplines(is, 1);
    std::vector<std::complex<double>> up(numAO, 0), dn(numAO, 0);
    std::vector<std::complex<double>> upkp(numAO, 0), dnkp(numAO, 0);
    if (is.tellg() > endspinors)
      break;
    while (std::getline(is,aline))
    {
      if (aline.size() == 0)
        break;

      char sp[4];
      char info[12];
      double up_r,up_i,dn_r,dn_i;
      int bidx;
      std::sscanf(aline.c_str(), "%3c%5d%2c%12c%2c%lf%lf%lf%lf",sp,&bidx,sp,info,sp,&up_r,&up_i,&dn_r,&dn_i);
      bidx -= 1; //count from 0
      assert(bidx >= 0);
      parsewords(info, currentWords);

      double norm              = 1.0;
      std::string label = currentWords[3];
      normMapType::iterator it = normMap.find(label);
      if (it != normMap.end())
        norm = it->second;
      else
      {
        std::cerr << "Unknown basis function type. Aborting" << std::endl;
        abort();
      }

      std::complex<double> a = {up_r, up_i};
      std::complex<double> b = {dn_r, dn_i};
      std::complex<double> c = {-dn_r, dn_i};
      std::complex<double> d = {up_r, -up_i};
      up[bidx] = a * norm;
      dn[bidx] = b * norm;
      upkp[bidx] = c * norm;
      dnkp[bidx] = d * norm;

    }

    upcoeff.push_back(up);
    upcoeff.push_back(upkp);
    dncoeff.push_back(dn);
    dncoeff.push_back(dnkp);

  }

  std::cout << "Found " << upcoeff.size() << " spinors" << std::endl;

}

void DiracParser::dumpHDF5(const std::string& fname)
{
  std::string str;

  h5file = fname + ".orbs.h5";
  hdf_archive hout;
  hout.create(h5file.c_str(), H5F_ACC_TRUNC);

  hout.push("application", true);
  str = "dirac";
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
  SpeciesSet& species(IonSystem.getSpeciesSet());
  for (int isp = 0; isp < species.size(); isp++)
  {
    str = "species_" + std::to_string(isp);
    hout.push(str, true);
    int at, core, zeff;
    at   = species(AtomicNumberIndex, isp);
    zeff = species(IonChargeIndex, isp);
    core = at - zeff;
    hout.write(at, "atomic_number");
    hout.write(zeff, "charge");
    hout.write(core, "core");
    hout.write(species.speciesName[isp], "name");
    hout.pop();
  }
  std::vector<int> ids(NumberOfAtoms);
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    ids[i] = IonSystem.GroupID[i];
  }
  hout.write(ids, "species_ids");
  hout.pop();
  hout.pop();

  str = "basisset";
  hout.push(str, true);
  hout.write(NumberOfSpecies, "NbElements");
  str = "LCAOBSet";
  hout.write(str, "name");
  for (int nb = 0; nb < basisset.size(); nb++)
  {
    str = "atomicBasisSet" + std::to_string(nb);
    hout.push(str, true);
    atBasisSet& bs = basisset[nb];
    hout.write(bs.name, "name");
    hout.write(bs.normalized, "normalized");
    hout.write(bs.grid_type, "grid_type");
    hout.write(bs.expandYlm, "expandYlm");
    hout.write(bs.elementType, "elementType");
    hout.write(bs.angular, "angular");
    hout.write(bs.grid_ri, "grid_ri");
    hout.write(bs.grid_rf, "grid_rf");
    hout.write(bs.grid_npts, "grid_npts");
    int nbg = bs.basisGroups.size();
    hout.write(nbg, "NbBasisGroups");
    for (int ng = 0; ng < bs.basisGroups.size(); ng++)
    {
      str = "basisGroup" + std::to_string(ng);
      hout.push(str, true);
      basisGroup& bg = bs.basisGroups[ng];
      int nrf        = bg.radfuncs.size();
      hout.write(nrf, "NbRadFunc");
      hout.write(bg.rid, "rid");
      hout.write(bg.l, "l");
      hout.write(bg.n, "n");
      hout.write(bg.type, "type");
      hout.push("radfunctions", true);
      for (int np = 0; np < bg.radfuncs.size(); np++)
      {
        str = "DataRad" + std::to_string(np);
        hout.push(str, true);
        primBasis& p = bg.radfuncs[np];
        hout.write(p.first, "exponent");
        hout.write(p.second, "contraction");
        hout.pop();
      }
      hout.pop();

      hout.pop();
    }
    hout.pop();
  }
  hout.pop();

  str = "Nb_KPTS";
  hout.push(str, true);
  int nk = 1;
  hout.write(nk, "Nbkpts");
  hout.pop();

  numMO = upcoeff.size();
  Matrix<double> up_real(numMO, numAO);
  Matrix<double> up_imag(numMO, numAO);
  Matrix<double> dn_real(numMO, numAO);
  Matrix<double> dn_imag(numMO, numAO);
  for (int mo = 0; mo < numMO; mo++)
  {
    for (int ao = 0; ao < numAO; ao++)
    {
      up_real[mo][ao] = upcoeff[mo][ao].real();
      up_imag[mo][ao] = upcoeff[mo][ao].imag();
      dn_real[mo][ao] = dncoeff[mo][ao].real();
      dn_imag[mo][ao] = dncoeff[mo][ao].imag();
    }
  }
  str = "KPTS_0";
  hout.push(str, true);
  hout.write(up_real, "eigenset_0");
  hout.write(up_imag, "eigenset_0_imag");
  hout.write(dn_real, "eigenset_1");
  hout.write(dn_imag, "eigenset_1_imag");
  hout.pop();
  hout.close();
}

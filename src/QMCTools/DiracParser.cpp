//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCTools/DiracParser.h"
#include "io/hdf/hdf_archive.h"
#include <cstdio>
#include <algorithm>

using primExpCoeff = std::pair<double, double>;
using basisFunc    = std::vector<primExpCoeff>;

fermIrrep::fermIrrep(std::string label_in, int nSpinors, int numAO) : label(label_in)
{
  energies.resize(nSpinors, 0);
  orbtypes.resize(nSpinors, OrbType::VIRTUAL);
  spinor_mo_coeffs.resize(nSpinors);
  for (int mo = 0; mo < nSpinors; mo++)
  {
    spinor_mo_coeffs[mo].resize(numAO);
    for (int ao = 0; ao < numAO; ao++)
      spinor_mo_coeffs[mo][ao].resize(4, 0);
  }
}

fermIrrep fermIrrep::generate_kramers_pair()
{
  std::string new_label = label;
  assert(label[1] == '1');
  new_label[1] = '2';
  int nspinor  = get_num_spinors();
  int numAO    = get_num_ao();
  fermIrrep kp(new_label, nspinor, numAO);
  kp.energies = energies;
  for (int mo = 0; mo < nspinor; mo++)
  {
    for (int ao = 0; ao < numAO; ao++)
    {
      kp.spinor_mo_coeffs[mo][ao][0] = -spinor_mo_coeffs[mo][ao][2];
      kp.spinor_mo_coeffs[mo][ao][1] = spinor_mo_coeffs[mo][ao][3];
      kp.spinor_mo_coeffs[mo][ao][2] = spinor_mo_coeffs[mo][ao][0];
      kp.spinor_mo_coeffs[mo][ao][3] = -spinor_mo_coeffs[mo][ao][1];
    }
  }
  return kp;
}

cosciRep::cosciRep(std::string in_label, int nstates) : label(in_label) { states.resize(nstates); }

void cosciRep::printInfo(std::ostream& os, int& tot_state_count)
{
  os << "  Representation: " << label << " with " << states.size() << " states" << std::endl;
  os << "state#     Energies and Ndets: " << std::endl;
  for (int i = 0; i < states.size(); i++)
  {
    os << "      " << tot_state_count << " " << states[i].energy << " " << states[i].coeffs.size() << std::endl;
    tot_state_count++;
  }
}

DiracParser::DiracParser(int argc, char** argv) : QMCGaussianParserBase(argc, argv)
{
  ECP          = false;
  BohrUnit     = true;
  PBC          = false;
  isSpinor     = true;
  angular_type = "cartesian";
  expandYlm    = "Dirac";
  normMap      = {
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
  std::string dirac_out = fname;
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

  search(fin, "*SCF: Set-up for");
  skiplines(fin, 2);
  getwords(currentWords, fin);
  if (currentWords[1] == "Closed")
  {
    NumberOfBeta  = 0;
    NumberOfAlpha = std::stoi(currentWords[6]);
    NumberOfEls   = NumberOfAlpha;
  }
  else if (currentWords[1] == "Open")
  {
    NumberOfBeta = 0;
    search(fin, "* Shell specifications:");
    lookFor(fin, "Total", aline);
    parsewords(aline.c_str(), currentWords);
    NumberOfAlpha = std::stoi(currentWords[1]);
    NumberOfEls   = NumberOfAlpha;
  }
  else
  {
    std::cerr << "Error number of electrons" << std::endl;
    abort();
  }

  IonSystem.create({NumberOfAtoms});
  GroupName.resize(NumberOfAtoms);
  getGeometry(fin);
  getGaussianCenters(fin);
  getSpinors(fin);
  getWF(fin);
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
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
  SizeOfBasisSet = std::stoi(currentWords[4]);

  //now overwrite charge for pseudsized atoms
  is.clear();
  is.seekg(pivot_begin);
  while (lookFor(is, "Nuclear Gaussian exponent for atom", aline))
  {
    parsewords(aline.c_str(), currentWords);
    int z = std::stoi(currentWords[7]);
    std::getline(is, aline);
    if (aline.size() == 0)
      break;
    //found an ECP to replace
    ECP = true;
    getwords(currentWords, is);
    int zeff                         = std::stoi(currentWords[6]);
    zeffMap.find(IonName[z])->second = zeff;
  }
  is.clear();
  is.seekg(pivot_begin);

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
            bas.n = basisset[iat].basisGroups.size();
            bas.l = l;
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
            bg.n = basisset[iat].basisGroups.size();
            bg.l = l;
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

  //Store data to QMCGaussian output format
  std::map<std::string, int> basisDataMap;
  int nUniqAt = 0;
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    std::map<std::string, int>::iterator it(basisDataMap.find(GroupName[i]));
    if (it == basisDataMap.end())
    {
      basisDataMap[GroupName[i]] = nUniqAt++;
    }
  }
  gBound.resize(NumberOfAtoms + 1);
  gShell.clear();
  gNumber.clear();
  gExp.clear();
  gC0.clear();
  gC1.clear();
  int gtot = 0;
  for (int i = 0; i < NumberOfAtoms; i++)
  {
    std::map<std::string, int>::iterator it(basisDataMap.find(GroupName[i]));
    if (it == basisDataMap.end())
    {
      std::cerr << "Error in parser.\n";
      abort();
    }
    gBound[i] = gtot;
    int indx  = it->second;
    gtot += basisset[indx].basisGroups.size();
    for (int k = 0; k < basisset[indx].basisGroups.size(); k++)
    {
      int l  = basisset[indx].basisGroups[k].l;
      int sh = (l == 0) ? 1 : l + 2;
      gShell.push_back(sh);
    }
    for (int k = 0; k < basisset[indx].basisGroups.size(); k++)
      gNumber.push_back(basisset[indx].basisGroups[k].radfuncs.size());
    for (int k = 0; k < basisset[indx].basisGroups.size(); k++)
    {
      for (int c = 0; c < basisset[indx].basisGroups[k].radfuncs.size(); c++)
      {
        gExp.push_back(basisset[indx].basisGroups[k].radfuncs[c].first);
        gC0.push_back(basisset[indx].basisGroups[k].radfuncs[c].second);
      }
    }
  }
  gBound[NumberOfAtoms] = gtot;
}

void DiracParser::getSpinors(std::istream& is)
{
  std::cout << std::endl;
  std::cout << "Reading spinor info" << std::endl;
  std::cout << "========================================================================" << std::endl;
  is.clear();
  is.seekg(pivot_begin);

  search(is, "Output from DBLGRP");
  skiplines(is, 2);
  getwords(currentWords, is);
  //look for "* Nirrep fermion irreps: irrep1 irrep2 ..."
  int nirrep = currentWords.size() - 4;
  std::vector<std::string> labels(nirrep);
  for (int i = 0; i < nirrep; i++)
    labels[i] = currentWords[4 + i];

  for (int i = 0; i < nirrep; i++)
  {
    is.clear();
    is.seekg(pivot_begin);
    search(is, "VECINP: Vector print");
    std::string irstr = "- Orbitals in fermion ircop " + labels[i];
    lookFor(is, irstr, aline);
    parsewords(aline.c_str(), currentWords);
    if (currentWords.back() == ":1..oo")
    {
      search(is, "* Occupation of subblocks");
      lookFor(is, labels[i], aline);
      lookFor(is, "tot.num. of pos.erg shells:", aline);
      parsewords(aline.c_str(), currentWords);
      int nj       = currentWords.size() - 4;
      int nspinors = 0;
      for (int j = 0; j < nj; j++)
        nspinors += std::stoi(currentWords[4 + j]);
      irreps.push_back(fermIrrep(labels[i], nspinors, SizeOfBasisSet));
    }
    else
    {
      std::string splitby = ":.";
      aline               = currentWords.back();
      parsewords(aline.c_str(), currentWords, splitby);
      if (currentWords[0] != "1")
      {
        std::cerr << "Error: For vector printing, need to use 1..oo or 1..N" << std::endl;
        abort();
      }
      int nspinors = std::stoi(currentWords[1]);
      irreps.push_back(fermIrrep(labels[i], nspinors, SizeOfBasisSet));
    }
  }

  std::cout << "Found " << nirrep << " fermion irreps." << std::endl;
  for (int i = 0; i < nirrep; i++)
    std::cout << "  irrep " << irreps[i].get_label() << " with " << irreps[i].get_num_spinors() << " spinors and "
              << irreps[i].get_num_ao() << " AO coefficients." << std::endl;

  search(is, "Coefficients from DFCOEF", aline);
  std::streampos startspinors = is.tellg();

  for (int i = 0; i < nirrep; i++)
  {
    is.clear();
    is.seekg(startspinors);
    std::string start_irrep = "Fermion ircop " + irreps[i].get_label();
    search(is, start_irrep);
    for (int mo = 0; mo < irreps[i].get_num_spinors(); mo++)
    {
      lookFor(is, "Electronic eigenvalue no.", aline);
      parsewords(aline.c_str(), currentWords);
      irreps[i].energies[mo] = std::stod(currentWords.back());
      skiplines(is, 1);
      while (std::getline(is, aline))
      {
        if (aline.size() == 0)
          break;
        if (std::string(aline).find("*********") != std::string::npos)
        {
          std::cerr << "ERROR parsing line: " << std::endl;
          std::cerr << aline << std::endl;
          std::cerr << "One of the printed AO coefficients is outside the default DIRAC print format" << std::endl;
          std::cerr << "In order to continue, please change the following in DIRAC/src/dirac/dirout.F (around line 427)"
                    << std::endl;
          std::cerr << "     100  FORMAT(3X,I5,2X,A12,2X,4F14.10)" << std::endl;
          std::cerr << " to " << std::endl;
          std::cerr << "     100  FORMAT(3X,I5,2X,A12,2X,4F20.10)" << std::endl;
          std::cerr
              << " and recompile DIRAC. Then rerun your particular calulcation in order to get accurate AO coefficients"
              << std::endl;
          abort();
        }
        parsewords(aline.c_str(), currentWords);
        if (currentWords.size() != 9)
        {
          std::cerr << "ERROR parsing line: " << std::endl;
          std::cerr << aline << std::endl;
          std::cerr << "Expected line to be parsed into vector<string> of length 9" << std::endl;
          std::cerr
              << "Either recompile DIRAC in DIRAC/src/dirac/dirout.F (around line 427) to avoid this issue from now on"
              << std::endl;
          std::cerr << "     100  FORMAT(3X,I5,2X,A12,2X,4F14.10)" << std::endl;
          std::cerr << " to " << std::endl;
          std::cerr << "     100  FORMAT(3X,I5,2X,A12,2X,4F20.10)" << std::endl;
          std::cerr << " or just add a space appropriately to this line" << std::endl;
          abort();
        }
        int bidx = std::stoi(currentWords[0]) - 1;

        double norm              = 1.0;
        std::string label        = currentWords[4];
        normMapType::iterator it = normMap.find(label);
        if (it != normMap.end())
          norm = it->second;
        else
        {
          std::cerr << "Unknown basis function type. Aborting" << std::endl;
          abort();
        }

        double up_r = std::stod(currentWords[5]);
        double up_i = std::stod(currentWords[6]);
        double dn_r = std::stod(currentWords[7]);
        double dn_i = std::stod(currentWords[8]);

        irreps[i].spinor_mo_coeffs[mo][bidx][0] = up_r * norm;
        irreps[i].spinor_mo_coeffs[mo][bidx][1] = up_i * norm;
        irreps[i].spinor_mo_coeffs[mo][bidx][2] = dn_r * norm;
        irreps[i].spinor_mo_coeffs[mo][bidx][3] = dn_i * norm;
      }
    }

    std::cout << "Found coefficients for " << irreps[i].get_label() << std::endl;
    kp_irreps.push_back(irreps[i].generate_kramers_pair());
    std::cout << "Generated kramers pair with irrep " << kp_irreps[i].get_label() << std::endl;
  }

  std::cout << "Now we have the following spinors" << std::endl;
  for (int i = 0; i < irreps.size(); i++)
  {
    std::cout << "  irrep " << irreps[i].get_label() << " with " << irreps[i].get_num_spinors() << " spinors and "
              << irreps[i].get_num_ao() << " AO coefficients." << std::endl;
    std::cout << "  irrep " << kp_irreps[i].get_label() << " with " << kp_irreps[i].get_num_spinors() << " spinors and "
              << kp_irreps[i].get_num_ao() << " AO coefficients." << std::endl;
  }

  numMO = 0;
  for (int i = 0; i < irreps.size(); i++)
    numMO += 2 * irreps[i].get_num_spinors(); //irrep and kp
}

void DiracParser::getWF(std::istream& is)
{
  is.clear();
  is.seekg(pivot_begin);
  std::cout << std::endl;
  std::cout << "Parsing wave function info" << std::endl;
  std::cout << "========================================================================" << std::endl;

  if (lookFor(is, "Resolution of open-shell states", aline))
  {
    bool closed = lookFor(is, "No open-shell electrons", aline);
    if (closed)
    {
      std::cout << "Found single determinant wave function" << std::endl;
      getSingleDet(is);
    }
    else
    {
      std::cout << "Found Complete Open-Shell CI (COSCI) wave function" << std::endl;
      getCOSCI(is);
    }
  }
  else
  {
    std::cout << "Found single determinant wave function" << std::endl;
    getSingleDet(is);
  }
}

void DiracParser::getSingleDet(std::istream& is)
{
  //energy order spinors across fermion irreps
  std::cout << "Sorting spinors by energy" << std::endl;
  std::vector<std::pair<double, std::pair<int, int>>> idx;
  for (int ir = 0; ir < irreps.size(); ir++)
  {
    for (int mo = 0; mo < irreps[ir].get_num_spinors(); mo++)
    {
      std::pair<int, int> irmo(ir, mo);
      std::pair<double, std::pair<int, int>> enirmo(irreps[ir].energies[mo], irmo);
      idx.push_back(enirmo);
    }
  }
  std::sort(idx.begin(), idx.end());

  //store in QMCGaussian data format
  //storing in EigVec by up_real, dn_real, up_imag, dn_imag
  //only sort energy for irrep, since kramers pairs kp_irreps are degenerate
  EigVec.clear();
  EigVal_alpha.clear();
  EigVal_beta.clear();
  std::vector<int> spinor_component = {0, 2, 1, 3};
  for (int d = 0; d < 4; d++)
  {
    for (int i = 0; i < idx.size(); i++)
    {
      if (d == 0)
      {
        //store it twice, first spinor and its kramers pair
        EigVal_alpha.push_back(idx[i].first);
        EigVal_alpha.push_back(idx[i].first);
      }
      int ir = idx[i].second.first;
      int mo = idx[i].second.second;
      for (int ao = 0; ao < irreps[ir].get_num_ao(); ao++)
        EigVec.push_back(irreps[ir].spinor_mo_coeffs[mo][ao][spinor_component[d]]);
      for (int ao = 0; ao < irreps[ir].get_num_ao(); ao++)
        EigVec.push_back(kp_irreps[ir].spinor_mo_coeffs[mo][ao][spinor_component[d]]);
    }
  }
}

void DiracParser::parseCOSCIOrbInfo(std::istream& is, const int irrep_idx, OrbType type)
{
  std::string orbtype_str;
  if (type == OrbType::CORE)
    orbtype_str = "Core";
  else if (type == OrbType::ACTIVE)
    orbtype_str = "Active";
  else
  {
    std::cerr << "Orb type must be CORE or ACTIVE" << std::endl;
    abort();
  }

  is.seekg(pivot_begin);
  search(is, "Resolution of open-shell states");
  std::string tmp = "- " + orbtype_str + " orbitals";
  search(is, tmp);
  lookFor(is, irreps[irrep_idx].label, aline);
  parsewords(aline.c_str(), currentWords);
  assert(currentWords.back() == irreps[irrep_idx].label);
  skiplines(is, 1);
  getwords(currentWords, is);
  if (currentWords[0] == "Index")
  {
    int norb  = std::stoi(currentWords[2]);
    int count = 0;
    while (count < norb)
    {
      getline(is, aline);
      if (aline.size() == 0)
        continue;
      parsewords(aline.c_str(), currentWords);
      for (int i = 0; i < currentWords.size(); i++, count++)
      {
        int idx                            = std::stoi(currentWords[i]) - 1;
        irreps[irrep_idx].orbtypes[idx]    = type;
        kp_irreps[irrep_idx].orbtypes[idx] = type;
      }
      if (count > norb)
      {
        std::cerr << "Error: Read in more indices than exist for irrep " << irreps[irrep_idx].label << std::endl;
        abort();
      }
      else if (count == norb)
        break;
    }
  }
}

int DiracParser::sortAndStoreCOSCIOrbs(OrbType type, const int spinor_component)
{
  int total = 0;
  for (int ir = 0; ir < irreps.size(); ir++)
  {
    std::vector<std::pair<double, std::pair<int, int>>> idx;
    for (int mo = 0; mo < irreps[ir].get_num_spinors(); mo++)
    {
      if (irreps[ir].orbtypes[mo] == type)
      {
        std::pair<int, int> irmo(ir, mo);
        std::pair<double, std::pair<int, int>> enirmo(irreps[ir].energies[mo], irmo);
        idx.push_back(enirmo);
      }
    }

    for (int i = 0; i < idx.size(); i++)
    {
      if (spinor_component == 0)
      {
        EigVal_alpha.push_back(idx[i].first);
        total++;
      }
      int ir = idx[i].second.first;
      int mo = idx[i].second.second;
      for (int ao = 0; ao < irreps[ir].get_num_ao(); ao++)
        EigVec.push_back(irreps[ir].spinor_mo_coeffs[mo][ao][spinor_component]);
    }
    //now store KP for this irrep
    for (int i = 0; i < idx.size(); i++)
    {
      if (spinor_component == 0)
      {
        EigVal_alpha.push_back(idx[i].first);
        total++;
      }
      int ir = idx[i].second.first;
      int mo = idx[i].second.second;
      for (int ao = 0; ao < irreps[ir].get_num_ao(); ao++)
        EigVec.push_back(kp_irreps[ir].spinor_mo_coeffs[mo][ao][spinor_component]);
    }
  }

  return total;
}

void DiracParser::getCOSCI(std::istream& is)
{
  is.clear();
  multideterminant = true;

  //find closed info. Info is initialized as virtual.
  //Just need to update core and active
  for (int i = 0; i < irreps.size(); i++)
  {
    parseCOSCIOrbInfo(is, i, OrbType::CORE);
    parseCOSCIOrbInfo(is, i, OrbType::ACTIVE);
  }

  std::cout << std::endl;
  std::cout << "Orbital Info" << std::endl;
  std::cout << "------------------------------------" << std::endl;
  for (int i = 0; i < irreps.size(); i++)
  {
    std::cout << "irrep: " << irreps[i].label << std::endl;
    int closed = 0;
    int active = 0;
    int virt   = 0;
    for (int j = 0; j < irreps[i].orbtypes.size(); j++)
    {
      if (irreps[i].orbtypes[j] == OrbType::CORE)
        closed += 1;
      else if (irreps[i].orbtypes[j] == OrbType::ACTIVE)
        active += 1;
      else if (irreps[i].orbtypes[j] == OrbType::VIRTUAL)
        virt += 1;
    }
    std::cout << "  closed  : " << closed << std::endl;
    std::cout << "  active  : " << active << std::endl;
    std::cout << "  virtual : " << virt << std::endl;
    std::cout << "  total   : " << closed + active + virt << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Sorting spinors into DIRAC COSCI order" << std::endl;

  //store in QMCGaussian data format
  //storing in EigVec by up_real, dn_real, up_imag, dn_imag
  EigVec.clear();
  EigVal_alpha.clear();
  EigVal_beta.clear();
  //save real part of up component of spinor in EigVec first
  int total_core    = sortAndStoreCOSCIOrbs(OrbType::CORE, 0);
  int total_active  = sortAndStoreCOSCIOrbs(OrbType::ACTIVE, 0);
  int total_virtual = sortAndStoreCOSCIOrbs(OrbType::VIRTUAL, 0);
  //save real part of dn compeonent of spinor in Eigvec
  sortAndStoreCOSCIOrbs(OrbType::CORE, 2);
  sortAndStoreCOSCIOrbs(OrbType::ACTIVE, 2);
  sortAndStoreCOSCIOrbs(OrbType::VIRTUAL, 2);
  //save imag part of up compeonent of spinor in Eigvec
  sortAndStoreCOSCIOrbs(OrbType::CORE, 1);
  sortAndStoreCOSCIOrbs(OrbType::ACTIVE, 1);
  sortAndStoreCOSCIOrbs(OrbType::VIRTUAL, 1);
  //save imag part of dn compeonent of spinor in Eigvec
  sortAndStoreCOSCIOrbs(OrbType::CORE, 3);
  sortAndStoreCOSCIOrbs(OrbType::ACTIVE, 3);
  sortAndStoreCOSCIOrbs(OrbType::VIRTUAL, 3);


  //set occstrs for core and virtual
  std::string core_occstr;
  for (int i = 0; i < total_core; i++)
    core_occstr += "1";
  std::string virt_occstr;
  for (int i = 0; i < total_virtual; i++)
    virt_occstr += "0";
  //active occstr found from COSCI states below

  is.seekg(pivot_begin);
  search(is, "Resolution of open-shell states");
  search(is, "Orbital Representation");

  std::vector<std::streampos> rep_pos;
  while (lookFor(is, "Representation", aline))
    rep_pos.push_back(is.tellg());
  is.clear();

  for (int i = 0; i < rep_pos.size(); i++)
  {
    is.seekg(rep_pos[i]);
    lookFor(is, "Population analysis", aline);
    parsewords(aline.c_str(), currentWords);
    std::string label = currentWords.back();
    getwords(currentWords, is);
    int nstates = std::stoi(currentWords[2]);
    cosciRep rep(label, nstates);
    skiplines(is, 3);
    for (int j = 0; j < nstates; j++)
    {
      std::vector<double> ci_coeffs;
      std::vector<std::string> ci_occs;
      double energy = 0.0;
      while (getline(is, aline))
      {
        if (aline.size() == 0)
          break;
        parsewords(aline.c_str(), currentWords);
        if (currentWords.size() == 1)
          energy = std::stod(currentWords[0]);
        if (currentWords[0] == "Sum")
        {
          getline(is, aline);
          break;
        }
        if (currentWords.size() == 5 && currentWords[0] != "Sum")
        {
          ci_coeffs.push_back(std::stod(currentWords[3]));
          std::string tmp = core_occstr + currentWords[1] + virt_occstr;
          ci_occs.push_back(tmp);
        }
      }
      //finished reading CI coeffs for state
      ciState ci;
      ci.energy     = energy;
      ci.occstrings = ci_occs;
      ci.coeffs     = ci_coeffs;
      rep.states[j] = ci;
    }
    cosciReps.push_back(rep);
  }

  std::cout << std::endl;
  std::cout << "COSCI State Info" << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "Found " << cosciReps.size() << " representations" << std::endl;
  int state_count = 0;
  for (int i = 0; i < cosciReps.size(); i++)
    cosciReps[i].printInfo(std::cout, state_count);
  std::cout << "Saving wave function for target state " << target_state << std::endl;
  std::cout << "note: if you want another state run with --TargetState #_of_desired_state shown above" << std::endl;

  //store info for desired state in QMCGausianParserBase structures
  bool found  = false;
  state_count = 0;
  for (int i = 0; i < cosciReps.size(); i++)
  {
    for (int j = 0; j < cosciReps[i].states.size(); j++)
    {
      if (state_count == target_state)
      {
        CIcoeff    = cosciReps[i].states[j].coeffs;
        CIalpha    = cosciReps[i].states[j].occstrings;
        CIbeta     = CIalpha; //just store it. It isn't written to h5
        found      = true;
        ci_nstates = CIalpha[0].size();
        ci_size    = CIcoeff.size();
        ci_nca     = 0;
        ci_nea     = NumberOfEls;
      }
      state_count++;
    }
  }
  if (!found)
  {
    std::cerr << "Could not find requested state" << std::endl;
    abort();
  }
}

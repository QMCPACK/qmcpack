//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//
// File created by: Kevin Gasperich, kgasperich@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "RMGParser.h"
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>


RMGParser::RMGParser() { PBC = true; }

RMGParser::RMGParser(int argc, char** argv) : QMCGaussianParserBase(argc, argv) { PBC = true; }

void RMGParser::getCell(const std::string& fname)
{
  // from LCAOHDFParser (maybe put into base class?)
  X.resize(3);
  Y.resize(3);
  Z.resize(3);

  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }
  hin.push("supercell");
  Matrix<double> LatticeVec(3, 3);
  hin.read(LatticeVec, "primitive_vectors");
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

void RMGParser::dumpPBC(const std::string& psi_tag, const std::string& ion_tag)
{
  std::cout << " RMGParser::dumpPBC " << std::endl;
  if (!Structure)
  {
    xmlDocPtr doc_p      = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr qm_root_p = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    if (PBC)
      xmlAddChild(qm_root_p, createCell());
    xmlAddChild(qm_root_p, createIonSet());
    xmlAddChild(qm_root_p, createElectronSet(ion_tag));
    xmlDocSetRootElement(doc_p, qm_root_p);
    std::string fname = Title + ".structure.xml";
    xmlSaveFormatFile(fname.c_str(), doc_p, 1);
    xmlFreeDoc(doc_p);
    Structure = true;
  }
  xmlDocPtr doc      = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  {
    //wavefunction
    xmlNodePtr wfPtr = xmlNewNode(NULL, (const xmlChar*)"wavefunction");
    xmlNewProp(wfPtr, (const xmlChar*)"name", (const xmlChar*)psi_tag.c_str());
    xmlNewProp(wfPtr, (const xmlChar*)"target", (const xmlChar*)"e");
    {
      xmlNodePtr detsetPtr = xmlNewNode(NULL, (const xmlChar*)"determinantset");
      xmlNewProp(detsetPtr, (const xmlChar*)"type", (const xmlChar*)"bspline");
      xmlNewProp(detsetPtr, (const xmlChar*)"href", (const xmlChar*)h5file.c_str());
      xmlNewProp(detsetPtr, (const xmlChar*)"source", (const xmlChar*)ion_tag.c_str());
      // this just reproduces the behavior of pw2qmcpack
      xmlNewProp(detsetPtr, (const xmlChar*)"sort", (const xmlChar*)"1");
      xmlNewProp(detsetPtr, (const xmlChar*)"tilematrix", (const xmlChar*)"1 0 0 0 1 0 0 0 1");
      xmlNewProp(detsetPtr, (const xmlChar*)"twistnum", (const xmlChar*)"0");
      xmlNewProp(detsetPtr, (const xmlChar*)"version", (const xmlChar*)"0.10");
      {
        std::ostringstream up_size, down_size;
        up_size << NumberOfAlpha;
        down_size << NumberOfBeta;

        xmlNodePtr slaterdetPtr = xmlNewNode(NULL, (const xmlChar*)"slaterdeterminant");

        xmlNodePtr updetPtr = xmlNewNode(NULL, (const xmlChar*)"determinant");
        xmlNewProp(updetPtr, (const xmlChar*)"id", (const xmlChar*)"updet");
        xmlNewProp(updetPtr, (const xmlChar*)"size", (const xmlChar*)up_size.str().c_str());

        xmlNodePtr occUpPtr = xmlNewNode(NULL, (const xmlChar*)"occupation");
        xmlNewProp(occUpPtr, (const xmlChar*)"mode", (const xmlChar*)"ground");
        xmlNewProp(occUpPtr, (const xmlChar*)"spindataset", (const xmlChar*)"0");
        xmlAddChild(updetPtr, occUpPtr);


        xmlNodePtr downdetPtr = xmlNewNode(NULL, (const xmlChar*)"determinant");
        xmlNewProp(downdetPtr, (const xmlChar*)"id", (const xmlChar*)"downdet");
        xmlNewProp(downdetPtr, (const xmlChar*)"size", (const xmlChar*)down_size.str().c_str());

        xmlNodePtr occDownPtr = xmlNewNode(NULL, (const xmlChar*)"occupation");
        xmlNewProp(occDownPtr, (const xmlChar*)"mode", (const xmlChar*)"ground");
        if (NumberOfSpins == 2)
        {
          xmlNewProp(downdetPtr, (const xmlChar*)"ref", (const xmlChar*)"updet");
          xmlNewProp(occDownPtr, (const xmlChar*)"spindataset", (const xmlChar*)"1");
        }
        else if (NumberOfSpins == 1)
        {
          xmlNewProp(occDownPtr, (const xmlChar*)"spindataset", (const xmlChar*)"0");
        }
        else
        {
          std::cerr << "Error: Number of spins should be 1 or 2. (" << NumberOfSpins << ")" << std::endl;
          abort();
        }
        xmlAddChild(downdetPtr, occDownPtr);

        xmlAddChild(slaterdetPtr, updetPtr);
        xmlAddChild(slaterdetPtr, downdetPtr);

        xmlAddChild(detsetPtr, slaterdetPtr);
      }
      xmlAddChild(wfPtr, detsetPtr);
      if (addJastrow)
      {
        std::cout << "Adding Two-Body and One-Body jastrows with rcut=\"10\" and size=\"10\"" << std::endl;
        if (NumberOfEls > 1)
        {
          xmlAddChild(wfPtr, createJ2());
        }
        xmlAddChild(wfPtr, createJ1());
        if (NumberOfEls > 1)
        {
          std::cout << "Adding Three-Body jastrows with rcut=\"5\"" << std::endl;
          xmlAddChild(wfPtr, createJ3());
        }
      }
    }
    xmlAddChild(qm_root, wfPtr);
  }
  xmlDocSetRootElement(doc, qm_root);
  std::string fname = Title + ".wf" + WFS_name + ".xml";
  xmlSaveFormatFile(fname.c_str(), doc, 1);
  xmlFreeDoc(doc);
}

void RMGParser::parse(const std::string& fname)
{
  hdf_archive hin;

  if (!hin.open(fname.c_str(), H5F_ACC_RDONLY))
  {
    std::cerr << "Could not open H5 file" << std::endl;
    abort();
  }

  hin.push("application");
  hin.read(CodeName, "code");
  std::cout << "Converting Wavefunction from the " << CodeName << " Code" << std::endl;
  hin.pop();

  hin.push("electrons");
  std::vector<double> Nalpha_Nbeta(2);
  hin.read(Nalpha_Nbeta, "number_of_electrons");
  NumberOfAlpha = Nalpha_Nbeta[0];
  NumberOfBeta  = Nalpha_Nbeta[1];
  NumberOfEls   = NumberOfAlpha + NumberOfBeta;
  std::cout << "Number of (alpha,beta) electrons: (" << NumberOfAlpha << ", " << NumberOfBeta << ")" << std::endl;

  hin.read(NumberOfSpins, "number_of_spins");
  std::cout << "Number of spins: " << NumberOfSpins << std::endl;
  hin.pop();

  hin.push("atoms");
  hin.read(NumberOfAtoms, "number_of_atoms");
  std::cout << "Number of atoms: " << NumberOfAtoms << std::endl;
  getCell(fname);

  IonSystem.create(NumberOfAtoms);
  GroupName.resize(NumberOfAtoms);
  Matrix<double> IonPos(NumberOfAtoms, 3);
  //atomic numbers
  std::vector<int> atomic_number;
  std::vector<double> q;
  ECP_names.clear();
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
    std::string ecpName = "";
    hin.read(zint, "atomic_number");
    atomic_number.push_back(zint);
    hin.read(z, "valence_charge");
    q.push_back(z);
    hin.read(Name, "name");
    try
    {
      hin.read(ecpName, "pseudopotential");
      ECP_names.push_back(ecpName);
      ECP = true;
    }
    catch (...)
    {
      std::cerr << "WARNING: no ECP found for " << speciesName << ":" << Name << std::endl;
      ECP_names.push_back(ecpName);
    }
    hin.pop();
  }
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
  hin.pop(); // atoms


  // this is just a messy workaround to avoid modifying the behavior of convert4qmc
  // if (PBC) then convert4qmc will try to print the first 3 elements of STwist_Coord
  STwist_Coord.resize(3);
  try
  {
    hin.push("Super_Twist");
    Matrix<double> STVec(1, 3);
    hin.read(STVec, "Coord");
    hin.pop();
    STwist_Coord[0] = STVec[0][0];
    STwist_Coord[1] = STVec[0][1];
    STwist_Coord[2] = STVec[0][2];
  }
  catch (...)
  {
    std::cout << "Could not find Super_Twist, using [0,0,0] (not yet implemented in basic RMG interface)" << std::endl;
    STwist_Coord[0] = 0;
    STwist_Coord[1] = 0;
    STwist_Coord[2] = 0;
  }
}

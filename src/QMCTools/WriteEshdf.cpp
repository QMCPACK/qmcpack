/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////

#include "WriteEshdf.h"
#include "XmlRep.h"
#include "FftContainer.h"
#include "HdfHelpers.h"
#include <sstream>
#include <map>
using namespace std;
using namespace hdfhelper;

class KPoint 
{
public:
  double kx;
  double ky;
  double kz;
  KPoint() { kx = 0; ky = 0; kz = 0; }
  KPoint(double x, double y, double z) { kx = x; ky = y; kz = z; }
  KPoint(const KPoint& kp) { kx = kp.kx; ky = kp.ky; kz = kp.kz; }
  KPoint& operator=(const KPoint& kp) { kx = kp.kx; ky = kp.ky; kz = kp.kz; return *this; }
  bool operator==(const KPoint& kp) const 
  {
    if ((std::abs(kx-kp.kx) < 1e-7) && (std::abs(ky-kp.ky) < 1e-7) && (std::abs(kz-kp.kz) < 1e-7)) 
    {
      return true;
    }
    return false;
  }
  bool operator<(const KPoint& kp) const 
  {
    if (abs(kx-kp.kx) > 1e-7) 
    {
      return kx < kp.kx;
    } 
    else if (abs(ky-kp.ky) > 1e-7) 
    {
      return ky < kp.ky;
    } 
    else if (abs(kz-kp.kz) > 1e-7) 
    {
      return kz < kp.kz;
    }
    return false;
  }
};


EshdfFile::EshdfFile(const string& hdfFileName) 
{
  file_ = H5Fcreate(hdfFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

EshdfFile::~EshdfFile() 
{
  H5Fclose(file_);
}


int EshdfFile::wrapped(int i, int size) const 
{
  if (i < size/2) 
  {
    return i;
  } 
  else 
  {
    return wrapped(i-size, size);
  }
}

double EshdfFile::getOccupation(const XmlNode* nd) const 
{
  double result = 0.0;
  vector<double> occupancies;
  const XmlNode& densMat = nd->getChild("density_matrix");
  densMat.getValue(occupancies);
  for (int i = 0; i < occupancies.size(); i++) 
  {
    result += occupancies[i];
  }
  return result;
}

void EshdfFile::handleSpinGroup(const XmlNode* nd, hid_t groupLoc, double& nocc, FftContainer& cont) 
{
  nocc = getOccupation(nd);
  int stateCounter = 0;
  vector<double> eigvals;

  for (int chIdx = 0; chIdx < nd->getNumChildren(); chIdx++) 
  {
    if (nd->getChild(chIdx).getName() == "grid_function") 
    {
      //cout << "Working on state " << stateCounter << endl;
      stringstream statess;
      statess << "state_" << stateCounter;
      hid_t state_group = makeHDFGroup(statess.str(), groupLoc);
      
      // HACK_HACK_HACK!!!
      eigvals.push_back(-5000.0+stateCounter);

      const XmlNode& eigFcnNode = nd->getChild(chIdx);
      readInEigFcn(eigFcnNode, cont);
      // write eigfcn to proper place
      hsize_t psig_dims[]={static_cast<hsize_t>(cont.fullSize),2};

      vector<double> temp;
      for (int i = 0; i < cont.fullSize; i++) 
      {
	temp.push_back(cont.kspace[i][0]);
	temp.push_back(cont.kspace[i][1]);
      }

      writeNumsToHDF("psi_g", temp, state_group, 2, psig_dims);
      stateCounter++;
    }
  }
  writeNumsToHDF("number_of_states", stateCounter, groupLoc);
  writeNumsToHDF("eigenvalues", eigvals, groupLoc);
}

void EshdfFile::readInEigFcn(const XmlNode& nd, FftContainer& cont) 
{
  const string type = nd.getAttribute("type");
  const string encoding = nd.getAttribute("encoding");
  
  if (encoding != "text") 
  {
    cerr << "Don't yet know how to handle encoding of wavefunction values other than text" << endl;
    exit(1);
  }
  
  vector<double> values;
  nd.getValue(values);
  const double fixnorm = 1/std::sqrt(static_cast<double>(cont.fullSize));
  if (type == "complex") 
  {
    int index = 0;
    for (int ix = 0; ix < cont.getNx(); ix++) 
    {
      for (int iy = 0; iy < cont.getNy(); iy++) 
      {
	for (int iz = 0; iz < cont.getNz(); iz++) 
        {
	  const int qbx = cont.getQboxIndex(ix,iy,iz);
	  cont.rspace[index][0] = values[2*qbx]*fixnorm;
	  cont.rspace[index][1] = values[2*qbx+1]*fixnorm;
	  index++;
	}
      }
    }
  } 
  else if (type == "double") 
  {
    int index = 0;
    for (int ix = 0; ix < cont.getNx(); ix++) 
    {
      for (int iy = 0; iy < cont.getNy(); iy++) 
      {
	for (int iz = 0; iz < cont.getNz(); iz++) 
        {
	  const int qbx = cont.getQboxIndex(ix,iy,iz);
	  cont.rspace[index][0] = values[qbx] * fixnorm;
	  cont.rspace[index][1] = 0.0;
	  index++;
	}
      }
    }
  }
  //cout << "in readInEigFcn, before fft, real space L2 norm = " << cont.getL2NormRS() << endl;
  cont.executeFFT();
  cont.fixKsNorm(fixnorm);
  //cout << "in readInEigFcn, after fft, k space L2 norm = " << cont.getL2NormKS() << endl;
}

void EshdfFile::writeApplication(const string& appName, int major, int minor, int sub) 
{
  vector<int> version{major, minor, sub};

  hid_t h_app = makeHDFGroup("application", file_);
  { 
    writeStringToHDF("code", appName, h_app);
    writeNumsToHDF("version", version, h_app);
  }
}


void EshdfFile::writeVersion() 
{
  vector<int> version{2,1,0};
  writeNumsToHDF("version", version, file_);
}

void EshdfFile::writeCreator() 
{
  string nameStr("creator");
  vector<int> version{0,1,0};
  hid_t h_creator = makeHDFGroup("creator", file_);
  {
    writeStringToHDF("program_name", "qboxConverter", h_creator);
    writeNumsToHDF("version", version, h_creator);
  }
}

void EshdfFile::writeFormat() 
{
  writeStringToHDF("format", "ES-HDF", file_);
}

void EshdfFile::writeBoilerPlate(const string& appName, const XmlNode& qboxSample) 
{
  const XmlNode& description = qboxSample.getChild("description");
  string desString = description.getValue();
  int versionStart = desString.find("qbox-");
  versionStart += 5;
  string versionStr = desString.substr(versionStart);
  
  const int firstDotIdx = versionStr.find_first_of('.');
  const int secondDotIdx = versionStr.find_last_of('.');
  const int major = stoi(versionStr.substr(0,firstDotIdx));
  const int minor = stoi(versionStr.substr(firstDotIdx+1,secondDotIdx-firstDotIdx-1));
  const int sub = stoi(versionStr.substr(secondDotIdx+1));

  writeApplication(appName, major, minor, sub);
  writeVersion();
  writeCreator();
  writeFormat();
}

void EshdfFile::writeSupercell(const XmlNode& qboxSample) 
{
  // grab the primitive translation vectors from the atomset tag's attributes and put the entries in the vector ptvs
  const XmlNode& atomset = qboxSample.getChild("atomset");
  const XmlNode& unit_cell = atomset.getChild("unit_cell");

  stringstream ss;
  ss << unit_cell.getAttribute("a") << "  " << unit_cell.getAttribute("b") << "  " << unit_cell.getAttribute("c");

  vector<double> ptvs;
  double temp;
  while (ss >> temp) 
  {
    ptvs.push_back(temp);
  }

  // write the ptvs to the supercell group of the hdf file
  hid_t supercell_group = makeHDFGroup("supercell", file_);
  hsize_t dims[]={3,3};
  writeNumsToHDF("primitive_vectors", ptvs, supercell_group, 2, dims);  
}


void EshdfFile::writeAtoms(const XmlNode& qboxSample) 
{
  const XmlNode& atomset = qboxSample.getChild("atomset");

  //make group
  hid_t atoms_group = makeHDFGroup("atoms", file_);
  
  map<string, int> SpeciesNameToInt;  
  //go through each species, extract: 
  //   atomic_number, mass, name, pseudopotential and valence_charge
  //   then write to a file, also set up mapping between species name and number
  int speciesNum = 0;
  for (int i = 0; i < atomset.getNumChildren(); i++) 
  {
    if (atomset.getChild(i).getName() == "species") 
    {
      const XmlNode& species = atomset.getChild(i);

      string spName = species.getAttribute("name");
      SpeciesNameToInt[spName] = speciesNum;

      int atomic_number;
      species.getChild("atomic_number").getValue(atomic_number);
      double mass;
      species.getChild("mass").getValue(mass);
      string name = species.getChild("symbol").getValue();
      int val_charge;
      species.getChild("norm_conserving_pseudopotential").getChild("valence_charge").getValue(val_charge);

      stringstream gname;
      gname << "species_" << speciesNum;
      hid_t species_group = makeHDFGroup(gname.str(), atoms_group);
      writeNumsToHDF("atomic_number", atomic_number, species_group);
      writeNumsToHDF("mass", mass, species_group);
      writeNumsToHDF("valence_charge", val_charge, species_group);
      writeStringToHDF("name", name, species_group);
      writeStringToHDF("pseudopotential", "unknown", species_group);
      speciesNum++;
    }
  }
  writeNumsToHDF("number_of_species", speciesNum, atoms_group);
  
  // go through atoms and extract their position and type 
  std::vector<int> species_ids;
  std::vector<double> positions;
  hsize_t atNum = 0;
  for (int i = 0; i < atomset.getNumChildren(); i++) 
  {
    if (atomset.getChild(i).getName() == "atom") 
    {
      const XmlNode& atNode = atomset.getChild(i);
      species_ids.push_back(SpeciesNameToInt[atNode.getAttribute("species")]);
      atNode.getChild("position").getValue(positions);
      atNum++;
    }
  }
  hsize_t dims[]={atNum,3};  
  writeNumsToHDF("positions", positions, atoms_group, 2, dims);
  writeNumsToHDF("species_ids", species_ids, atoms_group);
  writeNumsToHDF("number_of_atoms", static_cast<int>(atNum), atoms_group);
}

void EshdfFile::writeElectrons(const XmlNode& qboxSample) 
{
  const XmlNode& wfnNode = qboxSample.getChild("wavefunction");
  int nspin, nel;
  wfnNode.getAttribute("nspin", nspin);
  wfnNode.getAttribute("nel", nel);
  const XmlNode& gridNode = wfnNode.getChild("grid");
  int nx, ny, nz;
  gridNode.getAttribute("nx", nx);
  gridNode.getAttribute("ny", ny);
  gridNode.getAttribute("nz", nz);

  FftContainer fftCont(nx,ny,nz);

  vector<KPoint> kpts;
  map<KPoint, const XmlNode*> kptToUpNode;
  map<KPoint, const XmlNode*> kptToDnNode;

  for (int i = 0; i < wfnNode.getNumChildren(); i++) 
  {
    if (wfnNode.getChild(i).getName() == "slater_determinant") 
    {
      const XmlNode& sdNode = wfnNode.getChild(i);
      int spinIdx = sdNode.getAttributeIndex("spin");
      string species;
      if (spinIdx >= 0) 
      {
	species = sdNode.getAttribute(spinIdx);
      } 
      else 
      {
	species = "up";
      }

      KPoint kpt;
      stringstream ss;
      ss.str(sdNode.getAttribute("kpoint"));
      ss >> kpt.kx;
      ss >> kpt.ky;
      ss >> kpt.kz;

      bool newKpt = true;
      for (int j = 0; j < kpts.size(); j++) 
      {
	if (kpts[j] == kpt) 
        {
	  newKpt = false;
	}
      }
      if (newKpt) 
      {
	kpts.push_back(kpt);
      }
      
      if (species == "up") 
      {
	kptToUpNode[kpt] = std::addressof(sdNode);
      } 
      else 
      {
	kptToDnNode[kpt] = std::addressof(sdNode);
      }
    }
  }
  hid_t electrons_group = makeHDFGroup("electrons", file_);
  if (kpts.size() > 1) 
  {
    std::cerr << "Warning: Due to limitations of the current tool, extreme care" << std::endl;
    std::cerr << "is required if tiling to a supercell from qbox calculations with" << std::endl;
    std::cerr << "multiple k-points.  Specifically spo eigenvalues are not properly" << std::endl;
    std::cerr << "included, so improper choice of orbitals may result." << std::endl;
  }
  writeNumsToHDF("number_of_kpoints", static_cast<int>(kpts.size()), electrons_group);
  writeNumsToHDF("number_of_spins", nspin, electrons_group);

  double avgNup = 0.0;
  double avgNdn = 0.0;
  // go through kpt by kpt and write 
  for (int i = 0; i < kpts.size(); i++) 
  {
    stringstream kptElemName;
    kptElemName << "kpoint_" << i;
    hid_t kpt_group = makeHDFGroup(kptElemName.str(), electrons_group);
    vector<double> kvector;
    kvector.push_back(kpts[i].kx);
    kvector.push_back(kpts[i].ky);
    kvector.push_back(kpts[i].kz);
    //writeNumsToHDF("numsym", 1, kpt_group);
    //writeNumsToHDF("symgroup", 1, kpt_group);
    writeNumsToHDF("reduced_k", kvector, kpt_group);
    writeNumsToHDF("weight", 1.0/static_cast<double>(kpts.size()), kpt_group); 

    if (i == 0) 
    {
      // figure out the order of the g-vectors
      // write them to an integer array and then put it in the element
      vector<int> gvectors;
      for (int ix = 0; ix < nx; ix++) 
      {
	for (int iy = 0; iy < ny; iy++) 
        {
	  for (int iz = 0; iz < nz; iz++) 
          {
	    gvectors.push_back(wrapped(ix,nx));
	    gvectors.push_back(wrapped(iy,ny));
	    gvectors.push_back(wrapped(iz,nz));
	  }
	}
      }    
      
      hsize_t dims[]={static_cast<hsize_t>(nx*ny*nz),3};
      writeNumsToHDF("gvectors", gvectors, kpt_group, 2, dims);
      writeNumsToHDF("number_of_gvectors", nx*ny*nz, kpt_group);
    }


    // here is where we will read in both the occupations for the 
    // kpoint (species depenent if needed)
    // and also the wavefunction (in real space)
    // fourier transform it and write it into the appropriate hdf element
    double nup = 0;
    double ndn = 0;

    hid_t up_spin_group = makeHDFGroup("spin_0", kpt_group);
    const XmlNode* upnd = kptToUpNode[kpts[i]];    
    handleSpinGroup(upnd, up_spin_group, nup, fftCont);


    if (nspin == 2) 
    {
      hid_t dn_spin_group = makeHDFGroup("spin_1", kpt_group);
      const XmlNode* dnnd = kptToDnNode[kpts[i]];
      handleSpinGroup(dnnd, dn_spin_group, ndn, fftCont);
    } 
    else 
    {
      ndn = nup;
    }

    avgNup += nup / static_cast<double>(kpts.size());
    avgNdn += ndn / static_cast<double>(kpts.size());
  }
  vector<int> nels;
  nels.push_back(static_cast<int>(std::floor(avgNup+0.1)));
  nels.push_back(static_cast<int>(std::floor(avgNdn+0.1)));
  writeNumsToHDF("number_of_electrons", nels, electrons_group);
}

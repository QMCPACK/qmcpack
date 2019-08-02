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
    writeStringToHDF("program_name", "convertpw4qmc", h_creator);
    writeNumsToHDF("version", version, h_creator);
  }
}

void EshdfFile::writeFormat() 
{
  writeStringToHDF("format", "ES-HDF", file_);
}

void EshdfFile::writeQboxBoilerPlate(const XmlNode& qboxSample) 
{
  const string appName = "qbox";
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

int EshdfFile::getIntsOnly(const string& str) const
{
  stringstream numerals;
  int temp;
  for (int i = 0; i < str.size(); ++i) 
  {
    string blah = str.substr(i,1);
    if (stringstream(blah) >> temp) {
      numerals << str[i];
    }
  }
  int result;
  numerals >> result;
  return result;
}
	 

void EshdfFile::writeQEBoilerPlate(const XmlNode& qeXml)
{
  const string appName = "espresso";
  const XmlNode& generalInfo = qeXml.getChild("general_info");
  const XmlNode& creatorNode = generalInfo.getChild("creator");
  const string versionStr = creatorNode.getAttribute("VERSION");
  int minor = 0;
  int sub = 0;
  const int firstDotIdx = versionStr.find_first_of('.');
  const int secondDotIdx = versionStr.find_last_of('.');
  const int major = getIntsOnly(versionStr.substr(0,firstDotIdx));
  if (firstDotIdx == secondDotIdx) // this means on subersion is provided
  {
    minor = getIntsOnly(versionStr.substr(firstDotIdx+1));
  } 
  else {
    minor = getIntsOnly(versionStr.substr(firstDotIdx+1,secondDotIdx-firstDotIdx-1));
    sub = getIntsOnly(versionStr.substr(secondDotIdx+1));
  }  
  writeApplication(appName, major, minor, sub);
  writeVersion();
  writeCreator();
  writeFormat();
}

void EshdfFile::writeQboxSupercell(const XmlNode& qboxSample) 
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

void EshdfFile::writeQESupercell(const XmlNode& qeXml) 
{
  // grab the primitive translation vectors the output tag's atomic structure part
  const XmlNode& output = qeXml.getChild("output");
  const XmlNode& atomStructure = output.getChild("atomic_structure");
  const XmlNode& cell = atomStructure.getChild("cell");
  const XmlNode& a = cell.getChild("a1");
  const XmlNode& b = cell.getChild("a2");
  const XmlNode& c = cell.getChild("a3");

  stringstream ss;
  ss << a.getValue() << "  " << b.getValue() << "  " << c.getValue();

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

void EshdfFile::writeQEAtoms(const XmlNode& qeXml) 
{
  const XmlNode& output_xml = qeXml.getChild("output");
  const XmlNode& atomic_species_xml = output_xml.getChild("atomic_species");

  //make group
  hid_t atoms_group = makeHDFGroup("atoms", file_);
  
  map<string, int> species_name_to_int;  
  //go through each species, extract:
  //   name, mass
  //   in future, it would be good to extract atomic_number, and valence_charge
  //   then write to a file, also set up mapping between species name and number

  int species_num = 0;
  for (int i = 0; i < atomic_species_xml.getNumChildren(); i++)
  {
    if (atomic_species_xml.getChild(i).getName() == "species")
    {
      const XmlNode& species_xml = atomic_species_xml.getChild(i);
      string sp_name = species_xml.getAttribute("name");
      species_name_to_int[sp_name] = species_num;
      double mass;
      species_xml.getChild("mass").getValue(mass);

      stringstream gname;
      gname << "species_" << species_num;
      hid_t species_group = makeHDFGroup(gname.str(), atoms_group);
      writeStringToHDF("name", sp_name, species_group);
      writeNumsToHDF("mass", mass, species_group);
      species_num++;
    }
  }      
  writeNumsToHDF("number_of_species", species_num, atoms_group);

  const XmlNode& atomic_structure_xml = output_xml.getChild("atomic_structure");
  const XmlNode& atomic_positions_xml = atomic_structure_xml.getChild("atomic_positions");
  // go through atoms and extract their position and type 
  std::vector<int> species_ids;
  std::vector<double> positions;
  hsize_t at_num = 0;
  for (int i = 0; i < atomic_positions_xml.getNumChildren(); i++) 
  {
    if (atomic_positions_xml.getChild(i).getName() == "atom") 
    {
      const XmlNode& at_node_xml = atomic_positions_xml.getChild(i);
      species_ids.push_back(species_name_to_int[at_node_xml.getAttribute("name")]);
      at_node_xml.getValue(positions); // this will append the three numbers to the verctor
      at_num++;
    }
  }
  hsize_t dims[]={at_num,3};  
  writeNumsToHDF("positions", positions, atoms_group, 2, dims);
  writeNumsToHDF("species_ids", species_ids, atoms_group);
  writeNumsToHDF("number_of_atoms", static_cast<int>(at_num), atoms_group);
}


void EshdfFile::writeQboxAtoms(const XmlNode& qboxSample) 
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

void EshdfFile::handleDensity(const XmlNode& qeXml, const string& dir_name, int spinpol, hid_t el_group) {
  const XmlNode& output_xml = qeXml.getChild("output");
  // read in the grid to put in the charge density 
  const XmlNode& basis_set_xml = output_xml.getChild("basis_set");
  const XmlNode& fft_grid = basis_set_xml.getChild("fft_grid");
  int nr1;
  int nr2;
  int nr3;
  fft_grid.getAttribute("nr1", nr1);
  fft_grid.getAttribute("nr2", nr2);
  fft_grid.getAttribute("nr3", nr3);
  // also figure out how many gvectors there will be
  const XmlNode& ngm_xml = basis_set_xml.getChild("ngm");
  int num_dens_gvecs;
  ngm_xml.getValue(num_dens_gvecs);

  // now open the hdf file and read the necessary quantities
  herr_t status;
  hid_t file, dset;
  const string dens_fname = dir_name + "charge-density.hdf5";
  if (!file_exists(dens_fname)) {
    cout << "could not find " << dens_fname << ", make sure espresso is" << endl;
    cout << "compiled with hdf support and you do not explicitly set" << endl;
    cout << "wf_collect=.false. in your input" << endl;
    exit(1);
  }
  file = H5Fopen(dens_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  const string gvectors_str ="MillerIndices";
  int* gvecs = new int[3*num_dens_gvecs];
  readNumsFromHDF(gvectors_str, gvecs, file);

  const string totdens_str = "rhotot_g";
  double* dens = new double[2*num_dens_gvecs];
  readNumsFromHDF(totdens_str, dens, file);
  
  double* diffdens = new double[2*num_dens_gvecs];
  if (spinpol == 1) {
    const string diffdens_str = "rhodiff_g";
    readNumsFromHDF(diffdens_str, diffdens, file);
  }

  // now need to write everything out
  hid_t density_group = makeHDFGroup("density", el_group);

  hsize_t dims[]={static_cast<hsize_t>(num_dens_gvecs),3};
  writeNumsToHDF("gvectors", gvecs, density_group, 2, dims);
  delete[] gvecs;
  
  vector<int> grid;
  grid.push_back(nr1);
  grid.push_back(nr2);
  grid.push_back(nr3);
  writeNumsToHDF("mesh", grid, density_group);
  writeNumsToHDF("number_of_gvectors", num_dens_gvecs, density_group);
  hsize_t dims_dens[]={static_cast<hsize_t>(num_dens_gvecs),2};
  if (spinpol == 0) {
    hid_t spin_0_group = makeHDFGroup("spin_0", density_group);
    writeNumsToHDF("density_g", dens, spin_0_group, 2, dims_dens);
  } else {
    // do spin up
    double* working = new double[2*num_dens_gvecs];
    for (int i = 0; i < num_dens_gvecs*2; i++) {
      working[i] = (dens[i]+diffdens[i])/2.0;
    }
    hid_t spin_0_group = makeHDFGroup("spin_0", density_group);
    writeNumsToHDF("density_g", working, spin_0_group, 2, dims_dens);

    // do spin down
    for (int i = 0; i < num_dens_gvecs*2; i++) {
      working[i] = (dens[i]-diffdens[i])/2.0;
    }
    hid_t spin_1_group = makeHDFGroup("spin_1", density_group);
    writeNumsToHDF("density_g", working, spin_1_group, 2, dims_dens);
    delete[] working;
  }
  delete[] dens;
  delete[] diffdens;
}

void EshdfFile::processKPts(const XmlNode& band_structure_xml,
			    vector<vector<double > >& eigenvals,
			    vector<vector<double > >& occupations,
			    vector<KPoint>& kpts,
			    vector<double>& weights,
			    vector<int>& ngvecs) 
{
  int nks = 0;
  for (int i = 0; i < band_structure_xml.getNumChildren(); i++) 
  {
    const XmlNode& child_xml = band_structure_xml.getChild(i);
    if (child_xml.getName() == "ks_energies")
    {
      vector<double> eigs;
      vector<double> occs;
      vector<double> kpt_vec;
      KPoint kpt;
      double weight;
      int ngvec;
      
      nks++;
      const XmlNode& k_point_xml = child_xml.getChild("k_point");
      k_point_xml.getAttribute("weight", weight);
      k_point_xml.getValue(kpt_vec);
      kpt.kx=kpt_vec[0];
      kpt.ky=kpt_vec[1];
      kpt.kz=kpt_vec[2];
      const XmlNode& npw_xml = child_xml.getChild("npw");
      npw_xml.getValue(ngvec);
      const XmlNode& eig_xml = child_xml.getChild("eigenvalues");
      eig_xml.getValue(eigs);
      const XmlNode& occs_xml = child_xml.getChild("occupations");
      occs_xml.getValue(occs);
      
      eigenvals.push_back(eigs);
      occupations.push_back(occs);
      kpts.push_back(kpt);
      weights.push_back(weight);
      ngvecs.push_back(ngvec);
    }
  }
}

void EshdfFile::getNumElectrons(vector<vector<double> >& occupations, 
				vector<double>& weights, int& nup, int& ndn, 
				int spinpol, int noncol) 
{
  nup = 0;
  ndn = 0;

  double nup_flt = 0.0;
  double ndn_flt = 0.0;

  for (int i = 0; i < weights.size(); i++) 
  {
    if (noncol == 1) 
    {
      double sum = 0.0;
      for (int j = 0; j < occupations[i].size(); j++) 
      {
	sum += occupations[i][j];
      }
      nup_flt += sum*weights[i];
    }
    else if (spinpol == 0) 
    {
      double sum = 0.0;
      for (int j = 0; j < occupations[i].size(); j++) 
      {
	sum += occupations[i][j];
      }
      nup_flt += sum*weights[i]/2.0;
      ndn_flt += sum*weights[i]/2.0;
    }
    else 
    {
      double sum_up = 0.0;
      double sum_dn = 0.0;
      for (int j = 0; j < occupations[i].size()/2; j++) 
      {
	sum_up += occupations[i][j];
	sum_dn += occupations[i][j+occupations[i].size()/2];
      }
      nup_flt += sum_up*weights[i];
      ndn_flt += sum_dn*weights[i];
    }
  }
  nup = static_cast<int>(round(nup_flt));
  ndn = static_cast<int>(round(ndn_flt));
}

void EshdfFile::handleKpt(int kpt_num, const std::string& dir_name, KPoint& kpt, 
			  const std::vector<double>& eigenvalues, int ngvec, int maxgvec_idx,
			  double weight, int spinpol, int noncol, hid_t electrons_group)
{
  stringstream ss;
  ss << "kpoint_" << kpt_num;
  hid_t kpt_group = makeHDFGroup(ss.str(), electrons_group);
  writeNumsToHDF("weight", weight, kpt_group);
  vector<double> kpt_vector;
  kpt_vector.push_back(kpt.kx);
  kpt_vector.push_back(kpt.ky);
  kpt_vector.push_back(kpt.kz);
  writeNumsToHDF("reduced_k", kpt_vector, kpt_group);

  if (kpt_num == 0) 
  {
    herr_t status;
    hid_t file, dset;
    string fname;
    if (spinpol == 1) 
    {
      stringstream ss;
      ss << dir_name << "wfcup" << (kpt_num+1) << ".hdf5";
      fname = ss.str();
    }
    else 
    {
      stringstream ss;
      ss << dir_name << "wfc" << (kpt_num+1) << ".hdf5";
      fname = ss.str();
    }
    if (!file_exists(fname)) {
      cout << "could not find " << fname << ", make sure espresso is" << endl;
      cout << "compiled with hdf support and you do not explicitly set" << endl;
      cout << "wf_collect=.false. in your input" << endl;
      exit(1);
    }
    file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    int* gvecs = new int[3*ngvec];
    const string gvectors_str ="MillerIndices";
    readNumsFromHDF(gvectors_str, gvecs, file);
    hsize_t dims[]={static_cast<hsize_t>(ngvec),3};
    writeNumsToHDF("gvectors", gvecs, kpt_group, 2, dims);
    writeNumsToHDF("number_of_gvectors", ngvec, kpt_group);

    delete[] gvecs;
    H5Fclose(file);
  }

  //dfile_id = H5Fopen (FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  //dataset_id = H5Dopen2 (file_id, DATASETNAME, H5P_DEFAULT);
  
  /* Specify size and shape of subset to write. */
  
  //offset[0] = 1;
  //offset[1] = 2;
  
  //count[0]  = DIM0_SUB;  
  //count[1]  = DIM1_SUB;
  
  //stride[0] = 1;
  //stride[1] = 1;
  
  //block[0] = 1;
  //block[1] = 1;
  
  /* Create memory space with size of subset. Get file dataspace 
     and select subset from file dataspace. */
  
  //dimsm[0] = DIM0_SUB;
  //dimsm[1] = DIM1_SUB;
  //memspace_id = H5Screate_simple (RANK, dimsm, NULL); 
  
  //dataspace_id = H5Dget_space (dataset_id);
  //status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset,
  //                              stride, count, block);

    
}


void EshdfFile::writeQEElectrons(const XmlNode& qeXml, const string& dir_name) {
  // make electrons group in hdf file
  hid_t electrons_group = makeHDFGroup("electrons", file_);
  const XmlNode& output_xml = qeXml.getChild("output");
  // need to figure out if this is spin polarized, or if it is noncolinear
  int spinpol = 0;
  int noncol = 0;
  const XmlNode& band_structure_xml = output_xml.getChild("band_structure");
  const XmlNode& lsda_xml = band_structure_xml.getChild("lsda");
  string lsda_string_bool = lsda_xml.getValue();
  const XmlNode& noncolin_xml = band_structure_xml.getChild("noncolin");
  string noncolin_string_bool = noncolin_xml.getValue();
  if (lsda_string_bool == "true") spinpol = 1;
  if (noncolin_string_bool == "true") noncol = 1;

  // scrape xml file and associated hdf for density and write out
  handleDensity(qeXml, dir_name, spinpol, electrons_group);

  // read in information about kpts from the xml
  vector<vector<double> > eigenvals;
  vector<vector<double> > occupations;
  vector<KPoint> kpts;
  vector<double> weights;
  vector<int> ngvecs;
  processKPts(band_structure_xml, eigenvals, occupations, kpts, weights, ngvecs);

  // write number of kpoints, number of spins and spinors if appropriate
  writeNumsToHDF("number_of_kpoints", static_cast<int>(kpts.size()), electrons_group);
  writeNumsToHDF("has_spinors", noncol, electrons_group);
  if (noncol == 0) 
  {
    int nspins = 1;
    if (spinpol == 1) nspins = 2;
    writeNumsToHDF("number_of_spins", nspins, electrons_group);
  }

  // figure out how many electrons of each spin and write to file
  int nup;
  int ndn;  
  getNumElectrons(occupations, weights, nup, ndn, spinpol, noncol);
  vector<int> nels;
  nels.push_back(nup);
  nels.push_back(ndn);
  writeNumsToHDF("number_of_electrons", nels, electrons_group);


  int maxgvecs = 0;
  int maxgvecs_index = -1;
  for (int i = 0; i < ngvecs.size(); i++) 
  {
    if (ngvecs[i] > maxgvecs) {
      maxgvecs = ngvecs[i];
      maxgvecs_index = i;
    }
  }

  // now read and write the appropriate coefficients to the kpoint_xxx group
  for (int i = 0; i < kpts.size(); i++) 
  {
    handleKpt(i, dir_name, kpts[i], eigenvals[i], maxgvecs, maxgvecs_index, weights[i], 
	      spinpol, noncol, electrons_group);
  }
}
  


void EshdfFile::writeQboxElectrons(const XmlNode& qboxSample) 
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

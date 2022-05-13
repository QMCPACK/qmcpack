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
#include <sstream>
#include <map>
using namespace std;
//using namespace hdfhelper;
using namespace qmcplusplus;

class KPoint
{
public:
  double kx;
  double ky;
  double kz;
  KPoint()
  {
    kx = 0;
    ky = 0;
    kz = 0;
  }
  KPoint(double x, double y, double z)
  {
    kx = x;
    ky = y;
    kz = z;
  }
  KPoint(const KPoint& kp)
  {
    kx = kp.kx;
    ky = kp.ky;
    kz = kp.kz;
  }
  KPoint& operator=(const KPoint& kp)
  {
    kx = kp.kx;
    ky = kp.ky;
    kz = kp.kz;
    return *this;
  }
  bool operator==(const KPoint& kp) const
  {
    if ((std::abs(kx - kp.kx) < 1e-7) && (std::abs(ky - kp.ky) < 1e-7) && (std::abs(kz - kp.kz) < 1e-7))
    {
      return true;
    }
    return false;
  }
  bool operator<(const KPoint& kp) const
  {
    if (abs(kx - kp.kx) > 1e-7)
    {
      return kx < kp.kx;
    }
    else if (abs(ky - kp.ky) > 1e-7)
    {
      return ky < kp.ky;
    }
    else if (abs(kz - kp.kz) > 1e-7)
    {
      return kz < kp.kz;
    }
    return false;
  }
};


EshdfFile::EshdfFile(const string& hdfFileName) { outfile_.create(hdfFileName); }

EshdfFile::~EshdfFile() { outfile_.close(); }


int EshdfFile::wrapped(int i, int size) const
{
  if (i < size / 2)
  {
    return i;
  }
  else
  {
    return wrapped(i - size, size);
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

// before you enter here, better be in the a spin group of the hdf file
void EshdfFile::handleSpinGroup(const XmlNode* nd, double& nocc, FftContainer& cont)
{
  nocc             = getOccupation(nd);
  int stateCounter = 0;
  vector<double> eigvals;

  for (int chIdx = 0; chIdx < nd->getNumChildren(); chIdx++)
  {
    if (nd->getChild(chIdx).getName() == "grid_function")
    {
      //cout << "Working on state " << stateCounter << endl;
      stringstream statess;
      statess << "state_" << stateCounter;
      outfile_.push(statess.str());

      // HACK_HACK_HACK!!!
      // QBOX does not write out the eigenvalues for the states in the
      // sample file, so just make sure they are treated in ascending order
      eigvals.push_back(-5000.0 + stateCounter);

      const XmlNode& eigFcnNode = nd->getChild(chIdx);
      readInEigFcn(eigFcnNode, cont);
      // write eigfcn to proper place
      array<int, 2> psig_dims{cont.fullSize, 2};

      vector<double> temp;
      for (int i = 0; i < cont.fullSize; i++)
      {
        temp.push_back(cont.kspace[i][0]);
        temp.push_back(cont.kspace[i][1]);
      }

      outfile_.writeSlabReshaped(temp, psig_dims, "psi_g");
      stateCounter++;
      outfile_.pop();
    }
  }

  outfile_.write(stateCounter, "number_of_states");
  outfile_.write(eigvals, "eigenvalues");
}

void EshdfFile::readInEigFcn(const XmlNode& nd, FftContainer& cont)
{
  const string type     = nd.getAttribute("type");
  const string encoding = nd.getAttribute("encoding");

  if (encoding != "text")
  {
    cerr << "Don't yet know how to handle encoding of wavefunction values other than text" << endl;
    exit(1);
  }

  vector<double> values;
  nd.getValue(values);
  const double fixnorm = 1 / std::sqrt(static_cast<double>(cont.fullSize));
  if (type == "complex")
  {
    int index = 0;
    for (int ix = 0; ix < cont.getNx(); ix++)
    {
      for (int iy = 0; iy < cont.getNy(); iy++)
      {
        for (int iz = 0; iz < cont.getNz(); iz++)
        {
          const int qbx         = cont.getQboxIndex(ix, iy, iz);
          cont.rspace[index][0] = values[2 * qbx] * fixnorm;
          cont.rspace[index][1] = values[2 * qbx + 1] * fixnorm;
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
          const int qbx         = cont.getQboxIndex(ix, iy, iz);
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

  outfile_.push("application");
  string str = appName;
  outfile_.write(str, "code");
  outfile_.write(version, "version");
  outfile_.pop();
}


void EshdfFile::writeVersion()
{
  vector<int> version{2, 1, 0};
  outfile_.write(version, "version");
}

void EshdfFile::writeCreator()
{
  vector<int> version{0, 1, 0};
  outfile_.push("creator");
  string tmp_str = "convertpw4qmc";
  outfile_.write(tmp_str, "program_name");
  outfile_.write(version, "version");
  outfile_.pop();
}

void EshdfFile::writeFormat()
{
  string tmp_str = "ES-HDF";
  outfile_.write(tmp_str, "format");
}

void EshdfFile::writeQboxBoilerPlate(const XmlNode& qboxSample)
{
  const string appName       = "qbox";
  const XmlNode& description = qboxSample.getChild("description");
  string desString           = description.getValue();
  int versionStart           = desString.find("qbox-");
  versionStart += 5;
  string versionStr = desString.substr(versionStart);

  const int firstDotIdx  = versionStr.find_first_of('.');
  const int secondDotIdx = versionStr.find_last_of('.');
  const int major        = stoi(versionStr.substr(0, firstDotIdx));
  const int minor        = stoi(versionStr.substr(firstDotIdx + 1, secondDotIdx - firstDotIdx - 1));
  const int sub          = stoi(versionStr.substr(secondDotIdx + 1));

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
    string blah = str.substr(i, 1);
    if (stringstream(blah) >> temp)
    {
      numerals << str[i];
    }
  }
  int result;
  numerals >> result;
  return result;
}


void EshdfFile::writeQEBoilerPlate(const XmlNode& qeXml)
{
  const string appName       = "espresso";
  const XmlNode& creatorNode = qeXml.getChild("general_info").getChild("creator");
  const string versionStr    = creatorNode.getAttribute("VERSION");
  int minor                  = 0;
  int sub                    = 0;
  const int firstDotIdx      = versionStr.find_first_of('.');
  const int secondDotIdx     = versionStr.find_last_of('.');
  const int major            = getIntsOnly(versionStr.substr(0, firstDotIdx));
  if (firstDotIdx == secondDotIdx) // this means on subersion is provided
  {
    minor = getIntsOnly(versionStr.substr(firstDotIdx + 1));
  }
  else
  {
    minor = getIntsOnly(versionStr.substr(firstDotIdx + 1, secondDotIdx - firstDotIdx - 1));
    sub   = getIntsOnly(versionStr.substr(secondDotIdx + 1));
  }
  writeApplication(appName, major, minor, sub);
  writeVersion();
  writeCreator();
  writeFormat();
}

void EshdfFile::writeQboxSupercell(const XmlNode& qboxSample)
{
  // grab the primitive translation vectors from the atomset tag's attributes and put the entries in the vector ptvs
  const XmlNode& atomset   = qboxSample.getChild("atomset");
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
  array<int, 2> dims{3, 3};
  outfile_.push("supercell");
  outfile_.writeSlabReshaped(ptvs, dims, "primitive_vectors");
  outfile_.pop();
}

void EshdfFile::writeQESupercell(const XmlNode& qeXml)
{
  // grab the primitive translation vectors the output tag's atomic structure part
  const XmlNode& cell = qeXml.getChild("output").getChild("atomic_structure").getChild("cell");
  const XmlNode& a    = cell.getChild("a1");
  const XmlNode& b    = cell.getChild("a2");
  const XmlNode& c    = cell.getChild("a3");

  stringstream ss;
  ss << a.getValue() << "  " << b.getValue() << "  " << c.getValue();

  vector<double> ptvs;
  double temp;
  while (ss >> temp)
  {
    ptvs.push_back(temp);
  }

  // write the ptvs to the supercell group of the hdf file
  array<int, 2> dims{3, 3};
  outfile_.push("supercell");
  outfile_.writeSlabReshaped(ptvs, dims, "primitive_vectors");
  outfile_.pop();
}

void EshdfFile::writeQEAtoms(const XmlNode& qeXml)
{
  const XmlNode& atomic_species_xml = qeXml.getChild("output").getChild("atomic_species");

  //make group
  outfile_.push("atoms");

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
      const XmlNode& species_xml   = atomic_species_xml.getChild(i);
      string sp_name               = species_xml.getAttribute("name");
      species_name_to_int[sp_name] = species_num;
      double mass;
      species_xml.getChild("mass").getValue(mass);

      stringstream gname;
      gname << "species_" << species_num;
      outfile_.push(gname.str());
      outfile_.write(sp_name, "name");
      outfile_.write(mass, "mass");
      species_num++;
      outfile_.pop();
    }
  }
  outfile_.write(species_num, "number_of_species");

  const XmlNode& atomic_positions_xml =
      qeXml.getChild("output").getChild("atomic_structure").getChild("atomic_positions");
  // go through atoms and extract their position and type
  std::vector<int> species_ids;
  std::vector<double> positions;
  int at_num = 0;
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
  array<int, 2> dims{at_num, 3};
  outfile_.writeSlabReshaped(positions, dims, "positions");
  outfile_.write(species_ids, "species_ids");
  outfile_.write(at_num, "number_of_atoms");

  outfile_.pop();
}


void EshdfFile::writeQboxAtoms(const XmlNode& qboxSample)
{
  const XmlNode& atomset = qboxSample.getChild("atomset");

  //make group
  outfile_.push("atoms");

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

      string spName            = species.getAttribute("name");
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
      outfile_.push(gname.str());
      outfile_.write(atomic_number, "atomic_number");
      outfile_.write(mass, "mass");
      outfile_.write(val_charge, "valence_charge");
      outfile_.write(name, "name");
      string tmp_str = "unknown";
      outfile_.write(tmp_str, "pseudopotential");
      speciesNum++;
      outfile_.pop();
    }
  }
  outfile_.write(speciesNum, "number_of_species");

  // go through atoms and extract their position and type
  std::vector<int> species_ids;
  std::vector<double> positions;
  int at_num = 0;
  for (int i = 0; i < atomset.getNumChildren(); i++)
  {
    if (atomset.getChild(i).getName() == "atom")
    {
      const XmlNode& atNode = atomset.getChild(i);
      species_ids.push_back(SpeciesNameToInt[atNode.getAttribute("species")]);
      atNode.getChild("position").getValue(positions);
      at_num++;
    }
  }
  array<int, 2> dims{at_num, 3};
  outfile_.writeSlabReshaped(positions, dims, "positions");
  outfile_.write(species_ids, "species_ids");
  outfile_.write(at_num, "number_of_atoms");
  outfile_.pop();
}

hdf_archive EshdfFile::openHdfFileForRead(const string& fname)
{
  hdf_archive file;
  bool result = file.open(fname);
  if (!result)
  {
    cout << "could not find " << fname << ", make sure espresso is" << endl;
    cout << "compiled with hdf support and you do not explicitly set" << endl;
    cout << "wf_collect=.false. in your input" << endl;
    exit(1);
  }
  return file;
}

// need to be in the electrons group when entering this
void EshdfFile::handleDensity(const XmlNode& qeXml, const string& dir_name, int spinpol)
{
  const XmlNode& output_xml = qeXml.getChild("output");
  // read in the grid to put in the charge density
  const XmlNode& basis_set_xml = output_xml.getChild("basis_set");
  const XmlNode& fft_grid      = basis_set_xml.getChild("fft_grid");
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
  const string dens_fname = dir_name + "charge-density.hdf5";
  hdf_archive densfile    = openHdfFileForRead(dens_fname);

  vector<int> readDims;
  densfile.getShape<int>("MillerIndices", readDims);
  vector<int> gvecs;
  array<int, 2> shape{readDims[0], readDims[1]};
  densfile.readSlabReshaped(gvecs, shape, "MillerIndices");

  vector<double> dens;
  densfile.read(dens, "rhotot_g");

  vector<double> diffdens;
  if (spinpol == 1)
  {
    densfile.read(diffdens, "rhodiff_g");
  }

  // now need to write everything out
  outfile_.push("density");
  array<int, 2> dims{num_dens_gvecs, 3};
  outfile_.writeSlabReshaped(gvecs, dims, "gvectors");

  vector<int> grid{nr1, nr2, nr3};
  outfile_.write(grid, "mesh");
  outfile_.write(num_dens_gvecs, "number_of_gvectors");

  array<int, 2> dims_dens{num_dens_gvecs, 2};
  if (spinpol == 0)
  {
    outfile_.push("spin_0");
    outfile_.writeSlabReshaped(dens, dims_dens, "density_g");
    outfile_.pop();
  }
  else
  {
    // do spin up
    vector<double> working(2 * num_dens_gvecs);
    for (int i = 0; i < num_dens_gvecs * 2; i++)
    {
      working[i] = (dens[i] + diffdens[i]) / 2.0;
    }
    outfile_.push("spin_0");
    outfile_.writeSlabReshaped(working, dims_dens, "density_g");
    outfile_.pop();

    // do spin down
    for (int i = 0; i < num_dens_gvecs * 2; i++)
    {
      working[i] = (dens[i] - diffdens[i]) / 2.0;
    }
    outfile_.push("spin_1");
    outfile_.writeSlabReshaped(working, dims_dens, "density_g");
    outfile_.pop();
  }
  outfile_.pop(); // get out of density group
}

void EshdfFile::processKPts(const XmlNode& band_structure_xml,
                            const vector<double>& ptvs,
                            vector<vector<double>>& eigenvals,
                            vector<vector<double>>& occupations,
                            vector<KPoint>& kpts,
                            vector<double>& weights,
                            vector<int>& ngvecs)
{
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

      const XmlNode& k_point_xml = child_xml.getChild("k_point");
      k_point_xml.getAttribute("weight", weight);
      k_point_xml.getValue(kpt_vec);
      kpt.kx = kpt_vec[0] * ptvs[0] + kpt_vec[1] * ptvs[1] + kpt_vec[2] * ptvs[2];
      kpt.ky = kpt_vec[0] * ptvs[3] + kpt_vec[1] * ptvs[4] + kpt_vec[2] * ptvs[5];
      kpt.kz = kpt_vec[0] * ptvs[6] + kpt_vec[1] * ptvs[7] + kpt_vec[2] * ptvs[8];

      child_xml.getChild("npw").getValue(ngvec);
      child_xml.getChild("eigenvalues").getValue(eigs);
      child_xml.getChild("occupations").getValue(occs);

      eigenvals.push_back(eigs);
      occupations.push_back(occs);
      kpts.push_back(kpt);
      weights.push_back(weight);
      ngvecs.push_back(ngvec);
    }
  }
}

void EshdfFile::getNumElectrons(vector<vector<double>>& occupations,
                                vector<double>& weights,
                                int& nup,
                                int& ndn,
                                int spinpol,
                                int noncol)
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
      nup_flt += sum * weights[i];
      ndn_flt = -1.0;
    }
    else if (spinpol == 0)
    {
      double sum = 0.0;
      for (int j = 0; j < occupations[i].size(); j++)
      {
        sum += occupations[i][j];
      }
      nup_flt += sum * weights[i] / 2.0;
      ndn_flt += sum * weights[i] / 2.0;
    }
    else
    {
      double sum_up = 0.0;
      double sum_dn = 0.0;
      for (int j = 0; j < occupations[i].size() / 2; j++)
      {
        sum_up += occupations[i][j];
        sum_dn += occupations[i][j + occupations[i].size() / 2];
      }
      nup_flt += sum_up * weights[i];
      ndn_flt += sum_dn * weights[i];
    }
  }
  nup = static_cast<int>(round(nup_flt));
  ndn = static_cast<int>(round(ndn_flt));
}

vector<double> EshdfFile::getPtvs(const XmlNode& qeXml)
{
  const XmlNode& reciprocal_lattice_xml = qeXml.getChild("output").getChild("basis_set").getChild("reciprocal_lattice");
  vector<double> rlv; // reciprocal_lattice_vectors
  reciprocal_lattice_xml.getChild("b1").getValue(rlv);
  reciprocal_lattice_xml.getChild("b2").getValue(rlv);
  reciprocal_lattice_xml.getChild("b3").getValue(rlv);

  const double det = rlv[0] * (rlv[4] * rlv[8] - rlv[5] * rlv[7]) - rlv[1] * (rlv[3] * rlv[8] - rlv[5] * rlv[6]) +
      rlv[2] * (rlv[3] * rlv[7] - rlv[4] * rlv[6]);
  const double invdet = 1.0 / det;
  vector<double> ptv;
  ptv.push_back(invdet * (rlv[4] * rlv[8] - rlv[5] * rlv[7]));
  ptv.push_back(invdet * (rlv[5] * rlv[6] - rlv[3] * rlv[8]));
  ptv.push_back(invdet * (rlv[3] * rlv[7] - rlv[4] * rlv[6]));
  ptv.push_back(invdet * (rlv[2] * rlv[7] - rlv[1] * rlv[8]));
  ptv.push_back(invdet * (rlv[0] * rlv[8] - rlv[2] * rlv[6]));
  ptv.push_back(invdet * (rlv[1] * rlv[6] - rlv[0] * rlv[7]));
  ptv.push_back(invdet * (rlv[1] * rlv[5] - rlv[2] * rlv[4]));
  ptv.push_back(invdet * (rlv[2] * rlv[3] - rlv[0] * rlv[5]));
  ptv.push_back(invdet * (rlv[0] * rlv[4] - rlv[1] * rlv[3]));
  return ptv;
}

void EshdfFile::readKptGvecs(int kpt_num, const string& dir_name, int spinpol, momap_t& morefmap)

{
  stringstream ss;
  string fname;
  if (spinpol == 1)
    ss << dir_name << "wfcup" << (kpt_num + 1) << ".hdf5";
  else
    ss << dir_name << "wfc" << (kpt_num + 1) << ".hdf5";
  fname                = ss.str();
  hdf_archive wfc_file = openHdfFileForRead(fname);
  vector<int> readShape;
  wfc_file.getShape<int>("MillerIndices", readShape);
  vector<int> gvecs;
  array<int, 2> shape{readShape[0], readShape[1]};
  wfc_file.readSlabReshaped(gvecs, shape, "MillerIndices");
  wfc_file.close();

  for (int i = 0; i < readShape[0]; i++)
  {
    vector<int> gv;
    for (int d = 0; d < readShape[1]; d++)
      gv.push_back(gvecs[i * readShape[1] + d]);
    mopair_t p(gv, std::complex<double>(0, 0));
    morefmap.insert(p);
  }
}

// should already be in the electrons group when you call this
void EshdfFile::handleKpt(int kpt_num,
                          const std::string& dir_name,
                          KPoint& kpt,
                          const std::vector<double>& eigenvalues,
                          double weight,
                          int spinpol,
                          int noncol,
                          const momap_t& moref)
{
  stringstream ss;
  ss << "kpoint_" << kpt_num;
  outfile_.push(ss.str());
  outfile_.write(weight, "weight");

  vector<double> kpt_vector{kpt.kx, kpt.ky, kpt.kz};
  outfile_.write(kpt_vector, "reduced_k");

  //reread gvecs for this kpts to create a map to coeffs
  //then will merge with the maximal gvec map in moref
  //this will fill unused gvecs with zeros
  string fname;
  if (spinpol == 1)
  {
    stringstream ss;
    ss << dir_name << "wfcup" << (kpt_num + 1) << ".hdf5";
    fname = ss.str();
  }
  else
  {
    stringstream ss;
    ss << dir_name << "wfc" << (kpt_num + 1) << ".hdf5";
    fname = ss.str();
  }
  hdf_archive wfc_file = openHdfFileForRead(fname);
  vector<int> readShape;
  wfc_file.getShape<int>("MillerIndices", readShape);
  vector<int> gvecs;
  array<int, 2> shape{readShape[0], readShape[1]};
  wfc_file.readSlabReshaped(gvecs, shape, "MillerIndices");
  wfc_file.close();

  vector<vector<int>> gvs;
  for (int i = 0; i < readShape[0]; i++)
  {
    vector<int> g;
    for (int d = 0; d < readShape[1]; d++)
      g.push_back(gvecs[i * readShape[1] + d]);
    gvs.push_back(g);
  }

  // now do the non-gvectors related stuff
  hdf_archive spin_0_file;
  hdf_archive spin_1_file;

  if (spinpol == 1)
  {
    stringstream ss0;
    ss0 << dir_name << "wfcup" << (kpt_num + 1) << ".hdf5";
    spin_0_file.open(ss0.str());
    stringstream ss1;
    ss1 << dir_name << "wfcdw" << (kpt_num + 1) << ".hdf5";
    spin_1_file.open(ss1.str());
  }
  else
  {
    stringstream ss0;
    ss0 << dir_name << "wfc" << (kpt_num + 1) << ".hdf5";
    spin_0_file.open(ss0.str());
  }

  // set up storage space for the coefficients to be read from the file
  vector<double> upcoefs;
  vector<double> dncoefs;
  int ng = gvs.size();
  if (spinpol == 0)
  {
    if (noncol == 0)
      upcoefs.resize(ng * 2);
    if (noncol == 1)
    {
      upcoefs.resize(ng * 4);
      dncoefs.resize(ng * 2);
    }
  }
  else
  {
    upcoefs.resize(ng * 2);
    dncoefs.resize(ng * 2);
  }

  vector<momap_t> states_up;
  vector<momap_t> states_dn;
  int states_to_loop = eigenvalues.size();
  if (spinpol == 1)
    states_to_loop /= 2;
  for (int state = 0; state < states_to_loop; state++)
  {
    // set the elements of coefs to 0
    for (int i = 0; i < upcoefs.size(); i++)
      upcoefs[i] = 0.0;
    for (int i = 0; i < dncoefs.size(); i++)
      dncoefs[i] = 0.0;

    // do what is necessary to read only this band's coefficients
    array<int, 2> read_from{state, -1};
    if (spinpol == 0)
    {
      spin_0_file.readSlabSelection(upcoefs, read_from, "evc");
    }
    else
    {
      spin_0_file.readSlabSelection(upcoefs, read_from, "evc");
      spin_1_file.readSlabSelection(dncoefs, read_from, "evc");
    }

    if (spinpol == 0)
    {
      if (noncol == 0)
      {
        momap_t moup;
        for (int i = 0; i < ng; i++)
        {
          complex<double> c(upcoefs[i * 2], upcoefs[i * 2 + 1]);
          mopair_t p(gvs[i], c);
          moup.insert(p);
        }
        //now fill in rest of gvecs with zeros
        moup.insert(moref.begin(), moref.end());
        states_up.push_back(moup);
      }
      else
      {
        //dn part of spinor in second half of upcoefs
        for (int i = 0; i < ng * 2; i++)
          dncoefs[i] = upcoefs[i + 2 * ng];
        momap_t moup, modn;
        for (int i = 0; i < ng; i++)
        {
          complex<double> c(upcoefs[i * 2], upcoefs[i * 2 + 1]);
          mopair_t p(gvs[i], c);
          moup.insert(p);
          c = complex<double>(dncoefs[i * 2], dncoefs[i * 2 + 1]);
          p = mopair_t(gvs[i], c);
          modn.insert(p);
        }
        //now fill in rest of gvecs with zeros
        moup.insert(moref.begin(), moref.end());
        modn.insert(moref.begin(), moref.end());
        states_up.push_back(moup);
        states_dn.push_back(modn);
      }
    }
    else
    {
      momap_t moup, modn;
      for (int i = 0; i < ng; i++)
      {
        complex<double> c(upcoefs[i * 2], upcoefs[i * 2 + 1]);
        mopair_t p(gvs[i], c);
        moup.insert(p);
        c = complex<double>(dncoefs[i * 2], dncoefs[i * 2 + 1]);
        p = mopair_t(gvs[i], c);
        modn.insert(p);
      }
      moup.insert(moref.begin(), moref.end());
      modn.insert(moref.begin(), moref.end());
      states_up.push_back(moup);
      states_dn.push_back(modn);
    }
  }
  spin_0_file.close();
  spin_1_file.close();

  //now write to eshdf
  if (kpt_num == 0)
  {
    vector<int> allgvs;
    int nallgvs = moref.size();
    int dim     = moref.begin()->first.size();
    array<int, 2> shape{nallgvs, dim};
    for (auto& v : moref)
    {
      for (int d = 0; d < dim; d++)
        allgvs.push_back(v.first[d]);
    }
    outfile_.writeSlabReshaped(allgvs, shape, "gvectors");
    outfile_.write(nallgvs, "number_of_gvectors");
  }

  for (int state = 0; state < states_up.size(); state++)
  {
    stringstream ss;
    ss << "state_" << state;
    outfile_.push("spin_0");
    outfile_.push(ss.str());
    vector<double> c;
    for (auto& v : states_up[state])
    {
      c.push_back(v.second.real());
      c.push_back(v.second.imag());
    }
    array<int, 2> dims{static_cast<int>(states_up[state].size()), 2};
    outfile_.writeSlabReshaped(c, dims, "psi_g");
    outfile_.pop();
    outfile_.pop();
    if (noncol == 1 || spinpol == 1)
    {
      outfile_.push("spin_1");
      outfile_.push(ss.str());
      vector<double> c;
      for (auto& v : states_dn[state])
      {
        c.push_back(v.second.real());
        c.push_back(v.second.imag());
      }
      array<int, 2> dims{static_cast<int>(states_dn[state].size()), 2};
      outfile_.writeSlabReshaped(c, dims, "psi_g");
      outfile_.pop();
      outfile_.pop();
    }
  }

  // now all the states are written, so write out eigenvalues and number of states
  vector<double> eigval = eigenvalues;
  if (spinpol == 0)
  {
    outfile_.push("spin_0");
    int int_eig_sz = static_cast<int>(eigenvalues.size());
    outfile_.write(int_eig_sz, "number_of_states");
    outfile_.write(eigval, "eigenvalues");
    outfile_.pop();
    if (noncol == 1)
    {
      outfile_.push("spin_1");
      outfile_.write(int_eig_sz, "number_of_states");
      outfile_.write(eigval, "eigenvalues");
      outfile_.pop();
    }
  }
  else // spin polarized case
  {
    int totstates = eigenvalues.size();
    int upstates  = totstates / 2;
    int dnstates  = totstates / 2;

    vector<double> upeig;
    vector<double> dneig;
    for (int i = 0; i < upstates; i++)
    {
      upeig.push_back(eigenvalues[i]);
      dneig.push_back(eigenvalues[i + upstates]);
    }

    outfile_.push("spin_0");
    outfile_.write(upstates, "number_of_states");
    outfile_.write(upeig, "eigenvalues");
    outfile_.pop();

    outfile_.push("spin_1");
    outfile_.write(dnstates, "number_of_states");
    outfile_.write(dneig, "eigenvalues");
    outfile_.pop();
  }

  outfile_.pop(); // get out of the kpoint_ kpt_num group
}

void EshdfFile::writeQEElectrons(const XmlNode& qeXml, const string& dir_name)
{
  // make electrons group in hdf file
  outfile_.push("electrons");
  const XmlNode& output_xml = qeXml.getChild("output");

  // need to figure out if this is spin polarized, or if it is noncolinear
  int spinpol                       = 0;
  int noncol                        = 0;
  const XmlNode& band_structure_xml = output_xml.getChild("band_structure");

  const XmlNode& lsda_xml     = band_structure_xml.getChild("lsda");
  string lsda_string_bool     = lsda_xml.getValue();
  const XmlNode& noncolin_xml = band_structure_xml.getChild("noncolin");
  string noncolin_string_bool = noncolin_xml.getValue();
  if (lsda_string_bool == "true")
    spinpol = 1;
  if (noncolin_string_bool == "true")
    noncol = 1;

  // scrape xml file and associated hdf for density and write out
  handleDensity(qeXml, dir_name, spinpol);

  // read in information about kpts from the xml
  vector<double> ptv = getPtvs(qeXml);
  vector<vector<double>> eigenvals;
  vector<vector<double>> occupations;
  vector<KPoint> kpts;
  vector<double> weights;
  vector<int> ngvecs;
  processKPts(band_structure_xml, ptv, eigenvals, occupations, kpts, weights, ngvecs);

  // write number of kpoints, number of spins and spinors if appropriate
  int int_kpt_sz = static_cast<int>(kpts.size());
  outfile_.write(int_kpt_sz, "number_of_kpoints");
  outfile_.write(noncol, "has_spinors");

  if (noncol == 0)
  {
    int nspins = 1;
    if (spinpol == 1)
      nspins = 2;
    outfile_.write(nspins, "number_of_spins");
  }

  // figure out how many electrons of each spin and write to file
  int nup;
  int ndn;
  getNumElectrons(occupations, weights, nup, ndn, spinpol, noncol);
  vector<int> nels;
  nels.push_back(nup);
  nels.push_back(ndn);
  outfile_.write(nels, "number_of_electrons");

  //find maximal set of gvecs
  momap_t moref;
  for (int i = 0; i < kpts.size(); i++)
    readKptGvecs(i, dir_name, spinpol, moref);
  //moref now has maximal set of gvecs
  for (int i = 0; i < kpts.size(); i++)
    handleKpt(i, dir_name, kpts[i], eigenvals[i], weights[i], spinpol, noncol, moref);
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

  FftContainer fftCont(nx, ny, nz);

  vector<KPoint> kpts;
  map<KPoint, const XmlNode*> kptToUpNode;
  map<KPoint, const XmlNode*> kptToDnNode;

  for (int i = 0; i < wfnNode.getNumChildren(); i++)
  {
    if (wfnNode.getChild(i).getName() == "slater_determinant")
    {
      const XmlNode& sdNode = wfnNode.getChild(i);
      int spinIdx           = sdNode.getAttributeIndex("spin");
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
  outfile_.push("electrons");
  if (kpts.size() > 1)
  {
    std::cout << "Warning: Due to limitations of the current tool, extreme care" << std::endl;
    std::cout << "is required if tiling to a supercell from qbox calculations with" << std::endl;
    std::cout << "multiple k-points.  Specifically spo eigenvalues are not properly" << std::endl;
    std::cout << "included, so improper choice of orbitals may result." << std::endl;
  }
  const int int_kpts_sz = static_cast<int>(kpts.size());
  outfile_.write(int_kpts_sz, "number_of_kpoints");
  outfile_.write(nspin, "number_of_spins");

  double avgNup = 0.0;
  double avgNdn = 0.0;
  // go through kpt by kpt and write
  for (int i = 0; i < kpts.size(); i++)
  {
    stringstream kptElemName;
    kptElemName << "kpoint_" << i;
    outfile_.push(kptElemName.str());

    vector<double> kvector;
    kvector.push_back(kpts[i].kx);
    kvector.push_back(kpts[i].ky);
    kvector.push_back(kpts[i].kz);
    //writeNumsToHDF("numsym", 1, kpt_group);
    //writeNumsToHDF("symgroup", 1, kpt_group);
    outfile_.write(kvector, "reduced_k");
    const double dbl_weight = 1.0 / static_cast<double>(kpts.size());
    outfile_.write(dbl_weight, "weight");

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
            gvectors.push_back(wrapped(ix, nx));
            gvectors.push_back(wrapped(iy, ny));
            gvectors.push_back(wrapped(iz, nz));
          }
        }
      }

      array<int, 2> dims{nx * ny * nz, 3};
      outfile_.writeSlabReshaped(gvectors, dims, "gvectors");
      const int int_tot_num = nx * ny * nz;
      outfile_.write(int_tot_num, "number_of_gvectors");
    }


    // here is where we will read in both the occupations for the
    // kpoint (species depenent if needed)
    // and also the wavefunction (in real space)
    // fourier transform it and write it into the appropriate hdf element
    double nup = 0;
    double ndn = 0;

    outfile_.push("spin_0");
    const XmlNode* upnd = kptToUpNode[kpts[i]];
    handleSpinGroup(upnd, nup, fftCont);
    outfile_.pop();

    if (nspin == 2)
    {
      outfile_.push("spin_1");
      const XmlNode* dnnd = kptToDnNode[kpts[i]];
      handleSpinGroup(dnnd, ndn, fftCont);
      outfile_.pop();
    }
    else
    {
      ndn = nup;
    }

    avgNup += nup / static_cast<double>(kpts.size());
    avgNdn += ndn / static_cast<double>(kpts.size());
  }
  outfile_.pop();
  vector<int> nels;
  nels.push_back(static_cast<int>(std::floor(avgNup + 0.1)));
  nels.push_back(static_cast<int>(std::floor(avgNdn + 0.1)));
  outfile_.write(nels, "number_of_electrons");
  outfile_.pop();
}

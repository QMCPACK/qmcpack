#ifndef AFQMC_READWFN_HPP
#define AFQMC_READWFN_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctype.h>

#include "Utilities/SimpleParser.h"
#include "readWfn.h"

#include "AFQMC/config.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/csr_hdf5_readers.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"

namespace qmcplusplus
{
namespace
{
void read_header(std::ifstream& in, std::string& type, int& wfn_type, bool& fullMOMat, bool& Cstyle, int& ndet)
{
  std::vector<std::string> words;
  getwords(words, in);
  do
  {
    if (words.size() == 0)
      app_error() << "Format error in ASCII integral file. End of file in header. \n";
    for (std::vector<std::string>::iterator it = words.begin(); it != words.end(); it++)
    {
      if (*it == "&FCI")
      {
        // do nothing
      }
      else if (*it == "Type" || *it == "TYPE" || *it == "type")
      {
        if (it + 1 == words.end())
        {
          app_error() << "Format error in ASCII integral file. NAEB \n";
          APP_ABORT("Format error in ASCII integral file. NAEB \n");
        }
        type = *(it + 1);
        it++;
      }
      else if (*it == "NCI" || *it == "nci")
      {
        if (it + 1 == words.end())
        {
          app_error() << "Format error in ASCII integral file. NETOT \n";
          APP_ABORT("Format error in ASCII integral file. NETOT \n");
        }
        ndet = atoi((++it)->c_str());
      }
      else if (*it == "UHF" || *it == "GHF")
      {
        if (it + 1 == words.end())
        {
          app_error() << "Format error in ASCII integral file. UHF/GHF \n";
          APP_ABORT("Format error in ASCII integral file. UHF/GHF \n");
        }
        wfn_type = atoi((++it)->c_str());
        switch (wfn_type)
        {
        case 0: {
          app_log() << "Reading a RHF-type trial wave-function. \n";
          break;
        }
        case 1: {
          app_log() << "Reading a UHF-type trial wave-function. \n";
          break;
        }
        case 2: {
          app_log() << "Reading a GHF-type trial wave-function. \n";
          break;
        }
        default: {
          app_error() << "Unknown wave-function type in AFQMC/Utilities/readWfn.hpp: " << wfn_type << std::endl;
          APP_ABORT("Unknown wave-function type in AFQMC/Utilities/readWfn.hpp \n");
        }
        }
      }
      else if (*it == "FullMO" || *it == "FULLMO")
      {
        fullMOMat = true;
      }
      else if (*it == "CMajor")
      {
        Cstyle = false;
      }
    }
    getwords(words, in);
    if (words.size() == 0)
      app_error() << "Format error in ASCII integral file. End of file in header. \n";
  } while ((words[0].find(std::string("/")) == std::string::npos &&
            words[0].find(std::string("&END")) == std::string::npos));
}

void skip_determinant(std::ifstream& in, bool Cstyle, bool fullMOMat, int NMO, int NAEA)
{
  int nread;
  afqmc::ComplexType dummy;
  if (Cstyle)
  {
    nread = fullMOMat ? NMO : NAEA;
    for (int i = 0; i < NMO; i++)
      for (int j = 0; j < nread; j++)
      {
        in >> dummy;
        if (in.fail())
        {
          app_error() << "Problems reading ASCII file in readWfn.hpp.  \n";
          app_error() << i << " " << j << std::endl;
          in.close();
          APP_ABORT("Problems reading ASCII file in readWfn.hpp.  \n");
        }
      }
  }
  else
  {
    nread = fullMOMat ? NMO : NAEA;
    for (int j = 0; j < nread; j++)
      for (int i = 0; i < NMO; i++)
      {
        in >> dummy;
        if (in.fail())
        {
          app_error() << "Problems reading ASCII file in readWfn.hpp. \n";
          app_error() << i << " " << j << std::endl;
          in.close();
          APP_ABORT("Problems reading wfn in readWfn.hpp.\n");
        }
      }
  }
}

template<class Mat>
void read_mat(std::ifstream& in, Mat&& OrbMat, bool Cstyle, bool fullMOMat, int NMO, int NAEA)
{
  int nread;
  afqmc::ComplexType dummy;
  if (Cstyle)
  {
    nread = fullMOMat ? NMO : NAEA;
    for (int i = 0; i < NMO; i++)
      for (int j = 0; j < nread; j++)
      {
        in >> dummy;
        if (j < NAEA)
          OrbMat[i][j] = dummy;
        if (in.fail())
        {
          app_error() << "Problems reading ASCII file in readWfn.hpp.  \n";
          app_error() << i << " " << j << std::endl;
          in.close();
          APP_ABORT("Problems reading ASCII file in readWfn.hpp.  \n");
        }
      }
  }
  else
  {
    nread = fullMOMat ? NMO : NAEA;
    for (int j = 0; j < nread; j++)
      for (int i = 0; i < NMO; i++)
      {
        in >> dummy;
        if (j < NAEA)
          OrbMat[i][j] = dummy;
        if (in.fail())
        {
          app_error() << "Problems reading ASCII file in readWfn.hpp. \n";
          app_error() << i << " " << j << std::endl;
          in.close();
          APP_ABORT("Problems reading wfn in readWfn.hpp.\n");
        }
      }
  }
}

} // namespace

namespace afqmc
{
WALKER_TYPES getWalkerType(std::string filename)
{
  std::ifstream in;
  in.open(filename.c_str());
  if (in.fail())
  {
    app_error() << "Problems opening file:  " << filename << std::endl;
    APP_ABORT("Problems opening ASCII integral file. \n");
  }

  bool fullMOMat   = false;
  bool Cstyle      = true;
  int wfn_type     = 0;
  int ndet_in_file = -1;
  std::string type;

  read_header(in, type, wfn_type, fullMOMat, Cstyle, ndet_in_file);
  in.close();

  if (wfn_type == 0)
    return CLOSED;
  else if (wfn_type == 1)
    return COLLINEAR;
  else if (wfn_type == 2)
    return NONCOLLINEAR;
  else
    return UNDEFINED_WALKER_TYPE;
}

WALKER_TYPES getWalkerTypeHDF5(std::string filename, std::string type)
{
  hdf_archive dump;
  if (!dump.open(filename, H5F_ACC_RDONLY))
  {
    std::cerr << " Error opening wavefunction file in read_info_from_wfn. \n";
    APP_ABORT("");
  }
  if (!dump.push("Wavefunction", false))
  {
    std::cerr << " Error in getWalkerTypeHDF5: Group Wavefunction found. \n";
    APP_ABORT("");
  }
  if (!dump.push(type, false))
  {
    std::cerr << " Error in getWalkerTypeHDF5: Group " << type << " not found. \n";
    APP_ABORT("");
  }

  std::vector<int> Idata(5);
  if (!dump.readEntry(Idata, "dims"))
  {
    std::cerr << " Error in getWalkerTypeHDF5: Problems reading dims. \n";
    APP_ABORT("");
  }

  dump.pop();
  int wfn_type = Idata[3];

  if (wfn_type == 1)
    return CLOSED;
  else if (wfn_type == 2)
    return COLLINEAR;
  else if (wfn_type == 3)
    return NONCOLLINEAR;
  else
    return UNDEFINED_WALKER_TYPE;
}

std::string getWfnType(std::ifstream& in)
{
  in.clear();
  in.seekg(0, std::ios::beg);

  bool fullMOMat   = false;
  bool Cstyle      = true;
  int wfn_type     = 0;
  int ndet_in_file = -1;
  std::string type;
  read_header(in, type, wfn_type, fullMOMat, Cstyle, ndet_in_file);

  return type;
}

/*
 * Reads ndets from the ascii file. 
 */
void read_general_wavefunction(std::ifstream& in,
                               int& ndets,
                               WALKER_TYPES walker_type,
                               boost::mpi3::shared_communicator& comm,
                               int NMO,
                               int NAEA,
                               int NAEB,
                               std::vector<PsiT_Matrix>& PsiT,
                               std::vector<ComplexType>& ci)
{
  in.clear();
  in.seekg(0, std::ios::beg);

  assert(walker_type != UNDEFINED_WALKER_TYPE);
  std::string type;
  bool fullMOMat   = false;
  bool Cstyle      = true;
  int wfn_type     = 0;
  int ndet_in_file = -1;

  /*
   * type:
   *   - matrix: Slater matrix for all terms in the expansion
   */
  read_header(in, type, wfn_type, fullMOMat, Cstyle, ndet_in_file);

  std::transform(type.begin(), type.end(), type.begin(), (int (*)(int))tolower);

  /*
   * Expected order of inputs and tags:
   * Coefficients: 
   * Determinant: 
   */

  if (ndets <= 0)
    ndets = ndet_in_file;

  if (ndet_in_file < ndets)
    APP_ABORT("Error: Requesting too many determinants from wfn file.\n");
  if (type != "matrix")
    APP_ABORT("Error: Expects type=matrix in read_general_wavefunction.\n");
  if (wfn_type == 1 && walker_type == CLOSED)
    APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");
  if (wfn_type == 2 && (walker_type == CLOSED || walker_type == COLLINEAR))
    APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");

  if (walker_type == COLLINEAR && (NAEA == NAEB && wfn_type == 0))
    app_log() << "  MESSAGE: Using walker_type=colinear with a closed-shell wfn with NAEA==NAEB. \n"
              << "           Consider doing a closed shell calculation ( walker_type=closed in WalkerSet)\n";

  ci.reserve(ndets);

  std::string tag;
  in >> tag;
  if (tag != "Coefficients:")
    APP_ABORT(" Error: Expecting Coefficients: tag in wavefunction file. \n");
  ComplexType dummy;
  for (int i = 0; i < ndet_in_file; i++)
  {
    in >> dummy;
    if (i < ndets)
      ci.emplace_back(dummy);
  }

  if (type == "matrix")
  {
    PsiT.reserve((walker_type != COLLINEAR) ? ndets : 2 * ndets);

    if (wfn_type == 0)
    {
      boost::multi::array<ComplexType, 2> OrbMat({NMO, NAEA});
      for (int i = 0, q = 0; i < ndets; i++)
      {
        if (comm.rank() == 0)
        {
          in >> tag >> q;
          if (tag != "Determinant:" || q != i + 1)
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in, OrbMat, Cstyle, fullMOMat, NMO, NAEA);
        }
        if (walker_type == CLOSED)
        {
          PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
        }
        else if (walker_type == COLLINEAR)
        {
          PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
          PsiT.emplace_back(
              csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat(OrbMat.extension(0), {0, NAEB}), 1e-8,
                                                                       'H', comm));
        }
        else if (walker_type == NONCOLLINEAR)
        {
          APP_ABORT(" Error in readWfn: wfn_type==closed with walker_type=noncollinear.\n");
          /*
          boost::multi::array<ComplexType,2> Mat({2*NMO,NAEA+NAEB});
          if(comm.rank()==0) {
            std::fill_n(Mat.origin(),2*NMO*(NAEA+NAEB),ComplexType(0.0));
            Mat({0,NMO},{0,NAEA}) = OrbMat;
            Mat({NMO,2*NMO},{NAEA,NAEA+NAEB}) = OrbMat({0,NMO},{0,NAEB});
          }  
          PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(
                                        OrbMat,1e-8,'H',comm));
*/
        }
      }
    }
    else if (wfn_type == 1)
    {
      if (walker_type != COLLINEAR)
        APP_ABORT(" Error in readWfn: wfn_type==collinear with walker_type!=collinear.\n");
      boost::multi::array<ComplexType, 2> OrbMat({NMO, NAEA});
      for (int i = 0, q = 0; i < ndets; i++)
      {
        if (comm.rank() == 0)
        {
          in >> tag >> q;
          if (tag != "Determinant:" || q != i + 1)
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in, OrbMat, Cstyle, fullMOMat, NMO, NAEA);
        }

        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
        if (comm.rank() == 0)
          read_mat(in, OrbMat(OrbMat.extension(0), {0, NAEB}), Cstyle, fullMOMat, NMO, NAEB);
        PsiT.emplace_back(
            csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat(OrbMat.extension(0), {0, NAEB}), 1e-8, 'H',
                                                                     comm));
      }
    }
    else if (wfn_type == 2)
    {
      if (walker_type != NONCOLLINEAR)
        APP_ABORT(" Error in readWfn: wfn_type==collinear with walker_type!=collinear.\n");
      if (NAEB != 0)
        APP_ABORT(" Error in readWfn: walker_type==collinear with NAEB!=0.\n");
      boost::multi::array<ComplexType, 2> OrbMat({2 * NMO, NAEA});
      for (int i = 0, q = 0; i < ndets; i++)
      {
        if (comm.rank() == 0)
        {
          in >> tag >> q;
          if (tag != "Determinant:" || q != i + 1)
            APP_ABORT(" Error: Expecting Determinant: # tag in wavefunction file. \n");
          read_mat(in, OrbMat, Cstyle, fullMOMat, 2 * NMO, NAEA);
        }
        PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
      }
    } //type
  }
  else
  {
    APP_ABORT("Error: Unknown wavefunction type in file. Expected type matrix .\n");
  }
}

// only root reads
/*
 * The number of terms returned by this routine in PsiT depends on wfn_type, NOT on walker_type.
 */
ph_excitations<int, ComplexType> read_ph_wavefunction(std::ifstream& in,
                                                      int& ndets,
                                                      WALKER_TYPES walker_type,
                                                      boost::mpi3::shared_communicator& comm,
                                                      int NMO,
                                                      int NAEA,
                                                      int NAEB,
                                                      std::vector<PsiT_Matrix>& PsiT)
{
  if (comm.root())
  {
    in.clear();
    in.seekg(0, std::ios::beg);
  }

  assert(walker_type != UNDEFINED_WALKER_TYPE);
  bool fullMOMat   = false;
  bool Cstyle      = true;
  int wfn_type     = 0;
  int ndet_in_file = -1;
  bool mixed       = false;
  std::string type;

  /*
   * Expected order of inputs and tags:
   * Reference: 
   * Configurations: 
   */

  /*
   * type:
   *   - occ: All determinants are specified with occupation numbers 
   *   - mixed: mixed representation. Reference determinant in matrix form, 
   *            determinant list (including reference) in occupation numbers. 
   *
   * wfn_type:
   *   - 0: excitations out of a RHF reference  
   *          NOTE: Does not mean perfect pairing, means excitations from a single reference
   *   - 1: excitations out of a UHF reference  
   */
  if (comm.root())
  {
    read_header(in, type, wfn_type, fullMOMat, Cstyle, ndet_in_file);

    std::transform(type.begin(), type.end(), type.begin(), (int (*)(int))tolower);


    if (walker_type != COLLINEAR)
      APP_ABORT(" Error: walker_type!=COLLINEAR not yet implemented in read_ph_wavefunction.\n");

    if (ndet_in_file < ndets)
      APP_ABORT("Error: Requesting too many determinants from wfn file.\n");
    if (type != "occ" && type != "mixed")
      APP_ABORT("Error: Expect type=occ or type==mixed in read_ph_wavefunction.\n");
    if (wfn_type == 1 && walker_type == CLOSED)
      APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");
    if (wfn_type == 2 && (walker_type == CLOSED || walker_type == COLLINEAR))
      APP_ABORT("Error in read_wavefunction: walker_type < wfn_type. \n");

    if (type == "mixed")
      mixed = true;

    comm.broadcast_n(&mixed, 1, 0);
    comm.broadcast_n(&ndet_in_file, 1, 0);
    comm.broadcast_n(&wfn_type, 1, 0);
    comm.broadcast_n(&fullMOMat, 1, 0);
  }
  else
  {
    comm.broadcast_n(&mixed, 1, 0);
    comm.broadcast_n(&ndet_in_file, 1, 0);
    comm.broadcast_n(&wfn_type, 1, 0);
    comm.broadcast_n(&fullMOMat, 1, 0);
  }

  if (ndets <= 0)
    ndets = ndet_in_file;

  if (mixed)
  { // read reference
    int nmo_ = (walker_type == NONCOLLINEAR ? 2 * NMO : NMO);
    if (not comm.root())
      nmo_ = 0; // only root reads matrices
    if (not fullMOMat)
      APP_ABORT("Error: Wavefunction type mixed requires fullMOMat=true.\n");
    PsiT.reserve((wfn_type != 1) ? 1 : 2);

    boost::multi::array<ComplexType, 2> OrbMat({nmo_, nmo_});
    if (comm.root())
    {
      std::string tag;
      in >> tag;
      if (tag != "Reference:")
        APP_ABORT(" Error: Expecting Reference tag in wavefunction file. \n");
      read_mat(in, OrbMat, Cstyle, fullMOMat, nmo_, nmo_);
    }
    PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
    if (wfn_type == 1)
    {
      if (comm.root())
        read_mat(in, OrbMat, Cstyle, fullMOMat, NMO, NMO);
      PsiT.emplace_back(csr::shm::construct_csr_matrix_single_input<PsiT_Matrix>(OrbMat, 1e-8, 'H', comm));
    }
  }

  if (comm.root())
  {
    std::string tag;
    in >> tag;
    if (tag != "Configurations:")
    {
      app_error() << " tag: " << tag << std::endl;
      APP_ABORT(" Error: Expecting Configurations: tag in wavefunction file. \n");
    }
  }

  ComplexType ci;
  // count number of k-particle excitations
  // counts[0] has special meaning, it must be equal to NAEA+NAEB.
  std::vector<size_t> counts_alpha(NAEA + 1);
  std::vector<size_t> counts_beta(NAEB + 1);
  // ugly but need dynamic memory allocation
  std::vector<std::vector<int>> unique_alpha(NAEA + 1);
  std::vector<std::vector<int>> unique_beta(NAEB + 1);
  // reference configuration, taken as the first one right now
  std::vector<int> refa;
  std::vector<int> refb;
  // space to read configurations
  std::vector<int> confg;
  // space for excitation string identifying the current configuration
  std::vector<int> exct;
  // record file position to come back
  std::vector<int> Iwork; // work arrays for permutation calculation
  std::streampos start;
  if (comm.root())
  {
    confg.reserve(NAEA);
    Iwork.resize(2 * NAEA);
    exct.reserve(2 * NAEA);
    start = in.tellg();
    for (int i = 0; i < ndets; i++)
    {
      in >> ci;
      if (in.fail())
        APP_ABORT(" Error: Reading wfn file.\n");
      // alpha
      confg.clear();
      for (int k = 0, q = 0; k < NAEA; k++)
      {
        in >> q;
        if (in.fail())
          APP_ABORT(" Error: Reading wfn file.\n");
        if (q < 1 || q > NMO)
          APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
        confg.emplace_back(q - 1);
      }
      if (i == 0)
      {
        refa = confg;
      }
      else
      {
        int np = get_excitation_number(true, refa, confg, exct, ci, Iwork);
        push_excitation(exct, unique_alpha[np]);
      }
      // beta
      confg.clear();
      for (int k = 0, q = 0; k < NAEB; k++)
      {
        in >> q;
        if (in.fail())
          APP_ABORT(" Error: Reading wfn file.\n");
        if (q <= NMO || q > 2 * NMO)
          APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
        confg.emplace_back(q - 1);
      }
      if (i == 0)
      {
        refb = confg;
      }
      else
      {
        int np = get_excitation_number(true, refb, confg, exct, ci, Iwork);
        push_excitation(exct, unique_beta[np]);
      }
    }
    // now that we have all unique configurations, count
    for (int i = 1; i <= NAEA; i++)
      counts_alpha[i] = unique_alpha[i].size();
    for (int i = 1; i <= NAEB; i++)
      counts_beta[i] = unique_beta[i].size();
  }
  comm.broadcast_n(counts_alpha.begin(), counts_alpha.size());
  comm.broadcast_n(counts_beta.begin(), counts_beta.size());
  // using int for now, but should move to short later when everything works well
  // ph_struct stores the reference configuration on the index [0]
  ph_excitations<int, ComplexType> ph_struct(ndets, NAEA, NAEB, counts_alpha, counts_beta, shared_allocator<int>(comm));

  if (comm.root())
  {
    in.clear();
    in.seekg(start);
    std::map<int, int> refa2loc;
    for (int i = 0; i < NAEA; i++)
      refa2loc[refa[i]] = i;
    std::map<int, int> refb2loc;
    for (int i = 0; i < NAEB; i++)
      refb2loc[refb[i]] = i;
    // add reference
    ph_struct.add_reference(refa, refb);
    // add unique configurations
    // alpha
    for (int n = 1; n < unique_alpha.size(); n++)
      for (std::vector<int>::iterator it = unique_alpha[n].begin(); it < unique_alpha[n].end(); it += (2 * n))
        ph_struct.add_alpha(n, it);
    // beta
    for (int n = 1; n < unique_beta.size(); n++)
      for (std::vector<int>::iterator it = unique_beta[n].begin(); it < unique_beta[n].end(); it += (2 * n))
        ph_struct.add_beta(n, it);
    // read configurations
    int alpha_index;
    int beta_index;
    int np;
    for (int i = 0; i < ndets; i++)
    {
      in >> ci;
      if (in.fail())
        APP_ABORT(" Error: Reading wfn file.\n");
      confg.clear();
      for (int k = 0, q = 0; k < NAEA; k++)
      {
        in >> q;
        if (in.fail())
          APP_ABORT(" Error: Reading wfn file.\n");
        if (q < 1 || q > NMO)
          APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
        confg.emplace_back(q - 1);
      }
      np = get_excitation_number(true, refa, confg, exct, ci, Iwork);
      alpha_index =
          ((np == 0) ? (0)
                     : (find_excitation(exct, unique_alpha[np]) + ph_struct.number_of_unique_smaller_than(np)[0]));
      confg.clear();
      for (int k = 0, q = 0; k < NAEB; k++)
      {
        in >> q;
        if (in.fail())
          APP_ABORT(" Error: Reading wfn file.\n");
        if (q <= NMO || q > 2 * NMO)
          APP_ABORT("Error: Bad occupation number in wavefunction file. \n");
        confg.emplace_back(q - 1);
      }
      np = get_excitation_number(true, refb, confg, exct, ci, Iwork);
      beta_index =
          ((np == 0) ? (0) : (find_excitation(exct, unique_beta[np]) + ph_struct.number_of_unique_smaller_than(np)[1]));
      ph_struct.add_configuration(alpha_index, beta_index, ci);
    }
  }
  comm.barrier();
  return ph_struct;
}

void read_ph_wavefunction_hdf(hdf_archive& dump,
                              std::vector<ComplexType>& ci_coeff,
                              std::vector<int>& occs,
                              int& ndets,
                              WALKER_TYPES walker_type,
                              boost::mpi3::shared_communicator& comm,
                              int NMO,
                              int NAEA,
                              int NAEB,
                              std::vector<PsiT_Matrix>& PsiT,
                              std::string& type)
{
  using Alloc = shared_allocator<ComplexType>;
  assert(walker_type != UNDEFINED_WALKER_TYPE);
  int wfn_type = 0;
  bool mixed   = false;

  /*
   * Expected order of inputs and tags:
   * Reference:
   * Configurations:
   */

  /*
   * type:
   *   - occ: All determinants are specified with occupation numbers
   *
   * wfn_type:
   *   - 0: excitations out of a RHF reference
   *          NOTE: Does not mean perfect pairing, means excitations from a single reference
   *   - 1: excitations out of a UHF reference
   */
  WALKER_TYPES wtype;
  getCommonInput(dump, NMO, NAEA, NAEB, ndets, ci_coeff, wtype, comm.root());
  if (walker_type != wtype && NAEA != NAEB)
    APP_ABORT(" Error: When Different walker_type between hdf5 and xml inputs when NAEA!=NAEB (not allowed). \n");
  if (walker_type != COLLINEAR)
    APP_ABORT(" Error: walker_type!=COLLINEAR not yet implemented in read_ph_wavefunction.\n");

  int type_;
  if (!dump.readEntry(type_, "type"))
  {
    app_error() << " Error in WavefunctionFactory::fromHDF5(): Problems reading type. \n";
    APP_ABORT("");
  }
  if (type_ == 0)
    type = "occ";
  else
    type = "mixed";
  if (type == "mixed")
    mixed = true;

  if (mixed)
  {
    PsiT.reserve((wfn_type != 1) ? 1 : 2);

    if (!dump.push(std::string("PsiT_") + std::to_string(0), false))
    {
      app_error() << " Error in WavefunctionFactory: Group PsiT not found. \n";
      APP_ABORT("");
    }
    PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, Alloc>(dump, comm));
    dump.pop();
    if (wfn_type == 1)
    {
      if (wtype == CLOSED)
      {
        if (!dump.push(std::string("PsiT_") + std::to_string(0), false))
        {
          app_error() << " Error in WavefunctionFactory: Group PsiT not found. \n";
          APP_ABORT("");
        }
      }
      else if (wtype == COLLINEAR)
      {
        if (!dump.push(std::string("PsiT_") + std::to_string(1), false))
        {
          app_error() << " Error in WavefunctionFactory: Group PsiT not found. \n";
          APP_ABORT("");
        }
      }
      PsiT.emplace_back(csr_hdf5::HDF2CSR<PsiT_Matrix, Alloc>(dump, comm));
      dump.pop();
    }
  }
  if (!dump.readEntry(occs, "occs"))
    APP_ABORT("Error reading occs array.\n");
  comm.barrier();
}

ph_excitations<int, ComplexType> build_ph_struct(std::vector<ComplexType> ci_coeff,
                                                 boost::multi::array_ref<int, 2>& occs,
                                                 int ndets,
                                                 boost::mpi3::shared_communicator& comm,
                                                 int NMO,
                                                 int NAEA,
                                                 int NAEB)
{
  using Alloc = shared_allocator<ComplexType>;

  ComplexType ci;
  // count number of k-particle excitations
  // counts[0] has special meaning, it must be equal to NAEA+NAEB.
  std::vector<size_t> counts_alpha(NAEA + 1);
  std::vector<size_t> counts_beta(NAEB + 1);
  // ugly but need dynamic memory allocation
  std::vector<std::vector<int>> unique_alpha(NAEA + 1);
  std::vector<std::vector<int>> unique_beta(NAEB + 1);
  // reference configuration, taken as the first one right now
  std::vector<int> refa;
  std::vector<int> refb;
  // space to read configurations
  std::vector<int> confg;
  // space for excitation string identifying the current configuration
  std::vector<int> exct;
  // record file position to come back
  std::vector<int> Iwork; // work arrays for permutation calculation
  std::streampos start;
  if (comm.root())
  {
    confg.reserve(NAEA);
    Iwork.resize(2 * NAEA);
    exct.reserve(2 * NAEA);
    for (int i = 0; i < ndets; i++)
    {
      ci = ci_coeff[i];
      // alpha
      confg.clear();
      for (int k = 0, q = 0; k < NAEA; k++)
      {
        q = occs[i][k];
        if (q < 0 || q >= NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      if (i == 0)
      {
        refa = confg;
      }
      else
      {
        int np = get_excitation_number(true, refa, confg, exct, ci, Iwork);
        push_excitation(exct, unique_alpha[np]);
      }
      // beta
      confg.clear();
      for (int k = 0, q = 0; k < NAEB; k++)
      {
        q = occs[i][NAEA + k];
        if (q < NMO || q >= 2 * NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      if (i == 0)
      {
        refb = confg;
      }
      else
      {
        int np = get_excitation_number(true, refb, confg, exct, ci, Iwork);
        push_excitation(exct, unique_beta[np]);
      }
    }
    // now that we have all unique configurations, count
    for (int i = 1; i <= NAEA; i++)
      counts_alpha[i] = unique_alpha[i].size();
    for (int i = 1; i <= NAEB; i++)
      counts_beta[i] = unique_beta[i].size();
  }
  comm.broadcast_n(counts_alpha.begin(), counts_alpha.size());
  comm.broadcast_n(counts_beta.begin(), counts_beta.size());
  // using int for now, but should move to short later when everything works well
  // ph_struct stores the reference configuration on the index [0]
  ph_excitations<int, ComplexType> ph_struct(ndets, NAEA, NAEB, counts_alpha, counts_beta, shared_allocator<int>(comm));

  if (comm.root())
  {
    std::map<int, int> refa2loc;
    for (int i = 0; i < NAEA; i++)
      refa2loc[refa[i]] = i;
    std::map<int, int> refb2loc;
    for (int i = 0; i < NAEB; i++)
      refb2loc[refb[i]] = i;
    // add reference
    ph_struct.add_reference(refa, refb);
    // add unique configurations
    // alpha
    for (int n = 1; n < unique_alpha.size(); n++)
      for (std::vector<int>::iterator it = unique_alpha[n].begin(); it < unique_alpha[n].end(); it += (2 * n))
        ph_struct.add_alpha(n, it);
    // beta
    for (int n = 1; n < unique_beta.size(); n++)
      for (std::vector<int>::iterator it = unique_beta[n].begin(); it < unique_beta[n].end(); it += (2 * n))
        ph_struct.add_beta(n, it);
    // read configurations
    int alpha_index;
    int beta_index;
    int np;
    for (int i = 0; i < ndets; i++)
    {
      ci = ci_coeff[i];
      confg.clear();
      for (int k = 0, q = 0; k < NAEA; k++)
      {
        q = occs[i][k];
        if (q < 0 || q >= NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      np = get_excitation_number(true, refa, confg, exct, ci, Iwork);
      alpha_index =
          ((np == 0) ? (0)
                     : (find_excitation(exct, unique_alpha[np]) + ph_struct.number_of_unique_smaller_than(np)[0]));
      confg.clear();
      for (int k = 0, q = 0; k < NAEB; k++)
      {
        q = occs[i][NAEA + k];
        if (q < NMO || q >= 2 * NMO)
          APP_ABORT("Error: Bad occupation number " << q << " in determinant " << i << " in wavefunction file. \n");
        confg.emplace_back(q);
      }
      np = get_excitation_number(true, refb, confg, exct, ci, Iwork);
      beta_index =
          ((np == 0) ? (0) : (find_excitation(exct, unique_beta[np]) + ph_struct.number_of_unique_smaller_than(np)[1]));
      ph_struct.add_configuration(alpha_index, beta_index, ci);
    }
  }
  comm.barrier();
  return ph_struct;
}

/*
 * Read trial wavefunction information from file.
*/
void getCommonInput(hdf_archive& dump,
                    int NMO,
                    int NAEA,
                    int NAEB,
                    int& ndets_to_read,
                    std::vector<ComplexType>& ci,
                    WALKER_TYPES& walker_type,
                    bool root)
{
  // check for consistency in parameters
  std::vector<int> dims(5);
  if (!dump.readEntry(dims, "dims"))
  {
    app_error() << " Error in getCommonInput(): Problems reading dims. \n";
    APP_ABORT("");
  }
  if (NMO != dims[0])
  {
    app_error() << " Error in getCommonInput(): Inconsistent NMO . \n";
    APP_ABORT("");
  }
  if (NAEA != dims[1])
  {
    app_error() << " Error in getCommonInput(): Inconsistent  NAEA. \n";
    APP_ABORT("");
  }
  if (NAEB != dims[2])
  {
    app_error() << " Error in getCommonInput(): Inconsistent  NAEB. \n";
    APP_ABORT("");
  }
  walker_type = afqmc::initWALKER_TYPES(dims[3]);
  // just read walker_type, to allow flexibility
  //if(walker_type != dims[3]) {
  //  app_error()<<" Error in getCommonInput(): Inconsistent  walker_type. \n";
  //  APP_ABORT("");
  //}
  if (ndets_to_read < 1)
    ndets_to_read = dims[4];
  app_log() << " - Number of determinants in trial wavefunction: " << ndets_to_read << "\n";
  if (ndets_to_read > dims[4])
  {
    app_error() << " Error in getCommonInput(): Inconsistent  ndets_to_read. \n";
    APP_ABORT("");
  }
  ci.resize(ndets_to_read);
  if (!dump.readEntry(ci, "ci_coeffs"))
  {
    app_error() << " Error in getCommonInput(): Problems reading ci_coeffs. \n";
    APP_ABORT("");
  }
  app_log() << " - Coefficient of first determinant: " << ci[0] << "\n";
}


// modify for multideterminant case based on type
int readWfn(std::string fileName, boost::multi::array<ComplexType, 3>& OrbMat, int NMO, int NAEA, int NAEB, int det)
{
  std::ifstream in;
  in.open(fileName.c_str());
  if (in.fail())
  {
    app_error() << "Problems opening ASCII integral file:  " << fileName << std::endl;
    APP_ABORT("Problems opening ASCII integral file. \n");
  }

  bool fullMOMat = false;
  bool Cstyle    = true;
  int wfn_type;
  int ndet = 1;
  std::string type;

  read_header(in, type, wfn_type, fullMOMat, Cstyle, ndet);

  if (ndet != 1)
    APP_ABORT("Error: readWfn is for single determinant wave functions. \n");
  if (type != "matrix")
    APP_ABORT("Error: Only type=matrix accepted in readWfn. \n");

  ComplexType dummy;
  std::string tag;
  in >> tag >> dummy;
  if (tag != "Coefficients:")
    APP_ABORT(" Error: Expecting Coefficients: tag in wavefunction file. \n");

  int q;
  in >> tag >> q;
  if (tag != "Determinant:" || q != 1)
    APP_ABORT(" Error: Expecting Determinant: 1 tag in wavefunction file. \n");

  if (wfn_type == 0)
  {
    OrbMat.reextent({1, NMO, NAEA});
    read_mat(in, OrbMat[0], Cstyle, fullMOMat, NMO, NAEA);
  }
  else if (wfn_type == 1)
  {
    OrbMat.reextent({2, NMO, NAEA});
    read_mat(in, OrbMat[0], Cstyle, fullMOMat, NMO, NAEA);
    read_mat(in, OrbMat[1], Cstyle, fullMOMat, NMO, NAEB);
  }
  else if (wfn_type == 2)
  {
    OrbMat.reextent({1, 2 * NMO, NAEA + NAEB});
    read_mat(in, OrbMat[0], Cstyle, fullMOMat, 2 * NMO, NAEA + NAEB);
  } //type

  in.close();

  return wfn_type;
}

} // namespace afqmc

} // namespace qmcplusplus
#endif

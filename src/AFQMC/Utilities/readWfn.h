#ifndef AFQMC_READWFN_H
#define AFQMC_READWFN_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctype.h>

#include "Utilities/SimpleParser.h"

#include "hdf/hdf_archive.h"
#include "AFQMC/config.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/csr_matrix_construct.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"

namespace qmcplusplus
{
namespace afqmc
{
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
                               std::vector<ComplexType>& ci);

ph_excitations<int, ComplexType> read_ph_wavefunction(std::ifstream& in,
                                                      int& ndets,
                                                      WALKER_TYPES walker_type,
                                                      boost::mpi3::shared_communicator& comm,
                                                      int NMO,
                                                      int NAEA,
                                                      int NAEB,
                                                      std::vector<PsiT_Matrix>& PsiT);

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
                              std::string& type);

ph_excitations<int, ComplexType> build_ph_struct(std::vector<ComplexType> ci_coeff,
                                                 boost::multi::array_ref<int, 2>& occs,
                                                 int ndets,
                                                 boost::mpi3::shared_communicator& comm,
                                                 int NMO,
                                                 int NAEA,
                                                 int NAEB);

void getCommonInput(hdf_archive& dump,
                    int NMO,
                    int NAEA,
                    int NAEB,
                    int& ndets_to_read,
                    std::vector<ComplexType>& ci,
                    WALKER_TYPES& walker_type,
                    bool root);

WALKER_TYPES getWalkerType(std::string filename);
WALKER_TYPES getWalkerTypeHDF5(std::string filename, std::string type);

std::string getWfnType(std::ifstream& in);

// modify for multideterminant case based on type
int readWfn(std::string fileName,
            boost::multi::array<ComplexType, 3>& OrbMat,
            int NMO,
            int NAEA,
            int NAEB,
            int det = 0);

} // namespace afqmc

} // namespace qmcplusplus

#endif

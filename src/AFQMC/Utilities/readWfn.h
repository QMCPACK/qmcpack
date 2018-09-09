#ifndef AFQMC_READWFN_H
#define AFQMC_READWFN_H

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<ctype.h>

#include "Utilities/SimpleParser.h"

#include "AFQMC/config.h"
#include "boost/multi_array.hpp"
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
void read_general_wavefunction(std::string filename, int& ndets, 
        std::string& type, WALKER_TYPES walker_type,
        boost::mpi3::shared_communicator& comm, int NMO, int NAEA, int NAEB,
        std::vector<PsiT_Matrix>& PsiT, std::vector<ComplexType>& ci);

ph_excitations read_ph_wavefunction(std::string filename, int& ndets, 
        std::string& type, WALKER_TYPES walker_type,
        boost::mpi3::shared_communicator& comm, int NMO, int NAEA, int NAEB,
        std::vector<PsiT_Matrix>& PsiT);


WALKER_TYPES getWalkerType(std::string filename);

// modify for multideterminant case based on type
int readWfn( std::string fileName, boost::multi_array<ComplexType,3>& OrbMat, int NMO, int NAEA, int NAEB, int det = 0);

bool readWfn( std::string fileName, ComplexMatrix& OrbMat, int NMO, int NAEA, int NAEB);

}

}

#endif


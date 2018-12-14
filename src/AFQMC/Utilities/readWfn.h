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

namespace qmcplusplus
{

namespace afqmc
{

/*
 * Reads ndets from the ascii file. 
 * If pureSD == false, PsiT contains the Slater Matrices of all the terms in the expansion.
 *      For walker_type==1, PsiT contains 2*ndets terms including both Alpha/Beta components.
 * If pureSD == true, PsiT contains only the reference determinant and excitations contains
 *      the occupation strings of all the determinants in the expansion, including the reference.  
 */ 
void read_wavefunction(std::string filename, int& ndets, std::string& type, WALKER_TYPES walker_type,
        boost::mpi3::shared_communicator& comm, int NMO, int NAEA, int NAEB,
        std::vector<PsiT_Matrix>& PsiT, std::vector<ComplexType>& ci,
        std::vector<int>& excitations); 

WALKER_TYPES getWalkerType(std::string filename);

// modify for multideterminant case based on type
int readWfn( std::string fileName, boost::multi_array<ComplexType,3>& OrbMat, int NMO, int NAEA, int NAEB, int det = 0);

bool readWfn( std::string fileName, ComplexMatrix& OrbMat, int NMO, int NAEA, int NAEB);

}

}

#endif


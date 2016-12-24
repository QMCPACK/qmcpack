#ifndef AFQMC_COMPARELIBRARIES
#define AFQMC_COMPARELIBRARIES

#include <vector>
#include "AFQMC/config.h"

namespace qmcplusplus
{

void compare_libraries(int NMO, int NAEA, int NAEB, std::vector<s2D<ComplexType> >& Propg_H1, std::vector<IndexType>& Propg_H1_indx, std::vector<s2D<ValueType> >& Vuv, std::vector<s2D<ComplexType> >& vn, std::vector<IndexType>& vn_indx); 

}

#endif

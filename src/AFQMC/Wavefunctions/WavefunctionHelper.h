#ifndef QMCPLUSPLUS_AFQMC_WAVEFUNCTIONHELPER_H
#define QMCPLUSPLUS_AFQMC_WAVEFUNCTIONHELPER_H

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cassert>
#include<complex>
#include<map>

#include "AFQMC/config.0.h"

// Helper functions for slater determinant routines 

namespace qmcplusplus
{

// taking from FCIQMC code. I'm sure there's a nicer way to do this using STL
int cmpDets(int NAEA, int NAEB, int* n, double &sg, std::vector<IndexType>::iterator sdet1, std::vector<IndexType>::iterator sdet2, std::vector<IndexType>& work );

}
#endif

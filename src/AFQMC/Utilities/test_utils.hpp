//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_TEST_UTILS_HPP
#define QMCPLUSPLUS_AFQMC_TEST_UTILS_HPP

#include<complex>
#include "io/hdf_archive.h"

namespace qmcplusplus
{
namespace afqmc
{

template<typename T>
struct TEST_DATA
{
  int NMO,NAEA,NAEB;
  T E0,E1,E2;
  T Xsum,Vsum;
};

inline std::tuple<int,int,int> read_info_from_hdf(std::string fileName)
{
    hdf_archive dump;
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      std::cerr<<" Error opening integral file in SparseGeneralHamiltonian. \n";
      APP_ABORT("");
    }
    if(!dump.push("Hamiltonian",false)) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Group not Hamiltonian found. \n";
      APP_ABORT("");
    }

    std::vector<int> Idata(8);
    if(!dump.read(Idata,"dims")) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Problems reading dims. \n";
      APP_ABORT("");
    }

    dump.pop();

    return std::make_tuple(Idata[3],Idata[4],Idata[5]);
}

template<typename T>
TEST_DATA<T>  read_test_results_from_hdf(std::string fileName)
{
    hdf_archive dump;
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      std::cerr<<" Error opening integral file in SparseGeneralHamiltonian. \n";
      APP_ABORT("");
    }
    if(!dump.push("Hamiltonian",false)) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Group not Hamiltonian found. \n";
      APP_ABORT("");
    }

    std::vector<int> Idata(8);
    if(!dump.read(Idata,"dims")) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Problems reading dims. \n";
      APP_ABORT("");
    }
    dump.pop();

    T E0(0),E1(0),E2(0),Xsum(0),Vsum(0);

    if(dump.push("TEST_RESULTS",false)) {
      dump.read(E0,"E0");
      dump.read(E1,"E1");
      dump.read(E2,"E2");
      dump.read(Xsum,"Xsum");
      dump.read(Vsum,"Vsum");
      dump.pop();
    }

    return TEST_DATA<T>{Idata[3],Idata[4],Idata[5],E0,E1,E2,Xsum,Vsum};
}

}
}

#endif

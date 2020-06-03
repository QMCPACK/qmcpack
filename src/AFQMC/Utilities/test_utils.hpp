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
    if(!dump.readEntry(Idata,"dims")) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Problems reading dims. \n";
      APP_ABORT("");
    }

    dump.pop();

    return std::make_tuple(Idata[3],Idata[4],Idata[5]);
}

inline std::tuple<int,int,int> read_info_from_wfn(std::string fileName, std::string type)
{
    hdf_archive dump;
    if(!dump.open(fileName,H5F_ACC_RDONLY)) {
      std::cerr<<" Error opening wavefunction file in read_info_from_wfn. \n";
      APP_ABORT("");
    }
    if(!dump.push("Wavefunction",false)) {
      std::cerr<<" Error in read_info_from_wfn(): Group not Wavefunction found. \n";
      APP_ABORT("");
    }
    if(!dump.push(type,false)) {
      std::cerr<<" Error in read_info_from_wfn(): Group " << type << " not found. \n";
      APP_ABORT("");
    }

    std::vector<int> Idata(5);
    if(!dump.readEntry(Idata,"dims")) {
      std::cerr<<" Error in read_info_from_wfn: Problems reading dims. \n";
      APP_ABORT("");
    }

    dump.pop();

    return std::make_tuple(Idata[0],Idata[1],Idata[2]);
}

template<typename T>
TEST_DATA<T>  read_test_results_from_hdf(std::string fileName, std::string wfn_type="")
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
    if(!dump.readEntry(Idata,"dims")) {
      std::cerr<<" Error in HamiltonianFactory::fromHDF5(): Problems reading dims. \n";
      APP_ABORT("");
    }
    dump.pop();

    T E0(0),E1(0),E2(0),Xsum(0),Vsum(0);

    if(dump.push("TEST_RESULTS",false)) {
      dump.read(E0,wfn_type+"_E0");
      dump.read(E1,wfn_type+"_E1");
      dump.read(E2,wfn_type+"_E2");
      dump.read(Xsum,wfn_type+"_Xsum");
      dump.read(Vsum,wfn_type+"_Vsum");
      dump.pop();
    }

    return TEST_DATA<T>{Idata[3],Idata[4],Idata[5],E0,E1,E2,Xsum,Vsum};
}

// Create a fake output hdf5 filename for unit tests.
inline std::string create_test_hdf(std::string& wfn_file, std::string& hamil_file)
{
  std::size_t startw = wfn_file.find_last_of("\\/");
  std::size_t endw = wfn_file.find_last_of(".");
  std::string wfn_base = wfn_file.substr(startw+1,endw-startw-1);

  std::size_t starth = hamil_file.find_last_of("\\/");
  std::size_t endh = hamil_file.find_last_of(".");
  std::string ham_base = hamil_file.substr(starth+1,endh-starth-1);

  return wfn_base + "_" + ham_base + ".h5";
}


}
}

#endif

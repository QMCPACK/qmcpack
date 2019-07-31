//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_THCHAMILTONIAN_H
#define QMCPLUSPLUS_AFQMC_THCHAMILTONIAN_H

#include<iostream>
#include<vector>
#include<map>
#include<fstream>

#include "io/hdf_archive.h"

#include "AFQMC/config.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Numerics/ma_operations.hpp"

#include "AFQMC/Hamiltonians/OneBodyHamiltonian.hpp"
#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

class THCHamiltonian: public OneBodyHamiltonian
{

  public:

  THCHamiltonian(AFQMCInfo const& info, xmlNodePtr cur, boost::multi::array<ValueType,2>&& h,
                 TaskGroup_& tg_, ValueType nucE=0, ValueType fzcE=0):
                            OneBodyHamiltonian(info,std::move(h),nucE,fzcE),
                            TG(tg_),cutoff_cholesky(1e-6),fileName(""),
                            useHalfRotatedMuv(true)
  {
    std::string str("yes");
    ParameterSet m_param;
    m_param.add(cutoff_cholesky,"cutoff_cholesky","double");
    m_param.add(fileName,"filename","std::string");
    m_param.add(str,"useHalfRotatedMuv","std::string");
    m_param.put(cur);

    std::transform(str.begin(),str.end(),str.begin(),(int (*)(int))tolower);
    if(str == "no" || str == "false") useHalfRotatedMuv=false;
  }

  ~THCHamiltonian() {}

  THCHamiltonian(THCHamiltonian const& other) = default;
  THCHamiltonian(THCHamiltonian && other) = default;
  THCHamiltonian& operator=(THCHamiltonian const& other) = default;
  THCHamiltonian& operator=(THCHamiltonian && other) = default;

  ValueType getNuclearCoulombEnergy() const { return OneBodyHamiltonian::NuclearCoulombEnergy; }

  boost::multi::array<ValueType,2> getH1() const{ return OneBodyHamiltonian::getH1(); }

  HamiltonianOperations getHamiltonianOperations(bool pureSD, bool addCoulomb, WALKER_TYPES type,
            std::vector<PsiT_Matrix>& PsiT, double cutvn, double cutv2,
            TaskGroup_& TGprop, TaskGroup_& TGwfn, hdf_archive& dump);

  ValueType H(IndexType I, IndexType J) const
  {  return OneBodyHamiltonian::H(I,J); }

  ValueType H(IndexType I, IndexType J, IndexType K, IndexType L) const
  {
    APP_ABORT("Error: Calling H(I,J,K,L) in THCHamiltonian. \n");
    return ValueType(0.0);
  }

  protected:

  TaskGroup_& TG;

  std::string fileName;

  RealType cutoff_cholesky;

  bool useHalfRotatedMuv;

};

}

}

#endif

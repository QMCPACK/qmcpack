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

#ifndef QMCPLUSPLUS_AFQMC_HAMOPSIO_HPP
#define QMCPLUSPLUS_AFQMC_HAMOPSIO_HPP

#include <fstream>

#include "hdf/hdf_multi.h"
#include "hdf/hdf_archive.h"

#include "AFQMC/config.h"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/HamiltonianOperations/SparseTensorIO.hpp"
#include "AFQMC/HamiltonianOperations/THCOpsIO.hpp"
//#ifdef QMC_COMPLEX
//#include "AFQMC/HamiltonianOperations/KP3IndexFactorizationIO.hpp"
//#endif

namespace qmcplusplus
{
namespace afqmc
{
HamiltonianOperations loadHamOps(hdf_archive& dump,
                                 WALKER_TYPES type,
                                 int NMO,
                                 int NAEA,
                                 int NAEB,
                                 std::vector<PsiT_Matrix>& PsiT,
                                 TaskGroup_& TGprop,
                                 TaskGroup_& TGwfn,
                                 RealType cutvn,
                                 RealType cutv2)
{
  int hops_type = -1;
  if (TGwfn.Global().root())
  {
    dump.push("HamiltonianOperations", false);

    if (dump.is_group(std::string("THCOps")))
      hops_type = 1;
    else if (dump.is_group(std::string("KP3IndexFactorization")))
      hops_type = 3;
    else if (dump.is_group(std::string("SparseTensor")))
    {
      dump.push("SparseTensor", false);
      std::vector<int> type;
      if (!dump.readEntry(type, "type"))
      {
        app_error() << " Error in loadHamOps: Problems reading type dataset. \n";
        APP_ABORT("");
      }
      if (type[0] == 11)
        hops_type = 211;
      else if (type[0] == 12)
        hops_type = 212;
      else if (type[0] == 21)
        hops_type = 221;
      else if (type[0] == 22)
        hops_type = 222;
      else
      {
        app_error() << " Unknown SparseTensor/type: " << type[0] << std::endl;
        APP_ABORT("");
      }
      dump.pop();
    }
    else
    {
      app_error() << " Error in loadHamOps: Unknown hdf5 format. \n";
      APP_ABORT("");
    }
    dump.pop();
  }
  TGwfn.Global().broadcast_value(hops_type);

  if (hops_type == 1)
    return HamiltonianOperations(loadTHCOps(dump, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2));
  else if (hops_type == 211)
    return HamiltonianOperations(
        loadSparseTensor<ValueType, ValueType>(dump, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2));
  else if (hops_type == 212)
    return HamiltonianOperations(
        loadSparseTensor<ValueType, ComplexType>(dump, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2));
  else if (hops_type == 221)
    return HamiltonianOperations(
        loadSparseTensor<ComplexType, ValueType>(dump, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2));
  else if (hops_type == 222)
    return HamiltonianOperations(
        loadSparseTensor<ComplexType, ComplexType>(dump, type, NMO, NAEA, NAEB, PsiT, TGprop, TGwfn, cutvn, cutv2));
  //  else if(hops_type == 3)
  //    return  HamiltonianOperations(loadKP3IndexFactorization(dump,type,NMO,NAEA,NAEB,PsiT,TGprop,TGwfn,cutvn,cutv2));

  app_error() << " Error in loadHamOps: Unknown HOps type: " << hops_type << std::endl;
  ;
  APP_ABORT("");
  return HamiltonianOperations{};
}

} // namespace afqmc
} // namespace qmcplusplus

#endif

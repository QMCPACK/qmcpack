//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   initially refactored from DiracDeterminantCUDA.cpp
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantBatched.cpp
 * @brief implementation of DiracDeterminantBatched with a S(ingle)P(article)O(rbital)SetBase
 */

#include "DiracDeterminantBatched.h"
#include "Numerics/CUDA/cuda_inverse.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Numerics/CUDA/cuda_inverse.h"
#include "Numerics/DeterminantOperators.h"
#include <unistd.h>

namespace qmcplusplus
{
DiracDeterminantBatched::DiracDeterminantBatched(SPOSet* const &spos, int first) :
  DiracDeterminantBase(first), Phi(dynamic_cast<SPOSetBatched*>(spos))
{
  registerTimers();
}

DiracDeterminantBase* DiracDeterminantBatched::makeCopy(SPOSet* spo) const
{
  DiracDeterminantBatched* dclone= new DiracDeterminantBatched(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}

DiracDeterminantBatched* DiracDeterminantBatched::makeCopy(SPOSetBatched* spo) const
{
  DiracDeterminantBatched* dclone= new DiracDeterminantBatched(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}

  
/////////////////////////////////////
// Vectorized evaluation functions //
/////////////////////////////////////

void CheckAlign (void *p, std::string var)
{
  if ((unsigned long)p & 63)
    app_error() << "CUDA alignment error for variable " << var << "!\n";
}



}

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DIRACDETERMINANT_TCC
#define QMCPLUSPLUS_DIRACDETERMINANT_TCC

#include "DiracDeterminant.h"

namespace qmcplusplus
{

// // this is a special case in which we can't easily make this function private to
// // the class, but rather need to make it "local" to this compilation unit using
// // an anonymous namespace so it's not associate with the template class
// namespace
// {
// template<class T>
// T do_ratioT(ParticleSet& P, int iat)
// {
//   UpdateMode             = ORB_PBYP_RATIO;
//   const int WorkingIndex = iat - FirstIndex;
//   assert(WorkingIndex >= 0);
//   {
//     ScopedTimer local_timer(SPOVTimer);
//     Phi->evaluateValue(P, iat, psiV);
//   }
//   {
//     ScopedTimer local_timer(RatioTimer);
//     // This is an optimization.
//     // check invRow_id against WorkingIndex to see if getInvRow() has been called
//     // This is intended to save redundant compuation in TM1 and TM3
//     if (invRow_id != WorkingIndex)
//     {
//       invRow_id = WorkingIndex;
//       updateEng.getInvRow(psiM, WorkingIndex, invRow);
//     }
//     curRatio = simd::dot(invRow.data(), psiV.data(), invRow.size());
//   }
//   return curRatio;
// }
// } // namespace

} // namespace qmcplusplus


#endif
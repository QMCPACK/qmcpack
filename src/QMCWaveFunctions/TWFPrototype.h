//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


/**@file TWFPrototype.h
 *@brief Declaration of TWFPrototype
 */
#ifndef QMCPLUSPLUS_TWFPROTOTYPE_H
#define QMCPLUSPLUS_TWFPROTOTYPE_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

class TWFPrototype
{
  public:
    using ValueMatrix_t = SPOSet::ValueMatrix_t;

    TWFPrototype();
  private:
  
  std::vector<SPOSet*> spos;
  //std::vector<IndexType> groups;
  std::vector<ValueMatrix_t> psiM;
  std::vector<ValueMatrix_t> psiMinv;
  std::vector<WaveFunctionComponent*> J;  

};

/**@}*/
} // namespace qmcplusplus
#endif

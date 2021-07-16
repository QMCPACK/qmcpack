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
#include "Configuration.h"
#include "Particle/ParticleSet.h"
namespace qmcplusplus
{

class TWFPrototype
{
  public:
    using ValueMatrix_t = SPOSet::ValueMatrix_t;
    using GradMatrix_t = SPOSet::GradMatrix_t;
    using IndexType = QMCTraits::IndexType;
    using RealType = QMCTraits::RealType;

    TWFPrototype();
    void add_determinant(const IndexType groupid, SPOSet* spo);
    inline void add_jastrow(WaveFunctionComponent* j){J.push_back(j);};

    void get_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mmat);
    void get_egrad_elapl_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec,
                           std::vector<GradMatrix_t>& gmat, std::vector<ValueMatrix_t>& lmat);

    RealType evaluateLog(ParticleSet& P);
       
    
  private:
  
  std::vector<SPOSet*> spos;
  std::vector<IndexType> groups;
  std::vector<ValueMatrix_t> psiM;
  std::vector<ValueMatrix_t> psiMinv;
  std::vector<WaveFunctionComponent*> J;  

};

/**@}*/
} // namespace qmcplusplus
#endif

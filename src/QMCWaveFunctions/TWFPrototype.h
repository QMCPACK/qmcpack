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
    using HessMatrix_t = SPOSet::HessMatrix_t;
    using IndexType = QMCTraits::IndexType;
    using RealType = QMCTraits::RealType;

    using ValueVector_t = SPOSet::ValueVector_t;
    using GradVector_t = SPOSet::GradVector_t;

    TWFPrototype();
    void add_determinant(const ParticleSet& P, const IndexType groupid, SPOSet* spo);
    inline void add_jastrow(WaveFunctionComponent* j){J.push_back(j);};

    void get_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mmat);
    IndexType get_M_row(const ParticleSet& P, IndexType iel, ValueVector_t& val);
    IndexType get_igrad_row(const ParticleSet& P, const ParticleSet& source, IndexType iel, IndexType iat_source, std::vector<ValueVector_t>& dval);
    void get_egrad_elapl_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec,
                           std::vector<GradMatrix_t>& gmat, std::vector<ValueMatrix_t>& lmat);
    void get_igrad_M(const ParticleSet& P, const ParticleSet& source, int iat, 
                               std::vector<std::vector<ValueMatrix_t> >& dmvec);  
    void get_igrad_igradelapl_M(const ParticleSet& P, const ParticleSet& source, int iat, 
                               std::vector<std::vector<ValueMatrix_t> >& dmvec,  
                               std::vector<std::vector<ValueMatrix_t> >& dlmat);

    RealType evaluateLog(ParticleSet& P);
       
    //Whatever the group is labelled as in the particle set, sid corresponds to the group used to create determinant sid.  
    //  so sid=4 returns the number of particles and orbitals used for det #4.  Assuming a multispecies determinantal wfn like
    //  Prod_i Det(M_i).  
    inline IndexType num_orbitals(const IndexType sid){return num_orbs[sid];};  
    inline IndexType num_particles(const IndexType sid){return num_ptcls[sid];};
    //This takes a particle set group index, and returns the "determinant" this refers to. 
    IndexType get_group_index(const IndexType gid);
  private:
 
  std::vector<IndexType> num_ptcls;
  std::vector<IndexType> num_orbs; 
  std::vector<SPOSet*> spos;
  std::vector<IndexType> groups;
  std::vector<ValueMatrix_t> psiM;
  std::vector<ValueMatrix_t> psiMinv;
  std::vector<WaveFunctionComponent*> J;  

};

/**@}*/
} // namespace qmcplusplus
#endif

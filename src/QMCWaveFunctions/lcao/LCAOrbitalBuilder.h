//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_SOA_LCAO_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{

  /** BasisSetBuilder using new LCAOrbitalSet and Soa versions
   *
   * Reimplement MolecularBasisSetBuilder
   * - support both CartesianTensor and SphericalTensor
   */
  class LCAOrbitalBuilder: public BasisSetBuilder
  {
    public:
    typedef RealBasisSetBase<RealType> BasisSet_t;
    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, xmlNodePtr cur);
    ~LCAOrbitalBuilder();
    bool put(xmlNodePtr cur);
    bool putXML(xmlNodePtr cur);
    bool putH5();
    SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);

    private:

    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///source ParticleSet
    ParticleSet& sourcePtcl;
    ///localized basis set
    BasisSet_t* myBasisSet;
    ///apply cusp correction to molecular orbitals
    int radialOrbType;
    bool cuspCorr;
    std::string cuspInfo;
    ///Path to HDF5 Wavefunction
    std::string h5_path;

    /** create basis set
     *
     * Use ao_traits<T,I,J> to match (ROT)x(SH) combo
     */
    template<int I, int J> BasisSet_t* createBasisSet(xmlNodePtr cur);
    template<int I, int J> BasisSet_t* createBasisSetH5();

  };
}
#endif

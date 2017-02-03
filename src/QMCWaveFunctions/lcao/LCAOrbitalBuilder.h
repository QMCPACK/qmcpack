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
    
    
#ifndef QMCPLUSPLUS_LCAO_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_LCAO_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"
#include "QMCWaveFunctions/lcao/SoaSphericalTensor.h"
#include "QMCWaveFunctions/lcao/SoaSphericalBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"

namespace qmcplusplus
{

  /** BasisSetBuilder using new LCAOrbitalSet and Soa versions
   *
   * Reimplement MolecularBasisSetBuilder
   * - support both CartesianTensor and SphericalTensor
   */
  struct LCAOrbitalBuilder: public BasisSetBuilder
  {
    //for now, use the same cubic spline: use BsplineFunctor later
    typedef OneDimCubicSpline<real_type> RadFuncT;
    typedef SoaSphericalBasisSet<RadFuncT,SoaCartesianTensor<RealType> > XYZCOT;
    typedef SoaSphericalBasisSet<RadFuncT,SoaSphericalTensor<RealType> > YlmCOT;
    typedef SoaLocalizedBasisSet<XYZCOT> XYZBasisT;
    typedef SoaLocalizedBasisSet<YlmCOT> YlmBasisT;

    typedef BasisSetBase<ValueType> BasisSet_t;
    /** constructor
     * \param els reference to the electrons
     * \param ions reference to the ions
     */
    LCAOrbitalBuilder(ParticleSet& els, ParticleSet& ions, bool cusp=false, std::string cusp_info=""):
      targetPtcl(els), sourcePtcl(ions), thisBasisSet(0),cuspCorr(cusp),cuspInfo(cusp_info)
      {
        ClassName="MolecularBasisBuilder";
      }

    inline bool is_same(const xmlChar* a, const char* b)
    {
      return !strcmp((const char*)a,b);
    }

    bool put(xmlNodePtr cur);

    SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);

    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///source ParticleSet
    ParticleSet& sourcePtcl;
    ///MO basis with CartesianTensor
    XYZBasisT *xyzBasisSet;
    ///MO basis with SphericalTensor
    YlmBasisT *ylmBasisSet;
    ///save AtomiBasisBuilder<RFB>*
    std::map<std::string,BasisSetBuilder*> aoBuilders;
    ///apply cusp correction to molecular orbitals
    bool cuspCorr;
    std::string cuspInfo;
  };
}
#endif

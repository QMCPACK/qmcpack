//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H

#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"
#include "CPU/math.hpp"


namespace qmcplusplus
{
struct EGOSet : public SPOSet
{
  int KptMax;
  std::vector<PosType> K;
  std::vector<RealType> mK2;

  EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2);
  EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2, const std::vector<int>& d);

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<EGOSet>(*this); }

  void resetParameters(const opt_variables_type& optVariables) override {}
  void setOrbitalSetSize(int norbs) override {}

  inline void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override
  {
    const PosType& r = P.activeR(iat);
    RealType sinkr, coskr;
    for (int ik = 0; ik < KptMax; ik++)
    {
      qmcplusplus::sincos(dot(K[ik], r), &sinkr, &coskr);
      psi[ik] = ValueType(coskr, sinkr);
    }
  }

  /** generic inline function to handle a row
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi value row
   * @param dpsi gradient row
   * @param d2psi laplacian row
   */
  inline void evaluateVGL(const ParticleSet& P,
                          int iat,
                          ValueVector& psi,
                          GradVector& dpsi,
                          ValueVector& d2psi) override
  {
    const PosType& r = P.activeR(iat);
    RealType sinkr, coskr;
    for (int ik = 0; ik < KptMax; ik++)
    {
      qmcplusplus::sincos(dot(K[ik], r), &sinkr, &coskr);
      psi[ik]   = ValueType(coskr, sinkr);
      dpsi[ik]  = ValueType(-sinkr, coskr) * K[ik];
      d2psi[ik] = ValueType(mK2[ik] * coskr, mK2[ik] * sinkr);
    }
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override
  {
    for (int iat = first, i = 0; iat < last; ++iat, ++i)
    {
      ValueVector v(logdet[i], OrbitalSetSize);
      GradVector g(dlogdet[i], OrbitalSetSize);
      ValueVector l(d2logdet[i], OrbitalSetSize);
      evaluateVGL(P, iat, v, g, l);
    }
  }
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) override
  {
    APP_ABORT(
        "Incomplete implementation EGOSet::evaluate(P,first,last,lodget,dlodet,grad_grad_logdet,grad_grad_grad_logdet");
  }
};


/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasComplexOrbitalBuilder : public WaveFunctionComponentBuilder
{
public:
  ///constructor
  ElectronGasComplexOrbitalBuilder(Communicate* comm, ParticleSet& els);

  ///implement vritual function
  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;
};

/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasSPOBuilder : public SPOSetBuilder
{
protected:
  bool has_twist;
  PosType unique_twist;
  HEGGrid<RealType> egGrid;
  xmlNodePtr spo_node;

public:
  ///constructor
  ElectronGasSPOBuilder(ParticleSet& p, Communicate* comm, xmlNodePtr cur);

  /** initialize the Antisymmetric wave function for electrons
  *@param cur the current xml node
  */
  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;
  std::unique_ptr<SPOSet> createSPOSetFromIndices(indices_t& indices);
};
} // namespace qmcplusplus
#endif

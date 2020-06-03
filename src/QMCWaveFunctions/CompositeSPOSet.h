//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COMPOSITE_SPOSET_H
#define QMCPLUSPLUS_COMPOSITE_SPOSET_H

#include <QMCWaveFunctions/SPOSet.h>
#include <QMCWaveFunctions/BasisSetBase.h>
#include <QMCWaveFunctions/SPOSetBuilder.h>

namespace qmcplusplus
{
class CompositeSPOSet : public SPOSet
{
public:
  ///component SPOSets
  std::vector<SPOSet*> components;
  ///temporary storage for values
  std::vector<ValueVector_t*> component_values;
  ///temporary storage for gradients
  std::vector<GradVector_t*> component_gradients;
  ///temporary storage for laplacians
  std::vector<ValueVector_t*> component_laplacians;
  ///store the precomputed offsets
  std::vector<int> component_offsets;

  CompositeSPOSet();
  ~CompositeSPOSet();

  ///add a sposet component to this composite sposet
  void add(SPOSet* component);

  ///print out component info
  void report();

  //SPOSet interface methods
  ///size is determined by component sposets and nothing else
  inline void setOrbitalSetSize(int norbs) {}

  void resetTargetParticleSet(ParticleSet& P);

  SPOSet* makeClone() const;

  /** add sposet clones from another Composite SPOSet
     *   should only be used in makeClone functions following shallow copy
     */
  void clone_from(const CompositeSPOSet& master);

  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi);

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  ///unimplemented functions call this to abort
  inline void not_implemented(const std::string& method)
  {
    APP_ABORT("CompositeSPOSet::" + method + " has not been implemented");
  }


  //methods to be implemented in the future (possibly)
  void resetParameters(const opt_variables_type& optVariables);
  void evaluate(const ParticleSet& P, PosType& r, ValueVector_t& psi);
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& ddlogdet);
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            HessMatrix_t& ddlogdet,
                            GGGMatrix_t& dddlogdet);
};

struct CompositeSPOSetBuilder : public SPOSetBuilder
{
  CompositeSPOSetBuilder(Communicate* comm) : SPOSetBuilder(comm) {}

  //SPOSetBuilder interface
  SPOSet* createSPOSetFromXML(xmlNodePtr cur);

  SPOSet* createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input);
};
} // namespace qmcplusplus

#endif

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


#include "CompositeSPOSet.h"
#include "Utilities/IteratorUtility.h"
#include <algorithm>
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
namespace MatrixOperators
{
/** copy a small matrix (N, M1) to a big matrix (N, M2), M2>M1
     * @param small input matrix
     * @param big outout matrix
     * @param offset_c column offset
     *
     * @todo smater and more efficient matrix, move up for others
     * The columns [0,M1) are inserted into [offset_c,offset_c+M1).
     */
template<typename MAT1, typename MAT2>
inline void insert_columns(const MAT1& small, MAT2& big, int offset_c)
{
  const int c = small.cols();
  for (int i = 0; i < small.rows(); ++i)
    std::copy(small[i], small[i] + c, big[i] + offset_c);
}
} // namespace MatrixOperators

CompositeSPOSet::CompositeSPOSet()
{
  className      = "CompositeSPOSet";
  OrbitalSetSize = 0;
  component_offsets.reserve(4);
}

CompositeSPOSet::CompositeSPOSet(const CompositeSPOSet& other) : SPOSet(other)
{
  for (auto& element : other.components)
  {
    this->add(element->makeClone());
  }
}

CompositeSPOSet::~CompositeSPOSet() = default;

void CompositeSPOSet::add(std::unique_ptr<SPOSet> component)
{
  if (components.empty())
    component_offsets.push_back(0); //add 0

  int norbs = component->size();
  components.push_back(std::move(component));
  component_values.emplace_back(norbs);
  component_gradients.emplace_back(norbs);
  component_laplacians.emplace_back(norbs);

  OrbitalSetSize += norbs;
  component_offsets.push_back(OrbitalSetSize);
}

void CompositeSPOSet::report()
{
  app_log() << "CompositeSPOSet" << std::endl;
  app_log() << "  ncomponents = " << components.size() << std::endl;
  app_log() << "  components" << std::endl;
  for (int i = 0; i < components.size(); ++i)
  {
    app_log() << "    " << i << std::endl;
    components[i]->basic_report("      ");
  }
}

std::unique_ptr<SPOSet> CompositeSPOSet::makeClone() const { return std::make_unique<CompositeSPOSet>(*this); }

void CompositeSPOSet::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  int n = 0;
  for (int c = 0; c < components.size(); ++c)
  {
    SPOSet& component   = *components[c];
    ValueVector& values = component_values[c];
    component.evaluateValue(P, iat, values);
    std::copy(values.begin(), values.end(), psi.begin() + n);
    n += component.size();
  }
}

void CompositeSPOSet::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  int n = 0;
  for (int c = 0; c < components.size(); ++c)
  {
    SPOSet& component       = *components[c];
    ValueVector& values     = component_values[c];
    GradVector& gradients   = component_gradients[c];
    ValueVector& laplacians = component_laplacians[c];
    component.evaluateVGL(P, iat, values, gradients, laplacians);
    std::copy(values.begin(), values.end(), psi.begin() + n);
    std::copy(gradients.begin(), gradients.end(), dpsi.begin() + n);
    std::copy(laplacians.begin(), laplacians.end(), d2psi.begin() + n);
    n += component.size();
  }
}

#ifdef QMC_CUDA
void CompositeSPOSet::evaluate(const ParticleSet& P, PosType& r, ValueVector& psi)
{
  not_implemented("evaluate(P,r,psi)");
}
#endif

//methods to be implemented later
void CompositeSPOSet::resetParameters(const opt_variables_type& optVariables)
{
  for (int c = 0; c < components.size(); ++c)
    components[c]->resetParameters(optVariables);
}

void CompositeSPOSet::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& logdet,
                                           GradMatrix& dlogdet,
                                           ValueMatrix& d2logdet)
{
  const int nat = last - first;
  for (int c = 0; c < components.size(); ++c)
  {
    int norb = components[c]->size();
    ValueMatrix v(nat, norb);
    GradMatrix g(nat, norb);
    ValueMatrix l(nat, norb);
    components[c]->evaluate_notranspose(P, first, last, v, g, l);
    int n = component_offsets[c];
    MatrixOperators::insert_columns(v, logdet, n);
    MatrixOperators::insert_columns(g, dlogdet, n);
    MatrixOperators::insert_columns(l, d2logdet, n);
  }
}

void CompositeSPOSet::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& logdet,
                                           GradMatrix& dlogdet,
                                           HessMatrix& grad_grad_logdet)
{
  const int nat = last - first;
  for (int c = 0; c < components.size(); ++c)
  {
    int norb = components[c]->size();
    ValueMatrix v(nat, norb);
    GradMatrix g(nat, norb);
    HessMatrix h(nat, norb);
    components[c]->evaluate_notranspose(P, first, last, v, g, h);
    int n = component_offsets[c];
    MatrixOperators::insert_columns(v, logdet, n);
    MatrixOperators::insert_columns(g, dlogdet, n);
    MatrixOperators::insert_columns(h, grad_grad_logdet, n);
  }
}

void CompositeSPOSet::evaluate_notranspose(const ParticleSet& P,
                                           int first,
                                           int last,
                                           ValueMatrix& logdet,
                                           GradMatrix& dlogdet,
                                           HessMatrix& grad_grad_logdet,
                                           GGGMatrix& grad_grad_grad_logdet)
{
  not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
}


std::unique_ptr<SPOSet> CompositeSPOSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  std::vector<std::string> spolist;
  putContent(spolist, cur);
  if (spolist.empty())
  {
    return nullptr;
  }

  auto spo_now = std::make_unique<CompositeSPOSet>();
  for (int i = 0; i < spolist.size(); ++i)
  {
    const SPOSet* spo = sposet_builder_factory_.getSPOSet(spolist[i]);
    if (spo)
      spo_now->add(spo->makeClone());
  }
  return (spo_now->size()) ? std::unique_ptr<SPOSet>{std::move(spo_now)} : nullptr;
}

std::unique_ptr<SPOSet> CompositeSPOSetBuilder::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input)
{
  return createSPOSetFromXML(cur);
}

} // namespace qmcplusplus

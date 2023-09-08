//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National
// Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign Mark A. Berrill,
//                    berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/CompositeSPOSetT.h"

#include "OhmmsData/AttributeSet.h"
#include "Utilities/IteratorUtility.h"

#include <algorithm>

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
template <typename MAT1, typename MAT2>
inline void
insert_columns(const MAT1& small, MAT2& big, int offset_c)
{
    const int c = small.cols();
    for (int i = 0; i < small.rows(); ++i)
        std::copy(small[i], small[i] + c, big[i] + offset_c);
}
} // namespace MatrixOperators

template <typename T>
CompositeSPOSetT<T>::CompositeSPOSetT(const std::string& my_name) :
    SPOSetT<T>(my_name)
{
    this->OrbitalSetSize = 0;
    component_offsets.reserve(4);
}

template <typename T>
CompositeSPOSetT<T>::CompositeSPOSetT(const CompositeSPOSetT<T>& other) :
    SPOSetT<T>(other)
{
    for (auto& element : other.components) {
        this->add(element->makeClone());
    }
}

template <typename T>
CompositeSPOSetT<T>::~CompositeSPOSetT() = default;

template <typename T>
void
CompositeSPOSetT<T>::add(std::unique_ptr<SPOSetT<T>> component)
{
    if (components.empty())
        component_offsets.push_back(0); // add 0

    int norbs = component->size();
    components.push_back(std::move(component));
    component_values.emplace_back(norbs);
    component_gradients.emplace_back(norbs);
    component_laplacians.emplace_back(norbs);

    this->OrbitalSetSize += norbs;
    component_offsets.push_back(this->OrbitalSetSize);
}

template <typename T>
void
CompositeSPOSetT<T>::report()
{
    app_log() << "CompositeSPOSetT" << std::endl;
    app_log() << "  ncomponents = " << components.size() << std::endl;
    app_log() << "  components" << std::endl;
    for (int i = 0; i < components.size(); ++i) {
        app_log() << "    " << i << std::endl;
        components[i]->basic_report("      ");
    }
}

template <typename T>
std::unique_ptr<SPOSetT<T>>
CompositeSPOSetT<T>::makeClone() const
{
    return std::make_unique<CompositeSPOSetT<T>>(*this);
}

template <typename T>
void
CompositeSPOSetT<T>::evaluateValue(
    const ParticleSetT<T>& P, int iat, ValueVector& psi)
{
    int n = 0;
    for (int c = 0; c < components.size(); ++c) {
        SPOSetT<T>& component = *components[c];
        ValueVector& values = component_values[c];
        component.evaluateValue(P, iat, values);
        std::copy(values.begin(), values.end(), psi.begin() + n);
        n += component.size();
    }
}

template <typename T>
void
CompositeSPOSetT<T>::evaluateVGL(const ParticleSetT<T>& P, int iat,
    ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
    int n = 0;
    for (int c = 0; c < components.size(); ++c) {
        SPOSetT<T>& component = *components[c];
        ValueVector& values = component_values[c];
        GradVector& gradients = component_gradients[c];
        ValueVector& laplacians = component_laplacians[c];
        component.evaluateVGL(P, iat, values, gradients, laplacians);
        std::copy(values.begin(), values.end(), psi.begin() + n);
        std::copy(gradients.begin(), gradients.end(), dpsi.begin() + n);
        std::copy(laplacians.begin(), laplacians.end(), d2psi.begin() + n);
        n += component.size();
    }
}

template <typename T>
void
CompositeSPOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet, ValueMatrix& d2logdet)
{
    const int nat = last - first;
    for (int c = 0; c < components.size(); ++c) {
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

template <typename T>
void
CompositeSPOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet,
    HessMatrix& grad_grad_logdet)
{
    const int nat = last - first;
    for (int c = 0; c < components.size(); ++c) {
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

template <typename T>
void
CompositeSPOSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet,
    HessMatrix& grad_grad_logdet, GGGMatrix& grad_grad_grad_logdet)
{
    not_implemented(
        "evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
}

template <typename T>
std::unique_ptr<SPOSetT<T>>
CompositeSPOSetBuilderT<T>::createSPOSetFromXML(xmlNodePtr cur)
{
    std::vector<std::string> spolist;
    putContent(spolist, cur);
    if (spolist.empty()) {
        return nullptr;
    }

    auto spo_now = std::make_unique<CompositeSPOSetT<T>>(
        getXMLAttributeValue(cur, "name"));
    for (int i = 0; i < spolist.size(); ++i) {
        const SPOSetT<T>* spo = sposet_builder_factory_.getSPOSet(spolist[i]);
        if (spo)
            spo_now->add(spo->makeClone());
    }
    return (spo_now->size()) ? std::unique_ptr<SPOSetT<T>>{std::move(spo_now)} :
                               nullptr;
}

template <typename T>
std::unique_ptr<SPOSetT<T>>
CompositeSPOSetBuilderT<T>::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input)
{
    return createSPOSetFromXML(cur);
}

// Class concrete types from ValueType

#ifndef QMC_COMPLEX
#ifndef MIXED_PRECISION
template class CompositeSPOSetT<double>;
template class CompositeSPOSetBuilderT<double>;
#else
template class CompositeSPOSetT<float>;
template class CompositeSPOSetBuilderT<float>;
#endif
#else
#ifndef MIXED_PRECISION
template class CompositeSPOSetT<std::complex<double>>;
template class CompositeSPOSetBuilderT<std::complex<double>>;
#else
template class CompositeSPOSetT<std::complex<float>>;
template class CompositeSPOSetBuilderT<std::complex<float>>;
#endif
#endif

} // namespace qmcplusplus

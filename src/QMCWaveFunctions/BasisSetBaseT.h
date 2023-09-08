//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at
// Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore
//                    National Laboratory Jeremy McMinnis, jmcminis@gmail.com,
//                    University of Illinois at Urbana-Champaign Jaron T.
//                    Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign Mark A. Berrill,
//                    berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_BASISSETBASET_H
#define QMCPLUSPLUS_BASISSETBASET_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{
/** base class for a basis set
 *
 * Define a common storage for the derived classes and
 * provides  a minimal set of interfaces to get/set BasisSetSize.
 */
template <typename T>
struct BasisSetBaseT : public OrbitalSetTraits<T>
{
    enum
    {
        MAXINDEX = 2 + OHMMS_DIM
    };
    using RealType = typename OrbitalSetTraits<T>::RealType;
    using ValueType = typename OrbitalSetTraits<T>::ValueType;
    using IndexType = typename OrbitalSetTraits<T>::IndexType;
    using HessType = typename OrbitalSetTraits<T>::HessType;
    using IndexVector = typename OrbitalSetTraits<T>::IndexVector;
    using ValueVector = typename OrbitalSetTraits<T>::ValueVector;
    using ValueMatrix = typename OrbitalSetTraits<T>::ValueMatrix;
    using GradVector = typename OrbitalSetTraits<T>::GradVector;
    using GradMatrix = typename OrbitalSetTraits<T>::GradMatrix;
    using HessVector = typename OrbitalSetTraits<T>::HessVector;
    using HessMatrix = typename OrbitalSetTraits<T>::HessMatrix;
    using GGGType = TinyVector<HessType, OHMMS_DIM>;
    using GGGVector = Vector<GGGType>;
    using GGGMatrix = Matrix<GGGType>;

    /// size of the basis set
    IndexType BasisSetSize;
    /// index of the particle
    IndexType ActivePtcl;
    /// counter to keep track
    unsigned long Counter;
    /// phi[i] the value of the i-th basis set
    ValueVector Phi;
    /// dphi[i] the gradient of the i-th basis set
    GradVector dPhi;
    /// d2phi[i] the laplacian of the i-th basis set
    ValueVector d2Phi;
    /// grad_grad_Phi[i] the full hessian of the i-th basis set
    HessVector grad_grad_Phi;
    /// grad_grad_grad_Phi the full hessian of the i-th basis set
    GGGVector grad_grad_grad_Phi;
    /// container to store value, laplacian and gradient
    ValueMatrix Temp;

    ValueMatrix Y;
    GradMatrix dY;
    ValueMatrix d2Y;

    /// default constructor
    BasisSetBaseT() : BasisSetSize(0), ActivePtcl(-1), Counter(0)
    {
    }
    /// virtual destructor
    virtual ~BasisSetBaseT()
    {
    }
    /** resize the container */
    void
    resize(int ntargets)
    {
        if (BasisSetSize) {
            Phi.resize(BasisSetSize);
            dPhi.resize(BasisSetSize);
            d2Phi.resize(BasisSetSize);
            grad_grad_Phi.resize(BasisSetSize);
            grad_grad_grad_Phi.resize(BasisSetSize);
            Temp.resize(BasisSetSize, MAXINDEX);
            Y.resize(ntargets, BasisSetSize);
            dY.resize(ntargets, BasisSetSize);
            d2Y.resize(ntargets, BasisSetSize);
        }
        else {
            app_error() << "  BasisSetBase::BasisSetSize == 0" << std::endl;
        }
    }

    /// clone the basis set
    virtual BasisSetBaseT*
    makeClone() const = 0;
    /** return the basis set size */
    inline IndexType
    getBasisSetSize() const
    {
        return BasisSetSize;
    }

    /// resize the basis set
    virtual void
    setBasisSetSize(int nbs) = 0;

    virtual void
    evaluateWithHessian(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateWithThirdDeriv(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateThirdDerivOnly(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateForWalkerMove(const ParticleSetT<T>& P) = 0;
    virtual void
    evaluateForWalkerMove(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateForPtclMove(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateAllForPtclMove(const ParticleSetT<T>& P, int iat) = 0;
    virtual void
    evaluateForPtclMoveWithHessian(const ParticleSetT<T>& P, int iat) = 0;
};

/** Base for real basis set
 *
 * Equivalent to BasisSetBase with minimum requirements
 * Used by LCAO
 */
template <typename T>
struct SoaBasisSetBaseT
{
    using value_type = T;
    using vgl_type = VectorSoaContainer<T, OHMMS_DIM + 2>;
    using vgh_type = VectorSoaContainer<T, 10>;
    using vghgh_type = VectorSoaContainer<T, 20>;
    using OffloadMWVGLArray =
        Array<T, 3, OffloadPinnedAllocator<T>>; // [VGL, walker, Orbs]
    using OffloadMWVArray =
        Array<T, 2, OffloadPinnedAllocator<T>>; // [walker, Orbs]

    /// size of the basis set
    int BasisSetSize;

    virtual ~SoaBasisSetBaseT() = default;
    inline int
    getBasisSetSize()
    {
        return BasisSetSize;
    }

    virtual SoaBasisSetBaseT<T>*
    makeClone() const = 0;
    virtual void
    setBasisSetSize(int nbs) = 0;

    // Evaluates value, gradient, and laplacian for electron "iat".  Parks them
    // into a temporary data structure "vgl".
    virtual void
    evaluateVGL(const ParticleSetT<T>& P, int iat, vgl_type& vgl) = 0;
    // Evaluates value, gradient, and laplacian for electron "iat".  places them
    // in a offload array for batched code.
    virtual void
    mw_evaluateVGL(const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        OffloadMWVGLArray& vgl) = 0;
    // Evaluates value for electron "iat".  places it in a offload array for
    // batched code.
    virtual void
    mw_evaluateValue(const RefVectorWithLeader<ParticleSetT<T>>& P_list,
        int iat, OffloadMWVArray& v) = 0;
    // Evaluates value, gradient, and Hessian for electron "iat".  Parks them
    // into a temporary data structure "vgh".
    virtual void
    evaluateVGH(const ParticleSetT<T>& P, int iat, vgh_type& vgh) = 0;
    // Evaluates value, gradient, and Hessian, and Gradient Hessian for electron
    // "iat".  Parks them into a temporary data structure "vghgh".
    virtual void
    evaluateVGHGH(const ParticleSetT<T>& P, int iat, vghgh_type& vghgh) = 0;
    // Evaluates the x,y, and z components of ionic gradient associated with
    // "jion" of value.  Parks the raw data into "vgl" container.
    virtual void
    evaluateGradSourceV(const ParticleSetT<T>& P, int iat,
        const ParticleSetT<T>& ions, int jion, vgl_type& vgl) = 0;
    // Evaluates the x,y, and z components of ionic gradient associated with
    // "jion" value, gradient, and laplacian.
    //     Parks the raw data into "vghgh" container.
    virtual void
    evaluateGradSourceVGL(const ParticleSetT<T>& P, int iat,
        const ParticleSetT<T>& ions, int jion, vghgh_type& vghgh) = 0;
    virtual void
    evaluateV(const ParticleSetT<T>& P, int iat, value_type* restrict vals) = 0;
    virtual bool
    is_S_orbital(int mo_idx, int ao_idx)
    {
        return false;
    }

    /// Determine which orbitals are S-type.  Used for cusp correction.
    virtual void
    queryOrbitalsForSType(const std::vector<bool>& corrCenter,
        std::vector<bool>& is_s_orbital) const
    {
    }

    /** initialize a shared resource and hand it to collection
     */
    virtual void createResource(ResourceCollection& collection) const {}

    /** acquire a shared resource from collection
     */
    virtual void acquireResource(ResourceCollection& collection,
                                 const RefVectorWithLeader<SoaBasisSetBaseT>& bset_list) const
    {}

    /** return a shared resource to collection
     */
    virtual void releaseResource(ResourceCollection& collection,
                                 const RefVectorWithLeader<SoaBasisSetBaseT>& bset_list) const
    {}
};

} // namespace qmcplusplus
#endif

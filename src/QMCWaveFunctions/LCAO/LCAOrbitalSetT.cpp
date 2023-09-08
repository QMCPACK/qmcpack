//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

#include "LCAOrbitalSetT.h"

#include "CPU/BLAS.hpp"
#include "Numerics/MatrixOperators.h"
#include <ResourceCollection.h>

namespace qmcplusplus
{

template <class T>
struct LCAOrbitalSetT<T>::LCAOMultiWalkerMem : public Resource
{
    LCAOMultiWalkerMem() : Resource("LCAOrbitalSetT")
    {
    }
    LCAOMultiWalkerMem(const LCAOMultiWalkerMem&) : LCAOMultiWalkerMem()
    {
    }

    std::unique_ptr<Resource>
    makeClone() const override
    {
        return std::make_unique<LCAOMultiWalkerMem>(*this);
    }

    OffloadMWVGLArray phi_vgl_v; // [5][NW][NumMO]
    OffloadMWVGLArray basis_mw; // [5][NW][NumAO]
    OffloadMWVArray phi_v; // [NW][NumMO]
    OffloadMWVArray basis_v_mw; // [NW][NumMO]
};

template <class T>
LCAOrbitalSetT<T>::LCAOrbitalSetT(
    const std::string& my_name, std::unique_ptr<basis_type>&& bs) :
    SPOSetT<T>(my_name),
    BasisSetSize(bs ? bs->getBasisSetSize() : 0),
    Identity(true),
    basis_timer_(createGlobalTimer("LCAOrbitalSetT::Basis", timer_level_fine)),
    mo_timer_(createGlobalTimer("LCAOrbitalSetT::MO", timer_level_fine))
{
    if (!bs)
        throw std::runtime_error(
            "LCAOrbitalSetT cannot take nullptr as its  basis set!");
    myBasisSet = std::move(bs);
    Temp.resize(BasisSetSize);
    Temph.resize(BasisSetSize);
    Tempgh.resize(BasisSetSize);
    this->OrbitalSetSize = BasisSetSize;
    LCAOrbitalSetT<T>::checkObject();
}

template <class T>
LCAOrbitalSetT<T>::LCAOrbitalSetT(const LCAOrbitalSetT<T>& in) :
    SPOSetT<T>(in),
    myBasisSet(in.myBasisSet->makeClone()),
    C(in.C),
    BasisSetSize(in.BasisSetSize),
    C_copy(in.C_copy),
    Identity(in.Identity),
    basis_timer_(in.basis_timer_),
    mo_timer_(in.mo_timer_)
{
    Temp.resize(BasisSetSize);
    Temph.resize(BasisSetSize);
    Tempgh.resize(BasisSetSize);
    if (!in.Identity) {
        Tempv.resize(this->OrbitalSetSize);
        Temphv.resize(this->OrbitalSetSize);
        Tempghv.resize(this->OrbitalSetSize);
    }
    LCAOrbitalSetT<T>::checkObject();
}

template <class T>
void
LCAOrbitalSetT<T>::setOrbitalSetSize(int norbs)
{
    if (C)
        throw std::runtime_error("LCAOrbitalSetT::setOrbitalSetSize cannot "
                                 "reset existing MO coefficients");

    Identity = false;
    this->OrbitalSetSize = norbs;
    C = std::make_shared<ValueMatrix>(this->OrbitalSetSize, BasisSetSize);
    Tempv.resize(this->OrbitalSetSize);
    Temphv.resize(this->OrbitalSetSize);
    Tempghv.resize(this->OrbitalSetSize);
    LCAOrbitalSetT<T>::checkObject();
}

template <class T>
void
LCAOrbitalSetT<T>::checkObject() const
{
    if (Identity) {
        if (this->OrbitalSetSize != BasisSetSize)
            throw std::runtime_error(
                "LCAOrbitalSetT::checkObject OrbitalSetSize and BasisSetSize "
                "must be equal if Identity = true!");
        if (C)
            throw std::runtime_error("LCAOrbitalSetT::checkObject C should be "
                                     "nullptr if Identity = true!");
    }
    else {
        if (!C)
            throw std::runtime_error("LCAOrbitalSetT::checkObject C should not "
                                     "be nullptr if Identity = false!");
        if (this->OrbitalSetSize != C->rows())
            throw std::runtime_error("LCAOrbitalSetT::checkObject C rows "
                                     "doesn't match OrbitalSetSize.");
        if (BasisSetSize != C->cols())
            throw std::runtime_error("LCAOrbitalSetT::checkObject C columns "
                                     "doesn't match BasisSetSize.");
    }
}

template <class T>
void
LCAOrbitalSetT<T>::createResource(ResourceCollection& collection) const
{
    myBasisSet->createResource(collection);
    auto resource_index =
        collection.addResource(std::make_unique<LCAOMultiWalkerMem>());
}

template <class T>
void
LCAOrbitalSetT<T>::acquireResource(ResourceCollection& collection,
    const RefVectorWithLeader<SPOSetT<T>>& spo_list) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    spo_leader.myBasisSet->acquireResource(collection, extractBasisRefList(spo_list));
    spo_leader.mw_mem_handle_ = collection.lendResource<LCAOMultiWalkerMem>();
}

template <class T>
void
LCAOrbitalSetT<T>::releaseResource(ResourceCollection& collection,
    const RefVectorWithLeader<SPOSetT<T>>& spo_list) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    spo_leader.myBasisSet->releaseResource(collection, extractBasisRefList(spo_list));
    collection.takebackResource(spo_leader.mw_mem_handle_);
}

template <class T>
RefVectorWithLeader<typename LCAOrbitalSetT<T>::basis_type> LCAOrbitalSetT<T>::extractBasisRefList(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list) const
{
    RefVectorWithLeader<basis_type> basis_list(*spo_list.template getCastedLeader<LCAOrbitalSetT<T>>().myBasisSet);
    basis_list.reserve(spo_list.size());
    for (size_t iw = 0; iw < spo_list.size(); iw++)
      basis_list.push_back(*spo_list.template getCastedElement<LCAOrbitalSetT<T>>(iw).myBasisSet);
    return basis_list;
}

template <class T>
std::unique_ptr<SPOSetT<T>>
LCAOrbitalSetT<T>::makeClone() const
{
    return std::make_unique<LCAOrbitalSetT<T>>(*this);
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateValue(
    const ParticleSetT<T>& P, int iat, ValueVector& psi)
{
    if (Identity) { // PAY ATTENTION TO COMPLEX
        myBasisSet->evaluateV(P, iat, psi.data());
    }
    else {
        Vector<T> vTemp(Temp.data(0), BasisSetSize);
        this->myBasisSet->evaluateV(P, iat, vTemp.data());
        assert(psi.size() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), psi.size(), BasisSetSize);
        MatrixOperators::product(C_partial_view, vTemp, psi);
    }
}

/** Find a better place for other user classes, Matrix should be padded as well
 */
template <typename T, unsigned D>
static void
Product_ABt(const VectorSoaContainer<T, D>& A, const Matrix<T>& B,
    VectorSoaContainer<T, D>& C)
{
    constexpr char transa = 't';
    constexpr char transb = 'n';
    constexpr T zone(1);
    constexpr T zero(0);
    BLAS::gemm(transa, transb, B.rows(), D, B.cols(), zone, B.data(), B.cols(),
        A.data(), A.capacity(), zero, C.data(), C.capacity());
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_vgl_impl(const vgl_type& temp, ValueVector& psi,
    GradVector& dpsi, ValueVector& d2psi) const
{
    const size_t output_size = psi.size();
    std::copy_n(temp.data(0), output_size, psi.data());
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    for (size_t j = 0; j < output_size; j++) {
        dpsi[j][0] = gx[j];
        dpsi[j][1] = gy[j];
        dpsi[j][2] = gz[j];
    }
    std::copy_n(temp.data(4), output_size, d2psi.data());
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_vgh_impl(const vgh_type& temp, ValueVector& psi,
    GradVector& dpsi, HessVector& d2psi) const
{
    const size_t output_size = psi.size();
    std::copy_n(temp.data(0), output_size, psi.data());
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    const T* restrict hxx = temp.data(4);
    const T* restrict hxy = temp.data(5);
    const T* restrict hxz = temp.data(6);
    const T* restrict hyy = temp.data(7);
    const T* restrict hyz = temp.data(8);
    const T* restrict hzz = temp.data(9);

    for (size_t j = 0; j < output_size; j++) {
        dpsi[j][0] = gx[j];
        dpsi[j][1] = gy[j];
        dpsi[j][2] = gz[j];

        d2psi[j](0, 0) = hxx[j];
        d2psi[j](0, 1) = d2psi[j](1, 0) = hxy[j];
        d2psi[j](0, 2) = d2psi[j](2, 0) = hxz[j];
        d2psi[j](1, 1) = hyy[j];
        d2psi[j](2, 1) = d2psi[j](1, 2) = hyz[j];
        d2psi[j](2, 2) = hzz[j];
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_vghgh_impl(const vghgh_type& temp, int i,
    ValueMatrix& psi, GradMatrix& dpsi, HessMatrix& d2psi,
    GGGMatrix& dghpsi) const
{
    const size_t output_size = psi.cols();
    std::copy_n(temp.data(0), output_size, psi[i]);
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    const T* restrict hxx = temp.data(4);
    const T* restrict hxy = temp.data(5);
    const T* restrict hxz = temp.data(6);
    const T* restrict hyy = temp.data(7);
    const T* restrict hyz = temp.data(8);
    const T* restrict hzz = temp.data(9);
    const T* restrict gh_xxx = temp.data(10);
    const T* restrict gh_xxy = temp.data(11);
    const T* restrict gh_xxz = temp.data(12);
    const T* restrict gh_xyy = temp.data(13);
    const T* restrict gh_xyz = temp.data(14);
    const T* restrict gh_xzz = temp.data(15);
    const T* restrict gh_yyy = temp.data(16);
    const T* restrict gh_yyz = temp.data(17);
    const T* restrict gh_yzz = temp.data(18);
    const T* restrict gh_zzz = temp.data(19);

    for (size_t j = 0; j < output_size; j++) {
        dpsi[i][j][0] = gx[j];
        dpsi[i][j][1] = gy[j];
        dpsi[i][j][2] = gz[j];

        d2psi[i][j](0, 0) = hxx[j];
        d2psi[i][j](0, 1) = d2psi[i][j](1, 0) = hxy[j];
        d2psi[i][j](0, 2) = d2psi[i][j](2, 0) = hxz[j];
        d2psi[i][j](1, 1) = hyy[j];
        d2psi[i][j](2, 1) = d2psi[i][j](1, 2) = hyz[j];
        d2psi[i][j](2, 2) = hzz[j];

        dghpsi[i][j][0](0, 0) = gh_xxx[j]; // x|xx
        dghpsi[i][j][0](0, 1) = gh_xxy[j]; // x|xy
        dghpsi[i][j][0](0, 2) = gh_xxz[j]; // x|xz
        dghpsi[i][j][0](1, 0) = gh_xxy[j]; // x|yx = xxy
        dghpsi[i][j][0](1, 1) = gh_xyy[j]; // x|yy
        dghpsi[i][j][0](1, 2) = gh_xyz[j]; // x|yz
        dghpsi[i][j][0](2, 0) = gh_xxz[j]; // x|zx = xxz
        dghpsi[i][j][0](2, 1) = gh_xyz[j]; // x|zy = xyz
        dghpsi[i][j][0](2, 2) = gh_xzz[j]; // x|zz

        dghpsi[i][j][1](0, 0) = gh_xxy[j]; // y|xx = xxy
        dghpsi[i][j][1](0, 1) = gh_xyy[j]; // y|xy = xyy
        dghpsi[i][j][1](0, 2) = gh_xyz[j]; // y|xz = xyz
        dghpsi[i][j][1](1, 0) = gh_xyy[j]; // y|yx = xyy
        dghpsi[i][j][1](1, 1) = gh_yyy[j]; // y|yy
        dghpsi[i][j][1](1, 2) = gh_yyz[j]; // y|yz
        dghpsi[i][j][1](2, 0) = gh_xyz[j]; // y|zx = xyz
        dghpsi[i][j][1](2, 1) = gh_yyz[j]; // y|zy = yyz
        dghpsi[i][j][1](2, 2) = gh_yzz[j]; // y|zz

        dghpsi[i][j][2](0, 0) = gh_xxz[j]; // z|xx = xxz
        dghpsi[i][j][2](0, 1) = gh_xyz[j]; // z|xy = xyz
        dghpsi[i][j][2](0, 2) = gh_xzz[j]; // z|xz = xzz
        dghpsi[i][j][2](1, 0) = gh_xyz[j]; // z|yx = xyz
        dghpsi[i][j][2](1, 1) = gh_yyz[j]; // z|yy = yyz
        dghpsi[i][j][2](1, 2) = gh_yzz[j]; // z|yz = yzz
        dghpsi[i][j][2](2, 0) = gh_xzz[j]; // z|zx = xzz
        dghpsi[i][j][2](2, 1) = gh_yzz[j]; // z|zy = yzz
        dghpsi[i][j][2](2, 2) = gh_zzz[j]; // z|zz
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_vghgh_impl(const vghgh_type& temp, ValueVector& psi,
    GradVector& dpsi, HessVector& d2psi, GGGVector& dghpsi) const
{
    const size_t output_size = psi.size();
    std::copy_n(temp.data(0), output_size, psi.data());
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    const T* restrict hxx = temp.data(4);
    const T* restrict hxy = temp.data(5);
    const T* restrict hxz = temp.data(6);
    const T* restrict hyy = temp.data(7);
    const T* restrict hyz = temp.data(8);
    const T* restrict hzz = temp.data(9);
    const T* restrict gh_xxx = temp.data(10);
    const T* restrict gh_xxy = temp.data(11);
    const T* restrict gh_xxz = temp.data(12);
    const T* restrict gh_xyy = temp.data(13);
    const T* restrict gh_xyz = temp.data(14);
    const T* restrict gh_xzz = temp.data(15);
    const T* restrict gh_yyy = temp.data(16);
    const T* restrict gh_yyz = temp.data(17);
    const T* restrict gh_yzz = temp.data(18);
    const T* restrict gh_zzz = temp.data(19);

    for (size_t j = 0; j < output_size; j++) {
        dpsi[j][0] = gx[j];
        dpsi[j][1] = gy[j];
        dpsi[j][2] = gz[j];

        d2psi[j](0, 0) = hxx[j];
        d2psi[j](0, 1) = d2psi[j](1, 0) = hxy[j];
        d2psi[j](0, 2) = d2psi[j](2, 0) = hxz[j];
        d2psi[j](1, 1) = hyy[j];
        d2psi[j](2, 1) = d2psi[j](1, 2) = hyz[j];
        d2psi[j](2, 2) = hzz[j];

        dghpsi[j][0](0, 0) = gh_xxx[j]; // x|xx
        dghpsi[j][0](0, 1) = gh_xxy[j]; // x|xy
        dghpsi[j][0](0, 2) = gh_xxz[j]; // x|xz
        dghpsi[j][0](1, 0) = gh_xxy[j]; // x|yx = xxy
        dghpsi[j][0](1, 1) = gh_xyy[j]; // x|yy
        dghpsi[j][0](1, 2) = gh_xyz[j]; // x|yz
        dghpsi[j][0](2, 0) = gh_xxz[j]; // x|zx = xxz
        dghpsi[j][0](2, 1) = gh_xyz[j]; // x|zy = xyz
        dghpsi[j][0](2, 2) = gh_xzz[j]; // x|zz

        dghpsi[j][1](0, 0) = gh_xxy[j]; // y|xx = xxy
        dghpsi[j][1](0, 1) = gh_xyy[j]; // y|xy = xyy
        dghpsi[j][1](0, 2) = gh_xyz[j]; // y|xz = xyz
        dghpsi[j][1](1, 0) = gh_xyy[j]; // y|yx = xyy
        dghpsi[j][1](1, 1) = gh_yyy[j]; // y|yy
        dghpsi[j][1](1, 2) = gh_yyz[j]; // y|yz
        dghpsi[j][1](2, 0) = gh_xyz[j]; // y|zx = xyz
        dghpsi[j][1](2, 1) = gh_xyy[j]; // y|xy = xyy
        dghpsi[j][1](2, 2) = gh_yzz[j]; // y|zz

        dghpsi[j][2](0, 0) = gh_xzz[j]; // z|xx = xzz
        dghpsi[j][2](0, 1) = gh_xyz[j]; // z|xy = xyz
        dghpsi[j][2](0, 2) = gh_xzz[j]; // z|xz = xzz
        dghpsi[j][2](1, 0) = gh_xyz[j]; // z|yx = xyz
        dghpsi[j][2](1, 1) = gh_yyz[j]; // z|yy = yyz
        dghpsi[j][2](1, 2) = gh_yzz[j]; // z|yz = yzz
        dghpsi[j][2](2, 0) = gh_xzz[j]; // z|zx = xzz
        dghpsi[j][2](2, 1) = gh_yzz[j]; // z|zy = yzz
        dghpsi[j][2](2, 2) = gh_zzz[j]; // z|zz
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_ionderiv_v_row_impl(
    const vgl_type& temp, GradVector& dpsi) const
{
    const size_t output_size = dpsi.size();
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);

    for (size_t j = 0; j < output_size; j++) {
        // As mentioned in SoaLocalizedBasisSet, LCAO's have a nice property
        // that
        //  for an atomic center, the ion gradient is the negative of the
        //  elecron gradient. Hence minus signs for each of these.
        dpsi[j][0] = -gx[j];
        dpsi[j][1] = -gy[j];
        dpsi[j][2] = -gz[j];
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateVGL(const ParticleSetT<T>& P, int iat,
    ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
    // TAKE CARE OF IDENTITY
    {
        ScopedTimer local(basis_timer_);
        myBasisSet->evaluateVGL(P, iat, Temp);
    }

    if (Identity)
        evaluate_vgl_impl(Temp, psi, dpsi, d2psi);
    else {
        assert(psi.size() <= this->OrbitalSetSize);
        {
            ScopedTimer local(mo_timer_);
            ValueMatrix C_partial_view(C->data(), psi.size(), BasisSetSize);
            Product_ABt(Temp, C_partial_view, Tempv);
        }
        evaluate_vgl_impl(Tempv, psi, dpsi, d2psi);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateVGL(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
    const RefVector<ValueVector>& psi_v_list,
    const RefVector<GradVector>& dpsi_v_list,
    const RefVector<ValueVector>& d2psi_v_list) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    auto& phi_vgl_v = spo_leader.mw_mem_handle_.getResource().phi_vgl_v;

    phi_vgl_v.resize(QMCTraits::DIM_VGL, spo_list.size(), this->OrbitalSetSize);
    mw_evaluateVGLImplGEMM(spo_list, P_list, iat, phi_vgl_v);

    const size_t nw = phi_vgl_v.size(1);

    // TODO: make this cleaner?
    for (int iw = 0; iw < nw; iw++) {
        const size_t output_size = psi_v_list[iw].get().size();
        std::copy_n(phi_vgl_v.data_at(0, iw, 0), output_size,
            psi_v_list[iw].get().data());
        std::copy_n(phi_vgl_v.data_at(4, iw, 0), output_size,
            d2psi_v_list[iw].get().data());
        // grads are [dim, walker, orb] in phi_vgl_v
        //           [walker][orb, dim] in dpsi_v_list
        for (size_t idim = 0; idim < QMCTraits::DIM; idim++)
            BLAS::copy(output_size, phi_vgl_v.data_at(idim + 1, iw, 0), 1,
                &dpsi_v_list[iw].get().data()[0][idim], QMCTraits::DIM);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateVGLImplGEMM(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
    OffloadMWVGLArray& phi_vgl_v) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    auto& basis_mw = spo_leader.mw_mem_handle_.getResource().basis_mw;
    basis_mw.resize(QMCTraits::DIM_VGL, spo_list.size(), BasisSetSize);

    {
        ScopedTimer local(basis_timer_);
        myBasisSet->mw_evaluateVGL(P_list, iat, basis_mw);
    }

    if (Identity) {
        // output_size can be smaller than BasisSetSize
        const size_t output_size = phi_vgl_v.size(2);
        const size_t nw = phi_vgl_v.size(1);

        for (size_t idim = 0; idim < QMCTraits::DIM_VGL; idim++)
            for (int iw = 0; iw < nw; iw++)
                std::copy_n(basis_mw.data_at(idim, iw, 0), output_size,
                    phi_vgl_v.data_at(idim, iw, 0));
    }
    else {
        const size_t requested_orb_size = phi_vgl_v.size(2);
        assert(requested_orb_size <= this->OrbitalSetSize);
        {
            ScopedTimer local(mo_timer_);
            ValueMatrix C_partial_view(
                C->data(), requested_orb_size, BasisSetSize);
            // TODO: make class for general blas interface in Platforms
            // have instance of that class as member of LCAOrbitalSetT, call
            // gemm through that
            BLAS::gemm('T', 'N',
                requested_orb_size, // MOs
                spo_list.size() * QMCTraits::DIM_VGL, // walkers * DIM_VGL
                BasisSetSize, // AOs
                1, C_partial_view.data(), BasisSetSize, basis_mw.data(),
                BasisSetSize, 0, phi_vgl_v.data(), requested_orb_size);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateValue(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
    const RefVector<ValueVector>& psi_v_list) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    auto& phi_v = spo_leader.mw_mem_handle_.getResource().phi_v;
    phi_v.resize(spo_list.size(), this->OrbitalSetSize);
    mw_evaluateValueImplGEMM(spo_list, P_list, iat, phi_v);

    const size_t output_size = phi_v.size(1);
    const size_t nw = phi_v.size(0);

    for (int iw = 0; iw < nw; iw++)
        std::copy_n(
            phi_v.data_at(iw, 0), output_size, psi_v_list[iw].get().data());
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateValueImplGEMM(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
    OffloadMWVArray& phi_v) const
{
    assert(this == &spo_list.getLeader());
    auto& spo_leader = spo_list.template getCastedLeader<LCAOrbitalSetT<T>>();
    const size_t nw = spo_list.size();
    auto& basis_v_mw = spo_leader.mw_mem_handle_.getResource().basis_v_mw;
    basis_v_mw.resize(nw, BasisSetSize);

    myBasisSet->mw_evaluateValue(P_list, iat, basis_v_mw);

    if (Identity) {
        std::copy_n(basis_v_mw.data_at(0, 0), this->OrbitalSetSize * nw,
            phi_v.data_at(0, 0));
    }
    else {
        const size_t requested_orb_size = phi_v.size(1);
        assert(requested_orb_size <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), requested_orb_size, BasisSetSize);
        BLAS::gemm('T', 'N',
            requested_orb_size, // MOs
            spo_list.size(), // walkers
            BasisSetSize, // AOs
            1, C_partial_view.data(), BasisSetSize, basis_v_mw.data(),
            BasisSetSize, 0, phi_v.data(), requested_orb_size);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateDetRatios(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<const VirtualParticleSetT<T>>& vp_list,
    const RefVector<ValueVector>& psi_list,
    const std::vector<const T*>& invRow_ptr_list,
    std::vector<std::vector<T>>& ratios_list) const
{
    const size_t nw = spo_list.size();
    for (size_t iw = 0; iw < nw; iw++) {
        for (size_t iat = 0; iat < vp_list[iw].getTotalNum(); iat++) {
            spo_list[iw].evaluateValue(vp_list[iw], iat, psi_list[iw]);
            ratios_list[iw][iat] = simd::dot(psi_list[iw].get().data(),
                invRow_ptr_list[iw], psi_list[iw].get().size());
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateDetRatios(const VirtualParticleSetT<T>& VP,
    ValueVector& psi, const ValueVector& psiinv, std::vector<T>& ratios)
{
    Vector<T> vTemp(Temp.data(0), BasisSetSize);
    Vector<T> invTemp(Temp.data(1), BasisSetSize);

    {
        ScopedTimer local(mo_timer_);
        // when only a subset of orbitals is used, extract limited rows of C.
        Matrix<T> C_occupied(C->data(), psiinv.size(), BasisSetSize);
        MatrixOperators::product_Atx(C_occupied, psiinv, invTemp);
    }

    for (size_t j = 0; j < VP.getTotalNum(); j++) {
        {
            ScopedTimer local(basis_timer_);
            myBasisSet->evaluateV(VP, j, vTemp.data());
        }
        ratios[j] = simd::dot(vTemp.data(), invTemp.data(), BasisSetSize);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::mw_evaluateVGLandDetRatioGrads(
    const RefVectorWithLeader<SPOSetT<T>>& spo_list,
    const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
    const std::vector<const T*>& invRow_ptr_list, OffloadMWVGLArray& phi_vgl_v,
    std::vector<T>& ratios, std::vector<GradType>& grads) const
{
    assert(this == &spo_list.getLeader());
    assert(phi_vgl_v.size(0) == QMCTraits::DIM_VGL);
    assert(phi_vgl_v.size(1) == spo_list.size());

    mw_evaluateVGLImplGEMM(spo_list, P_list, iat, phi_vgl_v);
    // Device data of phi_vgl_v must be up-to-date upon return
    phi_vgl_v.updateTo();

    const size_t nw = spo_list.size();
    const size_t norb_requested = phi_vgl_v.size(2);
    for (int iw = 0; iw < nw; iw++) {
        ratios[iw] = simd::dot(
            invRow_ptr_list[iw], phi_vgl_v.data_at(0, iw, 0), norb_requested);
        GradType dphi;
        for (size_t idim = 0; idim < QMCTraits::DIM; idim++)
            dphi[idim] =
                simd::dot(invRow_ptr_list[iw],
                    phi_vgl_v.data_at(idim + 1, iw, 0), norb_requested) /
                ratios[iw];
        grads[iw] = dphi;
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateVGH(const ParticleSetT<T>& P, int iat,
    ValueVector& psi, GradVector& dpsi, HessVector& dhpsi)
{
    // TAKE CARE OF IDENTITY
    myBasisSet->evaluateVGH(P, iat, Temph);
    if (Identity)
        evaluate_vgh_impl(Temph, psi, dpsi, dhpsi);
    else {
        assert(psi.size() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), psi.size(), BasisSetSize);
        Product_ABt(Temph, C_partial_view, Temphv);
        evaluate_vgh_impl(Temphv, psi, dpsi, dhpsi);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateVGHGH(const ParticleSetT<T>& P, int iat,
    ValueVector& psi, GradVector& dpsi, HessVector& dhpsi, GGGVector& dghpsi)
{
    // APP_ABORT("LCAORbitalSet::evaluate(psi,gpsi,hpsi,ghpsi) not
    // implemented\n");

    // TAKE CARE OF IDENTITY
    myBasisSet->evaluateVGHGH(P, iat, Tempgh);
    if (Identity)
        evaluate_vghgh_impl(Tempgh, psi, dpsi, dhpsi, dghpsi);
    else {
        assert(psi.size() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), psi.size(), BasisSetSize);
        Product_ABt(Tempgh, C_partial_view, Tempghv);
        evaluate_vghgh_impl(Tempghv, psi, dpsi, dhpsi, dghpsi);
    }
}

/* implement using gemm algorithm */
template <class T>
inline void
LCAOrbitalSetT<T>::evaluate_vgl_impl(const vgl_type& temp, int i,
    ValueMatrix& logdet, GradMatrix& dlogdet, ValueMatrix& d2logdet) const
{
    const size_t output_size = logdet.cols();
    std::copy_n(temp.data(0), output_size, logdet[i]);
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    for (size_t j = 0; j < output_size; j++) {
        dlogdet[i][j][0] = gx[j];
        dlogdet[i][j][1] = gy[j];
        dlogdet[i][j][2] = gz[j];
    }
    std::copy_n(temp.data(4), output_size, d2logdet[i]);
}
template <class T>
void
LCAOrbitalSetT<T>::evaluate_vgh_impl(const vgh_type& temp, int i,
    ValueMatrix& psi, GradMatrix& dpsi, HessMatrix& d2psi) const
{
    const size_t output_size = psi.cols();
    std::copy_n(temp.data(0), output_size, psi[i]);
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    const T* restrict hxx = temp.data(4);
    const T* restrict hxy = temp.data(5);
    const T* restrict hxz = temp.data(6);
    const T* restrict hyy = temp.data(7);
    const T* restrict hyz = temp.data(8);
    const T* restrict hzz = temp.data(9);

    for (size_t j = 0; j < output_size; j++) {
        dpsi[i][j][0] = gx[j];
        dpsi[i][j][1] = gy[j];
        dpsi[i][j][2] = gz[j];

        d2psi[i][j](0, 0) = hxx[j];
        d2psi[i][j](0, 1) = d2psi[i][j](1, 0) = hxy[j];
        d2psi[i][j](0, 2) = d2psi[i][j](2, 0) = hxz[j];
        d2psi[i][j](1, 1) = hyy[j];
        d2psi[i][j](2, 1) = d2psi[i][j](1, 2) = hyz[j];
        d2psi[i][j](2, 2) = hzz[j];
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_ionderiv_v_impl(
    const vgl_type& temp, int i, GradMatrix& dpsi) const
{
    const size_t output_size = dpsi.cols();
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);

    for (size_t j = 0; j < output_size; j++) {
        // As mentioned in SoaLocalizedBasisSet, LCAO's have a nice property
        // that
        //  for an atomic center, the ion gradient is the negative of the
        //  elecron gradient. Hence minus signs for each of these.
        dpsi[i][j][0] = -gx[j];
        dpsi[i][j][1] = -gy[j];
        dpsi[i][j][2] = -gz[j];
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_ionderiv_vgl_impl(const vghgh_type& temp, int i,
    GradMatrix& dpsi, HessMatrix& dgpsi, GradMatrix& dlpsi) const
{
    const size_t output_size = dpsi.cols();
    const T* restrict gx = temp.data(1);
    const T* restrict gy = temp.data(2);
    const T* restrict gz = temp.data(3);
    const T* restrict hxx = temp.data(4);
    const T* restrict hxy = temp.data(5);
    const T* restrict hxz = temp.data(6);
    const T* restrict hyy = temp.data(7);
    const T* restrict hyz = temp.data(8);
    const T* restrict hzz = temp.data(9);
    const T* restrict gh_xxx = temp.data(10);
    const T* restrict gh_xxy = temp.data(11);
    const T* restrict gh_xxz = temp.data(12);
    const T* restrict gh_xyy = temp.data(13);
    const T* restrict gh_xzz = temp.data(15);
    const T* restrict gh_yyy = temp.data(16);
    const T* restrict gh_yyz = temp.data(17);
    const T* restrict gh_yzz = temp.data(18);
    const T* restrict gh_zzz = temp.data(19);

    for (size_t j = 0; j < output_size; j++) {
        // As mentioned in SoaLocalizedBasisSet, LCAO's have a nice property
        // that
        //  for an atomic center, the ion gradient is the negative of the
        //  elecron gradient. Hence minus signs for each of these.
        dpsi[i][j][0] = -gx[j];
        dpsi[i][j][1] = -gy[j];
        dpsi[i][j][2] = -gz[j];

        dgpsi[i][j](0, 0) = -hxx[j];
        dgpsi[i][j](0, 1) = dgpsi[i][j](1, 0) = -hxy[j];
        dgpsi[i][j](0, 2) = dgpsi[i][j](2, 0) = -hxz[j];
        dgpsi[i][j](1, 1) = -hyy[j];
        dgpsi[i][j](2, 1) = dgpsi[i][j](1, 2) = -hyz[j];
        dgpsi[i][j](2, 2) = -hzz[j];

        // Since this returns the ion gradient of the laplacian, we have to
        // trace the grad hessian vector.
        dlpsi[i][j][0] = -(gh_xxx[j] + gh_xyy[j] + gh_xzz[j]);
        dlpsi[i][j][1] = -(gh_xxy[j] + gh_yyy[j] + gh_yzz[j]);
        dlpsi[i][j][2] = -(gh_xxz[j] + gh_yyz[j] + gh_zzz[j]);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet, ValueMatrix& d2logdet)
{
    if (Identity) {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGL(P, iat, Temp);
            evaluate_vgl_impl(Temp, i, logdet, dlogdet, d2logdet);
        }
    }
    else {
        assert(logdet.cols() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), logdet.cols(), BasisSetSize);
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGL(P, iat, Temp);
            Product_ABt(Temp, C_partial_view, Tempv);
            evaluate_vgl_impl(Tempv, i, logdet, dlogdet, d2logdet);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet,
    HessMatrix& grad_grad_logdet)
{
    if (Identity) {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGH(P, iat, Temph);
            evaluate_vgh_impl(Temph, i, logdet, dlogdet, grad_grad_logdet);
        }
    }
    else {
        assert(logdet.cols() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), logdet.cols(), BasisSetSize);
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGH(P, iat, Temph);
            Product_ABt(Temph, C_partial_view, Temphv);
            evaluate_vgh_impl(Temphv, i, logdet, dlogdet, grad_grad_logdet);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluate_notranspose(const ParticleSetT<T>& P, int first,
    int last, ValueMatrix& logdet, GradMatrix& dlogdet,
    HessMatrix& grad_grad_logdet, GGGMatrix& grad_grad_grad_logdet)
{
    if (Identity) {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGHGH(P, iat, Tempgh);
            evaluate_vghgh_impl(Tempgh, i, logdet, dlogdet, grad_grad_logdet,
                grad_grad_grad_logdet);
        }
    }
    else {
        assert(logdet.cols() <= this->OrbitalSetSize);
        ValueMatrix C_partial_view(C->data(), logdet.cols(), BasisSetSize);
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateVGHGH(P, iat, this->Tempgh);
            Product_ABt(this->Tempgh, C_partial_view, this->Tempghv);
            evaluate_vghgh_impl(this->Tempghv, i, logdet, dlogdet,
                grad_grad_logdet, grad_grad_grad_logdet);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateGradSource(const ParticleSetT<T>& P, int first,
    int last, const ParticleSetT<T>& source, int iat_src, GradMatrix& gradphi)
{
    if (Identity) {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateGradSourceV(
                P, iat, source, iat_src, this->Temp);
            evaluate_ionderiv_v_impl(Temp, i, gradphi);
        }
    }
    else {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateGradSourceV(
                P, iat, source, iat_src, this->Temp);
            Product_ABt(this->Temp, *C, this->Tempv);
            evaluate_ionderiv_v_impl(this->Tempv, i, gradphi);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateGradSource(const ParticleSetT<T>& P, int first,
    int last, const ParticleSetT<T>& source, int iat_src, GradMatrix& grad_phi,
    HessMatrix& grad_grad_phi, GradMatrix& grad_lapl_phi)
{
    if (Identity) {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateGradSourceVGL(
                P, iat, source, iat_src, this->Tempgh);
            evaluate_ionderiv_vgl_impl(
                this->Tempgh, i, grad_phi, grad_grad_phi, grad_lapl_phi);
        }
    }
    else {
        for (size_t i = 0, iat = first; iat < last; i++, iat++) {
            myBasisSet->evaluateGradSourceVGL(
                P, iat, source, iat_src, this->Tempgh);
            Product_ABt(this->Tempgh, *C, this->Tempghv);
            evaluate_ionderiv_vgl_impl(
                this->Tempghv, i, grad_phi, grad_grad_phi, grad_lapl_phi);
        }
    }
}

template <class T>
void
LCAOrbitalSetT<T>::evaluateGradSourceRow(const ParticleSetT<T>& P, int iel,
    const ParticleSetT<T>& source, int iat_src, GradVector& gradphi)
{
    if (Identity) {
        myBasisSet->evaluateGradSourceV(P, iel, source, iat_src, this->Temp);
        evaluate_ionderiv_v_row_impl(this->Temp, gradphi);
    }
    else {
        myBasisSet->evaluateGradSourceV(P, iel, source, iat_src, this->Temp);
        Product_ABt(Temp, *C, this->Tempv);
        evaluate_ionderiv_v_row_impl(this->Tempv, gradphi);
    }
}

template <class T>
void
LCAOrbitalSetT<T>::applyRotation(
    const ValueMatrix& rot_mat, bool use_stored_copy)
{
    if (!use_stored_copy)
        *C_copy = *C;
    // gemm is out-of-place
    BLAS::gemm('N', 'T', BasisSetSize, this->OrbitalSetSize,
        this->OrbitalSetSize, RealType(1.0), C_copy->data(), BasisSetSize,
        rot_mat.data(), this->OrbitalSetSize, RealType(0.0), C->data(),
        BasisSetSize);

    /* debugging code
    app_log() << "PRINTING MO COEFFICIENTS AFTER ROTATION " << objectName <<
    std::endl; for (int j = 0; j < OrbitalSetSize; j++) for (int i = 0; i <
    BasisSetSize; i++)
      {
        app_log() << " " << std::right << std::fixed << std::setprecision(16) <<
    std::setw(23) << std::scientific
                  << *(C->data() + j * BasisSetSize + i);

        if ((j * BasisSetSize + i + 1) % 4 == 0)
          app_log() << std::endl;
      }
    */
}

// Class concrete types from ValueType
#ifndef QMC_COMPLEX
template class LCAOrbitalSetT<double>;
template class LCAOrbitalSetT<float>;
#else
template class LCAOrbitalSetT<std::complex<double>>;
template class LCAOrbitalSetT<std::complex<float>>;
#endif

} // namespace qmcplusplus

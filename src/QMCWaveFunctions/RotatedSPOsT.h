//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source
/// License. / See LICENSE file in top directory for details.
////
//// Copyright (c) QMCPACK developers.
////
//// File developed by: Sergio D. Pineda Flores,
/// sergio_pinedaflores@berkeley.edu, University of California, Berkeley / Eric
/// Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley /
/// Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////
//// File created by: Sergio D. Pineda Flores, sergio_pinedaflores@berkeley.edu,
/// University of California, Berkeley
////////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_ROTATEDSPOST_H
#define QMCPLUSPLUS_ROTATEDSPOST_H

#include "QMCWaveFunctions/SPOSetT.h"

namespace qmcplusplus
{
template <typename T>
class RotatedSPOsT;
namespace testing
{
OptVariablesType<double>&
getMyVarsFull(RotatedSPOsT<double>& rot);
OptVariablesType<float>&
getMyVarsFull(RotatedSPOsT<float>& rot);
std::vector<std::vector<double>>&
getHistoryParams(RotatedSPOsT<double>& rot);
std::vector<std::vector<float>>&
getHistoryParams(RotatedSPOsT<float>& rot);
} // namespace testing

template <class T>
class RotatedSPOsT : public SPOSetT<T>, public OptimizableObjectT<T>
{
public:
    using IndexType = typename SPOSetT<T>::IndexType;
    using RealType = typename SPOSetT<T>::RealType;
    using ValueType = typename SPOSetT<T>::ValueType;
    using GradType = typename SPOSetT<T>::GradType;
    using ComplexType = typename SPOSetT<T>::ComplexType;
    using FullRealType = typename SPOSetT<T>::FullRealType;
    using ValueVector = typename SPOSetT<T>::ValueVector;
    using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
    using GradVector = typename SPOSetT<T>::GradVector;
    using GradMatrix = typename SPOSetT<T>::GradMatrix;
    using HessVector = typename SPOSetT<T>::HessVector;
    using HessMatrix = typename SPOSetT<T>::HessMatrix;
    using GGGVector = typename SPOSetT<T>::GGGVector;
    using GGGMatrix = typename SPOSetT<T>::GGGMatrix;
    using OffloadMWVGLArray = typename SPOSetT<T>::OffloadMWVGLArray;
    template <typename DT>
    using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;

    // constructor
    RotatedSPOsT(
        const std::string& my_name, std::unique_ptr<SPOSetT<T>>&& spos);
    // destructor
    ~RotatedSPOsT() override;

    std::string
    getClassName() const override
    {
        return "RotatedSPOsT";
    }
    bool
    isOptimizable() const override
    {
        return true;
    }
    bool
    isOMPoffload() const override
    {
        return Phi->isOMPoffload();
    }
    bool
    hasIonDerivs() const override
    {
        return Phi->hasIonDerivs();
    }

    // Vector of rotation matrix indices
    using RotationIndices = std::vector<std::pair<int, int>>;

    // Active orbital rotation parameter indices
    RotationIndices m_act_rot_inds;

    // Full set of rotation values for global rotation
    RotationIndices m_full_rot_inds;

    // Construct a list of the matrix indices for non-zero rotation parameters.
    // (The structure for a sparse representation of the matrix)
    // Only core->active rotations are created.
    static void
    createRotationIndices(int nel, int nmo, RotationIndices& rot_indices);

    // Construct a list for all the matrix indices, including core->active,
    // core->core and active->active
    static void
    createRotationIndicesFull(int nel, int nmo, RotationIndices& rot_indices);

    // Fill in antisymmetric matrix from the list of rotation parameter indices
    // and a list of parameter values.
    // This function assumes rot_mat is properly sized upon input and is set to
    // zero.
    static void
    constructAntiSymmetricMatrix(const RotationIndices& rot_indices,
        const std::vector<RealType>& param, ValueMatrix& rot_mat);

    // Extract the list of rotation parameters from the entries in an
    // antisymmetric matrix This function expects rot_indices and param are the
    // same length.
    static void
    extractParamsFromAntiSymmetricMatrix(const RotationIndices& rot_indices,
        const ValueMatrix& rot_mat, std::vector<RealType>& param);

    // function to perform orbital rotations
    void
    apply_rotation(const std::vector<RealType>& param, bool use_stored_copy);

    // For global rotation, inputs are the old parameters and the delta
    // parameters. The corresponding rotation matrices are constructed,
    // multiplied together, and the new parameters extracted. The new rotation
    // is applied to the underlying SPO coefficients
    void
    applyDeltaRotation(const std::vector<RealType>& delta_param,
        const std::vector<RealType>& old_param,
        std::vector<RealType>& new_param);

    // Perform the construction of matrices and extraction of parameters for a
    // delta rotation. Split out and made static for testing.
    static void
    constructDeltaRotation(const std::vector<RealType>& delta_param,
        const std::vector<RealType>& old_param,
        const RotationIndices& act_rot_inds,
        const RotationIndices& full_rot_inds, std::vector<RealType>& new_param,
        ValueMatrix& new_rot_mat);

    // When initializing the rotation from VP files
    // This function applies the rotation history
    void
    applyRotationHistory();

    // This function applies the global rotation (similar to apply_rotation, but
    // for the full set of rotation parameters)
    void
    applyFullRotation(
        const std::vector<RealType>& full_param, bool use_stored_copy);

    // Compute matrix exponential of an antisymmetric matrix (result is rotation
    // matrix)
    static void
    exponentiate_antisym_matrix(ValueMatrix& mat);

    // Compute matrix log of rotation matrix to produce antisymmetric matrix
    static void
    log_antisym_matrix(const ValueMatrix& mat, ValueMatrix& output);

    // A particular SPOSet used for Orbitals
    std::unique_ptr<SPOSetT<T>> Phi;

    /// Set the rotation parameters (usually from input file)
    void
    setRotationParameters(const std::vector<RealType>& param_list);

    /// the number of electrons of the majority spin
    size_t nel_major_;

    std::unique_ptr<SPOSetT<T>>
    makeClone() const override;

    // myG_temp (myL_temp) is the Gradient (Laplacian) value of of the
    // Determinant part of the wfn myG_J is the Gradient of the all other parts
    // of the wavefunction (typically just the Jastrow).
    //       It represents \frac{\nabla\psi_{J}}{\psi_{J}}
    // myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}} . The
    // Laplacian portion IMPORTANT NOTE:  The value of P.L holds \nabla^2
    // ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J
    // will hold
    typename ParticleSetT<T>::ParticleGradient myG_temp, myG_J;
    typename ParticleSetT<T>::ParticleLaplacian myL_temp, myL_J;

    ValueMatrix Bbar;
    ValueMatrix psiM_inv;
    ValueMatrix psiM_all;
    GradMatrix dpsiM_all;
    ValueMatrix d2psiM_all;

    // Single Slater creation
    void
    buildOptVariables(size_t nel);

    // For the MSD case rotations must be created in MultiSlaterDetTableMethod
    // class
    void
    buildOptVariables(const RotationIndices& rotations,
        const RotationIndices& full_rotations);

    void
    evaluateDerivatives(ParticleSetT<T>& P, const OptVariablesType<T>& optvars,
        Vector<T>& dlogpsi, Vector<T>& dhpsioverpsi, const int& FirstIndex,
        const int& LastIndex) override;

    void
    evaluateDerivativesWF(ParticleSetT<T>& P,
        const OptVariablesType<T>& optvars, Vector<T>& dlogpsi, int FirstIndex,
        int LastIndex) override;

    void
    evaluateDerivatives(ParticleSetT<T>& P, const OptVariablesType<T>& optvars,
        Vector<T>& dlogpsi, Vector<T>& dhpsioverpsi, const T& psiCurrent,
        const std::vector<T>& Coeff, const std::vector<size_t>& C2node_up,
        const std::vector<size_t>& C2node_dn, const ValueVector& detValues_up,
        const ValueVector& detValues_dn, const GradMatrix& grads_up,
        const GradMatrix& grads_dn, const ValueMatrix& lapls_up,
        const ValueMatrix& lapls_dn, const ValueMatrix& M_up,
        const ValueMatrix& M_dn, const ValueMatrix& Minv_up,
        const ValueMatrix& Minv_dn, const GradMatrix& B_grad,
        const ValueMatrix& B_lapl, const std::vector<int>& detData_up,
        const size_t N1, const size_t N2, const size_t NP1, const size_t NP2,
        const std::vector<std::vector<int>>& lookup_tbl) override;

    void
    evaluateDerivativesWF(ParticleSetT<T>& P,
        const OptVariablesType<T>& optvars, Vector<ValueType>& dlogpsi,
        const ValueType& psiCurrent, const std::vector<ValueType>& Coeff,
        const std::vector<size_t>& C2node_up,
        const std::vector<size_t>& C2node_dn, const ValueVector& detValues_up,
        const ValueVector& detValues_dn, const ValueMatrix& M_up,
        const ValueMatrix& M_dn, const ValueMatrix& Minv_up,
        const ValueMatrix& Minv_dn, const std::vector<int>& detData_up,
        const std::vector<std::vector<int>>& lookup_tbl) override;

    // helper function to evaluatederivative; evaluate orbital rotation
    // parameter derivative using table method
    void
    table_method_eval(Vector<T>& dlogpsi, Vector<T>& dhpsioverpsi,
        const typename ParticleSetT<T>::ParticleLaplacian& myL_J,
        const typename ParticleSetT<T>::ParticleGradient& myG_J,
        const size_t nel, const size_t nmo, const T& psiCurrent,
        const std::vector<T>& Coeff, const std::vector<size_t>& C2node_up,
        const std::vector<size_t>& C2node_dn, const ValueVector& detValues_up,
        const ValueVector& detValues_dn, const GradMatrix& grads_up,
        const GradMatrix& grads_dn, const ValueMatrix& lapls_up,
        const ValueMatrix& lapls_dn, const ValueMatrix& M_up,
        const ValueMatrix& M_dn, const ValueMatrix& Minv_up,
        const ValueMatrix& Minv_dn, const GradMatrix& B_grad,
        const ValueMatrix& B_lapl, const std::vector<int>& detData_up,
        const size_t N1, const size_t N2, const size_t NP1, const size_t NP2,
        const std::vector<std::vector<int>>& lookup_tbl);

    void
    table_method_evalWF(Vector<T>& dlogpsi, const size_t nel, const size_t nmo,
        const T& psiCurrent, const std::vector<T>& Coeff,
        const std::vector<size_t>& C2node_up,
        const std::vector<size_t>& C2node_dn, const ValueVector& detValues_up,
        const ValueVector& detValues_dn, const ValueMatrix& M_up,
        const ValueMatrix& M_dn, const ValueMatrix& Minv_up,
        const ValueMatrix& Minv_dn, const std::vector<int>& detData_up,
        const std::vector<std::vector<int>>& lookup_tbl);

    void
    extractOptimizableObjectRefs(UniqueOptObjRefsT<T>& opt_obj_refs) override
    {
        opt_obj_refs.push_back(*this);
    }

    void
    checkInVariablesExclusive(OptVariablesType<T>& active) override
    {
        if (this->myVars.size())
            active.insertFrom(this->myVars);
    }

    void
    checkOutVariables(const OptVariablesType<T>& active) override
    {
        this->myVars.getIndex(active);
    }

    /// reset
    void
    resetParametersExclusive(const OptVariablesType<T>& active) override;

    void
    writeVariationalParameters(hdf_archive& hout) override;

    void
    readVariationalParameters(hdf_archive& hin) override;

    //*********************************************************************************
    // the following functions simply call Phi's corresponding functions
    void
    setOrbitalSetSize(int norbs) override
    {
        Phi->setOrbitalSetSize(norbs);
    }

    void
    checkObject() const override
    {
        Phi->checkObject();
    }

    void
    evaluateValue(const ParticleSetT<T>& P, int iat, ValueVector& psi) override
    {
        assert(psi.size() <= this->OrbitalSetSize);
        Phi->evaluateValue(P, iat, psi);
    }

    void
    evaluateVGL(const ParticleSetT<T>& P, int iat, ValueVector& psi,
        GradVector& dpsi, ValueVector& d2psi) override
    {
        assert(psi.size() <= this->OrbitalSetSize);
        Phi->evaluateVGL(P, iat, psi, dpsi, d2psi);
    }

    void
    evaluateDetRatios(const VirtualParticleSetT<T>& VP, ValueVector& psi,
        const ValueVector& psiinv, std::vector<T>& ratios) override
    {
        Phi->evaluateDetRatios(VP, psi, psiinv, ratios);
    }

    void
    evaluateDerivRatios(const VirtualParticleSetT<T>& VP,
        const OptVariablesType<T>& optvars, ValueVector& psi,
        const ValueVector& psiinv, std::vector<T>& ratios, Matrix<T>& dratios,
        int FirstIndex, int LastIndex) override;

    void
    evaluateVGH(const ParticleSetT<T>& P, int iat, ValueVector& psi,
        GradVector& dpsi, HessVector& grad_grad_psi) override
    {
        assert(psi.size() <= this->OrbitalSetSize);
        Phi->evaluateVGH(P, iat, psi, dpsi, grad_grad_psi);
    }

    void
    evaluateVGHGH(const ParticleSetT<T>& P, int iat, ValueVector& psi,
        GradVector& dpsi, HessVector& grad_grad_psi,
        GGGVector& grad_grad_grad_psi) override
    {
        Phi->evaluateVGHGH(
            P, iat, psi, dpsi, grad_grad_psi, grad_grad_grad_psi);
    }

    void
    evaluate_notranspose(const ParticleSetT<T>& P, int first, int last,
        ValueMatrix& logdet, GradMatrix& dlogdet,
        ValueMatrix& d2logdet) override
    {
        Phi->evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
    }

    void
    evaluate_notranspose(const ParticleSetT<T>& P, int first, int last,
        ValueMatrix& logdet, GradMatrix& dlogdet,
        HessMatrix& grad_grad_logdet) override
    {
        Phi->evaluate_notranspose(
            P, first, last, logdet, dlogdet, grad_grad_logdet);
    }

    void
    evaluate_notranspose(const ParticleSetT<T>& P, int first, int last,
        ValueMatrix& logdet, GradMatrix& dlogdet, HessMatrix& grad_grad_logdet,
        GGGMatrix& grad_grad_grad_logdet) override
    {
        Phi->evaluate_notranspose(P, first, last, logdet, dlogdet,
            grad_grad_logdet, grad_grad_grad_logdet);
    }

    void
    evaluateGradSource(const ParticleSetT<T>& P, int first, int last,
        const ParticleSetT<T>& source, int iat_src,
        GradMatrix& grad_phi) override
    {
        Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi);
    }

    void
    evaluateGradSource(const ParticleSetT<T>& P, int first, int last,
        const ParticleSetT<T>& source, int iat_src, GradMatrix& grad_phi,
        HessMatrix& grad_grad_phi, GradMatrix& grad_lapl_phi) override
    {
        Phi->evaluateGradSource(P, first, last, source, iat_src, grad_phi,
            grad_grad_phi, grad_lapl_phi);
    }

    //  void evaluateThirdDeriv(const ParticleSet& P, int first, int last,
    //  GGGMatrix& grad_grad_grad_logdet) {Phi->evaluateThridDeriv(P, first,
    //  last, grad_grad_grad_logdet); }

    /// Use history list (false) or global rotation (true)
    void
    set_use_global_rotation(bool use_global_rotation)
    {
        use_global_rot_ = use_global_rotation;
    }

    void
    mw_evaluateDetRatios(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<const VirtualParticleSetT<T>>& vp_list,
        const RefVector<ValueVector>& psi_list,
        const std::vector<const ValueType*>& invRow_ptr_list,
        std::vector<std::vector<ValueType>>& ratios_list) const override;

    void
    mw_evaluateValue(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        const RefVector<ValueVector>& psi_v_list) const override;

    void
    mw_evaluateVGL(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        const RefVector<ValueVector>& psi_v_list,
        const RefVector<GradVector>& dpsi_v_list,
        const RefVector<ValueVector>& d2psi_v_list) const override;

    void
    mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        const RefVector<ValueVector>& psi_v_list,
        const RefVector<GradVector>& dpsi_v_list,
        const RefVector<ValueVector>& d2psi_v_list,
        OffloadMatrix<ComplexType>& mw_dspin) const override;

    void
    mw_evaluateVGLandDetRatioGrads(
        const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        const std::vector<const ValueType*>& invRow_ptr_list,
        OffloadMWVGLArray& phi_vgl_v, std::vector<ValueType>& ratios,
        std::vector<GradType>& grads) const override;

    void
    mw_evaluateVGLandDetRatioGradsWithSpin(
        const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int iat,
        const std::vector<const ValueType*>& invRow_ptr_list,
        OffloadMWVGLArray& phi_vgl_v, std::vector<ValueType>& ratios,
        std::vector<GradType>& grads,
        std::vector<ValueType>& spingrads) const override;

    void
    mw_evaluate_notranspose(const RefVectorWithLeader<SPOSetT<T>>& spo_list,
        const RefVectorWithLeader<ParticleSetT<T>>& P_list, int first, int last,
        const RefVector<ValueMatrix>& logdet_list,
        const RefVector<GradMatrix>& dlogdet_list,
        const RefVector<ValueMatrix>& d2logdet_list) const override;

    void
    createResource(ResourceCollection& collection) const override;

    void
    acquireResource(ResourceCollection& collection,
        const RefVectorWithLeader<SPOSetT<T>>& spo_list) const override;

    void
    releaseResource(ResourceCollection& collection,
        const RefVectorWithLeader<SPOSetT<T>>& spo_list) const override;

private:
    /// true if SPO parameters (orbital rotation parameters) have been supplied
    /// by input
    bool params_supplied;
    /// list of supplied orbital rotation parameters
    std::vector<RealType> params;

    /// Full set of rotation matrix parameters for use in global rotation method
    OptVariablesType<T> myVarsFull;

    /// List of previously applied parameters
    std::vector<std::vector<RealType>> history_params_;

    /// Use global rotation or history list
    bool use_global_rot_ = true;

    static RefVectorWithLeader<SPOSetT<T>>
    extractPhiRefList(const RefVectorWithLeader<SPOSetT<T>>& spo_list);

    friend OptVariablesType<double>&
    testing::getMyVarsFull(RotatedSPOsT<double>& rot);
    friend OptVariablesType<float>&
    testing::getMyVarsFull(RotatedSPOsT<float>& rot);
    friend std::vector<std::vector<double>>&
    testing::getHistoryParams(RotatedSPOsT<double>& rot);
    friend std::vector<std::vector<float>>&
    testing::getHistoryParams(RotatedSPOsT<float>& rot);
};

} // namespace qmcplusplus

#endif

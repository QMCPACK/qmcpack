//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Nick Blunt, nicksblunt@gmail.com, University of Cambridge
//
// File created by: Nick Blunt, nicksblunt@gmail.com, University of Cambridge
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FDLRWFN_H
#define QMCPLUSPLUS_FDLRWFN_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include <QMCWaveFunctions/OrbitalBase.h>
#include <OhmmsPETE/TinyVector.h>

namespace qmcplusplus {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  A class for representing a linear response wave function as a finite difference
  ///         between two wave functions of the same ansatz but with different variables.
  ///
  ///    Psi_fdlr = Psi(x1 + d1, x2 + d2, ...) - Psi(x1 - d1, x2 - d2, ...)
  ///
  ///  which is (with an error that is 3rd order in d) equal to the linear response expansion
  ///
  ///    sum_i  2 * di * (dPsi/dxi)
  ///
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class FDLRWfn : public OrbitalBase {

  // protected data members
  private:

    // List of variables for the "x" parameters
    opt_variables_type x_vars;

    // List of variables for the delta or "d" parameters
    opt_variables_type d_vars;

    // List of variables for the "x" parameters
    // The difference between this and x_vars is that the names for the
    // variables in x_vars_driver have "x_" appended to the start. This is how
    // they are stored in the full list of optimizable parameters held by the
    // QMCDriver objects. It is helpful to have x_vars  and d_vars without this
    // naming difference, for the sake of adding "x" and "d" variables together
    // to get "x+d" and "x-d" objects.
    opt_variables_type x_vars_driver;

    // List of variables for the delta or "d" parameters
    // See comments for x_vars_driver for the difference between this
    // and the d_vars object.
    opt_variables_type d_vars_driver;

    // If true then optimize the "x" variables
    bool opt_x_vars;

    // If true then optimize the "d" variables
    bool opt_d_vars;

    // If true then enforce singlet symmetry on the FDLR wave function.
    bool singlet;
    // If true then enforce triplet symmetry on the FDLR wave function.
    bool triplet;

    // Pointer to the "x+d" wavefunction object
    TrialWaveFunction * m_wfn_xpd;

    // Pointer to the "x-d" wavefunction object
    TrialWaveFunction * m_wfn_xmd;

    // Vector to hold the derivatives of the logarithm of the "x+d" part
    // of the FDLR wave function with respect to the optimizable parameters
    // (i.e., w.r.t. z=x+d, *not* with respect to x and d separately). This
    // derivative is therefore evaluated at x+d.
    std::vector<RealType> dlogpsi_xpd;
    // Same as above, but for the derivatives of \frac{H \psi_+}{\psi_+} w.r.t
    // the same optimizable parameters.
    std::vector<RealType> dhpsioverpsi_xpd;

    // Vector to hold the derivatives of the logarithm of the "x-d" part
    // of the FDLR wave function with respect to the optimizable parameters
    // (i.e., w.r.t. z=x-d, *not* with respect to x and d separately). This
    // is the same as dlogpsi_xpd, but now evaluated at x-d.
    std::vector<RealType> dlogpsi_xmd;
    // Same as above, but for the derivatives of \frac{H \psi_-}{\psi_-} w.r.t
    // the same optimizable parameters.
    std::vector<RealType> dhpsioverpsi_xmd;

    // Vector to hold the derivatives of the gradient of the logarithm of
    // the "x+d" part of the FDLR wave function with respect to the
    // optimizable parameters. These objects are actually used to hold the dot
    // product of the above gradient with (P.G - G_FDLR), where P.G is the
    // total gradient of the wave function, and G_FDLR is the gradient of the
    // FDLR wave function specifically. This makes it smaller to store, and is
    // what is needed to calculate the dhpsioverpsi_fdlr vectors.
    std::vector<RealType> dgradlogpsi_xpd;
    // Same as above but for the "x-d" part of the FDLR wave function.
    std::vector<RealType> dgradlogpsi_xmd;

    // As above for dgradlogpsi_xpd and dgradlogpsi_xmd, but these vectors
    // hold the equivalent derivatives of the total FDLR wave function with
    // respect to the "x" and "d" parameters, respectively.
    std::vector<RealType> dgradlogpsi_fdlr_d;
    std::vector<RealType> dgradlogpsi_fdlr_x;

    // These two vectors are only used if x parameters are being optimized.
    // Vector to hold the derivatives of the logarithm of the FDLR wave
    // function w.r.t the "x" FDLR parameters.
    std::vector<RealType> dlogpsi_fdlr_x;
    // Same as above, but for the derivatives of
    // \frac{H \Psi_{FDLR}}{\{Psi_{FDLR}}.
    std::vector<RealType> dhpsioverpsi_fdlr_x;

    // These two vectors are the same the above two, but derivatives are taken
    // w.r.t the "d" FDLR parameters. Only used if d parameters are being
    // optimized.
    std::vector<RealType> dlogpsi_fdlr_d;
    std::vector<RealType> dhpsioverpsi_fdlr_d;

    ValueType curRatio;

    // gradients of FDLR wavefunction
    ParticleSet::ParticleGradient_t G_FDLR;
    // laplacians of FDLR wavefunction
    ParticleSet::ParticleLaplacian_t L_FDLR;


  public:

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  constructor
    ///
    /// \param[in]   wfn_1        Pointer to the TrialWaveFunction object which will hold the "x+d"
    ///                           wave function. IMPORTANT: On input, its parameters should equal
    ///                           the intended initial values for "d", *NOT* "x+d".
    /// \param[in]   wfn_2        Pointer to the TrialWaveFunction object which will hold the "x-d"
    ///                           wave function. IMPORTANT: On input, its parameters should equal
    ///                           the intended initial values for "x", *NOT* "x-d".
    /// \param[in]   P            The particle set.
    /// \param[in]   opt_x        If true then optimize the "x" parameters, otherwise don't.
    /// \param[in]   opt_d        If true then optimize the "d" parameters, otherwise don't.
    /// \param[in]   singlet_loc  If true, enforce singlet symmetry on WF parameters.
    /// \param[in]   triplet_loc  If true, enforce triplet symmetry on WF parameters.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////
    FDLRWfn(TrialWaveFunction * const wfn_1, TrialWaveFunction * const wfn_2, ParticleSet& P,
            const bool opt_x = true, const bool opt_d = true,
            const bool singlet_loc = false, const bool triplet_loc = false);

    ~FDLRWfn();

    void init_driver_vars_singlet_or_triplet();

    OrbitalBasePtr makeClone(ParticleSet& P) const;

    void checkInVariables(opt_variables_type& active);

    void checkOutVariables(const opt_variables_type& active);

    void extract_x_and_d_vars(const opt_variables_type& active);

    void extract_xpd_and_xmd_vars(const opt_variables_type& active, opt_variables_type& xpd_vars, opt_variables_type& xmd_vars);

    void resetParameters(const opt_variables_type& active);

    void reportStatus(std::ostream& os);

    void resetTargetParticleSet(ParticleSet& P);

    ValueType evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

    RealType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

    RealType evaluateLogFDLR(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L,
                             const RealType& logpsi_plus,                    const RealType logpsi_minus,
                             const RealType& phasevalue_plus,                const RealType phasevalue_minus,
                             const ParticleSet::ParticleGradient_t& G_plus,  const ParticleSet::ParticleGradient_t& G_minus,
                             const ParticleSet::ParticleLaplacian_t& L_plus, const ParticleSet::ParticleLaplacian_t& L_minus);

    GradType evalGrad(ParticleSet& P, int iat);

    ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

    void registerData(ParticleSet& P, WFBufferType& buf);

    RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

    void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

    ValueType ratio(ParticleSet& P, int iat);

    void acceptMove(ParticleSet& P, int iat);

    void restore(int iat);

    void evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars,
        std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi);

    void copy_recompute_vector(const opt_variables_type& vars_full, opt_variables_type& vars_part);

    void evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
        std::vector<ValueType>& ratios, Matrix<ValueType>& dratios);

    void resetPhaseDiff();
};

} // end namespace qmcplusplus

#endif

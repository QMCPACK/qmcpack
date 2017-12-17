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

#include "QMCWaveFunctions/FDLRWfn.h"

namespace qmcplusplus {

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
  FDLRWfn::FDLRWfn(TrialWaveFunction * const wfn_1, TrialWaveFunction * const wfn_2, ParticleSet& P,
          const bool opt_x, const bool opt_d, const bool singlet_loc, const bool triplet_loc) :
    m_wfn_xpd(wfn_1), m_wfn_xmd(wfn_2), opt_x_vars(opt_x), opt_d_vars(opt_d)
  {
    // Temporary hack: I get very odd results when turning off Optimizable
    // and when using the new optimizer (although OK with the old one),
    // which is prsumably related to CopyFromDerivativeBuffer being called
    // when Optimizable = false. For now, just set it to true, and don't
    // update anything in evaluateDerivatives.
    //if (opt_x || opt_d)
    //  Optimizable = true;
    //else
    //  Optimizable = false;
    Optimizable = true;

    OrbitalName = "FDLR";

    // Copy variables "x", which are currently encoded in the m_wfn_xmd wave
    // function, to the x_vars object, and similarly for "d" variables,
    // which are encoded in the m_wfn_xpd object.
    m_wfn_xmd->checkInVariables(x_vars);
    x_vars.resetIndex();
    m_wfn_xmd->checkOutVariables(x_vars);

    m_wfn_xpd->checkInVariables(d_vars);
    d_vars.resetIndex();
    m_wfn_xpd->checkOutVariables(d_vars);

    if (singlet_loc)
      singlet = true;
    else
      singlet = false;

    if (triplet_loc)
      triplet = true;
    else
      triplet = false;

    // Insert the variables "x" and "d" into the "driver" objects, with x
    // variables first, then d ones. Here, we have to give the x and d
    // variables unique names, so attach "x_" or "d_" to the start of the
    // variable names, appropriately.
    if (singlet_loc || triplet_loc)
    {
      init_driver_vars_singlet_or_triplet();
    }
    else
    {
      for (int i=0; i < x_vars.size(); i++) {
        std::stringstream sstr;
        sstr << "x_" << x_vars.NameAndValue[i].first;
        x_vars_driver.insert(sstr.str(), x_vars[i], true);
      }
      for (int i=0; i < d_vars.size(); i++) {
        std::stringstream sstr;
        sstr << "d_" << d_vars.NameAndValue[i].first;
        d_vars_driver.insert(sstr.str(), d_vars[i], true);
      }
    }

    opt_variables_type xpd_vars;
    xpd_vars.insertFromSum(x_vars,d_vars);
    // Make it so that m_wfn_xpd holds the parameters
    // x+d, where as before it just held d.
    m_wfn_xpd->resetParameters(xpd_vars);

    opt_variables_type xmd_vars;
    xmd_vars.insertFromDiff(x_vars,d_vars);
    // Make it so that m_wfn_xmd holds the parameters
    // x-d, where as before it just held x.
    m_wfn_xmd->resetParameters(xmd_vars);

    // Initialize G and L objects for both wave functions.
    m_wfn_xpd->G.resize(P.G.size());
    m_wfn_xmd->G.resize(P.G.size());
    m_wfn_xpd->L.resize(P.L.size());
    m_wfn_xmd->L.resize(P.L.size());
    G_FDLR.resize(P.G.size());
    L_FDLR.resize(P.L.size());

    dlogpsi_xpd.resize(xpd_vars.size(), 0.0);
    dlogpsi_xmd.resize(xmd_vars.size(), 0.0);
    dhpsioverpsi_xpd.resize(xpd_vars.size(), 0.0);
    dhpsioverpsi_xmd.resize(xpd_vars.size(), 0.0);
    dgradlogpsi_xpd.resize(xpd_vars.size(), 0.0);
    dgradlogpsi_xmd.resize(xmd_vars.size(), 0.0);

    if (opt_x_vars) {
      dlogpsi_fdlr_x.resize(x_vars_driver.size(), 0.0);
      dhpsioverpsi_fdlr_x.resize(x_vars_driver.size(), 0.0);
      dgradlogpsi_fdlr_x.resize(x_vars_driver.size(), 0.0);
    }
    if (opt_d_vars) {
      dlogpsi_fdlr_d.resize(d_vars_driver.size(), 0.0);
      dhpsioverpsi_fdlr_d.resize(d_vars_driver.size(), 0.0);
      dgradlogpsi_fdlr_d.resize(d_vars_driver.size(), 0.0);
    }
  };

  void FDLRWfn::init_driver_vars_singlet_or_triplet() {
    // The number of "up" parameters should be the same as the number of
    // "down ones, so halve the total number.
    int nvars_up = x_vars.size()/2;

    for (int i=0; i < nvars_up; i++) {
      // Take the optimizable orbital parameter name, and remove the first
      // part of the name, so that the only bit left is "orb_rot_0000_0001"
      // part. We don't want the first bit that refers to the name of the
      // "up" or "down" determinants.
      std::string str_full;
      str_full = x_vars.NameAndValue[i].first;
      str_full.erase(str_full.begin(), str_full.end()-17);

      // Add "x_" to the start, to separate from the "d" parameters below.
      std::stringstream sstr;
      sstr << "x_" << str_full;
      x_vars_driver.insert(sstr.str(), x_vars[i], true);
    }
    // Same as above, but now for the "d" variables.
    for (int i=0; i < nvars_up; i++) {
      std::string str_full;
      str_full = d_vars.NameAndValue[i].first;
      str_full.erase(str_full.begin(), str_full.end()-17);

      std::stringstream sstr;
      sstr << "d_" << str_full;
      d_vars_driver.insert(sstr.str(), d_vars[i], true);
    }
  }

  FDLRWfn::~FDLRWfn()
  {
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  return a clone of this object
  ///
  /// \param[in]      P  Reference to the particle set to be used by the clone.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  OrbitalBasePtr FDLRWfn::makeClone(ParticleSet& P) const {
    TrialWaveFunction* wfn_xpd_clone = m_wfn_xpd->makeClone(P);
    TrialWaveFunction* wfn_xmd_clone = m_wfn_xmd->makeClone(P);
    FDLRWfn* fdlr_clone = new FDLRWfn( wfn_xpd_clone, wfn_xmd_clone, P );
    return OrbitalBasePtr( fdlr_clone );
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Check this object's optimizable variables into the supplied overall list of
  ///         optimizable variables.
  ///
  /// \param[in,out]  active        the overall list of optimizable variables
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::checkInVariables(opt_variables_type& active) {
    if (opt_x_vars)
      active.insertFrom(x_vars_driver);
    if (opt_d_vars)
      active.insertFrom(d_vars_driver);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Records the positions of this object's optimizable variables in the supplied
  ///         overall list of optimizable variables.
  ///
  /// \param[in]      active        the overall list of optimizable variables
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::checkOutVariables(const opt_variables_type& active) {
    if (opt_x_vars) {
      x_vars_driver.num_active_vars = 0;
      x_vars_driver.getIndex(active);
    }
    if (opt_d_vars) {
      d_vars_driver.num_active_vars = 0;
      d_vars_driver.getIndex(active);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Reads the optimizable FDLR variables into the x_vars_driver and d_vars_driver
  ///         objects, and then subsequently copies these into the x_vars and d_vars objects,
  ///         ready for x_vars and d_vars to be added together and subtracted as needed.
  ///
  /// \param[in]      active        the overall list of optimizable variables
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::extract_x_and_d_vars(const opt_variables_type& active)
  {
    // Copy active variables to the x_vars_driver object.
    if (opt_x_vars) {
      for (int i=0; i < x_vars_driver.size(); i++)
      {
        int loc = x_vars_driver.where(i);
        if (loc >= 0)
          x_vars_driver[i] = active[loc];
      }
    }
    // Copy active variables to the d_vars_driver object.
    if (opt_d_vars) {
      for (int i=0; i < d_vars_driver.size(); i++)
      {
        int loc = d_vars_driver.where(i);
        if (loc >= 0)
          d_vars_driver[i] = active[loc];
      }
    }

    // Copy the "x" and "d" parameters from the 'driver' objects to the
    // x_vars and d_vars objects.
    for (int i=0; i < x_vars_driver.size(); i++)
      x_vars[i] = x_vars_driver[i];
    for (int i=0; i < d_vars_driver.size(); i++)
      d_vars[i] = d_vars_driver[i];

    // If singlet/triplet symmetry, then the variables are only stored once
    // (for both up- and down-electron determinants, which both have the
    // same paramaters with singlet symmetry) in the driver objects, whereas
    // for the x_vars and d_vars objects, they are stored twice, repeated,
    // because of the need for both the up and down determinants to have
    // access to the variables, without modifying those classes to deal with
    // the singlet special case. So we need to copy them across again.
    if (singlet)
    {
      for (int i=0; i < x_vars_driver.size(); i++)
        x_vars[i + x_vars_driver.size()] = x_vars_driver[i];
      for (int i=0; i < d_vars_driver.size(); i++)
        d_vars[i + d_vars_driver.size()] = d_vars_driver[i];
    }
    // For triplet symmetry, the two sets of parameters need to have opposite
    // signs, but only for the d_* variables.
    if (triplet)
    {
      for (int i=0; i < x_vars_driver.size(); i++)
        x_vars[i + x_vars_driver.size()] = x_vars_driver[i];
      for (int i=0; i < d_vars_driver.size(); i++)
        d_vars[i + d_vars_driver.size()] = -d_vars_driver[i];
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Firstly calls extract_x_and_d_vars to copy the x and d variables to x_vars and
  ///         d_vars, and then adds and subtracts these objects to create the xpd_vars and
  ///         xmd_vars objects, as required for constructing the FDLR wave function.
  ///
  /// \param[in]      active        the overall list of optimizable variables
  /// \param[out]     xpd_vars      the sum of "x" and "d" FDLR parameters
  /// \param[out]     xmd_vars      the difference of "x" and "d" FDLR parameters
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::extract_xpd_and_xmd_vars(const opt_variables_type& active, opt_variables_type& xpd_vars, opt_variables_type& xmd_vars)
  {
    // Set x_vars and d_vars using the input active.
    extract_x_and_d_vars(active);

    // Construct the "x+d" and "x-d" variable objects.
    xpd_vars.insertFromSum(x_vars, d_vars);
    xmd_vars.insertFromDiff(x_vars, d_vars);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Update the "x+d" and "x-d" TrialWaveFunction objects which define the FDLR
  ///         wave function, so that their parameters are consistent with those to be read
  ///         in from the updated total parameter list, as passed in through active.
  ///
  /// \param[in]      active         variable set to read new "x" and "d" parameters from.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::resetParameters(const opt_variables_type& active)
  {
    // Set xpd_vars and xmd_vars using the input active.
    opt_variables_type xpd_vars, xmd_vars;
    extract_xpd_and_xmd_vars(active, xpd_vars, xmd_vars);

    // Apply the summed "x+d" parameters to the m_wfn_xpd wave function.
    m_wfn_xpd->resetParameters(xpd_vars);
    // Apply the subtracted "x-d" parameters to the m_wfn_xmd wave function.
    m_wfn_xmd->resetParameters(xmd_vars);
  }

  void FDLRWfn::reportStatus(std::ostream& os) {
    throw std::runtime_error("FDLRWfn::reportStatus not yet implemented");
  }

  void FDLRWfn::resetTargetParticleSet(ParticleSet& P) {
    throw std::runtime_error("FDLRWfn::resetTargetParticleSet not yet implemented");
  }

  FDLRWfn::ValueType FDLRWfn::evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {
    throw std::runtime_error("FDLRWfn::evaluate not yet implemented");
    return 0.0;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Evaluates the log of the FDLR wave function, and adds the gradient and laplacian
  ///         of the log of the FDLR wave function to the running total gradient and laplacian.
  ///
  ///         Also stores the log and phase values for this wave function in this object's
  ///         LogValue and PhaseValue variables, and also stores the log value, phase value,
  ///         gradient and laplacian of the "x+d" and "x-d" wave functions in the corresponding
  ///         TrialWaveFunction objects.
  ///
  /// \param[in]          P      the particle set
  /// \param[in,out]      G      gradient to add to
  /// \param[in,out]      L      laplacian to add to
  ///
  /// \return  the log of the FDLR wave function
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  FDLRWfn::RealType FDLRWfn::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    // "x+d"
    m_wfn_xpd->evaluateLog(P);
    m_wfn_xpd->G = P.G;
    m_wfn_xpd->L = P.L;

    // "x-d"
    m_wfn_xmd->evaluateLog(P);
    m_wfn_xmd->G = P.G;
    m_wfn_xmd->L = P.L;

    // Update the log value, gradient and laplacian for the FDLR wave
    // function, given these objects for the individual "x+d" and "x-d"
    // wave functions.
    LogValue = evaluateLogFDLR(P, G, L,
                               m_wfn_xpd->getLogPsi(), m_wfn_xmd->getLogPsi(),
                               m_wfn_xpd->getPhase(), m_wfn_xmd->getPhase(),
                               m_wfn_xpd->G, m_wfn_xmd->G,
                               m_wfn_xpd->L, m_wfn_xmd->L);

    return LogValue;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Evaluates the log of the FDLR wave function, and adds the gradient and laplacian
  ///         of the log of the FDLR wave function to the running total gradient and laplacian.
  ///         This is done given the log value, gradient and laplacian for the individual
  ///         "x+d" and "x-d" wave functions, which are passed in as parameters.
  ///
  ///         Also stores the log and phase values for the FDLR wave function in this object's
  ///         LogValue and PhaseValue variables.
  ///
  /// \param[in]          P              the particle set
  /// \param[in,out]      G              gradient to add to
  /// \param[in,out]      L              laplacian to add to
  /// \param[in]          logpsi_plus    the log of the "x+d" wave function
  /// \param[in]          logpsi_minus   the log of the "x-d" wave function
  /// \param[in]          G_plus         the gradient of the log of the "x+d" wave function
  /// \param[in]          G_minus        the gradient of the log of the "x-d" wave function
  /// \param[in]          L_plus         the laplacian of the log of the "x+d" wave function
  /// \param[in]          L_minus        the laplacian of the log of the "x-d" wave function
  ///
  /// \return  the log of the FDLR wave function
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  FDLRWfn::RealType FDLRWfn::evaluateLogFDLR(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L,
                           const FDLRWfn::RealType& logpsi_plus,           const FDLRWfn::RealType logpsi_minus,
                           const FDLRWfn::RealType& phasevalue_plus,       const FDLRWfn::RealType phasevalue_minus,
                           const ParticleSet::ParticleGradient_t& G_plus,  const ParticleSet::ParticleGradient_t& G_minus,
                           const ParticleSet::ParticleLaplacian_t& L_plus, const ParticleSet::ParticleLaplacian_t& L_minus)
  {
    FDLRWfn::ValueType logpsi(0.0), psi(0.0), psi_plus(0.0), psi_minus(0.0);
    FDLRWfn::ValueType scaling_fac_1, scaling_fac_2;

    // Temporary space needed for calculating the gradient and the laplacian
    // of the FDLR wfn.
    ParticleSet::ParticleLaplacian_t G_FDLR_mag;
    ParticleSet::ParticleLaplacian_t L_temp_1;
    ParticleSet::ParticleLaplacian_t L_temp_2;
    G_FDLR_mag.resize(G.size());
    L_temp_1.resize(L.size());
    L_temp_2.resize(L.size());

    // ----Calculating the log of the FDLR wave function--------------------

    psi_plus = std::exp(logpsi_plus)*std::cos(phasevalue_plus);
    psi_minus = std::exp(logpsi_minus)*std::cos(phasevalue_minus);
    psi = psi_plus - psi_minus;
    //psi = std::exp(logpsi_plus)*std::cos(phasevalue_plus) - std::exp(logpsi_minus)*std::cos(phasevalue_minus);
    logpsi = evaluateLogAndPhase(psi, PhaseValue);

    // ----Calculating the gradient of the log of the FDLR wave function----

    scaling_fac_1 = 1/(1 - psi_minus/psi_plus);
    scaling_fac_2 = 1/(psi_plus/psi_minus - 1);

    G = scaling_fac_1 * G_plus - scaling_fac_2 * G_minus;

    // ----Calculating the laplacian of the log of the FDLR wave function---

    // Calculate \frac{ \del psi \cdot \del \psi}{ psi^2 }, for
    // \psi=\psi_+ and then for \psi=\psi_-.
    for (int i=0; i < G.size(); i++)
      G_FDLR_mag[i] = dot(G[i], G[i]);
    for (int i=0; i < G.size(); i++)
      L_temp_1[i] = dot(G_plus[i],  G_plus[i]);
    for (int i=0; i < G.size(); i++)
      L_temp_2[i] = dot(G_minus[i], G_minus[i]);

    // Calculate \frac{\del^2 \psi}{\psi} for \psi=\psi_+ and for \psi=\psi_-.
    L_temp_1 = L_plus  + L_temp_1;
    L_temp_2 = L_minus + L_temp_2;

    L = scaling_fac_1*L_temp_1 - scaling_fac_2*L_temp_2 - G_FDLR_mag;

    // ---------------------------------------------------------------------

    convert(logpsi, LogValue);

    return LogValue;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Evaluates and returns the gradient of the log of the FDLR wave function w.r.t. a
  ///         specified particle's position.
  ///
  /// \param[in]    P       the particle set
  /// \param[in]    iat     index of the particle in question
  ///
  /// \return  the one particle gradient of the log of the FDLR wave function.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  FDLRWfn::GradType FDLRWfn::evalGrad(ParticleSet& P, int iat)
  {
    FDLRWfn::RealType logpsi_plus = m_wfn_xpd->getLogPsi();
    FDLRWfn::RealType logpsi_minus = m_wfn_xmd->getLogPsi();
    FDLRWfn::RealType phasevalue_plus = m_wfn_xpd->getPhase();
    FDLRWfn::RealType phasevalue_minus = m_wfn_xmd->getPhase();

    FDLRWfn::ValueType psi_plus(0.0), psi_minus(0.0);
    psi_plus = std::exp(logpsi_plus)*std::cos(phasevalue_plus);
    psi_minus = std::exp(logpsi_minus)*std::cos(phasevalue_minus);

    FDLRWfn::ValueType scaling_fac_1 = 1/(1 - psi_minus/psi_plus);
    FDLRWfn::ValueType scaling_fac_2 = 1/(psi_plus/psi_minus - 1);

    FDLRWfn::GradType G_plus = m_wfn_xpd->evalGrad(P, iat);
    FDLRWfn::GradType G_minus = m_wfn_xmd->evalGrad(P, iat);
    FDLRWfn::GradType grad = scaling_fac_1*G_plus - scaling_fac_2*G_minus;

    return grad;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  For a given one particle move, this function evaluates the ratio of the new and
  ///         old FDLR wave function values.
  ///
  ///         Also, it calculates and adds in this FDLR wave function's contribution to grad_iat,
  ///         which is the gradient (evaluated at the new position) of the log of the overall
  ///         trial function w.r.t. the moved particle's coordinates.
  ///
  /// \param[in]      P              the particle set
  /// \param[in]      iat            the id number of the moved particle
  /// \param[in,out]  grad_iat       cumulative total of the gradient of the overall trial
  ///                                function's log w.r.t. the moved particle's coordinates
  ///
  /// \return  the ratio of new and old FDLR wave function values.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  FDLRWfn::ValueType FDLRWfn::ratioGrad(ParticleSet& P, int iat, FDLRWfn::GradType& grad_iat)
  {
    FDLRWfn::RealType logpsi_plus = m_wfn_xpd->getLogPsi();
    FDLRWfn::RealType logpsi_minus = m_wfn_xmd->getLogPsi();
    FDLRWfn::RealType phasevalue_plus = m_wfn_xpd->getPhase();
    FDLRWfn::RealType phasevalue_minus = m_wfn_xmd->getPhase();

    FDLRWfn::GradType G_plus, G_minus;

    // On output G_plus holds the gradient of the log of \psi_+.
    FDLRWfn::ValueType rat_plus = m_wfn_xpd->ratioGrad(P, iat, G_plus);
    // On output G_minus holds the gradient of the log of \psi_-.
    FDLRWfn::ValueType rat_minus = m_wfn_xmd->ratioGrad(P, iat, G_minus);

    FDLRWfn::ValueType psi_plus = std::exp(logpsi_plus)*std::cos(phasevalue_plus);
    FDLRWfn::ValueType psi_minus = std::exp(logpsi_minus)*std::cos(phasevalue_minus);

    FDLRWfn::ValueType psi_plus_new = rat_plus*psi_plus;
    FDLRWfn::ValueType psi_minus_new = rat_minus*psi_minus;

    FDLRWfn::ValueType scaling_fac_1 = 1/(1 - psi_minus/psi_plus);
    FDLRWfn::ValueType scaling_fac_2 = 1/(psi_plus/psi_minus - 1);

    FDLRWfn::ValueType scaling_fac_1_new = 1/(1 - psi_minus_new/psi_plus_new);
    FDLRWfn::ValueType scaling_fac_2_new = 1/(psi_plus_new/psi_minus_new - 1);

    // The gradient of the log of the FDLR wave function, evaluated at the
    // particles' new coordinates.
    FDLRWfn::GradType grad = scaling_fac_1_new * G_plus - scaling_fac_2_new * G_minus;

    // Calculate the ratio of new and old FDLR wave functions:
    // curRatio = \frac{ \psi_+(R_new) - \psi_-(R_new) }{ \psi_+(R_old) - \psi_-(R_old) }
    // This can be arranged to the following:
    curRatio = scaling_fac_1 * rat_plus - scaling_fac_2 * rat_minus;

    // Also sum the gradient of the log of the FDLR wave function at the
    // new coordinates into the running total.
    grad_iat += grad;

    return curRatio;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Initialize the object for the given particle positions and add its essential internal
  ///         data to the supplied buffer.
  ///
  ///         Also add in the FDLR wave function's gradient and laplacian to the running totals
  ///         being accumulated in P.G and P.L.
  ///
  /// \param[in]      P              the particle set
  /// \param[in]      buf            the buffer to add essential data to
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::registerData(ParticleSet& P, WFBufferType& buf)
  {
    m_wfn_xpd->registerData(P, buf);
    buf.add(m_wfn_xpd->G.first_address(),m_wfn_xpd->G.last_address());
    buf.add(m_wfn_xpd->L.first_address(),m_wfn_xpd->L.last_address());

    m_wfn_xmd->registerData(P, buf);
    buf.add(m_wfn_xmd->G.first_address(),m_wfn_xmd->G.last_address());
    buf.add(m_wfn_xmd->L.first_address(),m_wfn_xmd->L.last_address());

    buf.add(LogValue);
    buf.add(PhaseValue);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Put essential internal data in the supplied buffer.
  ///
  ///         Also add in the FDLR wave function's gradient and laplacian to the running totals
  ///         being accumulated in P.G and P.L.
  ///
  /// \param[in]      P              the particle set
  /// \param[in]      buf            the buffer to save essential data in
  /// \param[in]      fromscratch    ??? - I assume this specifies whether data is all re-evaluated?
  ///
  /// \return  the log of the FDLR wave function.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  FDLRWfn::RealType FDLRWfn::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
  {
    // Update buffer with data from the "x+d" part of the wave function.
    FDLRWfn::RealType logpsi_plus = m_wfn_xpd->updateBuffer(P, buf, fromscratch);
    m_wfn_xpd->G = P.G;
    m_wfn_xpd->L = P.L;
    buf.put(m_wfn_xpd->G.first_address(),m_wfn_xpd->G.last_address());
    buf.put(m_wfn_xpd->L.first_address(),m_wfn_xpd->L.last_address());

    // Update buffer with data from the "x-d" part of the wave function.
    FDLRWfn::RealType logpsi_minus = m_wfn_xmd->updateBuffer(P, buf, fromscratch);
    m_wfn_xmd->G = P.G;
    m_wfn_xmd->L = P.L;
    buf.put(m_wfn_xmd->G.first_address(),m_wfn_xmd->G.last_address());
    buf.put(m_wfn_xmd->L.first_address(),m_wfn_xmd->L.last_address());

    // Calculate data for the overall FDLR wave function.
    LogValue = evaluateLogFDLR(P, P.G, P.L,
                               logpsi_plus, logpsi_minus,
                               m_wfn_xpd->getPhase(), m_wfn_xmd->getPhase(),
                               m_wfn_xpd->G, m_wfn_xmd->G,
                               m_wfn_xpd->L, m_wfn_xmd->L);
    buf.put(LogValue);
    buf.put(PhaseValue);

    return LogValue;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Copy data that was stored in the internal buffer into the FDLR object and into the
  ///         OrbitalBase objects within the TrialWaveFunction objects.
  ///
  /// \param[in]    P         the particle set
  /// \param[in]    buf       the buffer to read from
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    m_wfn_xpd->copyFromBuffer(P, buf);
    buf.get(m_wfn_xpd->G.first_address(),m_wfn_xpd->G.last_address());
    buf.get(m_wfn_xpd->L.first_address(),m_wfn_xpd->L.last_address());

    m_wfn_xmd->copyFromBuffer(P, buf);
    buf.get(m_wfn_xmd->G.first_address(),m_wfn_xmd->G.last_address());
    buf.get(m_wfn_xmd->L.first_address(),m_wfn_xmd->L.last_address());

    LogValue = evaluateLogFDLR(P, P.G, P.L,
                               m_wfn_xpd->getLogPsi(), m_wfn_xmd->getLogPsi(),
                               m_wfn_xpd->getPhase(), m_wfn_xmd->getPhase(),
                               m_wfn_xpd->G, m_wfn_xmd->G,
                               m_wfn_xpd->L, m_wfn_xmd->L);
    buf.get(LogValue);
    buf.get(PhaseValue);
  }

  FDLRWfn::ValueType FDLRWfn::ratio(ParticleSet& P, int iat)
  {
    FDLRWfn::RealType logpsi_plus = m_wfn_xpd->getLogPsi();
    FDLRWfn::RealType logpsi_minus = m_wfn_xmd->getLogPsi();
    FDLRWfn::RealType phasevalue_plus = m_wfn_xpd->getPhase();
    FDLRWfn::RealType phasevalue_minus = m_wfn_xmd->getPhase();

    FDLRWfn::ValueType rat_plus = m_wfn_xpd->ratio(P, iat);
    FDLRWfn::ValueType rat_minus = m_wfn_xmd->ratio(P, iat);

    FDLRWfn::ValueType psi_plus = std::exp(logpsi_plus)*std::cos(phasevalue_plus);
    FDLRWfn::ValueType psi_minus = std::exp(logpsi_minus)*std::cos(phasevalue_minus);

    FDLRWfn::ValueType psi_plus_new = rat_plus*psi_plus;
    FDLRWfn::ValueType psi_minus_new = rat_minus*psi_minus;

    FDLRWfn::ValueType scaling_fac_1 = 1/(1 - psi_minus/psi_plus);
    FDLRWfn::ValueType scaling_fac_2 = 1/(psi_plus/psi_minus - 1);

    FDLRWfn::ValueType scaling_fac_1_new = 1/(1 - psi_minus_new/psi_plus_new);
    FDLRWfn::ValueType scaling_fac_2_new = 1/(psi_plus_new/psi_minus_new - 1);

    // Calculate the ratio of new and old FDLR wave functions:
    // curRatio = \frac{ \psi_+(R_new) - \psi_-(R_new) }{ \psi_+(R_old) - \psi_-(R_old) }
    // This can be arranged to the following:
    curRatio = scaling_fac_1 * rat_plus - scaling_fac_2 * rat_minus;

    return curRatio;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Updates data for the FDLR object and its internal OrbitalBase objects, to account
  ///         for an accepted move of a single particle.
  ///
  ///         Recalculate and store the log of the FDLR wave function.
  ///
  /// \param[in]     P        the particle set
  /// \param[in]     iat      the id number of the moved particle
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::acceptMove(ParticleSet& P, int iat)
  {
    m_wfn_xpd->acceptMove(P, iat);
    m_wfn_xmd->acceptMove(P, iat);

    PhaseValue += evaluatePhase(curRatio);
    LogValue += std::log(std::abs(curRatio));

    curRatio = 1.0;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Perform any necessary operations related to the fact that a proposed single
  ///         particle move has been rejected.
  ///
  /// \param[in]    iat     the id number of the moved particle
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::restore(int iat)
  {
    m_wfn_xpd->rejectMove(iat);
    m_wfn_xmd->rejectMove(iat);

    curRatio = 1.0;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief  Evaluate the FDLR wave function's contribution to the derivatives of log of the
  ///         overall trial function and the local energy w.r.t. the optimizable FDLR parameters.
  ///
  /// \param[in]      P              the particle set
  /// \param[in]      optvars        the current values of the optimizable parameters
  /// \param[in,out]  dlogpsi        the derivative of the log of the trial wave function w.r.t
  ///                                optimizable parameters
  /// \param[in,out]  dhpsioverpsi   the derivative of the trial wave function's local energy
  ///                                w.r.t optimizable parameters
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars,
      std::vector<FDLRWfn::RealType>& dlogpsi, std::vector<FDLRWfn::RealType>& dhpsioverpsi)
  {

    // The number of "x" parameters (which is also the number of "d"
    // parameters, hence the name - it is *not* the total number of both
    // sets of parameters together).
    int nvars_x_or_d;

    // Set xpd_vars and xmd_vars using the input optvars variables.
    opt_variables_type xpd_vars, xmd_vars;
    extract_xpd_and_xmd_vars(optvars, xpd_vars, xmd_vars);

    if (opt_d_vars || opt_x_vars) {
      copy_recompute_vector(optvars, xpd_vars);
      copy_recompute_vector(optvars, xmd_vars);
    }

    // Zero all vectors where the derivatives will be accumulated.
    std::fill(dlogpsi_xpd.begin(), dlogpsi_xpd.begin()+xpd_vars.size(), 0.0);
    std::fill(dlogpsi_xmd.begin(), dlogpsi_xmd.begin()+xmd_vars.size(), 0.0);
    std::fill(dhpsioverpsi_xpd.begin(), dhpsioverpsi_xpd.begin()+xpd_vars.size(), 0.0);
    std::fill(dhpsioverpsi_xmd.begin(), dhpsioverpsi_xmd.begin()+xmd_vars.size(), 0.0);

    std::fill(dgradlogpsi_xpd.begin(), dgradlogpsi_xpd.begin()+xpd_vars.size(), 0.0);
    std::fill(dgradlogpsi_xmd.begin(), dgradlogpsi_xmd.begin()+xmd_vars.size(), 0.0);

    // Difference between the total FDLR wave function gradients.
    ParticleSet::ParticleGradient_t G_diff;
    G_diff.resize(P.G.size());

    // Store the total FDLR wave function's G and L values in a temporary
    // particle set, because we will need to set P's G and L to equal those
    // of the "xpd" and "xmd" wave functions individually for the following
    // evaluateDerivatives call. evaluateDerivatives uses P.G for certain
    // wave function components.
    G_FDLR = P.G;
    L_FDLR = P.L;

    P.G = m_wfn_xpd->G;
    P.L = m_wfn_xpd->L;
    m_wfn_xpd->evaluateDerivatives(P, xpd_vars, dlogpsi_xpd, dhpsioverpsi_xpd);

    P.G = m_wfn_xmd->G;
    P.L = m_wfn_xmd->L;
    m_wfn_xmd->evaluateDerivatives(P, xmd_vars, dlogpsi_xmd, dhpsioverpsi_xmd);

    // Return G and L to their original values for the entire FDLR wave
    // function.
    P.G = G_FDLR;
    P.L = L_FDLR;

    // Calculate the log of the \psi_+ and \psi_- wave functions, and use
    // these values to calculate the required scaling factors for the
    // various FDLR wave function derivatives we're about to calculate.
    FDLRWfn::ValueType logpsi_plus = m_wfn_xpd->getLogPsi();
    FDLRWfn::ValueType logpsi_minus = m_wfn_xmd->getLogPsi();
    FDLRWfn::RealType phasevalue_plus = m_wfn_xpd->getPhase();
    FDLRWfn::RealType phasevalue_minus = m_wfn_xmd->getPhase();

    // Calculate the wave function values for "x+d" and "x-d" parts.
    FDLRWfn::ValueType psi_plus = std::exp(logpsi_plus)*std::cos(phasevalue_plus);
    FDLRWfn::ValueType psi_minus = std::exp(logpsi_minus)*std::cos(phasevalue_minus);

    // Equal to \frac{\psi_+}{\psi_-}.
    FDLRWfn::ValueType plus = psi_plus/psi_minus;
    //FDLRWfn::ValueType plus = std::exp(logpsi_plus - logpsi_minus);

    // Equal to \frac{\psi_-}{\psi_+}.
    FDLRWfn::ValueType minus = psi_minus/psi_plus;
    //FDLRWfn::ValueType minus = std::exp(logpsi_minus - logpsi_plus);

    FDLRWfn::ValueType scaling_fac_1 = 1/(1 - minus);
    FDLRWfn::ValueType scaling_fac_2 = 1/(plus - 1);

    G_FDLR = scaling_fac_1*m_wfn_xpd->G - scaling_fac_2*m_wfn_xmd->G;
    G_diff = P.G - G_FDLR;

    m_wfn_xpd->evaluateGradDerivatives(G_diff, dgradlogpsi_xpd);
    m_wfn_xmd->evaluateGradDerivatives(G_diff, dgradlogpsi_xmd);

    // To enforce singlet symmetry, we need to add together the derivatives
    // coming from both the up and down determinants, because now the same
    // optimizable parameters appear in both determinants, so the derivative
    // w.r.t those parameters is a sum, by the product rule.
    if (singlet)
    {
      nvars_x_or_d = xpd_vars.size()/2;

      for (int i=0; i<nvars_x_or_d; i++)
      {
        dlogpsi_xpd[i] += dlogpsi_xpd[i+nvars_x_or_d];
        dlogpsi_xmd[i] += dlogpsi_xmd[i+nvars_x_or_d];
        dhpsioverpsi_xpd[i] += dhpsioverpsi_xpd[i+nvars_x_or_d];
        dhpsioverpsi_xmd[i] += dhpsioverpsi_xmd[i+nvars_x_or_d];
        dgradlogpsi_xpd[i] += dgradlogpsi_xpd[i+nvars_x_or_d];
        dgradlogpsi_xmd[i] += dgradlogpsi_xmd[i+nvars_x_or_d];
      }
    } else {
      nvars_x_or_d = xpd_vars.size();
    }

    FDLRWfn::ValueType grad_dot_plus = 0.0;
    FDLRWfn::ValueType grad_dot_minus = 0.0;
    FDLRWfn::ValueType grad_dot_FDLR = 0.0;

    for (int i=0; i < m_wfn_xpd->G.size(); i++)
      grad_dot_plus += dot(G_diff[i], m_wfn_xpd->G[i]);
    for (int i=0; i < m_wfn_xmd->G.size(); i++)
      grad_dot_minus += dot(G_diff[i], m_wfn_xmd->G[i]);

    grad_dot_FDLR = scaling_fac_1*grad_dot_plus - scaling_fac_2*grad_dot_minus;

    for (int i=0; i<nvars_x_or_d; i++)
      dgradlogpsi_xpd[i] += grad_dot_plus * dlogpsi_xpd[i];
    for (int i=0; i<nvars_x_or_d; i++)
      dgradlogpsi_xmd[i] += grad_dot_minus * dlogpsi_xmd[i];

    // Calculate the kinetic energy of the "x+d", "x-d" and total FDLR
    // wave functions.

    ParticleSet::ParticleLaplacian_t G_plus_mag;
    ParticleSet::ParticleLaplacian_t G_minus_mag;
    G_plus_mag.resize(m_wfn_xpd->G.size());
    G_minus_mag.resize(m_wfn_xmd->G.size());
    G_plus_mag = 0.0;
    G_minus_mag = 0.0;

    // Calculate the magnitude of the gradient ratio vectors.
    for (int i=0; i < m_wfn_xpd->G.size(); i++)
      G_plus_mag[i] = dot(m_wfn_xpd->G[i], m_wfn_xpd->G[i]);
    for (int i=0; i < m_wfn_xmd->G.size(); i++)
      G_minus_mag[i] = dot(m_wfn_xmd->G[i], m_wfn_xmd->G[i]);

    ParticleSet::ParticleLaplacian_t L_temp_1;
    ParticleSet::ParticleLaplacian_t L_temp_2;
    L_temp_1.resize(m_wfn_xpd->L.size());
    L_temp_2.resize(m_wfn_xmd->L.size());

    // m_wfn_xpd->L stores the laplacian divided by the wave function, minus
    // the magnitude squared of (the gradient for the "x+d" wave function
    // divded by the wave function value), i.e.:
    //
    // m_wfn_xpd->L[i] = \frac{\nabla_i^2 \psi_+}{\psi_+} - \frac{nabla_i \psi_+}{\psi_+} \cdot \frac{nabla_i \psi_+}{\psi_+}
    //
    // so add the gradient squared to get the laplacian:
    L_temp_1 = m_wfn_xpd->L + G_plus_mag;
    // Similarly as for above, but now for the "x-d" wave function.
    L_temp_2 = m_wfn_xmd->L + G_minus_mag;

    // Now calculate the kinetic energies themselves.
    FDLRWfn::ValueType kinetic_plus = 0.0;
    FDLRWfn::ValueType kinetic_minus = 0.0;

    // L_temp_1 is currently a vector of laplacians (divded by the "plus" wave
    // function value) for each of the particles co-ordinates. The kinetic
    // energy involves a sum over all laplacians (i.e. over all particles), so
    // perform that sum, for both wave functions.
    for (int i=0; i < m_wfn_xpd->L.size(); i++)
      kinetic_plus += L_temp_1[i];
    for (int i=0; i < m_wfn_xmd->L.size(); i++)
      kinetic_minus += L_temp_2[i];

    // Get the final kinetic energy values.
    kinetic_plus *= -0.5;
    kinetic_minus *= -0.5;

    // The local kinetic energy for the FDLR wave function.
    FDLRWfn::ValueType kinetic_FDLR = scaling_fac_1*kinetic_plus - scaling_fac_2*kinetic_minus;

    // Now finally calculate the derivatives of the log of the FDLR wave
    // function and of the FDLR local energy, or "x" and/or "d" parameters,
    // as requested by the user.
    if (opt_x_vars)
    {
      for (int i=0; i < x_vars_driver.size(); i++)
      {
        dgradlogpsi_fdlr_x[i] = scaling_fac_1*dgradlogpsi_xpd[i] - scaling_fac_2*dgradlogpsi_xmd[i];

        dlogpsi_fdlr_x[i] = scaling_fac_1 * dlogpsi_xpd[i] - scaling_fac_2 * dlogpsi_xmd[i];

        dhpsioverpsi_fdlr_x[i] = scaling_fac_1 * (dhpsioverpsi_xpd[i] + kinetic_plus*dlogpsi_xpd[i])
                               - scaling_fac_2 * (dhpsioverpsi_xmd[i] + kinetic_minus*dlogpsi_xmd[i])
                               - kinetic_FDLR  * dlogpsi_fdlr_x[i]
                               - dgradlogpsi_fdlr_x[i] + grad_dot_FDLR * dlogpsi_fdlr_x[i];
      }
    }
    if (opt_d_vars)
    {
      for (int i=0; i < d_vars_driver.size(); i++)
      {
        dgradlogpsi_fdlr_d[i] = scaling_fac_1*dgradlogpsi_xpd[i] + scaling_fac_2*dgradlogpsi_xmd[i];

        dlogpsi_fdlr_d[i] = scaling_fac_1 * dlogpsi_xpd[i] + scaling_fac_2 * dlogpsi_xmd[i];

        dhpsioverpsi_fdlr_d[i] = scaling_fac_1 * (dhpsioverpsi_xpd[i] + kinetic_plus*dlogpsi_xpd[i])
                               + scaling_fac_2 * (dhpsioverpsi_xmd[i] + kinetic_minus*dlogpsi_xmd[i])
                               - kinetic_FDLR  * dlogpsi_fdlr_d[i]
                               - dgradlogpsi_fdlr_d[i] + grad_dot_FDLR * dlogpsi_fdlr_d[i];
      }
    }

    // Copy the FDLR derivatives to the FDLR section of the full optvars
    // vector, which contains the derivatives with respect to *all*
    // optimizable parameters in the wave function.
    if (opt_x_vars) {
      int x_loc = x_vars_driver.where(0);
      std::copy(dlogpsi_fdlr_x.begin(), dlogpsi_fdlr_x.begin()+x_vars_driver.size(), dlogpsi.begin()+x_loc);
      std::copy(dhpsioverpsi_fdlr_x.begin(), dhpsioverpsi_fdlr_x.begin()+x_vars_driver.size(), dhpsioverpsi.begin()+x_loc);
    }
    if (opt_d_vars) {
      int d_loc = d_vars_driver.where(0);
      std::copy(dlogpsi_fdlr_d.begin(), dlogpsi_fdlr_d.begin()+d_vars_driver.size(), dlogpsi.begin()+d_loc);
      std::copy(dhpsioverpsi_fdlr_d.begin(), dhpsioverpsi_fdlr_d.begin()+d_vars_driver.size(), dhpsioverpsi.begin()+d_loc);
    }

  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Copy the values of Recompute from vars_full to vars_part. Also perform some sanity
  ///        checking for safety. vars_full should be the full list of optimizable parameters,
  ///        while vars_part should contain parameters of the "x", "d", "x+d" or "x-d" type.
  ///
  /// \param[in]      vars_full  The full list of optimizable parameters
  /// \param[in,out]  vars_part  Partial list of parameters
  ///
  //////////////////////////////////////////////////////////////////////////////////////////////
  void FDLRWfn::copy_recompute_vector(const opt_variables_type& vars_full, opt_variables_type& vars_part)
  {
    int loc, nvars_x_or_d;

    // If enforcing singlet symmetry, then the "full" list of optimizable
    // parameters include each parameter only once. In the xpd_vars and
    // xmd_vars objects, passed in as vars_part, they appear twice, once
    // for the spin up determinant, and once for the spin down determinant.
    if (singlet)
      nvars_x_or_d = vars_part.size()/2;
    else
      nvars_x_or_d = vars_part.size();

    // Loop over the smaller input vector.
    for (int i=0; i < nvars_x_or_d; i++)
    {
      // Find the positions of the equivalent parameters in vars_full.
      if (opt_x_vars)
        loc = x_vars_driver.where(i);
      else if (opt_d_vars)
        loc = d_vars_driver.where(i);
      else
        throw std::runtime_error("No optimizable parameters for FDLR wave function, shouldn't be here.");

      if (loc >= 0)
      {
        vars_part.Recompute[i].second = vars_full.Recompute[loc].second;
        // For singlet, need to update the "up" and "down" determinant
        // parameters separately.
        if (singlet)
          vars_part.Recompute[i+nvars_x_or_d].second = vars_full.Recompute[loc].second;
      }
      else
      {
        throw std::runtime_error("Optimizable parameter not found in the list of all optimizable parameters.");
      }

      // Error checking.
      if (opt_x_vars && opt_d_vars)
      {
        // Make sure that, if both x and d variables exist in the full
        // optimizable parameters list, then their values of Recompute are
        // consistent. It would be very odd if not...
        if (vars_full.Recompute[i].second != vars_full.Recompute[i + nvars_x_or_d].second)
          throw std::runtime_error("Values of Recompute are not consistent in x and d variables.");
      }
    }
  }

  void FDLRWfn::evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
      std::vector<FDLRWfn::ValueType>& ratios, Matrix<FDLRWfn::ValueType>& dratios) {
    throw std::runtime_error("FDLRWfn::evaluateDerivRatios not yet implemented");
  }

  void FDLRWfn::resetPhaseDiff()
  {
    m_wfn_xpd->resetPhaseDiff();
    m_wfn_xmd->resetPhaseDiff();
  }

} // end namespace qmcplusplus

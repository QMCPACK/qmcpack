///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file src/QMCWaveFunctions/Fermion/SlaterDetOpt.cpp
///
/// \brief   Implementation file for a Slater determinant with optimizable orbitals.
///
/// \author  Eric Neuscamman
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/Fermion/SlaterDetOpt.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <QMCWaveFunctions/LCOrbitalSetOpt.h>
#include <Numerics/DeterminantOperators.h>
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus {

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Creates a slater determinant of the given spin with the given optimizable orbital set.
///
/// \param[in]      targetPtcl     the particle set
/// \param[in]      spo_ptr        pointer to the optimizable single particle orbital set
/// \param[in]      up_or_down     0 for up spin, 1 for down
/// \param[in]      nmo            number of optimizable molecular orbitals
///
///////////////////////////////////////////////////////////////////////////////////////////////////
SlaterDetOpt::SlaterDetOpt(ParticleSet & targetPtcl, SPOSetBase * spo_ptr, const int up_or_down)
  : m_spo(spo_ptr)
  , m_up_or_down(up_or_down)
  , m_nmo(spo_ptr->size())
{
  Optimizable=true;
  OrbitalName="SlaterDetOpt";
  this->resetTargetParticleSet(targetPtcl);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds the optimizable orbital set trial function component to the given trial function
///         if it wasn't already there.
///
/// \param[in]      twf            the trial wave function to add the component to
/// \param[in]      name           a name for the component
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::add_orbs_to_tf(TrialWaveFunction & twf, const std::string & name) {
  if ( std::find(twf.getOrbitals().begin(), twf.getOrbitals().end(), m_spo->tf_component()) == twf.getOrbitals().end() )
    twf.addOrbital(m_spo->tf_component(), name, false);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Makes a clone of the object that uses the supplied particle set.
///
/// \param[in]      tqp            the particle set the clone should use
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBasePtr SlaterDetOpt::makeClone(ParticleSet& tqp) const {
  return ( new SlaterDetOpt(tqp, m_spo->makeClone(), m_up_or_down) );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Destructor has nothing to do for now
///
///////////////////////////////////////////////////////////////////////////////////////////////////
SlaterDetOpt::~SlaterDetOpt() { }

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  reset which particle set we are using and initialize arrays accordingly
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::resetTargetParticleSet(ParticleSet& P) {

  // set which and how many particles we care about
  m_first = P.first(m_up_or_down);
  m_last = P.last(m_up_or_down);
  m_nel = m_last - m_first;

  // reset our optimizable orbitals object
  m_spo->resetTargetParticleSet(P);

  // resize matrices and arrays
  m_orb_val_mat_all.resize(m_nel, m_nmo);
  m_orb_der_mat_all.resize(m_nel, m_nmo);
  m_orb_lap_mat_all.resize(m_nel, m_nmo);
  m_orb_inv_mat.resize(m_nel, m_nel);
  m_orb_val_mat.resize(m_nel, m_nel);
  m_orb_der_mat.resize(m_nel, m_nel);
  m_orb_lap_mat.resize(m_nel, m_nel);
  m_orb_val_vec.resize(m_nel);
  m_orb_der_vec.resize(m_nel);
  m_orb_lap_vec.resize(m_nel);
  m_dp0.resize(  m_nel, m_nmo);
  m_dh0.resize(  m_nel, m_nmo);
  m_dh1.resize(3*m_nel, m_nmo);
  m_dh2.resize(  m_nel, m_nmo);
  m_work.resize( std::max( size_t(m_dh1.size()), size_t(std::max(m_nel, 10) * m_nel + m_nel * m_nel)) );
  m_pivot.resize(m_nel);

  // ensure some of the matrices are filled with zeros
  std::fill(m_dp0.begin(), m_dp0.end(), 0.0);
  std::fill(m_dh0.begin(), m_dh0.end(), 0.0);
  std::fill(m_dh1.begin(), m_dh1.end(), 0.0);
  std::fill(m_dh2.begin(), m_dh2.end(), 0.0);
  std::fill(m_work.begin(), m_work.end(), 0.0);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluates the determinant value and adds the gradient and laplacian of the
///         log of the determinant to the total gradient and laplacian
///
/// \param[in]      P              the particle set
/// \param[in,out]  G              gradient to add to
/// \param[in,out]  L              laplacian to add to
///
/// \return  the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::ValueType SlaterDetOpt::evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {
  RealType logval = evaluateLog(P, G, L);
  // Note that this class probably isn't implemented with complex wave
  // functions yet, but I'll leave this here anyway...
#if defined(QMC_COMPLEX)
  RealType magnitude = std::exp(logval);
  return std::complex<OHMMS_PRECISION>(std::cos(PhaseValue)*magnitude, std::sin(PhaseValue)*magnitude);
#else
  return std::cos(PhaseValue)*std::exp(logval);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluates the Slater determinant's matrix, its inverse, its first derivatives w.r.t.
///         particle positions, and its summed second derivatives w.r.t. particle positions.
///         Returns the log of the determinant value.
///
/// \param[in]      P              the particle set
/// \param[in]      all            whether to additionally evaluate data for all molecular orbitals,
///                                including those not present in the determinant
///
/// \return  the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::evaluate_matrices_from_scratch(ParticleSet& P, const bool all) {

  //app_log() << " EWN ENTERING SlaterDetOpt::evaluate_matrices_from_scratch(ParticleSet& P, const bool all)" << std::endl;

  // either evaluate all of the molecular orbitals
  if (all) {
    m_spo->evaluate_notranspose(P, m_first, m_last, m_orb_val_mat_all, m_orb_der_mat_all, m_orb_lap_mat_all);
    for (int i = 0; i < m_nel; i++)
    for (int j = 0; j < m_nel; j++) {
      m_orb_val_mat(i,j) = m_orb_val_mat_all(i,j);
      m_orb_der_mat(i,j) = m_orb_der_mat_all(i,j);
      m_orb_lap_mat(i,j) = m_orb_lap_mat_all(i,j);
    }

  // or just the molecular orbitals used by this determinant
  } else {
    m_spo->evaluate_notranspose(P, m_first, m_last, m_orb_val_mat, m_orb_der_mat, m_orb_lap_mat);
  }

  // m_orb_val_mat slow (i.e. first) index is now particles
  // m_orb_der_mat slow (i.e. first) index is now particles
  // m_orb_lap_mat slow (i.e. first) index is now particles

  // copy orbital values into inverse matrix, transposing so the slow index will end up as the one we want
  qmcplusplus::MatrixOperators::transpose(m_orb_val_mat, m_orb_inv_mat); // m_orb_inv_mat slow (i.e. first) index is now combinations

  // get the log of the determinant and the inverse of the orbital value matrix
  return InvertWithLog(m_orb_inv_mat.data(), m_nel, m_nel, &m_work.at(0), &m_pivot.at(0), PhaseValue); // m_orb_inv_mat slow (i.e. first) index is now particles

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluates the log of the determinant value and adds the gradient and laplacian of the
///         log of the determinant to the total gradient and laplacian
///
/// \param[in]      P              the particle set
/// \param[in,out]  G              gradient to add to
/// \param[in,out]  L              laplacian to add to
///
/// \return  the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {

  //throw std::runtime_error("SlaterDetOpt::evaluateLog (P, G, L) die");
  LogValue = this->evaluate_matrices_from_scratch(P, false);

  // compute gradient and laplacian parts
  const bool slow_loops = true;
  if ( slow_loops ) {
    // gradient components
    for (int a = 0, iat = m_first; a < m_nel; a++, iat++)
      for (int i = 0; i < m_nel; i++)
        for (int k = 0; k < 3; k++)
          G[iat][k] += m_orb_inv_mat(a,i) * m_orb_der_mat(a,i)[k];
    // laplacian components
    for (int a = 0, iat = m_first; a < m_nel; a++, iat++) {
      for (int i = 0; i < m_nel; i++)
        L[iat] += m_orb_inv_mat(a,i) * m_orb_lap_mat(a,i);
      for (int k = 0; k < 3; k++) {
        RealType Qdot_a_mu = 0.0;
        for (int i = 0; i < m_nel; i++)
          Qdot_a_mu += m_orb_inv_mat(a,i) * m_orb_der_mat(a,i)[k];
        L[iat] -= Qdot_a_mu * Qdot_a_mu;
      }
    }
  } else {
    // TO DO:  replace slow loops with BLAS calls
  }

  // return the log of the determinant
  return LogValue;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluates the log of the determinant value and adds the gradient and laplacian of the
///         log of the determinant to the total gradient and laplacian data.
///         Currently, this class computes everything from scratch and so does not use the buffer.
///
/// \param[in]      P              the particle set
/// \param[in,out]  G              gradient to add to
/// \param[in,out]  L              laplacian to add to
/// \param[in,out]  buf            buffer to load or store temporary data from
/// \param[in]      fillBuffer     whether to fill data into the buffer or read data from it
///
/// \return  the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::evaluateLog(ParticleSet& P,
                                                ParticleSet::ParticleGradient_t& G,
                                                ParticleSet::ParticleLaplacian_t& L,
                                                PooledData<RealType>& buf,
                                                bool fillBuffer ) {
  return this->evaluateLog(P,G,L);
  //APP_ABORT("ERROR:  SlaterDetOpt::evaluateLog(P,G,L,buff,bool) not implemented\n");
//  return 0.0;
}

// BEGIN EWN DEBUG 
//template<class T> void eval_grad_print_mat_opt(T & mat) {
//  for (int i = 0; i < mat.rows(); i++) {
//    for (int j = 0; j < mat.cols(); j++)
//      app_log() << "        " << mat(i,j);
//    app_log() << std::endl;
//  }
//}
// END EWN DEBUG 

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluates and returns the gradient of the log of the determinant w.r.t. a specified
///         particle's position.
///         This function assumes that m_orb_inv_mat and m_orb_der_mat have already been prepared.
///
/// \param[in]      P              the particle set
/// \param[in]      iat            index of the particle in question
///
/// \return  the one particle gradient of the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::GradType SlaterDetOpt::evalGrad(ParticleSet& P, int iat) {

  // TO DO:  replace the loop(s) with BLAS calls

  // compute and return gradient w.r.t. position of particle iat
  OrbitalBase::GradType g;
  g = 0.0;
  if ( iat >= m_first && iat < m_last ) {
    // BEGIN EWN DEBUG 
    //this->evaluate_matrices_from_scratch(P, true);
    //app_log() << "printing inverse matrix in evalGrad" << std::endl;
    //eval_grad_print_mat_opt(m_orb_inv_mat);
    //app_log() << "printing orbital grad matrix in evalGrad" << std::endl;
    //eval_grad_print_mat_opt(m_orb_der_mat);
    // END EWN DEBUG 
    const int a = iat - m_first;
    for (int i = 0; i < m_nel; i++)
      for (int k = 0; k < 3; k++)
        g[k] += m_orb_inv_mat(a,i) * m_orb_der_mat(a,i)[k];
  }
  return g;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  For a given one particle move, this function evaluates the ratio of the new and
///         old determinant values as well as the determinant's contribution to grad_iat, which is
///         the gradient (evaluated at the new position) of the log of the trial function w.r.t.
///         the moved particle's coordinates.
///
/// \param[in]      P              the particle set, which I think carries information about the particle's move
/// \param[in]      iat            the id number of the moved particle
/// \param[in,out]  grad_iat       cumulative total of the gradient of the trial function's log w.r.t. the moved particle's coordinates
///
/// \return  the ratio of new and old determinant values
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::ValueType SlaterDetOpt::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) {

  // no contribution if the particle is not part of this determinant
  if ( iat < m_first || iat >= m_last )
    return 1.0;

  // compute orbital values, gradients, and summed second derivatives for the particle's new position
  m_spo->evaluate(P, iat, m_orb_val_vec, m_orb_der_vec, m_orb_lap_vec);

  // compute the ratio of new to old determinant values
  curRatio = simd::dot(m_orb_inv_mat[iat-m_first], m_orb_val_vec.data(), m_nel);

  // compute the determinant's contribution to the gradient of the log of the trial function
  grad_iat += simd::dot(m_orb_inv_mat[iat-m_first], m_orb_der_vec.data(), m_nel) / curRatio;

  // return the ratio
  return curRatio;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
/// \return  ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::ValueType SlaterDetOpt::ratio(ParticleSet& P,
                                           int iat,
                                           ParticleSet::ParticleGradient_t& dG,
                                           ParticleSet::ParticleLaplacian_t& dL) {

  // no contribution if the particle is not part of this determinant
  if ( iat < m_first || iat >= m_last )
    return 1.0;

  // compute orbital values, gradients, and summed second derivatives for the particle's new position
  m_spo->evaluate(P, iat, m_orb_val_vec, m_orb_der_vec, m_orb_lap_vec);

  // compute the ratio of new to old determinant values
  curRatio = simd::dot(m_orb_inv_mat[iat-m_first], m_orb_val_vec.data(), m_nel);
  
  return curRatio;

//  throw std::runtime_error("SlaterDetOpt::ratio(P, iat, dG, dL) not implemented");
//  return 0.0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
/// \return  ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::ValueType SlaterDetOpt::ratio(ParticleSet& P, int iat) {
//  throw std::runtime_error("SlaterDetOpt::ratio(P, iat) not implemented");
//  return 0.0;
  // no contribution if the particle is not part of this determinant
  if ( iat < m_first || iat >= m_last )
    return 1.0;

  // compute orbital values, gradients, and summed second derivatives for the particle's new position
  m_spo->evaluate(P, iat, m_orb_val_vec, m_orb_der_vec, m_orb_lap_vec);

  // compute the ratio of new to old determinant values
  curRatio = simd::dot(m_orb_inv_mat[iat-m_first], m_orb_val_vec.data(), m_nel);
  return curRatio;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Updates the orbital matrix and its inverse for a single particle move.
///
/// \param[in]      P              the particle set, which I think carries information about the particle's move
/// \param[in]      iat            the id number of the moved particle
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::acceptMove(ParticleSet& P, int iat) {

  // nothing to do if particle does not belong to this determinant
  if ( iat < m_first || iat >= m_last )
    return;

  // Update Log and Phase values, then reset the ratio to 1.
  PhaseValue += evaluatePhase(curRatio);
  LogValue += std::log(std::abs(curRatio));
  curRatio = 1.0;

  // get row difference
  std::copy(m_orb_val_vec.begin(), m_orb_val_vec.end(), &m_work.at(0));
  BLAS::axpy(m_nel, -1.0, m_orb_val_mat[iat-m_first], &m_work.at(0));

  // get inverse times row difference (remembering that m_orb_inv_mat slow index is particles)
  BLAS::gemm('N', 'N', 1, m_nel, m_nel, 1.0, &m_work.at(0), 1, m_orb_inv_mat.data(), m_nel, 0.0, &m_work.at(m_nel), 1);

  // get inverse ratio
  const double ir = 1.0 / ( 1.0 + m_work.at(m_nel+iat-m_first) );

  // perform rank one update of the inverse matrix (keeping in mind that it is stored with particles as the slow index)
  std::copy(m_orb_inv_mat[iat-m_first], m_orb_inv_mat[iat-m_first] + m_nel, &m_work.at(0));
  BLAS::ger(m_nel, m_nel, -ir, &m_work.at(0), 1, &m_work.at(m_nel), 1, m_orb_inv_mat.data(), m_nel);

  // update the orbital matrix, orbital derivative matrix, and summed laplacian matrix
  std::copy(m_orb_val_vec.begin(), m_orb_val_vec.end(), m_orb_val_mat[iat-m_first]);
  std::copy(m_orb_der_vec.begin(), m_orb_der_vec.end(), m_orb_der_mat[iat-m_first]);
  std::copy(m_orb_lap_vec.begin(), m_orb_lap_vec.end(), m_orb_lap_mat[iat-m_first]);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Particle move rejected - simply reset the ratio of new to old wave functions
///
/// \param[in]      iat            the id number of the moved particle
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::restore(int iat)
{
  curRatio=1.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::update(ParticleSet& P,
                          ParticleSet::ParticleGradient_t& dG,
                          ParticleSet::ParticleLaplacian_t& dL,
                          int iat) {

  throw std::runtime_error("SlaterDetOpt::update(P, dG, dL, iat) not implemented");

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
/// \return  ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::evaluateLog(ParticleSet& P,BufferType& buf) {
  throw std::runtime_error("SlaterDetOpt::evaluateLog(P, buf) not implemented");
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Initialize the object for the given particle positions and add its essential internal
///         data to the supplied buffer
///
/// \param[in]      P              the particle set
/// \param[in]      buf            the buffer to add essential data to
///
/// \return  the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::registerData(ParticleSet& P, BufferType& buf) {

  //app_log() << " EWN ENTERING SlaterDetOpt::registerData(ParticleSet& P, BufferType& buf)" << std::endl;

  // evaluate the orbital matrix, its inverse, and the log of the determinant value
  LogValue = evaluateLog(P,P.G,P.L);

  // add the log of the wave function to the buffer
  buf.add(LogValue);

  // add the phase of the wave function to the buffer
  buf.add(PhaseValue);

  // add orbital matrix to the buffer
  buf.add(m_orb_val_mat.first_address(), m_orb_val_mat.last_address());

  // add the orbital derivative matrix to the buffer
  buf.add(&m_orb_der_mat(0,0)[0], &m_orb_der_mat(0,0)[0] + m_orb_der_mat.size() * m_orb_der_mat(0,0).size());

  // add orbital laplacian matrix to the buffer
  buf.add(m_orb_lap_mat.first_address(), m_orb_lap_mat.last_address());

  // add inverse matrix to the buffer
  buf.add(m_orb_inv_mat.first_address(),m_orb_inv_mat.last_address());

  // return log of determinant value
  return LogValue;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Stores data that will be needed later in computing derivative ratio information
///         for the linear method in the provided buffer.
///         Currently, this class recomputes everything on the fly and so stores no data here.
///
/// \param[in]      P              the particle set
/// \param[in]      buf            the buffer to save data in
/// \param[in]      storageType    a storage flag not used for this class
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::registerDataForDerivatives(ParticleSet& P, BufferType& buf, int storageType) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Initialize the object for the given particle positions and put its essential internal
///         data into the supplied buffer
///
/// \param[in]      P              the particle set
/// \param[in]      buf            the buffer to save essential data in
/// \param[in]      fromscratch    ???
///
/// \return  the log of the determinant value
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBase::RealType SlaterDetOpt::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch) {

  // BEGIN EWN DEBUG 
  //app_log() << " EWN ENTERING SlaterDetOpt::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch)" << std::endl;
  //app_log() << "printing P.R" << std::endl;
  //for (int i = 0; i < P.R.size(); i++)
  //  app_log() << P.R[i] << std::endl;
  // END EWN DEBUG 

  // !!! WE WILL EVENTUALLY WANT TO TAKE DIFFERENT ACTIONS DEPENDING ON THE VALUE OF fromscratch !!!

  // evaluate the orbital matrix, its inverse, and the log of the determinant value
  LogValue = evaluateLog(P,P.G,P.L);

  // add the log of the wave function to the buffer
  buf.put(LogValue);

  // add the phase of the wave function to the buffer
  buf.put(PhaseValue);

  // add orbital matrix to the buffer
  buf.put(m_orb_val_mat.first_address(), m_orb_val_mat.last_address());

  // add the orbital derivative matrix to the buffer
  buf.put(&m_orb_der_mat(0,0)[0], &m_orb_der_mat(0,0)[0] + m_orb_der_mat.size() * m_orb_der_mat(0,0).size());

  // add orbital laplacian matrix to the buffer
  buf.put(m_orb_lap_mat.first_address(), m_orb_lap_mat.last_address());

  // add inverse matrix to the buffer
  buf.put(m_orb_inv_mat.first_address(),m_orb_inv_mat.last_address());

  // return log of determinant value
  return LogValue;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Read in necessary internal data from the provided buffer
///
/// \param[in]      P              the particle set
/// \param[in]      buf            the buffer to read from
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::copyFromBuffer(ParticleSet& P, BufferType& buf) {

  // read in the log of the wave function
  buf.get(LogValue);

  // read in the phase of the wave function
  buf.get(PhaseValue);

  // read in the orbital matrix
  buf.get(m_orb_val_mat.first_address(), m_orb_val_mat.last_address());

  // read in the orbital derivative matrix
  buf.get(&m_orb_der_mat(0,0)[0], &m_orb_der_mat(0,0)[0] + m_orb_der_mat.size() * m_orb_der_mat(0,0).size());

  // read in the orbital laplacian matrix
  buf.get(m_orb_lap_mat.first_address(), m_orb_lap_mat.last_address());

  // read in the inverse matrix
  buf.get(m_orb_inv_mat.first_address(), m_orb_inv_mat.last_address());

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Nothing to check in, as the optimizable orbital variables are handled by the
///         orbital set's trial function component.
///
/// \param[in]      active         ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::checkInVariables(opt_variables_type& active) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Nothing to check out, as the optimizable orbital variables are handled by the
///         orbital set's trial function component.
///
/// \param[in]      active         ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::checkOutVariables(const opt_variables_type& active) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Nothing to reset, as the optimizable orbital variables are handled by the
///         orbital set's trial function component.
///
/// \param[in]      active         ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::resetParameters(const opt_variables_type& active) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::reportStatus(std::ostream& os) {
  throw std::runtime_error("SlaterDetOpt::reportStatus(os) not implemented");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluate the determinant's contribution to the derivatives of log of the trial function
///         and the local energy w.r.t. the molecular orbital coefficients, which will be used
///         by the orbital set's trial function component determine the corresponding derivatives
///         w.r.t. rotations between its molecular orbitals.
///
/// \param[in]      P              the particle set
/// \param[in]      optvars        unused here as we are just preparing what m_spo will need
/// \param[in,out]  dlogpsi        unused here as we are just preparing what m_spo will need
/// \param[in,out]  dhpsioverpsi   unused here as we are just preparing what m_spo will need
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::evaluateDerivatives(ParticleSet& P,
                                       const opt_variables_type& optvars,
                                       std::vector<RealType>& dlogpsi,
                                       std::vector<RealType>& dhpsioverpsi) {


  // Evaluate orbital data for all orbitals (not just those in this determinant).
  // Prepares:  m_orb_val_mat_all, m_orb_der_mat_all, m_orb_lap_mat_all, m_orb_val_mat, m_orb_der_mat, m_orb_lap_mat, and m_orb_inv_mat
  LogValue = this->evaluate_matrices_from_scratch(P, true);

  // print for debugging
  if ( false ) {

    std::vector<char> buff(1000, ' ');

    app_log() << "printing m_orb_val_mat_all" << std::endl;
    for (int p = 0; p < m_nmo; p++) { // loop over orbitals
      for (int a = 0; a < m_nel; a++) { // loop over particles
        const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_val_mat_all(a,p));
        for (int k = 0; k < len; k++)
          app_log() << buff[k];
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

    app_log() << "printing m_orb_der_mat_all" << std::endl;
    for (int p = 0; p < m_nmo; p++) { // loop over orbitals
      for (int a = 0; a < m_nel; a++) { // loop over particles
        for (int k = 0; k < 3; k++) { // loop over x,y,z directions
          const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_der_mat_all(a,p)[k]);
          for (int k = 0; k < len; k++)
            app_log() << buff[k];
        }
        app_log() << "        ";
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

    app_log() << "printing m_orb_lap_mat_all" << std::endl;
    for (int p = 0; p < m_nmo; p++) { // loop over orbitals
      for (int a = 0; a < m_nel; a++) { // loop over particles
        const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_lap_mat_all(a,p));
        for (int k = 0; k < len; k++)
          app_log() << buff[k];
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

    app_log() << "printing P.G" << std::endl;
    for (int a = 0; a < m_nel; a++) { // loop over particles
      for (int k = 0; k < 3; k++) { // loop over x,y,z directions
        const int len = std::sprintf(&buff[0], "  %12.6f", P.G[m_first+a][k]);
        for (int k = 0; k < len; k++)
          app_log() << buff[k];
      }
      app_log() << "        ";
    }
    app_log() << std::endl;
    app_log() << std::endl;

    app_log() << "printing m_orb_inv_mat" << std::endl;
    for (int b = 0; b < m_nel; b++) {
      for (int a = 0; a < m_nel; a++) {
        const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_inv_mat(a,b));
        for (int k = 0; k < len; k++)
          app_log() << buff[k];
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

    app_log() << "printing m_orb_der_mat" << std::endl;
    for (int p = 0; p < m_nel; p++) {
      for (int a = 0; a < m_nel; a++) {
        for (int k = 0; k < 3; k++) {
          const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_der_mat(a,p)[k]);
          for (int k = 0; k < len; k++)
            app_log() << buff[k];
        }
        app_log() << "        ";
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

    app_log() << "printing m_orb_lap_mat" << std::endl;
    for (int p = 0; p < m_nel; p++) {
      for (int a = 0; a < m_nel; a++) {
        const int len = std::sprintf(&buff[0], "  %12.6f", m_orb_lap_mat(a,p));
        for (int k = 0; k < len; k++)
          app_log() << buff[k];
      }
      app_log() << std::endl;
    }
    app_log() << std::endl;

  }

  // fill matrix of contributions to derivatives of log of det value w.r.t. molecular orbital values
  for (int a = 0; a < m_nel; a++) // loop over particles
  for (int p = 0; p < m_nel; p++) // loop over orbitals
    m_dp0(a,p) = m_orb_inv_mat(a,p);

  // construct temporary Y matrix
  RealType * const Ymat = &m_work.at(0);
  for (int b = 0; b < m_nel; b++) { // loop over particles
    const OrbitalBase::GradType g = P.G[m_first+b] - simd::dot(m_orb_inv_mat[b], m_orb_der_mat[b], m_nel);
    for (int q = 0; q < m_nel; q++) // loop over orbitals
      Ymat[q+b*m_nel] = 0.5 * m_orb_lap_mat(b,q) + qmcplusplus::dot(m_orb_der_mat(b,q), g);
  }

  // contract Y with inverse matrices to get contribution of local energy derivatives w.r.t. orbital values
  BLAS::gemm('N', 'T', m_nel, m_nel, m_nel, 1.0,    m_orb_inv_mat.data(), m_nel,                 Ymat, m_nel, 0.0, &m_work.at(m_nel*m_nel), m_nel);
  BLAS::gemm('N', 'N', m_nel, m_nel, m_nel, 1.0, &m_work.at(m_nel*m_nel), m_nel, m_orb_inv_mat.data(), m_nel, 0.0,           &m_work.at(0), m_nel);

  // fill result of contraction into top of local energy derivatives w.r.t. molecular orbital value matrix (derivatives w.r.t. virtual orbitals are zero and so we leave the bottom of the matrix alone)
  for (int b = 0; b < m_nel; b++) // loop over particles
  for (int q = 0; q < m_nel; q++) // loop over orbitals
    m_dh0(b,q) = m_work[ q + b * m_nel ];

  // fill matrices of contributions to local energy derivatives w.r.t. orbital first derivatives
  for (int a = 0; a < m_nel; a++) { // loop over particles
    const OrbitalBase::GradType g = simd::dot(m_orb_inv_mat[a], m_orb_der_mat[a], m_nel) - P.G[m_first+a];
    for (int v = 0; v < 3; v++) // loop over particle coordinates x,y,z
    for (int p = 0; p < m_nel; p++) // loop over orbitals
      m_dh1(a+v*m_nel,p) = m_orb_inv_mat(a,p) * g[v];
  }

  // fill matrix of contributions to local energy derivatives w.r.t. orbital second derivatives
  for (int a = 0; a < m_nel; a++) // loop over particles
  for (int p = 0; p < m_nel; p++) // loop over orbitals
      m_dh2(a,p) = -0.5 * m_orb_inv_mat(a,p);

  // use the work matrix to arrange the molecular orbital derivative data in the order needed
  for (int a = 0; a < m_nel; a++) // loop over particles
  for (int p = 0; p < m_nmo; p++) // loop over orbitals
  for (int k = 0; k < 3; k++) // loop over particle coordinates x,y,z
    m_work[ p + a*m_nmo + k*m_nmo*m_nel ] = m_orb_der_mat_all(a,p)[k];

  // add this determinant's contribution to the orbital linear combinations' derivatives
  this->tfc_ptr("SlaterDetOpt::evaluateDerivatives")->add_derivatives(m_nmo, m_nel, m_dp0.data(), m_dh0.data(), m_dh1.data(), m_dh2.data(),
                                                                      m_orb_val_mat_all.data(), &m_work.at(0), m_orb_lap_mat_all.data());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Return a pointer to the single particle orbital set's trial function component.
///
/// \param[in]      calling_func   name of calling function for use in error reporting
///
/// \return  a pointer to the single particle orbital set's trial function component.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
LCOrbitalSetOptTrialFunc * SlaterDetOpt::tfc_ptr(const std::string & calling_func) {

  // get pointer to the optimizable single particle orbital set's trial function component
  LCOrbitalSetOptTrialFunc * const ptr = dynamic_cast<LCOrbitalSetOptTrialFunc*>(m_spo->tf_component());

  // check that the pointer conversion was successful
  if ( !ptr ) {
    std::stringstream message;
    message << "dynamic_cast failure for trial function component pointer in " << calling_func;
    throw std::runtime_error(message.str());
  }

  // return the pointer
  return ptr;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Sets the single particle orbital set's optimizable rotations to be only those
///         between occupied and virtual molecular orbitals.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::set_spo_optimizable_rotations() {

  // add this determinant's contribution to the orbital linear combinations' derivatives
  this->tfc_ptr("SlaterDetOpt::set_spo_optimizable_rotations")->set_optimizable_rotation_ranges(0, m_nel, m_nel, m_nmo);

}

}

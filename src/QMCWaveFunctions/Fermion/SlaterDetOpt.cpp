//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//                    Nick Blunt, nicksblunt@gmail.com, University of Cambridge
//
// File created by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//////////////////////////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/Fermion/SlaterDetOpt.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <QMCWaveFunctions/LCOrbitalSetOpt.h>
#include <Numerics/DeterminantOperators.h>
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus {

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Creates a slater determinant of the given spin with the given optimizable orbital set.
///
/// \param[in]      ptcl           the particle set
/// \param[in]      spo_ptr        pointer to the optimizable single particle orbital set
/// \param[in]      up_or_down     0 for up spin, 1 for down
/// \param[in]      nmo            number of optimizable molecular orbitals
///
///////////////////////////////////////////////////////////////////////////////////////////////////
SlaterDetOpt::SlaterDetOpt(ParticleSet & ptcl, SPOSetBase * spo_ptr, const int up_or_down)
  : DiracDeterminantBase(spo_ptr, ptcl.first(up_or_down))
  , m_up_or_down(up_or_down)
  , m_nmo(spo_ptr->size())
  , m_first_var_pos(-1)
  , m_act_rot_inds()
{
  targetPtcl = &ptcl;

  Optimizable=true;
  OrbitalName="SlaterDetOpt";
  this->resetTargetParticleSet(*targetPtcl);

  m_nlc = Phi->OrbitalSetSize;
  m_nb = Phi->BasisSetSize;

  // make sure we didn't start with a bad m_nlc
  check_index_sanity();

  // by default set all rotations to be active
  m_act_rot_inds.resize( m_nlc * ( m_nlc - 1 ) / 2 );
  int rots_recorded = 0;
  for (int j = 1; j < m_nlc; j++)
  for (int i = 0; i < j; i++)
    m_act_rot_inds.at(rots_recorded++) = std::pair<int,int>(i,j);
  if ( m_act_rot_inds.size() != rots_recorded )
    throw std::runtime_error("wrong number of active rotations recorded in SlaterDetOpt constructor.");

  // make sure we didn't do something stupid
  check_index_sanity();

  // prepare matrices that will hold derivatives wrt orbital rotations
  this->initialize_matrices();

  // add this determinant's contribution to the orbital linear combinations' derivatives
  set_optimizable_rotation_ranges(0, m_nel, m_nel, m_nmo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Adds the optimizable orbital set trial function component to the given trial function
///         if it wasn't already there.
///
/// \param[in]      twf            the trial wave function to add the component to
/// \param[in]      name           a name for the component
///
///////////////////////////////////////////////////////////////////////////////////////////////////
//void SlaterDetOpt::add_orbs_to_tf(TrialWaveFunction & twf, const std::string & name) {
//  if ( std::find(twf.getOrbitals().begin(), twf.getOrbitals().end(), Phi->tf_component()) == twf.getOrbitals().end() )
//    twf.addOrbital(Phi->tf_component(), name, false);
//}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  verifies that the number of linear combinations and the list of active rotation
///         indices is sane
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::check_index_sanity() const {

  // ensure the number of linear combinations is not negative
  if ( m_nlc < 0 )
    throw std::runtime_error("check_index_sanity found a negative number of linear combinations");

  // throw an error if any active rotation index is unreasonable
  for (std::vector<std::pair<int,int> >::const_iterator it = m_act_rot_inds.begin(); it != m_act_rot_inds.end(); it++) {
    if ( it->first >= it->second ) {
      std::stringstream error_msg;
      error_msg << "check_index_sanity found an active rotation index pair ("
                << it->first << "," << it->second << ") in which the first index was not smaller than the second";
      throw std::runtime_error(error_msg.str());
    }
    if ( it->first < 0 || it->first >= m_nlc ) {
      std::stringstream error_msg;
      error_msg << it->first << " is an out of bounds first active rotation index in check_index_sanity";
      throw std::runtime_error(error_msg.str());
    }
    if ( it->second < 0 || it->second >= m_nlc ) {
      std::stringstream error_msg;
      error_msg << it->second << " is an out of bounds second active rotation index in check_index_sanity";
      throw std::runtime_error(error_msg.str());
    }
  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Makes a clone of the object that uses the supplied particle set.
///
/// \param[in]      tqp            the particle set the clone should use
///
///////////////////////////////////////////////////////////////////////////////////////////////////
OrbitalBasePtr SlaterDetOpt::makeClone(ParticleSet& tqp) const {
  SlaterDetOpt* clone = new SlaterDetOpt(tqp, Phi->makeClone(), m_up_or_down);

  clone->Optimizable=Optimizable;
  clone->myVars=myVars;
  clone->m_first_var_pos = m_first_var_pos;

  return clone;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Makes a clone (copy) of the object that uses the supplied single
///         particle orbital set.
///
/// \param[in]      spo       the single particle orbital set the copy should use
///
///////////////////////////////////////////////////////////////////////////////////////////////////
DiracDeterminantBase* SlaterDetOpt::makeCopy(SPOSetBasePtr spo) const
{
  SlaterDetOpt* copy = new SlaterDetOpt(*targetPtcl, spo, m_up_or_down);

  copy->myVars=myVars;
  copy->Optimizable=Optimizable;
  copy->m_first_var_pos = m_first_var_pos;

  return copy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Destructor has nothing to do for now
///
///////////////////////////////////////////////////////////////////////////////////////////////////
SlaterDetOpt::~SlaterDetOpt() { }

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  prepares the matrices that will hold derivatives w.r.t. orbital rotations by
///         ensuring they are the right size and that their elements are all zero
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::initialize_matrices() {

  // initialize matrix for Log(Psi) derivatives
  m_pder_mat.resize(m_nlc * m_nlc);
  std::fill(m_pder_mat.begin(), m_pder_mat.end(), 0.0);

  // initialize matrix for ( H Psi ) / Psi derivatives
  m_hder_mat.resize(m_nlc * m_nlc);
  std::fill(m_hder_mat.begin(), m_hder_mat.end(), 0.0);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Exponentiates a matrix
///
/// \param[in]      n              matrix dimensions
/// \param[in,out]  mat            On entry, the n by n matrix.
///                                On exit, the exponential of the matrix.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::exponentiate_matrix(const int n, RealType * const mat) {

  // save initial matrix and get some workspace
  std::vector<RealType> mat_a(mat, mat + n*n);
  std::vector<RealType> mat_b(mat, mat + n*n);
  std::vector<RealType> mat_c(mat, mat + n*n);
  RealType * ptr_a = &mat_a.at(0);
  RealType * ptr_b = &mat_b.at(0);
  RealType * ptr_c = &mat_c.at(0);

  // initialize output to identity matrix
  for (int j = 0; j < n; j++)
  for (int i = 0; i < n; i++)
    mat[i+n*j] = ( i == j ? 1.0 : 0.0 );

  // compute exponential of matrix
  for (int q = 1; q < 20; q++) {
    BLAS::axpy(n*n, RealType(1.0), ptr_b, mat);
    BLAS::gemm('N', 'N', n, n, n, RealType(1.0) / (q+1), ptr_a, n, ptr_b, n, RealType(0.0), ptr_c, n);
    std::swap(ptr_b, ptr_c);
    RealType max_elem = 0.0;
    for (int i = 0; i < n*n; i++)
      if ( std::abs(ptr_b[i]) > max_elem )
        max_elem = std::abs(ptr_b[i]);
    if ( max_elem < 1.0e-15 )
      break;
  }

}

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
  Phi->resetTargetParticleSet(P);
  targetPtcl = &P;

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
    Phi->evaluate_notranspose(P, m_first, m_last, m_orb_val_mat_all, m_orb_der_mat_all, m_orb_lap_mat_all);
    for (int i = 0; i < m_nel; i++)
    for (int j = 0; j < m_nel; j++) {
      m_orb_val_mat(i,j) = m_orb_val_mat_all(i,j);
      m_orb_der_mat(i,j) = m_orb_der_mat_all(i,j);
      m_orb_lap_mat(i,j) = m_orb_lap_mat_all(i,j);
    }

  // or just the molecular orbitals used by this determinant
  } else {
    Phi->evaluate_notranspose(P, m_first, m_last, m_orb_val_mat, m_orb_der_mat, m_orb_lap_mat);
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
  Phi->evaluate(P, iat, m_orb_val_vec, m_orb_der_vec, m_orb_lap_vec);

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
OrbitalBase::ValueType SlaterDetOpt::ratio(ParticleSet& P, int iat) {
//  throw std::runtime_error("SlaterDetOpt::ratio(P, iat) not implemented");
//  return 0.0;
  // no contribution if the particle is not part of this determinant
  if ( iat < m_first || iat >= m_last )
    return 1.0;

  // compute orbital values, gradients, and summed second derivatives for the particle's new position
  Phi->evaluate(P, iat, m_orb_val_vec, m_orb_der_vec, m_orb_lap_vec);

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
/// \brief  Initialize the object for the given particle positions and add its essential internal
///         data to the supplied buffer
///
/// \param[in]      P              the particle set
/// \param[in]      buf            the buffer to add essential data to
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::registerData(ParticleSet& P, WFBufferType& buf) {

  //app_log() << " EWN ENTERING SlaterDetOpt::registerData(ParticleSet& P, WFBufferType& buf)" << std::endl;

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

}

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
OrbitalBase::RealType SlaterDetOpt::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch) {

  // BEGIN EWN DEBUG 
  //app_log() << " EWN ENTERING SlaterDetOpt::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)" << std::endl;
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
void SlaterDetOpt::copyFromBuffer(ParticleSet& P, WFBufferType& buf) {

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
/// \brief  Check this object's optimizable variables into the supplied overall list of
///         optimizable variables.
///
/// \param[in,out]  active        the overall list of optimizable variables
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::checkInVariables(opt_variables_type& active) {

  // add these variables to the overall list of optimizable variables
  active.insertFrom(this->myVars);

  // reset my first variable's position to say that we don't know where it is yet
  m_first_var_pos = -1;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Records the positions of this object's optimizable variables in the supplied overall
///         list of optimizable variables, checks that the variables are stored contiguously,
///         and records the position of the first of this objects variables.
///
/// \param[in]      active        the overall list of optimizable variables
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::checkOutVariables(const opt_variables_type& active) {

  // record the positions of this object's optimizable variables within the overall list of optimizable variables
  myVars.getIndex(active);

  // ensure that this object's variables are stored contiguously
  for (int i = 0; i < myVars.size(); i++) {
    if ( myVars.where(i) - myVars.where(0) != i ) {
      std::stringstream error_msg;
      error_msg << "variable " << (i-1) << " was not contiguous with variable " << i << " in SlaterDetOpt::checkOutVariables";
      throw std::runtime_error(error_msg.str());
    }
  }

  // record the position of my first variable
  if ( myVars.size() > 0 )
    m_first_var_pos = myVars.where(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Nothing to reset, as the optimizable orbital variables are handled by the
///         orbital set's trial function component.
///
/// \param[in]      active         ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::resetParameters(const opt_variables_type& active)
{
  // read out the parameters that define the rotation into an antisymmetric matrix
  std::vector<RealType> rot_mat(m_nlc*m_nlc, 0.0);
  for (int i = 0; i < m_act_rot_inds.size(); i++) {
    const int p = m_act_rot_inds.at(i).first;
    const int q = m_act_rot_inds.at(i).second;
    //const RealType x = active[i + m_first_var_pos] - myVars[i];
    const RealType x = active[i + m_first_var_pos];
    rot_mat[p+q*m_nlc] =  x;
    rot_mat[q+p*m_nlc] = -x;
  }

  // exponentiate antisymmetric matrix to get the unitary rotation
  this->exponentiate_matrix(m_nlc, &rot_mat.at(0));

  // get the linear combination coefficients by applying the rotation to the old coefficients
  Phi->rotate_B(rot_mat);

  // Store the orbital rotations parameters internally in myVars
  for (int i = 0; i < m_act_rot_inds.size(); i++)
    myVars[i] = active[i + m_first_var_pos];

  //if (false)
  //  this->print_B();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Not yet implemented.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::reportStatus(std::ostream& os) {
  throw std::runtime_error("SlaterDetOpt::reportStatus(os) not implemented");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  add to the Log(Psi) and ( H Psi ) / Psi derivatives
///
/// \param[in]      nl         The number of molecular orbitals.
/// \param[in]      np         The number of particles over which to sum derivative contributions.
/// \param[in]      dp0        An nl by np column-major-ordered matrix of the derivatives
///                            of Log(Psi) with respect to the values of the
///                            molecular orbitals at each particle's position.
/// \param[in]      dh0        An nl by np column-major-ordered matrix of the derivatives
///                            of ( H Psi ) / Psi with respect to the values of the
///                            molecular orbitals at each particle's position.
/// \param[in]      dh1        Three nl by np column-major-ordered matrices (stored contiguously
///                            one after the other) of the derivatives of
///                            of ( H Psi ) / Psi with respect to the values of the
///                            molecular orbitals' first position derivatives (w.r.t. x,y,z)
///                            at each particle's position.
/// \param[in]      dh2        An nl by np column-major-ordered matrix of the derivatives
///                            of ( H Psi ) / Psi with respect to the values of the
///                            molecular orbitals' second position derivatives at each
///                            particle's position.  Note that we assume the derivatives of
///                            ( H Psi ) / Psi are the same for each of the three directions'
///                            (x,y,z) second derivatives and so dh2 is defined as the
///                            derivaties corresponding to the x coordinate's second derivative,
///                            NOT the sum of the derivatives for all three x, y, and z.
/// \param[in]      Bchi       An nl by np column-major-ordered matrix of the values of the
///                            molecular orbitals at each particle's position.
/// \param[in]      dBchi      Three nl by np column-major-ordered matrices (stored contiguously
///                            one after the other) of the first position derivatives of the
///                            molecular orbitals at each particle's position.
/// \param[in]      d2Bchi     An nl by np column-major-ordered matrix of the molecular orbitals'
///                            x-y-z summed second derivatives at each particle's position.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::add_derivatives(const int nl,
                                   const int np,
                                   const RealType * const dp0,
                                   const RealType * const dh0,
                                   const RealType * const dh1,
                                   const RealType * const dh2,
                                   const RealType * const Bchi,
                                   const RealType * const dBchi,
                                   const RealType * const d2Bchi) {

  // ensure the number of linear combinations is correct
  if ( nl != m_nlc ) {
    std::stringstream error_msg;
    error_msg << "supplied number of linear combinations (" << nl << ") does not match that held internally (" << m_nlc << ") in add_derivatives";
    throw std::runtime_error(error_msg.str());
  }

  // ensure orbital derivative matrices are the correct size
  if ( m_pder_mat.size() != nl * nl ) {
    std::stringstream error_msg;
    error_msg << "nl (" << nl << ") does not match size of m_pder_mat (" << m_pder_mat.size() << ") in add_derivatives";
    throw std::runtime_error(error_msg.str());
  }
  if ( m_hder_mat.size() != nl * nl ) {
    std::stringstream error_msg;
    error_msg << "nl (" << nl << ") does not match size of m_hder_mat (" << m_hder_mat.size() << ") in add_derivatives";
    throw std::runtime_error(error_msg.str());
  }

  // contract Log(Psi) derivatives with the current linear combination values and add result to the Log(Psi) derivatives w.r.t. the matrix C
  //for (int j = 0; j < nl; j++)
  //for (int i = 0; i < nl; i++)
  //  for (int k = 0; k < np; k++)
  //    m_pder_mat.at(i+j*nl) += dp0[i+k*nl] * Bchi[j+k*nl];
  BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dp0, nl, Bchi, nl, RealType(1.0), &m_pder_mat.at(0), nl);

  // compute products of ( H Psi ) / Psi derivatives with linear combination values and their derivatives and add results to energy derivatives w.r.t. the matrix C
  BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh0, nl, Bchi, nl, RealType(1.0), &m_hder_mat.at(0), nl);
  for (int i = 0; i < 3; i++)
    BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh1+i*nl*np, nl, dBchi+i*nl*np, nl, RealType(1.0), &m_hder_mat.at(0), nl);
  BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh2, nl, d2Bchi, nl, RealType(1.0), &m_hder_mat.at(0), nl);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  add to the \grad(\textrm{log}(Psi)) derivatives
///
/// \param[in]      nl         The number of molecular orbitals.
/// \param[in]      np         The number of particles over which to sum derivative contributions.
/// \param[in]      dh0        An nl by np column-major-ordered matrix of the derivatives
///                            of \grad(\textrm{log}(Psi)) with respect to the values of the
///                            molecular orbitals at each particle's position.
/// \param[in]      dh1        Three nl by np column-major-ordered matrices (stored contiguously
///                            one after the other) of the derivatives of
///                            of \grad(\textrm{log}(Psi)) with respect to the values of the
///                            molecular orbitals' first position derivatives (w.r.t. x,y,z)
///                            at each particle's position.
/// \param[in]      Bchi       An nl by np column-major-ordered matrix of the values of the
///                            molecular orbitals at each particle's position.
/// \param[in]      dBchi      Three nl by np column-major-ordered matrices (stored contiguously
///                            one after the other) of the first position derivatives of the
///                            molecular orbitals at each particle's position.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::add_grad_derivatives(const int nl,
                          const int np,
                          const RealType * const dh0,
                          const RealType * const dh1,
                          const RealType * const Bchi,
                          const RealType * const dBchi) {

  // ensure the number of linear combinations is correct
  if ( nl != m_nlc ) {
    std::stringstream error_msg;
    error_msg << "supplied number of linear combinations (" << nl << ") does not match that held internally (" << m_nlc << ") in LCOrbitalSetOptTrialFunc::add_grad_derivatives";
    throw std::runtime_error(error_msg.str());
  }

  // ensure orbital derivative matrices are the correct size
  if ( m_hder_mat.size() != nl * nl ) {
    std::stringstream error_msg;
    error_msg << "nl (" << nl << ") does not match size of m_hder_mat (" << m_hder_mat.size() << ") in LCOrbitalSetOptTrialFunc::add_grad_derivatives";
    throw std::runtime_error(error_msg.str());
  }

  BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh0, nl, Bchi, nl, RealType(1.0), &m_hder_mat.at(0), nl);
  for (int i = 0; i < 3; i++)
    BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh1+i*nl*np, nl, dBchi+i*nl*np, nl, RealType(1.0), &m_hder_mat.at(0), nl);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluate the determinant's contribution to the derivatives of log of the trial function
///         and the local energy w.r.t. the molecular orbital coefficients, which will be used
///         by the orbital set's trial function component determine the corresponding derivatives
///         w.r.t. rotations between its molecular orbitals.
///
/// \param[in]      P              the particle set
/// \param[in]      optvars        unused here as we are just preparing what Phi will need
/// \param[in,out]  dlogpsi        unused here as we are just preparing what Phi will need
/// \param[in,out]  dhpsioverpsi   unused here as we are just preparing what Phi will need
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
  add_derivatives(m_nmo, m_nel, m_dp0.data(), m_dh0.data(), m_dh1.data(), m_dh2.data(),
                  m_orb_val_mat_all.data(), &m_work.at(0), m_orb_lap_mat_all.data());



  // check that we have the position of the first of our variables in the overall list
  if ( myVars.size() > 0 && m_first_var_pos < 0 )
    throw std::runtime_error("position of first variable was not set on entry to SlaterDetOpt::evaluateDerivatives");

  // check that my number of variables is consistent with the number of active rotations
  if ( myVars.size() != m_act_rot_inds.size() ) {
    std::stringstream error_msg;
    error_msg << "mismatch between myVars.size() (" << myVars.size() << ") and m_act_rot_inds.size() (" << m_act_rot_inds.size() << ") in SlaterDetOpt::evaluateDerivatives";
    throw std::runtime_error(error_msg.str());
  }

  for (int i = 0; i < m_act_rot_inds.size(); i++) {
    const int p = m_act_rot_inds.at(i).first;
    const int q = m_act_rot_inds.at(i).second;
    dlogpsi.at(m_first_var_pos+i) += m_pder_mat.at(p+q*m_nlc) - m_pder_mat.at(q+p*m_nlc);
    dhpsioverpsi.at(m_first_var_pos+i) += m_hder_mat.at(p+q*m_nlc) - m_hder_mat.at(q+p*m_nlc);
    if ( false ) {
      std::vector<char> buff(1000, ' ');
      const int len = std::sprintf(&buff[0], " p = %4i   q = %4i     dlogpsi = %20.12f     dhpsioverpsi = %20.12f", p, q, dlogpsi.at(m_first_var_pos+i), dhpsioverpsi.at(m_first_var_pos+i));
      for (int k = 0; k < len; k++)
        app_log() << buff[k];
      app_log() << std::endl;
    }
  }

  // reset the internally stored derivatives to zero in preperation for the next sample
  this->initialize_matrices();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Evaluate the derivatives (with respect to optimizable parameters) of the gradient of
///         the logarithm of the optimizable determinant wave function. The dot product of this
///         object is then taken (with respect to the gradient vector components) with the
///         input G_in gradient vector. The resulting object is a vector in the optimizable
///         parameter components, which is returned by this function in dgradlogpsi.
///
/// \param[in]  G_in         Some gradient vector to be dotted with d(\grad(\textrm{log}(\psi)))
///                          where \psi is just the optimiziable determinant.
/// \param[out] dgradlogpsi  The dot product of G_in with d(\grad(\textrm{log}(\psi)))
///                          (not calculated here but in the evaluateGradDerivatives function in
///                          LCOrbitalSetOpt.h, for this particular wave function ansatz).
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                        std::vector<RealType>& dgradlogpsi) {

  // construct temporary Y matrix
  RealType * const Ymat = &m_work.at(0);

  for (int b = 0; b < m_nel; b++) { // loop over particles
    const OrbitalBase::GradType g = G_in[m_first+b];
    for (int q = 0; q < m_nel; q++) // loop over orbitals
      Ymat[q+b*m_nel] = - qmcplusplus::dot(m_orb_der_mat(b,q), g);
  }

  // contract Y with inverse matrices to get contribution of gradient
  // derivatives w.r.t. orbital values
  BLAS::gemm('N', 'T', m_nel, m_nel, m_nel, 1.0,    m_orb_inv_mat.data(), m_nel,                 Ymat, m_nel, 0.0, &m_work.at(m_nel*m_nel), m_nel);
  BLAS::gemm('N', 'N', m_nel, m_nel, m_nel, 1.0, &m_work.at(m_nel*m_nel), m_nel, m_orb_inv_mat.data(), m_nel, 0.0,           &m_work.at(0), m_nel);

  // fill result of contraction into top of gradient derivatives w.r.t.
  // molecular orbital value matrix (derivatives w.r.t. virtual orbitals are
  // zero and so we leave the bottom of the matrix alone)
  for (int b = 0; b < m_nel; b++) // loop over particles
  for (int q = 0; q < m_nel; q++) // loop over orbitals
    m_dh0(b,q) = m_work[ q + b * m_nel ];

  // fill matrices of contributions to gradient derivatives w.r.t. orbital
  // first derivatives
  for (int a = 0; a < m_nel; a++) { // loop over particles
    const OrbitalBase::GradType g = G_in[m_first+a];
    for (int v = 0; v < 3; v++) // loop over particle coordinates x,y,z
      for (int p = 0; p < m_nel; p++) // loop over orbitals
        m_dh1(a+v*m_nel,p) = m_orb_inv_mat(a,p) * g[v];
  }

  // use the work matrix to arrange the molecular orbital derivative data in
  // the order needed
  for (int a = 0; a < m_nel; a++) // loop over particles
    for (int p = 0; p < m_nmo; p++) // loop over orbitals
      for (int k = 0; k < 3; k++) // loop over particle coordinates x,y,z
        m_work[ p + a*m_nmo + k*m_nmo*m_nel ] = m_orb_der_mat_all(a,p)[k];

  // add this determinant's contribution to the orbital linear combinations' derivatives
  add_grad_derivatives(m_nmo, m_nel, m_dh0.data(), m_dh1.data(), m_orb_val_mat_all.data(), &m_work.at(0));

  // check that we have the position of the first of our variables in the
  // overall list.
  if ( myVars.size() > 0 && m_first_var_pos < 0 )
    throw std::runtime_error("position of first variable was not set on entry to "
                              "LCOrbitalSetOptTrialFunc::evaluateGradDerivatives");

  // check that my number of variables is consistent with the
  // number of active rotations.
  if ( myVars.size() != m_act_rot_inds.size() ) {
    std::stringstream error_msg;
    error_msg << "mismatch between myVars.size() (" << myVars.size()
              << ") and m_act_rot_inds.size() (" << m_act_rot_inds.size()
              << ") in LCOrbitalSetOptTrialFunc::evaluateGradDerivatives";
    throw std::runtime_error(error_msg.str());
  }

  // add derivatives to totals.
  for (int i = 0; i < m_act_rot_inds.size(); i++) {
    const int p = m_act_rot_inds.at(i).first;
    const int q = m_act_rot_inds.at(i).second;
    dgradlogpsi.at(m_first_var_pos+i) += m_hder_mat.at(p+q*m_nlc) - m_hder_mat.at(q+p*m_nlc);
  }

  // reset the internally stored derivatives to zero in preperation for the next sample
  this->initialize_matrices();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Set the optimizable rotations to be those between the two specified orbital sets.
///
/// \param[in]      istart         1st index in the 1st orbital set
/// \param[in]      iend           one past the last index in the 1st orbital set
/// \param[in]      jstart         1st index in the 2nd orbital set
/// \param[in]      jend           one past the last index in the 2nd orbital set
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::set_optimizable_rotation_ranges(const int istart, const int iend, const int jstart, const int jend) {

  // check sanity
  assert( istart >= 0 );
  assert( iend   >= 0 );
  assert( jstart >= 0 );
  assert( jend   >= 0 );
  assert( istart <= iend );
  assert( jstart <= jend );

  // remove any existing rotations
  m_act_rot_inds.clear();

  // add all rotations between the orbital sets [istart, iend) and [jstart, jend)
  for (int i = istart; i < iend; i++)
    for (int j = jstart; j < jend; j++)
      if ( i != j )
        m_act_rot_inds.push_back(std::pair<int,int>(std::min(i,j), std::max(i,j)));

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Build this object's list of its own optimizable variables, and optionally set them
///         to input values which may be provided in input_params. If input_params is empty on
///         input, then params_supplied should be false, and each parameter will be set to 0.
///         Then, also apply the initial rotation using the provided input parameters.
///
/// \param[in]    input_params     the input list of parameters - can be empty, if no parameters
///                                were supplied by the user
/// \param[in]    params_supplied  true if parameters are provided in input_params, false if
///                                input_params is empty
/// \param[in]    print_vars       if true, then print out the initialized values of the variables
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void SlaterDetOpt::buildOptVariables(std::vector<RealType>& input_params,
                       bool params_supplied, bool print_vars) {

  int p, q;
  int nparams_active = m_act_rot_inds.size();

  if (params_supplied) {
    int nparams_input = input_params.size();
    if (nparams_input != nparams_active)
      throw std::runtime_error("Number of parameters provided for orbital rotations "
                               "is not consistent with the expected number.");
  }

  for (int i=0; i< nparams_active; i++)
  {
    p = m_act_rot_inds.at(i).first;
    q = m_act_rot_inds.at(i).second;
    std::stringstream sstr;
    sstr << Phi->objectName
         << "_orb_rot_"
         << ( p <   10 ? "0" : "" )
         << ( p <  100 ? "0" : "" )
         << ( p < 1000 ? "0" : "" )
         << p
         << "_"
         << ( q <   10 ? "0" : "" )
         << ( q <  100 ? "0" : "" )
         << ( q < 1000 ? "0" : "" )
         << q;

    // If the user input parameteres, use those. Otherwise, initialize the
    // parameter to zero.
    if (params_supplied) {
      myVars.insert(sstr.str(), input_params[i]);
    } else {
      myVars.insert(sstr.str(), 0.0);
    }
  }

  if (print_vars) {
    // Print the current values of all the optimisable parameters,
    // hopefully with correct formatting.
    app_log() << std::string(16,' ') << "Parameter name" << std::string(15,' ') << "Value\n";
    myVars.print(app_log());
  }

  // The code below applies the initial rotation requested by the user.
  // This is basically doing the same as resetParameters, but that routine
  // is a bit too specialized for what we want to do here. So we just
  // rewrite the specific code we want, rather than calling that routine.

  // Read out the parameters that define the rotation into an antisymmetric matrix
  std::vector<RealType> rot_mat(m_nlc*m_nlc, 0.0);
  for (int i = 0; i < m_act_rot_inds.size(); i++) {
    const int p = m_act_rot_inds.at(i).first;
    const int q = m_act_rot_inds.at(i).second;
    rot_mat[p+q*m_nlc] =  myVars[i];
    rot_mat[q+p*m_nlc] = -myVars[i];
  }
  // Exponentiate antisymmetric matrix to get the unitary rotation.
  this->exponentiate_matrix(m_nlc, &rot_mat.at(0));
  // get the linear combination coefficients by applying the rotation to the old coefficients
  Phi->rotate_B(rot_mat);
}

}

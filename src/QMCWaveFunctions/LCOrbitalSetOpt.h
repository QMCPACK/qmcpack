///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file src/QMCWaveFunctions/LCOrbitalSetOpt.h
///
/// \brief   A class for an optimizable set of linear combinations of single particle orbitals.
///
/// \author  Eric Neuscamman
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_OPTIMIZABLE_H
#define QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_OPTIMIZABLE_H

#include <stdexcept>
#include <utility>
#include <iterator>
#include <cassert>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

#include <QMCWaveFunctions/SPOSetBase.h>
#include <QMCWaveFunctions/OrbitalBase.h>
#include <Numerics/MatrixOperators.h>
#include <Utilities/RandomGenerator.h>

namespace qmcplusplus {

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  A class for packaging the data for a set of optimizable linear combinations of the
///         single particle orbitals in a way that can be optimized as part of the trial wave
///         function.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
class LCOrbitalSetOptTrialFunc : public OrbitalBase {

  // protected data members
  protected:

    /// \brief  number of linear combinations of basis functions (i.e. molecular orbitals)
    int m_nlc;

    /// \brief  number of basis functions
    int m_nb;

    /// \brief  position of the first of this object's optimizable variables in the overall list of optimizable variables
    int m_first_var_pos;

    /// \brief  vector of active rotation indices, stored in pairs with the first element of the pair less than the second
    std::vector<std::pair<int,int> > m_act_rot_inds;

    /// \brief  The column-major-order m_nb by m_nlc matrix of orbital coefficients resulting from a rotation of the old coefficients.
    ///         Thus B = old_B * C, where C is a unitary orbital rotation matrix.
    std::vector<RealType> m_B;

    /// \brief  the column-major-order m_nb by m_nlc matrix of old orbital coefficients
    std::vector<RealType> m_old_B;

    /// \brief  the column-major-order m_nb by m_nlc initial orbital coefficients, from the start of the simulation
    std::vector<RealType> m_init_B;

    /// \brief  matrix of derivatives of Log(Psi) w.r.t. the m_nlc by m_nlc orbital rotation matrix C
    std::vector<RealType> m_pder_mat;

    /// \brief  matrix of derivatives of (H Psi) / Psi w.r.t. the m_nlc by m_nlc orbital rotation matrix C
    std::vector<RealType> m_hder_mat;

    /// \brief  name of the LCOrbitalSetOpt object that uses this trial function
    std::string m_name;

  // protected member functions
  protected:

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  verifies that the number of linear combinations and the list of active rotation
    ///         indices is sane
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void check_index_sanity() const {

      // ensure the number of linear combinations is not negative
      if ( m_nlc < 0 )
        throw std::runtime_error("LCOrbitalSetOptTrialFunc::check_index_sanity found a negative number of linear combinations");

      // throw an error if any active rotation index is unreasonable
      for (std::vector<std::pair<int,int> >::const_iterator it = m_act_rot_inds.begin(); it != m_act_rot_inds.end(); it++) {
        if ( it->first >= it->second ) {
          std::stringstream error_msg;
          error_msg << "LCOrbitalSetOptTrialFunc::check_index_sanity found an active rotation index pair ("
                    << it->first << "," << it->second << ") in which the first index was not smaller than the second";
          throw std::runtime_error(error_msg.str());
        }
        if ( it->first < 0 || it->first >= m_nlc ) {
          std::stringstream error_msg;
          error_msg << it->first << " is an out of bounds first active rotation index in LCOrbitalSetOptTrialFunc::check_index_sanity";
          throw std::runtime_error(error_msg.str());
        }
        if ( it->second < 0 || it->second >= m_nlc ) {
          std::stringstream error_msg;
          error_msg << it->second << " is an out of bounds second active rotation index in LCOrbitalSetOptTrialFunc::check_index_sanity";
          throw std::runtime_error(error_msg.str());
        }
      }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  prepares the matrices that will hold derivatives w.r.t. orbital rotations by
    ///         ensuring they are the right size and that their elements are all zero
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void initialize_matrices() {

      // initialize matrix for Log(Psi) derivatives
      m_pder_mat.resize(m_nlc * m_nlc);
      std::fill(m_pder_mat.begin(), m_pder_mat.end(), 0.0);

      // initialize matrix for ( H Psi ) / Psi derivatives
      m_hder_mat.resize(m_nlc * m_nlc);
      std::fill(m_hder_mat.begin(), m_hder_mat.end(), 0.0);

    }

  // public member functions
  public:

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  return the name of this class
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    static std::string name() { return "LCOrbitalSetOptTrialFunc"; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  return a pointer to the orbital coefficient matrix data
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    const RealType * ptr_to_B() const { return &(*m_B.begin()); }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  constructor
    ///
    /// \param[in]      nlc            the number of linear combinations (i.e. molecular orbitals)
    /// \param[in]      nb             the number of basis functions
    /// \param[in]      initial_B      nb by nlc column-major-ordered matrix of initial orbital coefficients
    /// \param[in]      name           name of the LCOrbitalSetOpt object that owns this trial function component
    /// \param[in]      mix_factor     factor controlling mixing of the initial orbitals
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    LCOrbitalSetOptTrialFunc(const int nlc, const int nb, const RealType * const initial_B, const std::string & name, double mix_factor)
      : m_nlc(nlc),
        m_nb(nb),
        m_first_var_pos(-1),
        m_act_rot_inds(),
        m_B(initial_B, initial_B + nlc*nb),
        m_old_B(initial_B, initial_B + nlc*nb),
        m_init_B(initial_B, initial_B + nlc*nb),
        m_name(name)
    {

      // set the name
      this->OrbitalName = LCOrbitalSetOptTrialFunc::name();

      // make sure we didn't start with a bad m_nlc
      this->check_index_sanity();

      // by default set all rotations to be active
      m_act_rot_inds.resize( m_nlc * ( m_nlc - 1 ) / 2 );
      int rots_recorded = 0;
      for (int j = 1; j < m_nlc; j++)
      for (int i = 0; i < j; i++)
        m_act_rot_inds.at(rots_recorded++) = std::pair<int,int>(i,j);
      if ( m_act_rot_inds.size() != rots_recorded )
        throw std::runtime_error("wrong number of active rotations recorded in LCOrbitalSetOptTrialFunc::LCOrbitalSetOptTrialFunc");

      // make sure we didn't do something stupid
      this->check_index_sanity();

      // prepare matrices that will hold derivatives wrt orbital rotations
      this->initialize_matrices();

      // if requested, mix the initial basis orbitals together
      if ( mix_factor != 0.0 ) {

        // mix
        for (int i = m_nb - 1; i >= 0; i--) {
          for (int j = 0; j < m_nlc; j++) {
            //if ( mix_factor > 0.5 )
            //  throw std::runtime_error("mix_factor grew too large.  Please choose a smaller value of orbital_mix_magnitude");
            m_B.at(i+j*m_nb) += mix_factor * 2.0 * ( Random() - 0.5 );
            //mix_factor *= 1.002;
            //m_B.at(i+k*m_nb) =  std::sqrt(1.0 - mix_factor) * m_old_B.at(i+k*m_nb) + std::sqrt(      mix_factor) * m_old_B.at(i+j*m_nb);
            //m_B.at(i+j*m_nb) = -std::sqrt(      mix_factor) * m_old_B.at(i+k*m_nb) + std::sqrt(1.0 - mix_factor) * m_old_B.at(i+j*m_nb);
          }
        }

        // re-orthonormalize
        for (int j = 0; j < m_nlc; j++) {
          const RealType norm = std::abs(std::sqrt(BLAS::dot(m_nb, &m_B.at(0+j*m_nb), &m_B.at(0+j*m_nb))));
          BLAS::scal(m_nb, 1.0 / norm, &m_B.at(0+j*m_nb));
          for (int k = j+1; k < m_nlc; k++) {
            const RealType x = BLAS::dot(m_nb, &m_B.at(0+j*m_nb), &m_B.at(0+k*m_nb));
            BLAS::axpy(m_nb, -x, &m_B.at(0+j*m_nb), 1, &m_B.at(0+k*m_nb), 1);
          }
        }

        // save the mixed orbitals
        m_old_B = m_B;
        m_init_B = m_B;

      }

      // print the orbitals
      this->print_B();

      //// print all orbitals for use in test input
      //app_log() << std::endl;
      //app_log() << "printing transposed molecular orbital coefficients for test" << std::endl;
      //for (int j = 0; j < m_nlc; j++) {
      //  for (int i = 0; i < m_nb; i++)
      //    app_log() << " " << std::right << std::scientific << std::setprecision(15) << std::setw(25) << m_B.at(i+j*m_nb);
      //  app_log() << std::endl;
      //}
      //app_log() << std::endl;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  print the molecular orbital coefficients, one MO per column
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void print_B() {
      app_log() << std::endl;
      app_log() << "printing molecular orbital coefficients" << std::endl;
      for (int i = 0; i < m_nb; i++) {
        for (int j = 0; j < m_nlc; j++)
          app_log() << " " << std::right << std::fixed << std::setprecision(16) << std::setw(22) << m_B.at(i+j*m_nb);
        app_log() << std::endl;
      }
      app_log() << std::endl;
    }

    //DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;

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
    void add_derivatives(const int nl,
                         const int np,
                         const RealType * const dp0,
                         const RealType * const dh0,
                         //const RealType * const dh1,
                         const RealType * const dh2,
                         const RealType * const Bchi,
                         //const RealType * const dBchi,
                         const RealType * const d2Bchi) {

      // ensure the number of linear combinations is correct
      if ( nl != m_nlc ) {
        std::stringstream error_msg;
        error_msg << "supplied number of linear combinations (" << nl << ") does not match that held internally (" << m_nlc << ") in LCOrbitalSetOptTrialFunc::add_derivatives";
        throw std::runtime_error(error_msg.str());
      }

      // ensure orbital derivative matrices are the correct size
      if ( m_pder_mat.size() != nl * nl ) {
        std::stringstream error_msg;
        error_msg << "nl (" << nl << ") does not match size of m_pder_mat (" << m_pder_mat.size() << ") in LCOrbitalSetOptTrialFunc::add_derivatives";
        throw std::runtime_error(error_msg.str());
      }
      if ( m_hder_mat.size() != nl * nl ) {
        std::stringstream error_msg;
        error_msg << "nl (" << nl << ") does not match size of m_hder_mat (" << m_hder_mat.size() << ") in LCOrbitalSetOptTrialFunc::add_derivatives";
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
      //for (int i = 0; i < 3; i++)
      //  BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh1+i*nl*np, nl, dBchi+i*nl*np, nl, RealType(1.0), &m_hder_mat.at(0), nl);
      BLAS::gemm('N', 'T', nl, nl, np, RealType(1.0), dh2, nl, d2Bchi, nl, RealType(1.0), &m_hder_mat.at(0), nl);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Exponentiates a matrix
    ///
    /// \param[in]      n              matrix dimensions
    /// \param[in,out]  mat            On entry, the n by n matrix.
    ///                                On exit, the exponential of the matrix.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void exponentiate_matrix(const int n, RealType * const mat) {

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

//    ///////////////////////////////////////////////////////////////////////////////////////////////////
//    /// \brief  don't know what this is supposed to do yet
//    ///
//    /// \param[in]      optvars        ???
//    ///
//    ///////////////////////////////////////////////////////////////////////////////////////////////////
//    void resetParameters(const opt_variables_type& optvars) {
//      //throw std::runtime_error("LCOrbitalSetOptTrialFunc::resetParameters not yet implemented");
//      app_log() << "WARNING: LCOrbitalSetOptTrialFunc::resetParameters is not doing anything" << std::endl;
//    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Reset the new orbital coefficients by reading out the orbital rotation from the
    ///         supplied variable set and then applying the orbital rotation to the old coefficients.
    ///
    ///         This is a specialization of the SPOSetBase class virtual funciton.
    ///
    /// \param[in]      active         variable set to read orbital rotations from
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void resetParameters(const opt_variables_type & active) {

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
      //BLAS::gemm('N', 'T', m_nb, m_nlc, m_nlc, RealType(1.0), &m_old_B.at(0), m_nb, &rot_mat.at(0), m_nlc, RealType(0.0), &m_B.at(0), m_nb);
      BLAS::gemm('N', 'T', m_nb, m_nlc, m_nlc, RealType(1.0), &m_init_B.at(0), m_nb, &rot_mat.at(0), m_nlc, RealType(0.0), &m_B.at(0), m_nb);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Copy this object's parameters from the supplied variable set, convert this object's
    ///         parameters to a standard form, and optionally copy the standard form parameters back
    ///         into the supplied variable set.
    ///
    /// \param[in,out]  active         the supplied variable set
    /// \param[in]      copy_back      whether to copy parameters back to the variable set at the end
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void putParametersInStandardForm(opt_variables_type & active, const bool copy_back) {

      // compute the updated orbital coefficients
      this->resetParameters(active);

      // record the updated orbital coefficients as the base coefficients from which we will make new rotations
      m_old_B = m_B;

      // Update the internally stored optimizable variables list to equal the values passed in. Then,
      // the next time this function is called, we'll be able to calculate by how much the parameters
      // changed, and therefore what orbital rotation to apply to update the orbitals appropriately.
      if ( copy_back )
        for (int i = 0; i < m_act_rot_inds.size(); i++)
          myVars[i] = active[i + m_first_var_pos];

      // print the new orbitals
      if ( copy_back )
        this->print_B();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  don't know what this is supposed to do yet
    ///
    /// \param[in,out]  os             ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void reportStatus(std::ostream& os) {
      throw std::runtime_error("LCOrbitalSetOptTrialFunc::reportStatus not yet implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  don't know what this is supposed to do yet
    ///
    /// \param[in,out]  P              ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void resetTargetParticleSet(ParticleSet& P) {
      throw std::runtime_error("LCOrbitalSetOptTrialFunc::resetTargetParticleSet not yet implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  don't know what this is supposed to do yet
    ///
    /// \param[in,out]  P              ???
    /// \param[in,out]  G              ???
    /// \param[in,out]  L              ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ValueType evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {
      throw std::runtime_error("LCOrbitalSetOptTrialFunc::evaluate not yet implemented");
      return 0.0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Returns ln(Phi), where Phi is the contribution from the LCOrbitalSetOptTrialFunc.  As this is a dummy
    //          function call, return ln(Phi)=0 and do nothing to gradients or laplacians.
    ///
    /// \param[in,out]  P              Particle set
    /// \param[in,out]  G              Gradient to be modified
    /// \param[in,out]  L              Laplacian to be modified
    //
    //  \return  0, as this doesn't contribute to the wavefunction
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    RealType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) {
      //throw std::runtime_error("LCOrbitalSetOptTrialFunc::evaluateLog not yet implemented");
      return 0.0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Returns ln(Phi), where Phi is the contribution from the LCOrbitalSetOptTrialFunc.  As this is a dummy
    //          function call, return ln(Phi)=0 and do nothing to gradients or laplacians.
    ///
    /// \param[in,out]  P              Particle set
    /// \param[in,out]  buf            walker buffer for storing orbital specific variables.  Do nothing to it.
    //
    //  \return  0, as this doesn't contribute to the wavefunction
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    RealType evaluateLog(ParticleSet& P, BufferType& buf) {
        throw std::runtime_error("LCOrbitalSetOpt::evaluateLog(P, buff)");
	return 0.0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  This class makes no direct contribution to the trial function value, gradient,
    ///         or laplacian, and so this function does nothing.
    ///
    /// \param[in]      P              the particle set
    /// \param[in,out]  G              gradient to add to
    /// \param[in,out]  L              laplacian to add to
    /// \param[in,out]  buf            buffer to load or store temporary data from
    /// \param[in]      fillBuffer     whether to fill data into the buffer or read data from it
    ///
    /// \return  zero (this class's contribution to the log of the trial function value)
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    RealType evaluateLog(ParticleSet& P,
                         ParticleSet::ParticleGradient_t& G,
                         ParticleSet::ParticleLaplacian_t& L,
                         PooledData<RealType>& buf,
                         bool fillBuffer ) {

//     throw std::runtime_error("LCOrbitalSetOptTrialFunc::evaluateLog(bunch stuff) not yet implemented"); 
//        throw std::runtime_error("LCOrbitalSetOpt::evaluateLog(P, G, L, buff, bool)");
     return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Initialize the object for the given particle positions and add any essential internal
    ///         data to the supplied buffer.
    ///
    /// \param[in]      P              the particle set
    /// \param[in]      buf            the buffer to add essential data to
    ///
    /// \return  zero, as this object does not contribute directly to the wave function value
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    RealType registerData(ParticleSet& P, BufferType& buf) { return 0; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Initialize the object for the given particle positions and put its essential internal
    ///         data into the supplied buffer
    ///
    /// \param[in]      P              the particle set
    /// \param[in]      buf            the buffer to save essential data in
    /// \param[in]      fromscratch    ???
    ///
    /// \return  zero, as this object does not contribute directly to the wave function value
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) { return 0; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Read in necessary internal data from the provided buffer
    ///
    /// \param[in]      P              the particle set
    /// \param[in]      buf            the buffer to read from
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void copyFromBuffer(ParticleSet& P, BufferType& buf) {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Dummy call from TrialWaveFunction to evaluate contribution to wavefunction ratio from LCOrbitalSetTrialFunc
    //		Does nothing by returning ratio=1.0 (rationew=ratioold*1.0)t
    ///
    /// \param[in,out]  P              the particle set
    /// \param[in]      iat            the atom id of particle being moved
    /// \param[in,out]  dG             the gradient to be added to
    /// \param[in,out]  dL             the laplacian to be added to
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ValueType ratio(ParticleSet& P, int iat, ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL) {
      return 1.0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Dummy call from TrialWaveFunction to evaluate contribution to wavefunction ratio from LCOrbitalSetTrialFunc
    //		Does nothing by returning ratio=1.0 (rationew=ratioold*1.0)t
    ///
    /// \param[in,out]  P              the particle set
    /// \param[in]      iat            the atom id of particle being moved
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ValueType ratio(ParticleSet& P, int iat) {
      return 1.0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  don't know what this is supposed to do yet
    ///
    /// \param[in,out]  P              ???
    /// \param[in,out]  dG             ???
    /// \param[in,out]  dL             ???
    /// \param[in]      iat            ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void update(ParticleSet& P, ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL, int iat) {
      throw std::runtime_error("LCOrbitalSetOptTrialFunc::update not yet implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  The orbital set does not need to do anything for a single particle move.
    ///
    /// \param[in]      P              the particle set, which I think carries information about the particle's move
    /// \param[in]      iat            the id number of the moved particle
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void acceptMove(ParticleSet& P, int iat) {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Nothing to do for a rejected single particle move
    ///
    /// \param[in]      iat            the id number of the moved particle
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void restore(int iat) {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  return a clone of this object
    ///
    /// \param[in,out]  tqp            ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    OrbitalBasePtr makeClone(ParticleSet& tqp) const {
      return OrbitalBasePtr( new LCOrbitalSetOptTrialFunc(*this) );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Add this trial function component's logpsi and hpsioverpsi derivatives w.r.t.
    ///         to optimizable parameters to the cumulative total of these derivatives and
    ///         then reset the internally held cumulative derivatives to zero.
    ///
    /// \param[in]      P             Object containing information on particle positions.
    /// \param[in]      active        not used here
    /// \param[in,out]  dlogpsi       vector of derivatives of  Log(Psi)  w.r.t optimizable parameters
    /// \param[in,out]  dhpsioverpsi  vector of derivatives of  ( H Psi ) / Psi  w.r.t optimizable parameters
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluateDerivatives(ParticleSet& P,
                             const opt_variables_type& active,
                             std::vector<RealType>& dlogpsi,
                             std::vector<RealType>& dhpsioverpsi) {

      // check that we have the position of the first of our variables in the overall list
      if ( myVars.size() > 0 && m_first_var_pos < 0 )
        throw std::runtime_error("position of first variable was not set on entry to LCOrbitalSetOptTrialFunc::evaluateDerivatives");

      // check that my number of variables is consistent with the number of active rotations
      if ( myVars.size() != m_act_rot_inds.size() ) {
        std::stringstream error_msg;
        error_msg << "mismatch between myVars.size() (" << myVars.size() << ") and m_act_rot_inds.size() (" << m_act_rot_inds.size() << ") in LCOrbitalSetOptTrialFunc::evaluateDerivatives";
        throw std::runtime_error(error_msg.str());
      }

      // add derivatives to totals
      ////app_log() << std::endl;
      ////for (int i = 0; i < m_act_rot_inds.size(); i++) {
      ////  const int p = m_act_rot_inds.at(i).first;
      ////  const int q = m_act_rot_inds.at(i).second;
      ////  if ( true ) {
      ////    std::vector<char> buff(1000, ' ');
      ////    const int len = std::sprintf(&buff[0], " p = %4i   q = %4i     dlogpsi = %20.12f     dhpsioverpsi = %20.12f", p, q, dlogpsi.at(m_first_var_pos+i), dhpsioverpsi.at(m_first_var_pos+i));
      ////    for (int k = 0; k < len; k++)
      ////      app_log() << buff[k];
      ////    app_log() << std::endl;
      ////  }
      ////}
      ////app_log() << std::endl;
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
      ////app_log() << std::endl;
      ////throw std::runtime_error("STOPPING HERE");

      // reset the internally stored derivatives to zero in preperation for the next sample
      this->initialize_matrices();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Build this object's list of its own optimizable variables, and optionally set them
    ///         to input values which may be provided in input_params. If input_params is empty on
    ///         input, then params_supplied should be false, and each parameter will be set to 0.
    ///         Then, also apply the initial rotation using the provided input parameters.
    ///
    /// \param[in]    input_params     the input list of parameters - can be empty, if no parameters
    ///                                were supplied by the user
    /// \param[in]    spo_name         name of the single particle basis set object
    /// \param[in]    params_supplied  true if parameters are provided in input_params, false if
    ///                                input_params is empty
    /// \param[in]    print_vars       if true, then print out the initialized values of the variables
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void buildOptVariables(std::vector<RealType>& input_params, const std::string & spo_name,
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
        sstr << spo_name
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
      // Get the linear combination coefficients by applying the rotation to
      // the old coefficients
      //BLAS::gemm('N', 'T', m_nb, m_nlc, m_nlc, RealType(1.0), &m_old_B.at(0), m_nb, &rot_mat.at(0), m_nlc, RealType(0.0), &m_B.at(0), m_nb);
      BLAS::gemm('N', 'T', m_nb, m_nlc, m_nlc, RealType(1.0), &m_init_B.at(0), m_nb, &rot_mat.at(0), m_nlc, RealType(0.0), &m_B.at(0), m_nb);
      // Record the updated orbital coefficients as the base coefficients from
      // which we will make new rotations
      m_old_B = m_B;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Check this object's optimizable variables into the supplied overall list of
    ///         optimizable variables.
    ///
    /// \param[in,out]  active        the overall list of optimizable variables
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void checkInVariables(opt_variables_type& active) {

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
    void checkOutVariables(const opt_variables_type& active) {

      // record the positions of this object's optimizable variables within the overall list of optimizable variables
      myVars.getIndex(active);

      // ensure that this object's variables are stored contiguously
      for (int i = 0; i < myVars.size(); i++) {
        if ( myVars.where(i) - myVars.where(0) != i ) {
          std::stringstream error_msg;
          error_msg << "variable " << (i-1) << " was not contiguous with variable " << i << " in LCOrbitalSetOptTrialFunc::checkOutVariables";
          throw std::runtime_error(error_msg.str());
        }
      }

      // record the position of my first variable
      if ( myVars.size() > 0 )
        m_first_var_pos = myVars.where(0);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  This object makes no contribution to the gradient directly.
    ///
    /// \param[in]      P              the particle set, which I think carries information about the particle's move
    /// \param[in]      iat            the id number of the moved particle
    ///
    /// \return  a GradType vector full of zeros
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    OrbitalBase::GradType evalGrad(ParticleSet& P, int iat) {

      // no gradient contribution directly, instead it should be handled by wave function components relying on the LCOrbitalSetOpt object
      OrbitalBase::GradType g;
      g = 0.0;
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
    OrbitalBase::ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) {
  
      // no gradient contribution directly, instead it should be handled by wave function components relying on the LCOrbitalSetOpt object
      return 1.0;
  
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
    void set_optimizable_rotation_ranges(const int istart, const int iend, const int jstart, const int jend) {

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

};

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  A class for a set of optimizable linear combinations of the single particle orbitals
///         provided by either another SPO set or by a basis set of the templated type BS.
///         We refer to the linear combinations of basis orbitals as molecular orbitals.
///         The set of molecular orbitals is optimized by rotating them amongst each other.
///
///         The molecular orbital coefficients themselves are stored in a subobject that is
///         used as part of the trial funciton, so that the necessary derivatives w.r.t.
///         rotations between the molecular orbitals can be collected during sampling.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class BS> class LCOrbitalSetOpt : public SPOSetBase {

  // protected data members
  protected:

    /// \brief  pointer that, if not null, will be used to evaluate the basis orbitals
    SPOSetBase * m_spo_set;

    /// \brief  pointer to the basis set that evaluates the basis orbitals if m_spo_set == 0
    BS * m_basis_set;

    /// \brief  pointer to a subobject that will be part of the trial function
    LCOrbitalSetOptTrialFunc * m_tf;

    /// \brief  the level of printing
    int m_report_level;

    /// \brief  workspace matrix 
    std::vector<ValueType> m_lc_coeffs;

    /// \brief  workspace matrix
    std::vector<ValueType> m_basis_vals;

    /// \brief  workspace matrix
    std::vector<ValueType> m_basis_der1;

    /// \brief  workspace matrix
    std::vector<ValueType> m_basis_der2;

    /// \brief  vector to put temporary orbital data in
    ValueVector_t m_temp_p;

    /// \brief  vector to put temporary gradient data in
    GradVector_t m_temp_g;

    /// \brief  vector to put temporary laplacian data in
    ValueVector_t m_temp_l;

    /// \brief  factor controlling how much to mix the initial orbitals
    double m_omixfac;

  // public data members
  public:

  // protected member functions
  protected:

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Ensures the supplied dimension is reasonable and throws an expection if not
    ///
    /// \param[in]      name           a name for the thing whose dimension in being checked
    /// \param[in]      caller         a name for the calling function
    /// \param[in]      n              the dimension
    /// \param[in]      s              the maximum allowed length for the dimension
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void check_input_dim(const std::string & name, const std::string & caller, const int n, const int s) const {

      // ensure dimension is nonzero
      if ( n <= 0 )
        throw std::runtime_error(name + " has a length less than one in " + caller);

      // ensure vector is not too long
      if ( n > s )
        throw std::runtime_error(name + " is too long in " + caller);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  sets the factor by which we will mix the initial orbitals
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void set_orbital_mixing_factor(const double factor) {
      m_omixfac = factor;
    }

  // public member functions
  public:

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  constructor from basis set and reporting level
    ///
    /// \param[in]      bs             pointer to the basis set to use
    /// \param[in]      rl             reporting level to use
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    LCOrbitalSetOpt(BS * const bs = 0, const int rl = 0) : m_spo_set(0), m_basis_set(0), m_tf(0), m_report_level(rl), m_omixfac(0) {

      // set the basis set
      if ( bs ) this->setBasisSet(bs);

      // initialize number of molecular orbitals as zero
      this->OrbitalSetSize = 0;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  constructor from SPO set and reporting level
    ///
    /// \param[in]      spo            pointer to the spo set to use
    /// \param[in]      rl             reporting level to use
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    LCOrbitalSetOpt(SPOSetBase * const spo, const int rl = 0) : m_spo_set(0), m_basis_set(0), m_tf(0), m_report_level(rl), m_omixfac(0) {

      // set the internal SPO set
      if ( spo ) this->setSPOSet(spo);

      // initialize number of molecular orbitals as zero
      this->OrbitalSetSize = 0;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  the destructor, which assumes deallocation of basis set is done elsewhere
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ~LCOrbitalSetOpt() {
      // if the trial function subobject exists, destroy it
      if ( m_tf ) delete m_tf;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  returns that this is indeed an LCOrbitalSetOpt object
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    bool is_of_type_LCOrbitalSetOpt() const { return true; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  creates this object's trial function component
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void prepare_tf_component() {
      if ( m_tf )
        throw std::runtime_error("an LCOrbitalSetOpt object tried to create its trial function component twice");
      app_log() << "creating LCOrbitalSetOpt trial function component with C.rows() = " << C.rows() << " and C.cols() = " << C.cols() << std::endl;
      m_tf = new LCOrbitalSetOptTrialFunc(this->OrbitalSetSize, this->BasisSetSize, C.data(), this->objectName, m_omixfac);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  returns this object's trial function component, creating it if necessary
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    OrbitalBase * tf_component() {
      if ( !m_tf )
        this->prepare_tf_component();
      return m_tf;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  set the basis set the object should use
    ///
    /// \param[in]      bs             pointer to the basis set to use
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void setBasisSet(BS * const bs) {

      // error if the pointer is empty
      if ( !bs )
        throw std::runtime_error("basis set pointer was empty in LCOrbitalSetOpt::setBasisSet");

      // remember the basis set
      m_basis_set = bs;

      // extract the number of single particle orbitals in the basis set
      this->BasisSetSize = m_basis_set->getBasisSetSize();

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  set the internal SPO set object to use
    ///
    /// \param[in]      spo            pointer to the SPO set to use
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void setSPOSet(SPOSetBase * const spo) {

      // error if the pointer is empty
      if ( !spo )
        throw std::runtime_error("spo set pointer was empty in LCOrbitalSetOpt::setSPOSet");

      // remember the basis set
      m_spo_set = spo;

      // extract the number of single particle orbitals in the basis set
      this->BasisSetSize = m_spo_set->OrbitalSetSize;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  clone the object (specialization of the base class virtual funciton)
    ///
    ///         This is a specialization of the SPOSetBase class virtual funciton.
    ///
    /// \return  a base class pointer to the clone
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    SPOSetBase * makeClone() const {

      // create a clone that contains a cloned spo set or basis set
      SPOSetBase * retval;
      if ( m_spo_set )
        retval = new LCOrbitalSetOpt(m_spo_set->makeClone(), m_report_level);
      else
        retval = new LCOrbitalSetOpt(m_basis_set->makeClone(), m_report_level);

      // set the number of molecular orbitals
      retval->setOrbitalSetSize(this->OrbitalSetSize);

      // return the clone
      return retval;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  reset the basis's target particleset
    ///
    ///         This is a specialization of the SPOSetBase class virtual funciton.
    ///
    /// \param[in,out]  P              ???
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void resetTargetParticleSet(ParticleSet & P) {

      if ( m_spo_set )
        m_spo_set->resetTargetParticleSet(P);
      else
        m_basis_set->resetTargetParticleSet(P);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Nothing to do for this as the variable rotations are handled by the trial function
    ///         component subobject.
    ///
    /// \param[in]      optvars        not used here
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void resetParameters(const opt_variables_type& optvars) {
      //app_log() << "WARNING: LCOrbitalSetOpt::resetParameters is not doing anything" << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  set the number of linear combinations of basis functions (i.e. molecular orbitals)
    ///
    ///         This is a specialization of the SPOSetBase class virtual funciton.
    ///
    /// \param[in]      norbs          how many linear combinations are desired
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void setOrbitalSetSize(int norbs) {

      // record the number of linear combinations (i.e. molecular orbitals)
      this->OrbitalSetSize = norbs;
      app_log() << "LCOrbitalSetOpt finished setOrbitalSetSize with norbs = " << norbs << std::endl;

    }

//    !!! this function does not appear to be callable via the base class pointer as it is not virtual in SPOSetBase
//    inline int getBasisSetSize() const
//    {
//      return (m_basis_set==0)? 0: m_basis_set->getBasisSetSize();
//    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Evaluates the values, x,y,z derivatives, and x-y-z-summed second derivatives of the
    ///         specified linear combinations of basis set orbitals at the specified particles'
    ///         positions.
    ///
    /// \param[in]      P            Object containing information on particle positions.
    /// \param[in]      mt           the move type: 'p' for particle move, 'w' for walker move
    /// \param[in]      ostart       Iterator for the start of the index range specifying which linear combinations of orbitals to evaluate.
    /// \param[in]      oend         Iterator for the end   of the index range specifying which linear combinations of orbitals to evaluate.
    /// \param[in]      pstart       Iterator for the start of the index range specifying which particles' positions to use.
    /// \param[in]      pend         Iterator for the end   of the index range specifying which particles' positions to use.
    /// \param[in,out]  vmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix
    ///                              of the values of the specified linear combinations for the specified particles' positions.
    /// \param[in,out]  gmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
    ///                              each element of which is a length 3 vector containing the x,y,z gradients of the values in vmat.
    /// \param[in,out]  lmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
    ///                              each element of which is the sum of x^2, y^2, and z^2 second derivatives of the values in vmat.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluate_notranspose_general(const ParticleSet& P,
                                      const char mt,
                                      std::vector<int>::const_iterator ostart,
                                      std::vector<int>::const_iterator oend,
                                      std::vector<int>::const_iterator pstart,
                                      std::vector<int>::const_iterator pend,
                                      ValueType * const vmat,
                                      GradType * const gmat,
                                      ValueType * const lmat) {

      // get the number of linear combinations
      const int no = std::distance(ostart, oend);

      // get the number of particles
      const int np = std::distance(pstart, pend);

      // get the number of basis orbitals
      const int nb = BasisSetSize;

      // check sanity of matrix dimensions
      assert( no > 0 );
      assert( np > 0 );
      assert( nb > 0 );
      assert( nb >= no );
      assert( nb >= np );

      // resize the temporary arrays if they are not big enough
      SPOSetBase::ensure_vector_is_big_enough(m_lc_coeffs, no * nb);
      SPOSetBase::ensure_vector_is_big_enough(m_basis_vals, np * nb);
      SPOSetBase::ensure_vector_is_big_enough(m_basis_der1, 3 * np * nb);
      SPOSetBase::ensure_vector_is_big_enough(m_basis_der2, np * nb);
      if ( m_temp_p.size() != nb ) m_temp_p.resize(nb);
      if ( m_temp_g.size() != nb ) m_temp_g.resize(nb);
      if ( m_temp_l.size() != nb ) m_temp_l.resize(nb);

      // get convenient name for iterator type
      typedef std::vector<int>::const_iterator Iter;

      // choose whether to use careful loops or BLAS copies for moving gradient data
      const bool careful_loops_for_grad = true;

      // Evaluate and store the basis values, derivatives, and second derivatives for each particle position.
      // We store these data in five column-major-ordered (# of basis states) by (# of particles) matrices,
      // ( 1 matrix in m_basis_vals, 1 matrix in m_basis_der2, and 3 matrices in m_basis_der1 )
      {
        int i = 0;
        for (Iter it = pstart; it != pend; it++, i++) {

          // evaluate basis set data using the internal spo set if we have one
          if ( m_spo_set ) m_spo_set->evaluate(P, *it, m_temp_p, m_temp_g, m_temp_l);

          // evaluate basis set data for a particle move
          else if ( mt == 'p' ) m_basis_set->evaluateAllForPtclMove(P, *it);

          // evaluate basis set data for a walker move
          else if ( mt == 'w' ) m_basis_set->evaluateForWalkerMove(P, *it);

          // error for no internal spo set and an unknown move type
          else throw std::runtime_error("unknown move type in LCOrbitalSetOpt::evaluate_notranspose_general");

          // sanity checks
          if ( m_basis_set ) {
            assert( m_basis_set->Phi.size() == nb );
            assert( m_basis_set->d2Phi.size() == nb );
            assert( m_basis_set->dPhi.size() == nb );
            assert( m_basis_set->dPhi[0].size() == 3 );
          }

          // get references to the basis set data
          ValueVector_t & data_p = ( m_spo_set ? m_temp_p : m_basis_set->Phi );
           GradVector_t & data_g = ( m_spo_set ? m_temp_g : m_basis_set->dPhi );
          ValueVector_t & data_l = ( m_spo_set ? m_temp_l : m_basis_set->d2Phi );

          // copy values into a column of the basis value matrix
          BLAS::copy(nb, &data_p[0], 1, &m_basis_vals[i*nb], 1);

          // copy summed 2nd derivatives into a column of the basis 2nd derivative matrix
          BLAS::copy(nb, &data_l[0], 1, &m_basis_der2[i*nb], 1);

          // copy 1st derivatives into columns of the three different basis 1st derivative matrices
          if ( careful_loops_for_grad ) {
            for (int p = 0; p < 3; p++)
            for (int j = 0; j < nb; j++)
              m_basis_der1[ j + i*nb + p*np*nb ] = data_g[j][p];
          } else {
            for (int p = 0; p < 3; p++)
              BLAS::copy(nb, &data_g[0][p], 3, &m_basis_der1[ i*nb + p*np*nb ], 1);
          }

        }
      }

      // Store the slice of the linear combination coefficient matrix that we need in a column-major-ordered
      // (# of linear combinations) by (# of basis states) matrix.
      {
        int i = 0;
        for (Iter it = ostart; it != oend; it++, i++)
          BLAS::copy(nb, m_tf->ptr_to_B() + (*it)*nb, 1, &m_lc_coeffs[i], no);
      }

      // print what is in C
      if ( false ) {
        app_log() << "printing C" << std::endl;
        std::vector<char> buff(1000, ' ');
        for (int i = 0; i < BasisSetSize; i++) {
          for (int j = 0; j < OrbitalSetSize; j++) {
            const int len = std::sprintf(&buff[0], "  %12.6f", m_tf->ptr_to_B()[i+j*BasisSetSize]);
            for (int k = 0; k < len; k++)
              app_log() << buff[k];
          }
          app_log() << std::endl;
        }
        app_log() << std::endl;
        throw std::runtime_error("done printing C");
      }

      // compute the matrix of linear combination values for each particle
      BLAS::gemm('N', 'N', no, np, nb, ValueType(1.0), &m_lc_coeffs[0], no, &m_basis_vals[0], nb, ValueType(0.0), vmat, no);

      // compute the matrix of summed 2nd derivatives of linear combinations for each particle
      BLAS::gemm('N', 'N', no, np, nb, ValueType(1.0), &m_lc_coeffs[0], no, &m_basis_der2[0], nb, ValueType(0.0), lmat, no);

      // compute the matrix of 1st derivatives of linear combinations for each particle (using m_basis_vals as temporary storage)
      for (int p = 0; p < 3; p++) {
        BLAS::gemm('N', 'N', no, np, nb, ValueType(1.0), &m_lc_coeffs[0], no, &m_basis_der1[p*np*nb], nb, ValueType(0.0), &m_basis_vals[0], no);
        if ( careful_loops_for_grad ) {
          for (int j = 0; j < np; j++)
          for (int i = 0; i < no; i++)
            gmat[i+j*no][p] = m_basis_vals[i+j*no];
        } else {
          BLAS::copy(no*np, &m_basis_vals[0], 1, &gmat[0][p], 3);
        }
      }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Evaluates the values, x,y,z derivatives, and x-y-z-summed second derivatives of the
    ///         linear combinations in the range [os,oe) at the particle positions in the range [ps,pe).
    ///
    /// \param[in]      P            Object containing information on particle positions.
    /// \param[in]      mt           the move type: 'p' for particle move, 'w' for walker move
    /// \param[in]      os           Beginning of the range specifying which linear combinations of orbitals to evaluate.
    /// \param[in]      oe           End of the range specifying which linear combinations of orbitals to evaluate.
    /// \param[in]      ps           Beginning of the range specifying which particles' positions to use.
    /// \param[in]      pe           End of the range specifying which particles' positions to use.
    /// \param[in,out]  vmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix
    ///                              of the values of the specified linear combinations for the specified particles' positions.
    /// \param[in,out]  gmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
    ///                              each element of which is a length 3 vector containing the x,y,z gradients of the values in vmat.
    /// \param[in,out]  lmat         On input, points to an array of length (# of linear combinations) * (# of particle).
    ///                              On exit, holds a column-major-ordered (# of linear combinations) by (# of particle) matrix,
    ///                              each element of which is the sum of x^2, y^2, and z^2 second derivatives of the values in vmat.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluate_notranspose_ranges(const ParticleSet& P,
                                     const char mt,
                                     const int os,
                                     const int oe,
                                     const int ps,
                                     const int pe,
                                     ValueType * const vmat,
                                     GradType * const gmat,
                                     ValueType * const lmat) {

      // check sanity
      if ( oe < os )
        throw std::runtime_error("orbitital end (oe) is less than start (os) in LCOrbitalSetOpt::evaluate_notranspose_ranges");
      if ( pe < ps )
        throw std::runtime_error("particle end (pe) is less than start (ps) in LCOrbitalSetOpt::evaluate_notranspose_ranges");

      // prepare orbital list
      std::vector<int>::const_iterator oend = SPOSetBase::prepare_index_vector_contiguous(os, oe, m_oidx);

      // prepare particle list
      std::vector<int>::const_iterator pend = SPOSetBase::prepare_index_vector_contiguous(ps, pe, m_pidx);

      // evaluate
      this->evaluate_notranspose_general(P, mt, m_oidx.begin(), oend, m_pidx.begin(), pend, vmat, gmat, lmat);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Evaluates the values, x,y,z derivatives, and x-y-z-summed second derivatives of the
    ///         linear combinations in the range [0,logdet.cols()) at the particle positions in the
    ///         range [first,last).
    ///
    /// \param[in]      P            Object containing information on particle positions.
    /// \param[in]      first        Beginning of the range specifying which particles' positions to use.
    /// \param[in]      last         End of the range specifying which particles' positions to use.
    /// \param[in,out]  logdet       On input, a row-major-ordered matrix of dimension (# of particle) by (# of linear combinations).
    ///                              On exit, holds the linear combinations' values for the specified particles' positions.
    /// \param[in,out]  dlogdet      On input, a row-major-ordered matrix of dimension (# of particle) by (# of linear combinations).
    ///                              On exit, holds the linear combinations' x,y,z gradients for the specified particles' positions.
    ///                              each element of which is a length 3 vector containing the x,y,z gradients of the values in vmat.
    /// \param[in,out]  d2logdet     On input, a row-major-ordered matrix of dimension (# of particle) by (# of linear combinations).
    ///                              On exit, each element is the sum of x^2, y^2, and z^2 second derivatives of the values in logdet.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluate_notranspose(const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) {

      //app_log() << "logdet.cols()        = " << logdet.cols() << std::endl;
      //app_log() << "logdet.rows()        = " << logdet.rows() << std::endl;
      //app_log() << "dlogdet.cols()       = " << dlogdet.cols() << std::endl;
      //app_log() << "dlogdet.rows()       = " << dlogdet.rows() << std::endl;
      //app_log() << "d2logdet.cols()      = " << d2logdet.cols() << std::endl;
      //app_log() << "d2logdet.rows()      = " << d2logdet.rows() << std::endl;
      //app_log() << "this->OrbitalSetSize = " << this->OrbitalSetSize << std::endl;
      // check sanity
      this->check_input_dim("logdet # of columns", "LCOrbitalSetOpt::evaluate_notranspose", logdet.cols(), this->OrbitalSetSize);
      if ( logdet.cols() != dlogdet.cols() || logdet.cols() != d2logdet.cols() )
        throw std::runtime_error("logdet, dlogdet, and d2logdet should have the same number of columns in LCOrbitalSetOpt::evaluate_notranspose");
      if ( logdet.rows() != dlogdet.rows() || logdet.rows() != d2logdet.rows() )
        throw std::runtime_error("logdet, dlogdet, and d2logdet should have the same number of rows in LCOrbitalSetOpt::evaluate_notranspose");

      // evaluate the first logdet.cols() orbitals for the particles in the range [first, last)
      this->evaluate_notranspose_ranges(P, 'w', 0, logdet.cols(), first, last, logdet.data(), dlogdet.data(), d2logdet.data());

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Evaluates the values of the linear combinations in the range [0,psi.size()) at a
    ///         particular particle's position.
    ///
    /// \param[in]      P            Object containing information on particle positions.
    /// \param[in]      iat          Index of the particle whose position will be used.
    /// \param[in,out]  psi          On input, a vector of dimension (# of linear combinations).
    ///                              On exit, holds the linear combinations' values for the specified particle position.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi) {

      // check input vector dimension
      this->check_input_dim("psi", "LCOrbitalSetOpt::evaluate", psi.size(), this->OrbitalSetSize);

      // resize temporary arrays if necessary
      if ( m_temp_g.size() != BasisSetSize ) m_temp_g.resize(BasisSetSize);
      if ( m_temp_l.size() != BasisSetSize ) m_temp_l.resize(BasisSetSize);

      // sanity check
      if ( m_temp_g.size() < psi.size() )
        throw std::runtime_error("unexpected too-small size of m_temp_g in LCOrbitalSetOpt::evaluate");

      // evaluate the first psi.size() orbitals for the particle with index iat
      this->evaluate_notranspose_ranges(P, 'p', 0, psi.size(), iat, iat+1, psi.data(), &m_temp_g[0], &m_temp_l[0]);

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Evaluates the values, x,y,z derivatives, and x-y-z-summed second derivatives of the
    ///         linear combinations in the range [0,psi.size()) at the position of particle iat.
    ///
    /// \param[in]      P            Object containing information on particle positions.
    /// \param[in]      iat          Index of the particle whose position will be used.
    /// \param[in,out]  psi          On input, a vector of dimension (# of linear combinations).
    ///                              On exit, holds the linear combinations' values for the specified particle position.
    /// \param[in,out]  dlogdet      On input, a vector of dimension (# of linear combinations).
    ///                              On exit, holds the linear combinations' x,y,z gradients for the specified particle position.
    /// \param[in,out]  d2logdet     On input, a vector of dimension (# of linear combinations).
    ///                              On exit, each element is the sum of x^2, y^2, and z^2 second derivatives for the specified particle position.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) {

      // check sanity
      this->check_input_dim("d2psi", "LCOrbitalSetOpt::evaluate", d2psi.size(), this->OrbitalSetSize);
      this->check_input_dim( "dpsi", "LCOrbitalSetOpt::evaluate",  dpsi.size(), this->OrbitalSetSize);
      this->check_input_dim(  "psi", "LCOrbitalSetOpt::evaluate",   psi.size(), this->OrbitalSetSize);
      if ( psi.size() != dpsi.size() || psi.size() != d2psi.size() )
        throw std::runtime_error("psi, dpsi, and d2psi vectors must be the same length in LCOrbitalSetOpt::evaluate");

      // evaluate the first psi.size() orbitals and derivatives for the particle with index iat
      this->evaluate_notranspose_ranges(P, 'p', 0, psi.size(), iat, iat+1, psi.data(), dpsi.data(), d2psi.data());

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  An evaluate function that has not yet been implemented.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi) {
      throw std::runtime_error("LCOrbitalSetOpt::evaluate(P, iat, psi, dpsi, grad_grad_psi) not implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  An evaluate function that has not yet been implemented.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluate_notranspose(const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet) {
      throw std::runtime_error("LCOrbitalSetOpt::evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet) not implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  An evaluate function that has not yet been implemented.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet) {
      throw std::runtime_error("LCOrbitalSetOpt::evaluate_notranspose(P, first, last, logdet, dlogdet, grad_grad_logdet, grad_grad_grad_logdet) not implemented");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  An evaluate function that has not yet been implemented.
    ///
    ///   ------------------------------------------------ Notes on old routine ------------------------------------------------
    ///
    ///      /** evaluate psiM for virtual moves
    ///      *
    ///      * For the i-th virtual move and the j-th orbital,
    ///      * \f$ psiM(i,j)= \sum_k phiM(i,k)*C(j,k) \f$
    ///      */
    ///
    ///      in this routine, the dgemm call multiplies a  (norbs) by (basis size)   by a   (basis size) by (nparticle)   matrix
    ///
    ///      this implies that the data inside C has all (basis size) orbital coefficients for each molecular orbital stored contiguously
    ///      i.e. if the C.data() array is viewed as a column major matrix, the coefficients for one MO would all be in one column
    ///
    ///      thus in its own row major ordering, psiM is an  (nparticle) by (norbs)  matrix
    ///
    ///      thus the untransposed matrix you use above, in your column major ordering, also needs to be (nparticle) by (norbs)
    ///
    ///   ----------------------------------------------------------------------------------------------------------------------
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluateValues(const ParticleSet& P, ValueMatrix_t& psiM) {
      throw std::runtime_error("LCOrbitalSetOpt::evaluateValues(P, psiM) not yet implemented");
      //ValueMatrix_t phiM(P.getTotalNum(),BasisSetSize);
      //m_basis_set->evaluateValues(P,phiM);
      //MatrixOperators::product_ABt(phiM,C,psiM);
      ////for(int i=0; i<psiM.rows(); ++i)
      ////  for(int j=0; j<psiM.cols(); ++j)
      ////    psiM(i,j)=simd::dot(C[j],phiM[i],BasisSetSize);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  An evaluate function that has not yet been implemented.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix_t& grad_grad_grad_logdet) {
      throw std::runtime_error("LCOrbitalSetOpt::evaluateThirdDeriv(P, first, last, grad_grad_grad_logdet) not yet implemented");
    }

};

} // end namespace qmcplusplus

#endif

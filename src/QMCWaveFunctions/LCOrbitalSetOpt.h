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
/// \brief  A class for a set of optimizable linear combinations of the single particle orbitals
///         provided by either another SPO set or by a basis set of the templated type BS.
///         We refer to the linear combinations of basis orbitals as molecular orbitals.
///         The set of molecular orbitals is optimized by rotating them amongst each other.
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class BS> class LCOrbitalSetOpt : public SPOSetBase {

  // protected data members
  protected:

    /// \brief  pointer that, if not null, will be used to evaluate the basis orbitals
    SPOSetBase * m_spo_set;

    /// \brief  pointer to the basis set that evaluates the basis orbitals if m_spo_set == 0
    BS * m_basis_set;

    /// \brief  number of linear combinations of basis functions (i.e. molecular orbitals)
    int m_nlc;

    /// \brief  number of basis functions
    int m_nb;

    /// \brief  the level of printing
    int m_report_level;

    /// For use by the LCOrbitalSetOpt class, derived from this:
    /// the column-major-order m_nb by m_nlc matrix of orbital coefficients
    /// resulting from a rotation of the old coefficients
    std::vector<RealType> m_B;

    /// the column-major-order m_nb by m_nlc initial orbital coefficients
    /// at the start of the simulation, from which rotations are performed
    std::vector<RealType> m_init_B;

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

    /// \brief  vector to hold orbital indices
    std::vector<int> m_oidx;

    /// \brief  vector to hold particle indices
    std::vector<int> m_pidx;

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Place the indices from a range of indices specified by iterators into a vector.
    ///         Note that the vector may be longer than the index range (it is not shrunk to fit it)
    ///         but that an iterator to the end of the range is returned.
    ///
    /// \param[in]      start       iterator for the start of the range (should dereference to int)
    /// \param[in]      end         iterator for the end   of the range (should dereference to int)
    /// \param[in]      vec         vector to store the range in
    ///
    /// \return  iterator to the end of the entered range
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    template<class IntIter> static std::vector<int>::iterator prepare_index_vector(IntIter start, IntIter end, std::vector<int> & vec) {

      // get the length
      int length = 0;
      for (IntIter s = start; s != end; s++)
        length++;

      // expand the vector if necessary
      ensure_vector_is_big_enough(vec, length);

      // put the values in the vector
      std::copy(start, end, vec.begin());

      // return an iterator to the end of the range inside the vector
      return ( vec.begin() + length );

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Place a range of contiguous indices [start,end) into a vector.
    ///         Note that the vector may be longer than the index range (it is not shrunk to fit it)
    ///         but an iterator to the end of the range returned.
    ///
    /// \param[in]      start       start of the range
    /// \param[in]      end         end of the range
    /// \param[in]      vec         vector to store the range in
    ///
    /// \return  iterator to the end of the entered range
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    static std::vector<int>::iterator prepare_index_vector_contiguous(const int start, const int end, std::vector<int> & vec) {

      // check sanity
      if ( end < start )
        throw std::runtime_error("end is less than start in prepare_index_vector_contiguous");

      // expand the vector if necessary
      ensure_vector_is_big_enough(vec, end - start);

      // put the range into the vector
      std::vector<int>::iterator it = vec.begin();
      for(int i = start; i < end; i++, it++)
        *it = i;

      // return the end of the range
      return it;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  Enlarges the supplied vector if it is not big enough
    ///
    /// \param[in,out]  v              the vector
    /// \param[in]      n              the minimum length we want the vector to have
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    template <class T> static void ensure_vector_is_big_enough(T & v, const size_t n) {
      if ( v.size() < n )
        v.resize(n);
    }

  // public member functions
  public:

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  initializes this object and its data. In particular, the number of atomic and
    ///         molecular orbitals, and the arrays to hold the rotated and unrotated orbitals
    ///         themselves. Also performs mixing of orbitals, if requested, and prints them.
    ///
    /// \param[in]     mix_factor     factor controlling mixing of the initial orbitals
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void init_LCOrbitalSetOpt(const double mix_factor) {
      m_omixfac = mix_factor;

      m_nlc = OrbitalSetSize;
      m_nb = BasisSetSize;

      m_B.resize(m_nlc*m_nb, 0.0);
      m_init_B.resize(m_nlc*m_nb, 0.0);

      std::copy( C->data(), C->data() + m_B.size(),      m_B.begin()      );
      std::copy( C->data(), C->data() + m_init_B.size(), m_init_B.begin() );

      // if requested, mix the initial basis orbitals together
      if ( mix_factor != 0.0 ) {

        // mix
        for (int i = m_nb - 1; i >= 0; i--) {
          for (int j = 0; j < m_nlc; j++) {
            m_B.at(i+j*m_nb) += mix_factor * 2.0 * ( Random() - 0.5 );
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

        m_init_B = m_B;
      }

      this->print_B();
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief rotate m_init_B to m_B
    /// \param[in]     rot_mat     rotation matrix
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void rotate_B(const std::vector<RealType> &rot_mat)
    {
      // get the linear combination coefficients by applying the rotation to the old coefficients
      BLAS::gemm('N', 'T', m_nb, m_nlc, m_nlc, RealType(1.0), m_init_B.data(),
                 m_nb, rot_mat.data(), m_nlc, RealType(0.0), m_B.data(), m_nb);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  constructor from basis set and reporting level
    ///
    /// \param[in]      bs             pointer to the basis set to use
    /// \param[in]      rl             reporting level to use
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    LCOrbitalSetOpt(BS * const bs = 0, const int rl = 0) : m_spo_set(0), m_basis_set(0), m_report_level(rl), m_omixfac(0) {

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
    LCOrbitalSetOpt(SPOSetBase * const spo, const int rl = 0) : m_spo_set(0), m_basis_set(0), m_report_level(rl), m_omixfac(0) {

      // set the internal SPO set
      if ( spo ) this->setSPOSet(spo);

      // initialize number of molecular orbitals as zero
      this->OrbitalSetSize = 0;

    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  the destructor, which assumes deallocation of basis set is done elsewhere
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ~LCOrbitalSetOpt() { }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief  returns that this is indeed an LCOrbitalSetOpt object
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    bool is_of_type_LCOrbitalSetOpt() const { return true; }

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
      LCOrbitalSetOpt * retval;

      if ( m_spo_set )
        retval = new LCOrbitalSetOpt(m_spo_set->makeClone(), m_report_level);
      else
        retval = new LCOrbitalSetOpt(m_basis_set->makeClone(), m_report_level);

      retval->C = C;
      retval->setOrbitalSetSize(this->OrbitalSetSize);
      retval->init_LCOrbitalSetOpt(0.0);

      retval->m_B = m_B;
      retval->m_init_B = m_init_B;

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
      ensure_vector_is_big_enough(m_lc_coeffs, no * nb);
      ensure_vector_is_big_enough(m_basis_vals, np * nb);
      ensure_vector_is_big_enough(m_basis_der1, 3 * np * nb);
      ensure_vector_is_big_enough(m_basis_der2, np * nb);
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
          BLAS::copy(nb, &(*m_B.begin()) + (*it)*nb, 1, &m_lc_coeffs[i], no);
      }

      // print what is in C
      if ( false ) {
        app_log() << "printing C" << std::endl;
        std::vector<char> buff(1000, ' ');
        for (int i = 0; i < BasisSetSize; i++) {
          for (int j = 0; j < OrbitalSetSize; j++) {
            const int len = std::sprintf(&buff[0], "  %12.6f", m_B[i+j*BasisSetSize]);
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
      std::vector<int>::const_iterator oend = prepare_index_vector_contiguous(os, oe, m_oidx);

      // prepare particle list
      std::vector<int>::const_iterator pend = prepare_index_vector_contiguous(ps, pe, m_pidx);

      // evaluate
      evaluate_notranspose_general(P, mt, m_oidx.begin(), oend, m_pidx.begin(), pend, vmat, gmat, lmat);

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
    ///      i.e. if the C->data() array is viewed as a column major matrix, the coefficients for one MO would all be in one column
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

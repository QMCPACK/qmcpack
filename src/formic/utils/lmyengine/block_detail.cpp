//////////////////////////////////////////////////////////////////////////////
/// \file  formic/utils/lmyengine/block_detai.cpp
///
/// \brief  implementation file for block-recursive algorithm functions
///
//////////////////////////////////////////////////////////////////////////////

#include<vector>
#include<list>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<sstream>
//#include<mpi.h>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/matrix.h>
#include<formic/utils/numeric.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lmyengine/block_detail.h>
#include<formic/utils/lmyengine/var_dependencies.h>
#include<formic/utils/lmyengine/eigen_solver.h>
#include<formic/utils/lmyengine/davidson_solver.h>
#include<formic/utils/lmyengine/spam_solver.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Get the begining and end indices and the length for each block of variables
///
/// \param[in]    nvar         the number of variables that will be divided into blocks
/// \param[in]    nblock       the number of blocks
/// \param[out]   block_beg    on exit, the length nblock vector of indices marking the begining of each block
/// \param[out]   block_end    on exit, the length nblock vector of block lengths
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::brlm_get_block_info(const int nvar, const int nblock, std::vector<int> & block_beg, std::vector<int> & block_end, std::vector<int> & block_len) {
  
  // initialize output vectors
  block_beg.assign(nblock, 0);
  block_end.assign(nblock, 0);
  block_len.assign(nblock, 0);

  // get the length of most blocks
  const int base_block_len = nvar / nblock;

  // get the residual
  const int base_block_rem = nvar % nblock;

  int beg = 0;

  // loop over blocks
  for (int i = 0; i < nblock; i++) {
    block_beg.at(i) = beg;
    block_len.at(i) = base_block_len + ( i < base_block_rem ? 1 : 0 );
    beg += block_len.at(i);
    block_end.at(i) = beg;
  }
  if ( beg != nvar ) 
    throw formic::Exception("after block index steup, beg = %i but should have been %i in cqmc::engin::brlm_get_block_info") % beg % nvar; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Computes the first so many important update directions in the basis of those that define 
///          the basis for the provided Hamiltinonian and overlap matrices
///
/// \param[in]    nkeep         the number of important update directions to keep
/// \param[in]    curr_e        the current energy
/// \param[in]    shift_i       magnitude of the identity shift
/// \param[in]    shift_s       magnitude of the overlap shift
/// \param[in]    ovl_thresh    threshold below which overlap eigenvalues will be truncated
/// \param[in]    hh            the Hamiltonian matrix in the basis of possible directions.
///                             the first of which is the current wavefunction
/// \param[in]    ss            the overlap     matrix in the basis of possible directions, 
///                             the first of which is the current wavefunction
/// \param[in]    dd            the identity shift matrix in the basis of possible directions,
///                             the first of which is the current wavefunction
/// \param[in]    ostream       the output stream
///
/// \return a matrix holding the best nkeep directions in its column, the first of which 
///         is always the current wavefunction
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
formic::Matrix<double> cqmc::engine::get_important_brlm_dirs(const int nkeep, 
                                                             const double curr_e,
                                                             const double shift_i,
                                                             const double shift_s,
                                                             const double ovl_thresh,
                                                             const formic::Matrix<double> & hh,
                                                             const formic::Matrix<double> & ss,
                                                             const formic::Matrix<double> & dd,
                                                             std::ostream & output) {

  ////////////////////////////////////////////////////////////////////////
  // Note in this function we deal with four sets of basis vectors
  //
  //  Basis 1: The initial basis, in which hh and ss are defined.
  //           The current wave function is the first basis vector.
  //
  //  Basis 2: Basis after projecting out the current wavefunction from 
  //           the second through through last basis vectors of basis 1.
  // 
  //  Basis 3: The basis consisting of the second through last vectors 
  //           of basis 2, i.e. everything except the current wfn.
  //
  //  Basis 4: The basis after projecting out all previous updates from 
  //           the vectors in basis 2 and then orthonormalizing, dropping
  //           vectors with pre-normalized norms below ovl_thresh.
  //           The first vector in this basis is still the current wfn.
  //
  /////////////////////////////////////////////////////////////////////////

  // get the number of vectors in basis 1
  const int nd = hh.rows();

  // sanity checks
  if ( hh.rows() != hh.cols() )
    throw formic::Exception("hh (%i by %i) should be square in cqmc::engine::get_important_brlm_dirs") % hh.rows() % hh.cols();
  if ( hh.rows() != ss.rows() || hh.cols() != ss.cols() ) 
    throw formic::Exception("ss dimension (%i by %i) should match hh dimension (%i by %i) in cqmc::engine::get_important_brlm_dirs")
                  % ss.rows() % ss.cols() % hh.rows() % hh.cols();
  //if ( std::abs( 1.0 - std::abs(ss.at(0,0)) ) > 1.0e-14 ) 
  //  throw formic::Exception("expected ss(0,0) to be 1.0 in cqmc::engine::get_important_brlm_dirs");

  // put the vectors defining basis 2 in the columns of a matrix 
  formic::Matrix<double> basis2 = formic::identity_matrix<double>(nd);
  basis2.at(0,0) = 1.0/std::abs(std::sqrt(ss(0,0)));
  for (int i = 1; i < nd; i++)
   basis2.at(0,i) -= (ss(i,0)/ss(0,0));
  
  // get H and S in basis 2, including the identity shift effect on H 
  formic::Matrix<double> hh2 = basis2.t() * ( hh + shift_i * dd ) * basis2;
  formic::Matrix<double> ss2 = basis2.t() * ss * basis2;
  //output << ss2.print("%12.6f", "ss2");

  //output << hh2.print("%12.4f", "hh2") << std::endl;

  // extract overlap matrix for basis 3
  formic::Matrix<double> ss3(nd-1,nd-1);
  for (int i = 0; i < nd-1; i++) {
    for (int j = 0; j < nd-1; j++) {
      ss3.at(i,j) = ss2.at(i+1,j+1);
    }
  }

  // initialize matrix of update direction coefficients in basis 2
  formic::Matrix<double> update_dirs(nd, nkeep, 1.0);

  // ietratively find the best nkeep directions
  for (int nkept = 0; nkept < nkeep; nkept++) {
    
    // make sure previous update directions are excluded by projecting them out
    formic::Matrix<double> proj_dirs = formic::identity_matrix<double>(nd-1);
    if ( nkept > 0 ) {
      
      // get the matrix of update coefficients that excludes the coeff on the currenr wfn
      formic::Matrix<double> bigV(nd-1, nkept);
      for (int j = 0; j < nkept; j++) {
        for (int i = 0; i < nd-1; i++) {
          bigV.at(i,j) = update_dirs.at(i+1, j);
        }
        formic::ColVec<double> current_col(nd-1, bigV.col_begin(j));
        current_col /= std::abs(std::sqrt(current_col.norm2()));
      }

      // get "bare" overlap between these vectors
      formic::Matrix<double> bares = bigV.t() * bigV;

      // check the smallest eigenvalue of the bare overlap
      formic::ColVec<double> bare_evals;
      formic::Matrix<double> bare_evecs;
      bares.sym_eig(bare_evals, bare_evecs);

      // project these directions out of the full set of directions
      proj_dirs -= bigV * bares.inv() * bigV.t();
    }

    // get overlap in these projected directions
    formic::Matrix<double> ps = proj_dirs.t() * ss3 * proj_dirs;
    //output << ps.print("%12.6f", "ps") << std::endl;

    // diagonalize the projected overlap
    formic::ColVec<double> ps_evals;
    formic::Matrix<double> ps_evecs;
    ps.sym_eig(ps_evals, ps_evecs);
    //output << ps_evecs.print("%12.6f", "ps_evecs") << std::endl;
    //output << ps_evals.print("%12.6f", "ps_evals") << std::endl;

    // count number of projected overlap eigenvalues above the threshold
    //const int nabove = std::accumulate(ps_evals.begin(), ps_evals.end(), int(0), [ovl_thresh] (int i, double x) { return i + ( x > ovl_thresh ? 1 : 0 ); });
    const int nabove = std::accumulate(ps_evals.begin(), ps_evals.end(), int(0), [ovl_thresh] (int i, double x) { return i + ( x > ovl_thresh ? 1 : 0 ); });

    const int nbelow = ps_evals.size() - nabove;

    // get coefficients for basis 4 vectors in terms of basis 2 vectors
    formic::Matrix<double> basis4(nd, nabove+1, 0.0);

    // initial wfn
    basis4.at(0,0) = 1.0;

    for (int i = 0; i < nabove; i++) {
      const double eval_sqrt = std::abs(std::sqrt(ps_evals.at(nbelow+i)));
      for (int j = 0; j < nd-1; j++) {
        basis4.at(j+1, i+1) = ps_evecs.at(j, i+nbelow) / eval_sqrt;
      }
    }
    //output << basis4.print("%12.6f", "basis4") << std::endl;

    // compute overlap in basis 4 and check that it is the identity matrix
    formic::Matrix<double> ss4 = basis4.t() * ss2 * basis4;
    //output << ss4.print("%12.6f", "ss4");
    for (int i = 0; i < ss4.rows(); i++) {
      for (int j = 0; j < ss4.cols(); j++) {
        if ( std::abs( ss4.at(i,j) - ( i == j ? 1.0 : 0.0) ) > 1.0e-6 ) 
          throw formic::Exception("failure to move to orthonormal basis in cqmc::engine::get_important_brlm_dirs, ss4(%i,%i) = %.8e")
                        % i % j % ss4.at(i,j);
       }
     }

     // compute Hamiltonian in basis 4
     formic::Matrix<double> hh4 = basis4.t() * hh2 * basis4;
     //output << hh4.print("%12.6f", "hh4") << std::endl;

     // apply overlap shift to the Hamiltonian
     for (int i = 1; i < hh4.rows(); i++) 
       hh4.at(i,i) += shift_s;

     // solve for the lowest energy in eigenvector in basis 4
     formic::ColVec<std::complex<double> > e_evals;
     formic::Matrix<std::complex<double> > e_evecs;
     hh4.nonsym_eig(e_evals, e_evecs);
     //output << e_evals.print("%12.4f", "hh4_evals") << std::endl;

     // find the lowest energy eigenvalue
     double lowest_eval = e_evals.at(0).real();
     int lowest_index = 0;
     for (int i = 1; i < e_evals.size(); i++) {
       if ( e_evals.at(i).real() < lowest_eval ) {
         lowest_eval = e_evals.at(i).real();
         lowest_index = i;
       }
     }

     output << boost::format("shift_s %.2e   vec %i   target = %20.12f + %16.2f i    delta = %17.12f") 
                   % shift_s % nkept % e_evals.at(lowest_index).real() % e_evals.at(lowest_index).imag() % (e_evals.at(lowest_index).real() - curr_e)
                << std::endl;

     // check that the eigenvalue is real
     if ( std::abs( e_evals.at(lowest_index).imag() ) >  1.0e-6 ) 
       throw formic::Exception("lowest_eval has an imaginary component of %.2e in cqmc::engine::get_lm_brnr_update") % e_evals.at(lowest_index).imag(); 

     // check that the eigenvector has weight on the initial wavefunction
     if ( std::abs( e_evecs.at(0, lowest_index) ) < 1.0e-6 ) 
       throw formic::Exception("lowest_evec has a small weight of %.2e on the initial wave function in cqmc::engine::get_lm_brnr_update") % std::abs(e_evecs.at(0,lowest_index));

     // attempt to get the eigenvector in terms of real numbers
     formic::ColVec<double> evec_real(e_evecs.rows(), 0.0);
     for (int i = 0; i < e_evecs.rows(); i++) {
       std::complex<double> val = e_evecs.at(i, lowest_index) / e_evecs.at(0, lowest_index);
       if ( std::abs(val.imag()) > 1.0e-6 ) 
         throw formic::Exception("lowest_evec element %i has an un-removed imaginary component of %.2e in cqmc::engine::get_lm_brnr_update") % i % val.imag();
       evec_real.at(i) = val.real();
     }

     // convert the eigenvector to basis 2
     ( basis4 * evec_real ) >>= update_dirs.col_begin(nkept);

     // check that this eigenvector still has weight on the initial wfn
     if ( std::abs( update_dirs.at(0, nkept) ) < 1.0e-6 )
       throw formic::Exception("basis 2 eigenvector has a small weight of %.2e on the initial wave function in cqmc::engine::get_lm_brnr_update")
                     % std::abs( update_dirs.at(0, nkept) );

     // scale the eigenvector in basis 2 so that the initial wavefunction has unit weight
     formic::ColVec<double>(update_dirs.rows(), update_dirs.col_begin(nkept)) /= *update_dirs.col_begin(nkept);

   }

   // convert the eigenvectors to basis 1
   update_dirs = basis2 * update_dirs;

   // scale the eigenvectors so that in basis 1 they have unit weight on the initial wave function 
   for (int nkept = 0; nkept < nkeep; nkept++) {
     if ( std::abs( *update_dirs.col_begin(nkept) ) < 1.0e-6 )
       throw formic::Exception("basis 1 eigenvector has an small weight of %.2e on the initial wave function in cqmc::engine::get_lm_brnr_update")
                     % std::abs( *update_dirs.col_begin(nkept) );
     formic::ColVec<double>(update_dirs.rows(), update_dirs.col_begin(nkept)) /= (*update_dirs.col_begin(nkept));
   }

   // clean up memory used by matrices and vectors
   formic::reusable_array_garbage_collect();

   // return the eigenvectors in basis 1
   return update_dirs;
 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   Computes the first so many important update directions in the basis of those that define 
///          the basis for the provided Hamiltinonian and overlap matrices
///
/// \param[in]    nkeep         the number of important update directions to keep
/// \param[in]    curr_e        the current energy
/// \param[in]    omega         energy shift used for harmonic davidson 
/// \param[in]    shift_i       magnitude of the identity shift
/// \param[in]    shift_s       magnitude of the overlap shift
/// \param[in]    ovl_thresh    threshold below which overlap eigenvalues will be truncated
/// \param[in]    hh            the Hamiltonian matrix in the basis of possible directions.
///                             the first of which is the current wavefunction
/// \param[in]    ss            the overlap     matrix in the basis of possible directions, 
///                             the first of which is the current wavefunction
/// \param[in]    dd            the identity shift matrix in the basis of possible directions,
///                             the first of which is the current wavefunction
/// \param[in]    ostream       the output stream
///
/// \return a matrix holding the best nkeep directions in its column, the first of which 
///         is always the current wavefunction
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
formic::Matrix<double> cqmc::engine::get_important_brlm_dirs_davidson(const formic::VarDeps * dep_ptr,
                                                                      const int nkeep, 
                                                                      const double omega,
                                                                      const double curr_cost,
                                                                      const double shift_i,
                                                                      const double shift_s,
                                                                      const double ovl_thresh,
                                                                      formic::Matrix<double> & hh,
                                                                      formic::Matrix<double> & ss,
                                                                      formic::Matrix<double> & dd,
                                                                      std::ostream & output) {


  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();

  // check to see whether nkeep is one and throw out an error if not
  if ( nkeep != 1 ) 
    throw formic::Exception("nkeep must be 1 in get_important_brlm_dirs_davidson!");

  // get the number of vectors in basis 1
  int nd = hh.rows();
  formic::mpi::bcast(&nd, 1);
  //MPI_Bcast(&nd, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // get the flag that whether we will use variable dependency system
  const bool use_var_deps = ( dep_ptr->n_ind() != dep_ptr->n_tot() );

  // get the flag that whether the calculation is ground or targeted excited 
  bool ground = true;
  if ( my_rank == 0 ) 
    ground = ( std::abs(ss.at(0,0) - 1.0) < 1.0e-6 );
  formic::mpi::bcast(&ground, 1);
  //MPI_Bcast(&ground, 1, MPI::BOOL, 0, MPI_COMM_WORLD);

  // initial energy/target fn value
  double init_cost = 0.0;
  if ( my_rank == 0 ) 
    init_cost = hh.at(0,0)/ss.at(0,0);
  formic::mpi::bcast(&init_cost, 1);
  //MPI_Bcast(&init_cost, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // sanity checks
  if ( my_rank == 0 ) { 
    if ( hh.rows() != hh.cols() ) 
      throw formic::Exception("hh (%i by %i) should be square in cqmc::engine::get_important_brlm_dirs") % hh.rows() % hh.cols();
    if ( hh.rows() != ss.rows() || hh.cols() != ss.cols() ) 
      throw formic::Exception("ss dimension (%i by %i) should match hh dimension (%i by %i) in cqmc::engine::get_important_brlm_dirs")
                    % ss.rows() % ss.cols() % hh.rows() % hh.cols();
  }

  // initialize matrix of update direction coefficients
  formic::Matrix<double> update_dirs(nd, nkeep, 1.0);

  // add identity shift to the Hamiltonian 
  formic::Matrix<double> m_hh = hh + shift_i * dd;

  // create eigensolver
  std::vector<double> vf_var;
  boost::shared_ptr< cqmc::engine::EigenSolver > eigensolver(new cqmc::engine::DavidsonLMHD(dep_ptr,
                                                                                            nd,
                                                                                            60,
                                                                                            1.0e-7,
                                                                                            ovl_thresh,
                                                                                            false, // var deps use?
                                                                                            true,  // chase lowest
                                                                                            false,  // chase closest
                                                                                            ground,  // ground state? 
                                                                                            false,  // variance correct
                                                                                            true, // build matrix
                                                                                            vf_var,
                                                                                            init_cost, // initial energy
                                                                                            0.0, // initial variance
                                                                                            omega, // omega
                                                                                            0.0, // var weight
                                                                                            20.0, // maximun energy change
                                                                                            0.0, // total weight
                                                                                            0.0, // total weight
                                                                                            hh,
                                                                                            hh,
                                                                                            m_hh,
                                                                                            ss,
                                                                                            ss));
  

  // iteratively find the best nkeep directions CURRENTLY ONLY WORK FOR NKEEP=1!!!!!
  for (int nkept = 0; nkept < nkeep; nkept++) {

    // reset the eigensolver
    eigensolver -> reset();

    // update shift
    eigensolver -> update_lm_shift(0.0, shift_s);

    // add first krylov vector
    { 
      formic::ColVec<double> temp(hh.cols());
      for (int j = 0; j < temp.size(); j++) 
        temp.at(j) = ( j == 0 ? 1.0 : 0.0);
      eigensolver -> add_krylov_vector(temp);
    }


    // solve the eigenvalue problem
    double davidson_eval;
    bool solve_shift = eigensolver -> solve(davidson_eval, output);

    // scale the eigenvector so that the initial wavefunction has the unit weight
    if ( my_rank == 0 ) {
      eigensolver -> convert_to_wf_coeff();
      formic::ColVec<double> evec_eigen = eigensolver -> wf_coeff();
    
      // compute the largest weight on the derivative vector
      double max_update_abs_value = std::abs(evec_eigen.at(1));

      // store this eigenvector
      evec_eigen >>= update_dirs.col_begin(nkept);

      for (int k = 0; k < nd; k++) {
        if ( k != 0 )
          max_update_abs_value = std::abs(evec_eigen.at(k)) > max_update_abs_value ? std::abs(evec_eigen.at(k)) : max_update_abs_value;
      }
      output << boost::format("The largest weight on the derivative vector for shift %.4e is %.6e") % shift_i % max_update_abs_value << std::endl << std::endl;
    }
  }

  // clean up memory used by matrices and vectors
  formic::reusable_array_garbage_collect();

  // return the eigenvectors in basis 1
  return update_dirs;
}


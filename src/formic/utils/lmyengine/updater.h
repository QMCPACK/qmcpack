/////////////////////////////////////////////////////////////////////////////////////////////////
// \file cqmc/chdengine/updater.h
//
//
// \brief header file for harmonic davidson method updater
//
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_UPDATER_HEADER
#define ENGINE_UPDATER_HEADER

#include<vector>

#include<boost/shared_ptr.hpp>

#include"formic/utils/zero_one.h"
#include"formic/utils/mpi_interface.h"
#include"formic/utils/lmyengine/matrix_builder.h"
#include"formic/utils/lmyengine/eigen_solver.h"
#include"formic/utils/lmyengine/davidson_solver.h"
#include"formic/utils/lmyengine/spam_solver.h"
#include"formic/utils/lmyengine/energy_target_accu.h"
#include"formic/utils/numeric.h"
#include"formic/utils/matrix.h"
#include"formic/utils/lmyengine/var_dependencies.h"

namespace cqmc {

  namespace engine {

    ////////////////////////////////////////////////////////////////////////////////////////////
    // \brief A class to perform harmonic linear method update of the current wave function 
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    template<typename S> class HDLinMethodUpdater {

      private:
        
      /// \brief harmonic Davidson shift(omega)
      double _omega;

      /// \brief weight of variance correction term
      double _var_weight;

      /// \brief bare derivative vectors 
      formic::Matrix<S> & _der_rat;

      /// \brief energy derivative vectors
      formic::Matrix<S> & _le_der;

      /// \brief |value/guiding|^2 list
      const std::vector<double> & _vgs;

      /// \brief weight list
      const std::vector<double> & _weight;

      /// \brief vector to store linear method shift scale
      const std::vector<double> & _shift_scale;

      /// \brief a flag to tell whether we are doing ground state calculation
      bool _ground_state;

      /// \brief a flag to tell whether to do variance correct calculation for excited states
      bool _variance_correct;

      /// \brief a flag to tell whether to build limear method matrix explicitly 
      bool _build_lm_matrix;

      public:
        
      HDLinMethodUpdater(formic::Matrix<S> & der_rat,
                         formic::Matrix<S> & le_der,
                         const std::vector<double> & vgs,
                         const std::vector<double> & weight,
                         const std::vector<double> & shift_scale,
                         const double omega,
                         const double var_weight,
                         const bool ground_state,
                         const bool variance_correct,
                         const bool build_lm_matrix)
      :_omega(omega),
      _var_weight(var_weight),
      _der_rat(der_rat),
      _le_der(le_der),
      _vgs(vgs),
      _weight(weight),
      _shift_scale(shift_scale),
      _ground_state(ground_state),
      _variance_correct(variance_correct),
      _build_lm_matrix(build_lm_matrix)
      {}

      void engine_update_build_matrix(const formic::VarDeps * dep_ptr,
                                      const int lm_krylov_iter,
                                      const int lm_spam_inner_iter,
                                      const double lm_eigen_thresh,
                                      const double lm_min_S_eval,
                                      const bool spam_use,
                                      const bool var_deps_use,
                                      const bool chase_lowest,
                                      const bool chase_closest,
                                      const bool print_matrix,
                                      const double init_cost,
                                      const double init_variance,
                                      const double lm_max_e_change,
                                      const double lm_ham_shift_i,
                                      const double lm_ham_shift_s,
                                      const double _omega,
                                      const double lm_max_update_abs,
                                      std::vector<S> & vf_var,
                                      std::vector<bool> & good_solve,
                                      std::vector<int> & shift_solved,
                                      std::ostream & output)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();

        // creates matrix builder 
        boost::shared_ptr< cqmc::engine::HamOvlpBuilderHD<S> > mbuilder( new cqmc::engine::HamOvlpBuilderHD<S>(_der_rat, 
                                                                                                               _le_der, 
                                                                                                               _le_der, 
                                                                                                               _vgs, 
                                                                                                               _weight, 
                                                                                                               _omega, 
                                                                                                               0,
                                                                                                               1, 
                                                                                                               false, 
                                                                                                               _ground_state, 
                                                                                                               _variance_correct,
                                                                                                               _build_lm_matrix,
                                                                                                               false, 
                                                                                                               print_matrix));
        
        // build the matrix 
        mbuilder -> MatrixBuild(output);

        // get the hamiltonian and overlap matrix 
        formic::Matrix<S> & hamiltonian = mbuilder -> ham();
        formic::Matrix<S> & overlap = mbuilder -> ovl();

        // if we want to correct the finite variance issue, get nomral linear method's overlap matrix 
        formic::Matrix<S> & LMoverlap = mbuilder -> lovl();

        const bool print_mats = false;

        boost::shared_ptr< cqmc::engine::EigenSolver<S> > eigensolver(new cqmc::engine::DavidsonLMHD<S>(dep_ptr,
                                                                                                        hamiltonian.cols(),
                                                                                                        lm_krylov_iter,
                                                                                                        lm_eigen_thresh,
                                                                                                        lm_min_S_eval,
                                                                                                        var_deps_use,
                                                                                                        chase_lowest,
                                                                                                        chase_closest,
                                                                                                        _ground_state,
                                                                                                        _variance_correct,
                                                                                                        _build_lm_matrix,
                                                                                                        vf_var,
                                                                                                        init_cost,
                                                                                                        init_variance,
                                                                                                        _omega,
                                                                                                        _var_weight,
                                                                                                        lm_max_e_change,
                                                                                                        mbuilder -> total_weight(),
                                                                                                        mbuilder -> vgsa(),
                                                                                                        _der_rat,
                                                                                                        _le_der,
                                                                                                        hamiltonian,
                                                                                                        overlap,
                                                                                                        LMoverlap));

        // set a flag to tell whether previous solve is good or not 
        bool previous_solve = false;

        // loop over different shift 
        for (int i = 0; i < _shift_scale.size(); i++) {

          // if this shift has been solved already, skip it
          if (shift_solved.at(i)) 
            continue;

          // get scaled shifts
          double x = _shift_scale.at(i) * lm_ham_shift_i;
          double y = _shift_scale.at(i) * lm_ham_shift_s;

          // print which shift's eigenpair is being solved
          if (my_rank == 0)
            output << boost::format("solving harmonic davidson linear method for identity shift %.4e, overlap shift %.4e and omega %.4e") % x % y % _omega << std::endl << std::endl;

          // reset the eigensolver for this shift
          eigensolver -> reset();
          eigensolver -> update_lm_shift(x, y);
          { 
            formic::ColVec<S> temp(hamiltonian.cols());
            for (int j = 0; j < temp.size(); j++) 
              temp.at(j) = ( j == 0 ? formic::unity(S()) : formic::zero(S()));
            eigensolver -> add_krylov_vector(temp);
          }

          // solve the eigenvalue problem
          double davidson_eval;
          bool solve_shift = eigensolver -> solve(davidson_eval, output);

          // converts the eigenvector to wave function coefficient and set the bad solve flag if the eigensolver found an eigenvector that is not dominated by the current state on root process 
          if (my_rank == 0) {
            eigensolver -> convert_to_wf_coeff();
            formic::ColVec<S> evec_eigen = eigensolver -> wf_coeff();
            double max_update_abs_value = std::abs(evec_eigen.at(1));
            // number of total variables + 1
            int n_tot = (var_deps_use ? (dep_ptr -> n_tot() + 1) : _der_rat.cols());

            for (int k = 0; k < n_tot; k++) {
              vf_var.at(i * n_tot + k) = evec_eigen.at(k);
              if ( k != 0 )
                max_update_abs_value = std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) > max_update_abs_value ? std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) : max_update_abs_value;
            }
            solve_shift = (solve_shift && (max_update_abs_value < lm_max_update_abs));
            output << boost::format("The largest weight on the derivative vector for shift %.4e is %.6e") % x % max_update_abs_value << std::endl << std::endl;
          }

          // record the result of this shift 
          good_solve.at(i) = solve_shift;
          
          // record the result of previous solve 
          previous_solve = solve_shift;

          // if the previous solve fails(imaginary energy gained), reset eigensolver 
          if ( !previous_solve && my_rank == 0 ) {
            output << boost::format("The largest weight on derivative vector is too large or it's a bad solve, and this update will not be used") << std::endl << std::endl;
            eigensolver -> reset();
          }

        }

        // save the wavefunction variables resulting from the update 
        formic::mpi::bcast(&vf_var.at(0), vf_var.size());

      }

      void engine_update_build_matrix(const formic::VarDeps * dep_ptr,
                                      const int lm_krylov_iter,
                                      const int lm_spam_inner_iter,
                                      const double lm_eigen_thresh,
                                      const double lm_min_S_eval,
                                      const bool spam_use,
                                      const bool var_deps_use,
                                      const bool chase_lowest,
                                      const bool chase_closest,
                                      const bool print_matrix,
                                      const double init_cost,
                                      const double init_variance,
                                      const double lm_max_e_change,
                                      const double lm_ham_shift_i,
                                      const double lm_ham_shift_s,
                                      const double _omega,
                                      const double lm_max_update_abs,
                                      formic::Matrix<S> & hh,
                                      formic::Matrix<S> & ss,
                                      formic::Matrix<S> & lss,
                                      std::vector<S> & vf_var,
                                      std::vector<bool> & good_solve,
                                      std::vector<int> & shift_solved,
                                      std::ostream & output)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();

        const bool print_mats = false;

        boost::shared_ptr< cqmc::engine::EigenSolver<S> > eigensolver(new cqmc::engine::DavidsonLMHD<S>(dep_ptr,
                                                                                                        hh.cols(),
                                                                                                        lm_krylov_iter,
                                                                                                        lm_eigen_thresh,
                                                                                                        lm_min_S_eval,
                                                                                                        var_deps_use,
                                                                                                        chase_lowest,
                                                                                                        chase_closest,
                                                                                                        _ground_state,
                                                                                                        _variance_correct,
                                                                                                        _build_lm_matrix,
                                                                                                        vf_var,
                                                                                                        init_cost,
                                                                                                        init_variance,
                                                                                                        _omega,
                                                                                                        _var_weight,
                                                                                                        lm_max_e_change,
                                                                                                        0.0,
                                                                                                        0.0,
                                                                                                        _der_rat,
                                                                                                        _le_der,
                                                                                                        hh,
                                                                                                        ss,
                                                                                                        lss));

        // set a flag to tell whether previous solve is good or not 
        bool previous_solve = false;

        // loop over different shift 
        for (int i = 0; i < _shift_scale.size(); i++) {

          // if this shift has been solved already, skip it
          if (shift_solved.at(i)) 
            continue;

          // get scaled shifts
          double x = _shift_scale.at(i) * lm_ham_shift_i;
          double y = _shift_scale.at(i) * lm_ham_shift_s;

          // print which shift's eigenpair is being solved
          if (my_rank == 0)
            output << boost::format("solving harmonic davidson linear method for identity shift %.4e, overlap shift %.4e and omega %.4e") % x % y % _omega << std::endl << std::endl;

          // reset the eigensolver for this shift
          eigensolver -> reset();
          eigensolver -> update_lm_shift(x, y);
          { 
            formic::ColVec<S> temp(hh.cols());
            for (int j = 0; j < temp.size(); j++) 
              temp.at(j) = ( j == 0 ? formic::unity(S()) : formic::zero(S()));
            eigensolver -> add_krylov_vector(temp);
          }

          // solve the eigenvalue problem
          double davidson_eval;
          bool solve_shift = eigensolver -> solve(davidson_eval, output);

          // converts the eigenvector to wave function coefficient and set the bad solve flag if the eigensolver found an eigenvector that is not dominated by the current state on root process 
          if (my_rank == 0) {
            eigensolver -> convert_to_wf_coeff();
            formic::ColVec<S> evec_eigen = eigensolver -> wf_coeff();
            double max_update_abs_value = std::abs(std::sqrt(formic::square_norm(evec_eigen.at(1))));
            // number of total variables + 1
            int n_tot = (var_deps_use ? (dep_ptr -> n_tot() + 1) : (dep_ptr -> n_tot() + 1));

            for (int k = 0; k < n_tot; k++) {
              vf_var.at(i * n_tot + k) = evec_eigen.at(k);
              if ( k != 0 )
                max_update_abs_value = std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) > max_update_abs_value ? std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) : max_update_abs_value;
            }
            solve_shift = (solve_shift && (max_update_abs_value < lm_max_update_abs));
            output << boost::format("The largest weight on the derivative vector for shift %.4e is %.6e") % x % max_update_abs_value << std::endl << std::endl;
          }

          // record the result of this shift 
          good_solve.at(i) = solve_shift;
          
          // record the result of previous solve 
          previous_solve = solve_shift;

          // if the previous solve fails(imaginary energy gained), reset eigensolver 
          if ( !previous_solve && my_rank == 0 ) {
            output << boost::format("The largest weight on derivative vector is too large or it's a bad solve, and this update will not be used") << std::endl << std::endl;
            eigensolver -> reset();
          }

        }

        // save the wavefunction variables resulting from the update 
        formic::mpi::bcast(&vf_var.at(0), vf_var.size());
      }

      void engine_update_no_matrix(const formic::VarDeps * dep_ptr,
                                   const int lm_krylov_iter,
                                   const int lm_spam_inner_iter,
                                   const double lm_eigen_thresh,
                                   const double lm_min_S_eval,
                                   const bool spam_use,
                                   const bool var_deps_use,
                                   const bool chase_lowest,
                                   const bool chase_closest,
                                   const bool print_matrix,
                                   const double init_cost,
                                   const double lm_max_e_change,
                                   const double lm_ham_shift_i,
                                   const double lm_ham_shift_s,
                                   const double _omega,
                                   const double lm_max_update_abs,
                                   std::vector<S> & vf_var,
                                   std::vector<bool> & good_solve,
                                   std::vector<int> & shift_solved,
                                   std::ostream & output)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();

        // get the number of total variables + 1
        const int Ntot = (var_deps_use ? (dep_ptr -> n_tot() + 1) : _der_rat.cols());

        // creates matrix builder 
        boost::shared_ptr< cqmc::engine::HamOvlpBuilderHD<S> > mbuilder( new cqmc::engine::HamOvlpBuilderHD<S>(_der_rat, 
                                                                                                               _le_der, 
                                                                                                               _le_der, 
                                                                                                               _vgs, 
                                                                                                               _weight, 
                                                                                                               _omega, 
                                                                                                               0,
                                                                                                               1, 
                                                                                                               false, 
                                                                                                               _ground_state, 
                                                                                                               false,
                                                                                                               _build_lm_matrix,
                                                                                                               false, 
                                                                                                               print_matrix));
        



        // since we don't build linear method matrix, modify derivative vector matrix and evaluate total weight and average of |value/guiding|^2 values 
        double prefactor = mbuilder -> MatrixAbsorb();     
          

        // create eigen solver
        boost::shared_ptr< cqmc::engine::EigenSolver<S> > eigensolver(new cqmc::engine::DavidsonLMHD<S>(dep_ptr,
                                                                                                        _der_rat.cols(),
                                                                                                        lm_krylov_iter,
                                                                                                        lm_eigen_thresh,
                                                                                                        lm_min_S_eval,
                                                                                                        var_deps_use,
                                                                                                        chase_lowest,
                                                                                                        chase_closest,
                                                                                                        _ground_state,
                                                                                                        false,
                                                                                                        _build_lm_matrix,
                                                                                                        vf_var,
                                                                                                        init_cost,
                                                                                                        0.0,
                                                                                                        _omega,
                                                                                                        0.0,
                                                                                                        lm_max_e_change,
                                                                                                        mbuilder -> total_weight(),
                                                                                                        mbuilder -> vgsa(),
                                                                                                        _der_rat,
                                                                                                        _le_der,
                                                                                                        _der_rat,
                                                                                                        _le_der,
                                                                                                        _le_der));


        // set a flag to tell whether previous solve is good or not 
        bool previous_solve = false;

        // loop over different shift 
        for (int i = 0; i < _shift_scale.size(); i++) {

          // if this shift has been solved already, skip it
          if (shift_solved.at(i)) 
            continue;

          // get scaled shifts
          double x = _shift_scale.at(i) * lm_ham_shift_i;
          double y = _shift_scale.at(i) * lm_ham_shift_s;

          // print which shift's eigenpair is being solved
          if (my_rank == 0)
            output << boost::format("solving harmonic davidson linear method for identity shift %.4e, overlap shift %.4e and omega %.4e with matrix-free Arnoldi") % x % y % _omega << std::endl << std::endl;

          // if previous solve succeeds, we only need to update hamiltonian vector product and hamiltonian projection 
          if ( previous_solve && i !=0 ) 
            eigensolver -> update_hvecs_sub(x, y); 

          // update shift
          eigensolver -> update_lm_shift(x, y);

          // if previous solve fails(imaginary energy) or this is the first shift, add initial guess to the solver(krylov subspace)
          if ( !previous_solve ) { 
            const int m = _der_rat.cols();
            formic::ColVec<S> temp(m);
            for (int j = 0; j < temp.size(); j++) 
              temp.at(j) = ( j == 0 ? formic::unity(S()) : formic::zero(S()));
            eigensolver -> add_krylov_vector(temp);
          }

          // solve the eigenvalue problem
          double davidson_eval;
          bool solve_shift = eigensolver -> iterative_solve(davidson_eval, output);

          // converts the eigenvector to wave function coefficient and set the bad solve flag if the eigensolver found an eigenvector that is not dominated by the current state on root process 
          if (my_rank == 0) {
            eigensolver -> convert_to_wf_coeff();
            formic::ColVec<S> evec_eigen = eigensolver -> wf_coeff();
            double max_update_abs_value = std::abs(std::sqrt(formic::square_norm(evec_eigen.at(1))));

            for (int k = 0; k < Ntot; k++) {
              vf_var.at(i * (Ntot) + k) = evec_eigen.at(k);
              if ( k != 0 )
                max_update_abs_value = std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) > max_update_abs_value ? std::abs(std::sqrt(formic::square_norm(evec_eigen.at(k)))) : max_update_abs_value;
            }
            solve_shift = (solve_shift && (max_update_abs_value < lm_max_update_abs));
            output << boost::format("The largest weight on the derivative vector for shift %.4e is %.6e") % x % max_update_abs_value << std::endl << std::endl;

            // record the result of this shift 
            good_solve.at(i) = solve_shift;
          
            // record the result of previous solve 
            previous_solve = solve_shift;

          }

          // broadcast solve results to all processes
          formic::mpi::bcast(&previous_solve, 1);

          // if the previous solve fails(imaginary energy gained), reset eigensolver 
          if ( !previous_solve && my_rank == 0 ) {
            output << boost::format("The largest weight on derivative vector is too large or it's a bad solve, and this update will not be used") << std::endl << std::endl;
            eigensolver -> reset();
          }

        }

        // save the wavefunction variables resulting from the update 
        formic::mpi::bcast(&vf_var.at(0), vf_var.size());

        // after solve for eigenvector, we recover the original derivative vectors
        mbuilder -> MatrixRecover(); 

      }
      
      void engine_update_spam(const formic::VarDeps * dep_ptr,
                              const int lm_krylov_iter,
                              const int lm_spam_inner_iter,
                              const int appro_degree,
                              const double lm_eigen_thresh,
                              const double lm_min_S_eval,
                              const bool spam_use,
                              const bool var_deps_use,
                              const bool chase_lowest,
                              const bool chase_closest,
                              const bool print_matrix,
                              const double init_cost,
                              const double lm_max_e_change,
                              const double lm_ham_shift_i,
                              const double lm_ham_shift_s,
                              const double _omega,
                              const double lm_max_update_abs,
                              std::vector<double> & vf_var,
                              std::vector<bool> & good_solve,
                              std::vector<int> & shift_solved,
                              std::ostream & output)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();

        // get the number of total variables + 1 
        const int Ntot = ( var_deps_use ? (dep_ptr -> n_tot() + 1) : _der_rat.cols());

        formic::Matrix<double> _der_rat_real(_der_rat.rows(), _der_rat.cols(), 0.0);
        formic::Matrix<double> _le_der_real(_le_der.rows(), _le_der.cols(), 0.0);
        for (int i=0; i<_der_rat.rows(); i++) {
          for (int j=0; j<_der_rat.cols(); j++) {
            _der_rat_real.at(i,j) = formic::real(_der_rat.at(i,j));
            _le_der_real.at(i,j)  = formic::real(_le_der.at(i,j));
          }
        }
        // creates matrix builder 
        boost::shared_ptr< cqmc::engine::HamOvlpBuilderHD<double> > mbuilder( new cqmc::engine::HamOvlpBuilderHD<double>(_der_rat_real, 
                                                                                                                         _le_der_real, 
                                                                                                                         _le_der_real, 
                                                                                                                         _vgs, 
                                                                                                                         _weight, 
                                                                                                                         _omega, 
                                                                                                                         0,
                                                                                                                         appro_degree, 
                                                                                                                         spam_use, 
                                                                                                                         _ground_state, 
                                                                                                                         false,
                                                                                                                         _build_lm_matrix,
                                                                                                                         false, 
                                                                                                                         print_matrix));
        



        // since we don't build linear method matrix, modify derivative vector matrix and evaluate total weight and average of |value/guiding|^2 values 
        double prefactor = mbuilder -> MatrixAbsorb();     
          

        // create eigensolver
        boost::shared_ptr< cqmc::engine::EigenSolver<double> > eigensolver(new cqmc::engine::SpamLMHD(dep_ptr, 
                                                                                                      _der_rat.cols(),
                                                                                                      lm_krylov_iter,
                                                                                                      lm_spam_inner_iter, 
                                                                                                      appro_degree,
                                                                                                      false,
                                                                                                      lm_eigen_thresh,
                                                                                                      lm_min_S_eval,
                                                                                                      prefactor,
                                                                                                      var_deps_use,
                                                                                                      chase_lowest,
                                                                                                      chase_closest,
                                                                                                      _ground_state,
                                                                                                      vf_var,
                                                                                                      init_cost,
                                                                                                      _omega,
                                                                                                      lm_max_e_change,
                                                                                                      mbuilder -> total_weight(),
                                                                                                      mbuilder -> vgsa(),
                                                                                                      _der_rat_real,
                                                                                                      _le_der_real,
                                                                                                      mbuilder -> approximate_der_vec(),
                                                                                                      mbuilder -> approximate_le_der()));



        // set a flag to tell whether previous solve is good or not 
        bool previous_solve = false;

        // loop over different shift 
        for (int i = 0; i < _shift_scale.size(); i++) {

          // if this shift has been solved already, skip it
          if (shift_solved.at(i)) 
            continue;

          // get scaled shifts
          double x = _shift_scale.at(i) * lm_ham_shift_i;
          double y = _shift_scale.at(i) * lm_ham_shift_s;

          // print which shift's eigenpair is being solved
          if (my_rank == 0)
            output << boost::format("solving harmonic davidson linear method for identity shift %.4e, overlap shift %.4e and omega %.4e with SPAM") % x % y % _omega << std::endl << std::endl;

          // if previous solve succeeds, we only need to update hamiltonian vector product and hamiltonian projection 
          if ( previous_solve && i !=0 ) 
            eigensolver -> update_hvecs_sub(x, y); 

          // update shift
          eigensolver -> update_lm_shift(x, y);

          // if previous solve fails(imaginary energy) or this is the first shift, add initial guess to the solver(krylov subspace)
          if ( !previous_solve ) { 
            const int m = _der_rat.cols();
            formic::Matrix<double> temp(m, 1);
            for (int j = 0; j < m; j++) 
              temp.at(j, 0) = ( j == 0 ? 1.0 : 0.0);
            eigensolver -> add_krylov_vectors_outer(temp);
          }

          // solve the eigenvalue problem
          double davidson_eval;
          bool solve_shift = eigensolver -> iterative_solve(davidson_eval, output);

          // converts the eigenvector to wave function coefficient and set the bad solve flag if the eigensolver found an eigenvector that is not dominated by the current state on root process 
          if (my_rank == 0) {
            eigensolver -> convert_to_wf_coeff();
            formic::ColVec<double> evec_eigen = eigensolver -> wf_coeff();
            double max_update_abs_value = std::abs(evec_eigen.at(1));

            for (int k = 0; k < Ntot; k++) {
              vf_var.at(i * (Ntot) + k) = evec_eigen.at(k);
              if ( k != 0 )
                max_update_abs_value = std::abs(evec_eigen.at(k)) > max_update_abs_value ? std::abs(evec_eigen.at(k)) : max_update_abs_value;
            }
            solve_shift = (solve_shift && (max_update_abs_value < lm_max_update_abs));
            output << boost::format("The largest weight on the derivative vector for shift %.4e is %.6e") % x % max_update_abs_value << std::endl << std::endl;

            // record the result of this shift 
            good_solve.at(i) = solve_shift;
          
            // record the result of previous solve 
            previous_solve = solve_shift;

          }

          
          // broadcast solve results to all processes
          formic::mpi::bcast(&previous_solve, 1);
          
          // if the previous solve fails(imaginary energy gained), reset eigensolver 
          if ( !previous_solve ) {
            if ( my_rank == 0 ) 
              output << boost::format("The largest weight on derivative vector is too large or it's a bad solve, and this update will not be used") << std::endl << std::endl;
            eigensolver -> reset();
          }

        }

        // save the wavefunction variables resulting from the update 
        formic::mpi::bcast(&vf_var.at(0), vf_var.size());

        // after solve for eigenvector, we recover the original derivative vectors
        mbuilder -> MatrixRecover(); 

      }
   }; // end of class definition
 } // end of engine
} // end of cqmc

#endif

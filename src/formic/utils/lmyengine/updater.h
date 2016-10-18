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

#include<formic/utils/matrix.h>
#include<formic/utils/lmyengine/var_dependencies.h>

namespace cqmc {

  namespace engine {

    ////////////////////////////////////////////////////////////////////////////////////////////
    // \brief A class to perform harmonic linear method update of the current wave function 
    //
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    class HDLinMethodUpdater {

      private:
        
      /// \brief harmonic Davidson shift(omega)
      double _omega;

      /// \brief weight of variance correction term
      double _var_weight;

      /// \brief bare derivative vectors 
      formic::Matrix<double> & _der_rat;

      /// \brief energy derivative vectors
      formic::Matrix<double> & _le_der;

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
        
      HDLinMethodUpdater(formic::Matrix<double> & der_rat,
                         formic::Matrix<double> & le_der,
                         const std::vector<double> & vgs,
                         const std::vector<double> & weight,
                         const std::vector<double> & shift_scale,
                         const double omega,
                         const double var_weight,
                         const bool ground_state,
                         const bool variance_correct,
                         const bool build_lm_matrix);

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
                                      std::vector<double> & vf_var,
                                      std::vector<bool> & good_solve,
                                      std::vector<int> & shift_solved,
                                      std::ostream & output);

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
                                      formic::Matrix<double> & hh,
                                      formic::Matrix<double> & ss,
                                      std::vector<double> & vf_var,
                                      std::vector<bool> & good_solve,
                                      std::vector<int> & shift_solved,
                                      std::ostream & output);

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
                                   std::vector<double> & vf_var,
                                   std::vector<bool> & good_solve,
                                   std::vector<int> & shift_solved,
                                   std::ostream & output);
      
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
                              std::ostream & output);


	
      };

    }

}

#endif

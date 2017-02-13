////////////////////////////////////////////////////////////////////////////////////////////////
// \brief matrix build class that builds harmonic davidson matrix
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EIGINE_MATRIX_BUILDER_HEADER
#define EIGINE_MATRIX_BUILDER_HEADER

#include<complex>
#include<vector>
#include<numeric>
#include<cassert>
#include<algorithm>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

//#include<mpi.h>

#include<formic/utils/lmyengine/var_dependencies.h>
#include<formic/utils/matrix.h>


namespace cqmc {

  namespace engine {

    class HamOvlpBuilderHD{

      private:

        /// \brief [out] harmonic davidson hamiltonian matrix 
        formic::Matrix<double> _hmat;
        std::vector<formic::Matrix<double> > _hmat_temp;

        /// \brief [out] harmonic davidson overlap matrix 
        formic::Matrix<double> _smat;
        std::vector<formic::Matrix<double> > _smat_temp;

        /// \brief [out] harmonic davidson approximate hamiltonian matrix used in SPAM algorithm 
        formic::Matrix<double> _hmat_appro;

        /// \brief [out] harmonic davidson approximate overlap matrix used in SPAM algorithm 
        formic::Matrix<double> _smat_appro;

        /// \brief [out] S^2 matrix 
        formic::Matrix<double> _ssmat;
        std::vector<formic::Matrix<double> > _ssmat_temp;

        /// \brief [out] normal linear method overlap matrix(only being built if the "variance correct" flag is set to be true
        formic::Matrix<double> _lsmat;

        /// \brief [in] the harmonic davidson shift 
        double _hd_shift;

        /// \brief [in] approximate degree
        int _appro_degree;

        /// \brief [in] number of optimizable parameters
        int _num_params;

        /// \brief [in] bare derivative ratios (<n|psi^x> / <n|psi>)
        formic::Matrix<double> & _der_rat;

        /// \brief [in] energy derivative ratio (<n|H|psi^x> / <n|psi>)
        formic::Matrix<double> & _le_der;

        /// \brief [in] S^2 detivative ration (<n|S^2|psi^x> / <n|psi>)
        formic::Matrix<double> & _ls_der;

        /// \brief [out] approximate derivative ratios 
        formic::Matrix<double> _der_rat_appro;

        /// \brief [out] approximate energy derivative ratio
        formic::Matrix<double> _le_der_appro;

        /// \brief [in] list of |value/guiding|^2 history
        const std::vector<double> & _vgs;

        /// \brief [out] average of |value/guiding|^2 value
        double _vgsa;

        /// \brief [in] list of weight history
        const std::vector<double> & _weight;

        /// \brief [out] total weight used in MPI reduce 
        double _total_weight;

        /// \brief flag to tell whether to use spam or not, this determines whether to build approximate matrix or not
        bool _spam_use;

        /// \brief flag to tell whether to do ground state calculation
        bool _ground_state;

        /// \brief flag to tell whether to do variance corrected calculation
        bool _variance_correct;

        /// \brief flag to tell whether to build matrix 
        bool _build_lm_matrix;

        /// \brief flag to tell whether to build S^2 matrix 
        bool _ss_build;

        /// \brief flag to tell whether to print the matrix after built
        bool _print_matrix;

      public:
        
      //////////////////////////////////////////////////////////////////////////////////////////////
      // \brief constructor that initializes the builder and revelent quantities
      //
      //
      //
      //////////////////////////////////////////////////////////////////////////////////////////////
      HamOvlpBuilderHD(formic::Matrix<double> & der_rat, 
                       formic::Matrix<double> & le_der,
                       formic::Matrix<double> & ls_der,
                       const std::vector<double> & vgs,
                       const std::vector<double> & weight,
                       const double hd_shift,
                       const int num_params,
                       const int appro_degree,
                       const bool spam_use,
                       const bool ground_state,
                       const bool variance_correct,
                       const bool build_lm_matrix,
                       const bool ss_build,
                       const bool print_matrix);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function that get parameters
      // 
      /////////////////////////////////////////////////////////////////////////////////////////////
      void get_param(const double hd_shift, 
                     const int num_params,
                     const int appro_degree, 
                     const bool spam_use, 
                     const bool ground_state, 
                     const bool variance_correct, 
                     const bool build_lm_matrix, 
                     const bool ss_build, 
                     const bool print_matrix);

      /////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief constructor used when we do not have derivative vectors stored in memory
      ///
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////
      //HamOvlpBuilderHD(const double hd_shift,
      //                 const bool ground_state,
      //                 const bool ss_build,
      //                 const bool print_matrix);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief build harmonic davidson matrix
      //
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void MatrixBuild(std::ostream & output);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief add contribution to matrix from this sample
      //
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void take_sample(std::vector<double> & der_rat_samp,
                       std::vector<double> & le_der_samp,
                       std::vector<double> & ls_der_samp,
                       double vgs_samp,
                       double weight_samp);

      /////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief finish sample by doing MPI communications 
      /// 
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////
      void finish_sample(const double total_weight);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief absorb |value/guiding|^2 value and weight value into derivative(used when 
      //        _build_lm_matrix flag to set to be false) 
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      double MatrixAbsorb();

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief revover the original derivative vectors(used after the calling of MatrixAbsorb 
      //        function)
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void MatrixRecover();

      /////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief reset this object 
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////
      void reset();
  
      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief  Converts matirx by combining derivatives for dependent variables into derivative
      //         for independent variables
      // \param[in]    deps    object decribing the variable dependencies
      // \param[out]   mat     the matrix to be converted
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void convert_to_ind_var_form(const formic::VarDeps * dep_ptr, formic::Matrix<double> & mat);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to project ground state wavefunction out of first derivative space
      // 
      //  after this function call, ground state wavefunction will be projected out of first derivative  
      //  space via gram-schmidt, thereby changing hamiltonian and overlap matrix 
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void proj_ground_out_fder();

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to calculate the pseudo inverse of overlap matrix 
      // 
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      formic::Matrix<double> ovlp_pseudo_inv(const double threshold, std::ostream & fout);

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to analyze each derivative vectors 
      //
      // 
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void derivative_analyze(std::ostream & fout);

      /// \brief returns total weight
      double total_weight();

      /// \brief returns average of |value/guiding|^2 value 
      double vgsa();

      /// \brief returns the approximate derivative vectors
      formic::Matrix<double> & approximate_der_vec();

      /// \brief returns the approximate energy derivatives
      formic::Matrix<double> & approximate_le_der();

      /// \brief returns the hamiltonian matrix 
      formic::Matrix<double> & ham();

      /// \brief returns the overlap matrix 
      formic::Matrix<double> & ovl();

      /// \brief return the S^2 matrix 
      formic::Matrix<double> & ssquare();

      /// \brief return the nomral linear method overlap matrix 
      formic::Matrix<double> & lovl();

      ///////////////////////////////////////////////////////////////////////////////////
      // \brief do D^(-1/2) transfrom on hamiltonian and overlap matrix 
      //
      //
      //
      ////////////////////////////////////////////////////////////////////////////////////

      void d_inv_half_trans();

    };
  }
}

#endif


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
#include"formic/utils/openmp.h"

#include"formic/utils/zero_one.h"
#include"formic/utils/lmyengine/var_dependencies.h"
#include"formic/utils/matrix.h"
#include"formic/utils/exception.h"
#include"formic/utils/lapack_interface.h"
#include"formic/utils/mpi_interface.h"


namespace cqmc {

  namespace engine {

    template<typename S> class HamOvlpBuilderHD{

      private:

        /// \brief [out] harmonic davidson hamiltonian matrix 
        formic::Matrix<S> _hmat;
        std::vector<formic::Matrix<S> > _hmat_temp;

        /// \brief [out] harmonic davidson overlap matrix 
        formic::Matrix<S> _smat;
        std::vector<formic::Matrix<S> > _smat_temp;

        /// \brief [out] harmonic davidson approximate hamiltonian matrix used in SPAM algorithm 
        formic::Matrix<S> _hmat_appro;

        /// \brief [out] harmonic davidson approximate overlap matrix used in SPAM algorithm 
        formic::Matrix<S> _smat_appro;

        /// \brief [out] S^2 matrix 
        formic::Matrix<S> _ssmat;
        std::vector<formic::Matrix<S> > _ssmat_temp;

        /// \brief [out] normal linear method overlap matrix(only being built if the "variance correct" flag is set to be true
        formic::Matrix<S> _lsmat;

        /// \brief [in] the harmonic davidson shift 
        double _hd_shift;

        /// \brief [in] approximate degree
        int _appro_degree;

        /// \brief [in] number of optimizable parameters
        int _num_params;

        /// \brief [in] bare derivative ratios (<n|psi^x> / <n|psi>)
        formic::Matrix<S> & _der_rat;

        /// \brief [in] energy derivative ratio (<n|H|psi^x> / <n|psi>)
        formic::Matrix<S> & _le_der;

        /// \brief [in] S^2 detivative ration (<n|S^2|psi^x> / <n|psi>)
        formic::Matrix<S> & _ls_der;

        /// \brief [out] approximate derivative ratios 
        formic::Matrix<S> _der_rat_appro;

        /// \brief [out] approximate energy derivative ratio
        formic::Matrix<S> _le_der_appro;

        /// \brief [out] one-body reduced density matrix
        std::vector<formic::Matrix<S> > _one_rdm_temp;
        formic::Matrix<S> _one_rdm;

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

        /// \brief flag to tell whether to compute 1rdm 
        bool _compute_rdm;

      public:
        
      //////////////////////////////////////////////////////////////////////////////////////////////
      // \brief constructor that initializes the builder and revelent quantities
      //
      //
      //
      //////////////////////////////////////////////////////////////////////////////////////////////
      HamOvlpBuilderHD(formic::Matrix<S> & der_rat, 
                       formic::Matrix<S> & le_der,
                       formic::Matrix<S> & ls_der,
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
                       const bool print_matrix)
      :_hd_shift(hd_shift),
      _appro_degree(appro_degree),
      _num_params(num_params),
      _der_rat(der_rat),
      _le_der(le_der),
      _ls_der(ls_der),
      _vgs(vgs),
      _weight(weight),
      _spam_use(spam_use),
      _ground_state(ground_state),
      _variance_correct(variance_correct),
      _build_lm_matrix(build_lm_matrix),
      _ss_build(ss_build),
      _print_matrix(print_matrix)
      {
        
        // number of threads
        int NumThreads = omp_get_max_threads();

        // thread number 
        int myThread = omp_get_thread_num();

        // check if matrix vector size matched the number of threads 
        _hmat_temp.resize(NumThreads);
        _smat_temp.resize(NumThreads);
        _ssmat_temp.resize(NumThreads);
        _one_rdm_temp.resize(NumThreads);

        // size the matrix correctly
        int ndim = _num_params + 1;
        for (int ip = 0; ip < NumThreads; ip++) {
          _hmat_temp[ip].reset(ndim, ndim, formic::zero(S()));
          _smat_temp[ip].reset(ndim, ndim, formic::zero(S()));
          _ssmat_temp[ip].reset(ndim, ndim, formic::zero(S()));
        }

      }

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
                     const bool print_matrix, 
                     const bool compute_rdm=false) {

        // get parameters
        _hd_shift = hd_shift; 
        _num_params = num_params;
        _appro_degree = appro_degree; 
        _spam_use = spam_use; 
        _ground_state = ground_state; 
        _variance_correct = variance_correct; 
        _build_lm_matrix = build_lm_matrix;
        _ss_build = ss_build; 
        _print_matrix = print_matrix;
        _compute_rdm = compute_rdm;

        // number of threads
        int NumThreads = omp_get_max_threads();

        // size the matrix correctly
        int ndim = _num_params + 1;
        for (int ip = 0; ip < NumThreads; ip++) {
          _hmat_temp[ip].reset(ndim, ndim, formic::zero(S()));
          _smat_temp[ip].reset(ndim, ndim, formic::zero(S()));
          _ssmat_temp[ip].reset(ndim, ndim, formic::zero(S()));
        }

      }

    void resetParamNumber(int new_num)
     {
     _num_params = new_num;
     int NumThreads = omp_get_max_threads();
     int my_rank = formic::mpi::rank();

     if(my_rank == 0)
     {    
        std::cout << "Changing size of _hmat_temp and _smat_temp inside matrix_builder" << std::endl;
     }

       // size the matrix correctly
       int ndim = _num_params + 1; 
       for (int ip = 0; ip < NumThreads; ip++) {
         _hmat_temp[ip].reset(ndim, ndim, 0.0);
         _smat_temp[ip].reset(ndim, ndim, 0.0);
         _ssmat_temp[ip].reset(ndim, ndim, 0.0);
       }

     }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief build harmonic davidson matrix
      //
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void MatrixBuild(std::ostream & output)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank(); 
        int num_rank = formic::mpi::size();
        //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
        //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

        // before computing the matrix, first combine bare derivative ratio and local energy derivative with respect to harmonic davidson shift 
        formic::Matrix<S> _hle_der;
        if ( !_ground_state ) {
          _hle_der = _hd_shift * _der_rat - _le_der;
        }

        // get the number of samples 
        const int N = _der_rat.rows();

        // absorb |value/guiding|^2 value and weight value into derivative vectors
        if ( !_ground_state ) {
          for (int i = 0; i < N; i ++) {
            _der_rat.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
            _hle_der.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
          }
        }

        if ( _ground_state ) {
          for (int i = 0; i < N; i ++) {
            _der_rat.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
            _le_der.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));

            // if requested, also absorb |value/guiding|^2 value and weight value into S^2 derivative vectors 
            if ( _ss_build ) 
              _ls_der.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
          }
        }

        // get the size of matrix 
        const int M = _der_rat.cols();

        // turn percentage to integer
        // first check percentage number is reasonable
        //assert( _percentage < 100);

        //const int incre = round(100 / _percentage);

        // initialize zero matrices
        _hmat.reset(M, M, formic::zero(S()));
        _smat.reset(M, M, formic::zero(S()));
        
        // if doing variance corrected calculation, need to build the normal linear method overlap matrix 
        if ( _variance_correct )
          _lsmat.reset(M, M, formic::zero(S()));

        // if requested, also initialize the S^2 matrix 
        if ( _ss_build ) 
          _ssmat.reset(M, M, formic::zero(S()));

        formic::Matrix<S> _hmat_temp(M, M, formic::zero(S()));
        formic::Matrix<S> _smat_temp(M, M, formic::zero(S()));
        formic::Matrix<S> _ssmat_temp(M, M, formic::zero(S()));
        formic::Matrix<S> _lsmat_temp(M, M, formic::zero(S()));

        // call blas level 3 routine to build the matrix
        if ( !_ground_state ) {
          // build the matrix via dgemm
          formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_der_rat.at(0,0), N, &_hle_der.at(0,0), N, 0.0, &_hmat_temp.at(0,0), M);
          formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_hle_der.at(0,0), N, &_hle_der.at(0,0), N, 0.0, &_smat_temp.at(0,0), M);

          // if doing variance correct calculation, build normal linear method's overlap matrix 
          if ( _variance_correct ) 
            formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_der_rat.at(0,0), N, &_der_rat.at(0,0), N, 0.0, &_lsmat_temp.at(0,0), M);

          // undo the changes to derivative vectors 
          for (int i = 0; i < N; i++) {
            _der_rat.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
          }
        }
        
        else {
          // build the matrix via dgemm
          formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_der_rat.at(0,0), N, &_le_der.at(0,0), N, 0.0, &_hmat_temp.at(0,0), M);
          formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_der_rat.at(0,0), N, &_der_rat.at(0,0), N, 0.0, &_smat_temp.at(0,0), M);

          // if requested, also build the S^2 matrix 
          if ( _ss_build )
            formic::xgemm('T', 'N', M, M, N, formic::unity(S()), &_der_rat.at(0,0), N, &_ls_der.at(0,0), N, 0.0, &_ssmat_temp.at(0,0), M);

          // undo the changes to derivative vectors 
          for (int i = 0; i < N; i++) {
            _der_rat.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
            _le_der.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));

            // if requested, also undo S^2 derivative vectors 
            if ( _ss_build ) 
              _ls_der.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
          }
        }

        // finalize the computing process(this may include conmunications between different processes)
        // first compute the total weight and average of |value/guiding|^2
        _total_weight = std::accumulate(_weight.begin(), _weight.end(), 0.0);
        _vgsa = 0.0;

        for (int i = 0; i < _vgs.size(); i++) {
          _vgsa += _weight.at(i) * _vgs.at(i);
        }

        // set up space for mpi reduction 
        const int nred = 2;
        double x[nred];
        double y[nred];

        // prepare for recduce 
        x[0] = _total_weight;
        x[1] = _vgsa;

        // all reduce 
        formic::mpi::allreduce(x, y, nred, MPI_SUM);

        // record results
        _total_weight = y[0];
        _vgsa = y[1] / _total_weight;

        // now do mpi reduce 
        formic::mpi::reduce(&_hmat_temp.at(0,0), &_hmat.at(0,0), _hmat_temp.size(), MPI_SUM);
        formic::mpi::reduce(&_smat_temp.at(0,0), &_smat.at(0,0), _smat_temp.size(), MPI_SUM);

        // if doing variance correct calculation, reduce nomral linear method overlap matrix 
        if ( _variance_correct ) 
          formic::mpi::reduce(&_lsmat_temp.at(0,0), &_lsmat.at(0,0), _lsmat_temp.size(), MPI_SUM);

        // if requested, also reduce S^2 matrix 
        if ( _ss_build ) 
          formic::mpi::reduce(&_ssmat_temp.at(0,0), &_ssmat.at(0,0), _ssmat_temp.size(), MPI_SUM);

        // finally
        _hmat /= (_total_weight * _vgsa);

        if ( _ss_build ) {
          _ssmat /= (_total_weight * _vgsa);
        }

        _smat /= (_total_weight * _vgsa);

        if ( _variance_correct ) 
          _lsmat /= (_total_weight * _vgsa);

        // print the matrix if requested
        if ( _print_matrix && my_rank == 0 ) {

          // hamiltonian 
          output << _hmat.print("%12.6f", "hamiltonian");

          // overlap
          output << _smat.print("%12.6f", "overlap");
        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief add contribution to matrix from this sample
      //
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void take_sample(std::vector<S> & der_rat_samp,
                       std::vector<S> & le_der_samp,
                       std::vector<S> & ls_der_samp,
                       double vgs_samp,
                       double weight_samp)
      {
        
        // dimension of matrix 
        const int ndim = der_rat_samp.size();

        // number of threads
        int NumThreads = omp_get_max_threads();

        // thread number 
        int myThread = omp_get_thread_num();
        //std::cout << boost::format("entering take_sample function in matrix build1") << std::endl;

        // check whether matrices are of the correct size and resize it if not 
        if ( _hmat_temp.at(myThread).rows() != _hmat_temp.at(myThread).cols() || _hmat_temp.at(myThread).rows() != der_rat_samp.size() ) 
          _hmat_temp.at(myThread).reset(ndim, ndim, formic::zero(S()));
        if ( _smat_temp.at(myThread).rows() != _smat_temp.at(myThread).cols() || _smat_temp.at(myThread).rows() != der_rat_samp.size() ) 
          _smat_temp.at(myThread).reset(ndim, ndim, formic::zero(S()));
        if ( _ssmat_temp.at(myThread).rows() != _ssmat_temp.at(myThread).cols() || _ssmat_temp.at(myThread).rows() != der_rat_samp.size() ) 
          _ssmat_temp.at(myThread).reset(ndim, ndim, formic::zero(S()));

        //std::cout << boost::format("entering take_sample function in matrix build3") << std::endl;

        // include the value to guiding square ratio in the weight 
        const double ww = weight_samp * vgs_samp;

        // take the complex conjugation of der_rat
        std::vector<S> der_rat_samp_ct(der_rat_samp.size(), formic::zero(S()));
        for (int i = 0; i < der_rat_samp.size(); i++)
          der_rat_samp_ct.at(i) = formic::conj(der_rat_samp.at(i));

        // for ground state calculations 
        if ( _ground_state ) {
           
          // add contribution to hamiltonian matrix 
          formic::xgemm('N', 'N', ndim, ndim, 1, ww * formic::unity(S()), &der_rat_samp_ct.at(0), ndim, &le_der_samp.at(0), 1, formic::unity(S()), &(_hmat_temp.at(myThread).at(0,0)), ndim);

          // add contribution to overlap matrix 
          formic::xgemm('N', 'N', ndim, ndim, 1, ww * formic::unity(S()), &der_rat_samp_ct.at(0), ndim, &der_rat_samp.at(0), 1, formic::unity(S()), &(_smat_temp.at(myThread).at(0,0)), ndim);

          // add contribution to spin matrix if requested 
          if ( _ss_build ) 
            formic::xgemm('N', 'N', ndim, ndim, 1, ww * formic::unity(S()), &der_rat_samp_ct.at(0), ndim, &ls_der_samp.at(0), 1, formic::unity(S()), &(_ssmat_temp.at(myThread).at(0,0)), ndim);

        }
        
        // for excited state calculations 
        else {
        
          // combine bare and energy derivative ratio with respect to harmonic davidson shift 
          std::vector<S> hle_der_samp(ndim, formic::zero(S()));
          for ( int i = 0; i < ndim; i++) 
            hle_der_samp.at(i) = _hd_shift * der_rat_samp.at(i) - le_der_samp.at(i);

          // take the complex conjugation of hle_der_samp
          std::vector<S> hle_der_samp_ct(ndim, formic::zero(S()));
          for (int i = 0; i < ndim; i++) 
            hle_der_samp_ct.at(i) = formic::conj(hle_der_samp.at(i));

          // add contribution to hailtonian matrix 
          formic::xgemm('N', 'N', ndim, ndim, 1, ww*formic::unity(S()), &der_rat_samp_ct.at(0), ndim, &hle_der_samp.at(0), 1, formic::unity(S()), &(_hmat_temp.at(myThread).at(0,0)), ndim);

          // add contribution to overlap matrix
          formic::xgemm('N', 'N', ndim, ndim, 1, ww*formic::unity(S()), &hle_der_samp_ct.at(0), ndim, &hle_der_samp.at(0), 1, formic::unity(S()), &(_smat_temp.at(myThread).at(0,0)), ndim);

          // add contribution to spin matrix if requested 
          formic::xgemm('N', 'N', ndim, ndim, 1, ww * formic::unity(S()), &der_rat_samp_ct.at(0), ndim, &der_rat_samp.at(0), 1, formic::unity(S()), &(_ssmat_temp.at(myThread).at(0,0)), ndim);

        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief add contribution to matrix from this sample
      //
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void take_sample(int nbasis, std::vector<S> _one_rdm_samp)
      {
        
        // number of threads
        int NumThreads = omp_get_max_threads();

        // thread number 
        int myThread = omp_get_thread_num();
        //std::cout << boost::format("entering take_sample function in matrix build1") << std::endl;

        // check whether matrices are of the correct size and resize it if not 
        if ( _one_rdm_temp.at(myThread).rows() != nbasis || _one_rdm_temp.at(myThread).cols() != nbasis ) 
          _one_rdm_temp.at(myThread).reset(nbasis, nbasis, formic::zero(S()));

        //std::cout << boost::format("entering take_sample function in matrix build3") << std::endl;
        for (int i = 0; i < nbasis; i++)
          for (int j = 0; j < nbasis; j++)
            _one_rdm_temp.at(myThread).at(i, j) += _one_rdm_samp.at(i*nbasis+j);

      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief finish sample by doing MPI communications 
      /// 
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////
      void finish_sample(const double total_weight)
      {

        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank(); 
        int num_rank = formic::mpi::size();

        //std::cout << "netering matrix build finish_sample1" << std::endl;

        // get the number of threads 
        int NumThreads = omp_get_max_threads();

        // get thread number 
        int myThread = omp_get_thread_num();

        // sum over threads
        for (int ip = 1; ip < NumThreads; ip++) {
          _hmat_temp[0] += _hmat_temp[ip];
          _smat_temp[0] += _smat_temp[ip];
          _ssmat_temp[0] += _ssmat_temp[ip];
          if ( _compute_rdm )
            _one_rdm_temp[0] += _one_rdm_temp[ip];
        }

        //std::cout << "netering matrix build finish_sample2" << std::endl;

        // temporary matrices used in mpi reduce
        _hmat.reset(_hmat_temp[0].rows(), _hmat_temp[0].cols(), formic::zero(S()));
        _smat.reset(_smat_temp[0].rows(), _smat_temp[0].cols(), formic::zero(S()));
        _ssmat.reset(_ssmat_temp[0].rows(), _ssmat_temp[0].cols(), formic::zero(S()));
        if ( _compute_rdm )
          _one_rdm.reset(_one_rdm_temp[0].rows(), _one_rdm_temp[0].cols(), formic::zero(S()));
        //std::cout << "netering matrix build finish_sample3" << std::endl;

        // mpi reduce 
        formic::mpi::reduce(&_hmat_temp[0].at(0,0), &_hmat.at(0,0), _hmat_temp[0].size(), MPI_SUM);
        formic::mpi::reduce(&_smat_temp[0].at(0,0), &_smat.at(0,0), _smat_temp[0].size(), MPI_SUM);
        formic::mpi::reduce(&_ssmat_temp[0].at(0,0), &_ssmat.at(0,0), _ssmat_temp[0].size(), MPI_SUM);   
        if ( _compute_rdm )
          formic::mpi::reduce(&_one_rdm_temp[0].at(0,0), &_one_rdm.at(0,0), _one_rdm_temp[0].size(), MPI_SUM);

        //std::cout << "netering matrix build finish_sample4" << std::endl;
       
        // compute the average 
        _hmat /= total_weight; 
        _smat /= total_weight;
        if ( _compute_rdm )
          _one_rdm /= total_weight;
        //std::cout << "total weight is " << total_weight << std::endl;

        if ( !_ground_state ) 
          _ssmat /= total_weight;
       
        //std::cout << "netering matrix build finish_sample4.5" << std::endl;
        // clear temporary matrices 
        for (int ip = 0; ip < NumThreads; ip++) {
          _hmat_temp[ip].reset(_hmat.rows(), _hmat.cols(), formic::zero(S()));
          _smat_temp[ip].reset(_smat.rows(), _smat.cols(), formic::zero(S()));
          _ssmat_temp[ip].reset(_ssmat.rows(), _ssmat.cols(), formic::zero(S()));
          if ( _compute_rdm )
            _one_rdm_temp[ip].reset(_one_rdm.rows(), _one_rdm.cols(), formic::zero(S()));
        }

        //std::cout << "netering matrix build finish_sample5" << std::endl;

        // print the matrix if requested
        if ( _print_matrix && my_rank == 0 ) {

          // hamiltonian 
          //output << _hmat.print("%12.6f", "hamiltonian");

          // overlap
          //output << _smat.print("%12.6f", "overlap");
        }

        if (my_rank == 0) {
          
          // discard the imaginary part of the diagonal elements of Hamiltonian matrix 
          //for (int i=0; i < _hmat.rows(); i++) {
          //  _hmat.at(i,i) = formic::convert(_hmat.at(i,i));
          //}

          // hamiltonian 
          std::cout << _hmat.print("%12.6f", "hamiltonian");

          // overlap
          std::cout << _smat.print("%12.6f", "overlap");

          std::cout << _ssmat.print("%12.6f", "LM overlap");

          if ( _compute_rdm )
            std::cout << _one_rdm.print("%12.6f", "one rdm");

        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief absorb |value/guiding|^2 value and weight value into derivative(used when 
      //        _build_lm_matrix flag to set to be false) 
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      double MatrixAbsorb()
      {
        // get the number of total samples on each process 
        const int Ns = _der_rat.rows();

        // get the number of samples that will be used to build approximate matrix 
        const int Ns_appro = Ns / _appro_degree;

        // absorb |value/guiding|^2 value and weight value into derivative vectors
        for (int i = 0; i < Ns; i++) {
          _der_rat.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
          _le_der.scale_row_by(i, std::sqrt(_vgs.at(i) * _weight.at(i)));
        }

        // evaluate average of |value/guiding|^2 value and total weight 
        _total_weight = std::accumulate(_weight.begin(), _weight.end(), 0.0);
        double _total_weight_appro = 0.0;
        for (int i = 0; i < Ns_appro; i++) 
          _total_weight_appro += _weight.at(i);

        _vgsa = 0.0;
        double _vgsa_appro = 0.0;

        for (int i = 0; i < _vgs.size(); i++) {
          _vgsa += _weight.at(i) * _vgs.at(i);
          if ( i < Ns_appro) 
            _vgsa_appro += _weight.at(i) * _vgs.at(i);
        }

        // set up space for mpi reduction 
        const int nred = 2;
        double x[nred];
        double y[nred];
        double x_appro[nred];
        double y_appro[nred];

        // prepare for reduce 
        x[0] = _total_weight;
        x[1] = _vgsa;
        x_appro[0] = _total_weight_appro;
        x_appro[1] = _vgsa_appro;

        // all reduce 
        formic::mpi::allreduce(x, y, nred, MPI_SUM);
        formic::mpi::allreduce(x_appro, y_appro, nred, MPI_SUM);

        // record results
        _total_weight = y[0];
        _total_weight_appro = y_appro[0];
        _vgsa = y[1] / _total_weight;
        _vgsa_appro = y_appro[1] / _total_weight_appro;

        _der_rat /= std::sqrt(_total_weight * _vgsa);
        _le_der /= std::sqrt(_total_weight * _vgsa);

        // if we are using SPAM solver, we create approximate derivative vectors
        if ( _spam_use ) {
          _der_rat_appro.reset(Ns_appro, _der_rat.cols());
          _le_der_appro.reset(Ns_appro, _le_der.cols());
          std::vector<int> row_inds, col_inds;
          for (int i = 0; i < _der_rat.cols(); i++) {
            col_inds.push_back(i);
          }
          for (int j = 0; j < Ns_appro; j++) {
            row_inds.push_back(j);
          }
          _der_rat_appro = _der_rat.slice(row_inds, col_inds);
          _le_der_appro = _le_der.slice(row_inds, col_inds);
        }

        return (_total_weight * _vgsa) / (_total_weight_appro * _vgsa_appro);

      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief revover the original derivative vectors(used after the calling of MatrixAbsorb 
      //        function)
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void MatrixRecover()
      {
        // get the number of samples 
        const int Ns = _der_rat.rows();

        // divide |value/guiding|^2 and weight value out of derivative vectors 
        for (int i = 0; i < Ns; i++) {
          _der_rat.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
          _le_der.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
        }

        // consider |value/guiding|^2 function 
        _der_rat *= std::sqrt(_total_weight * _vgsa);
        _le_der *= std::sqrt(_total_weight * _vgsa);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief reset this object 
      ///
      /////////////////////////////////////////////////////////////////////////////////////////////
      void reset()
      {
        // clear matrices
        _hmat.reset(0, 0, formic::zero(S()));
        _smat.reset(0, 0, formic::zero(S()));
        _ssmat.reset(0, 0, formic::zero(S()));
        _one_rdm.reset(0, 0, formic::zero(S()));
      }
  
      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief  Converts matirx by combining derivatives for dependent variables into derivative
      //         for independent variables
      // \param[in]    deps    object describing the variable dependencies
      // \param[out]   mat     the matrix to be converted
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void convert_to_ind_var_form(const formic::VarDeps * dep_ptr, formic::Matrix<S> & mat)
      {
        // get expected initial matrix dimension
        const int n = 1 + dep_ptr -> n_tot();

        // ensure the matirx is the correct size
        if ( mat.rows() != n || mat.cols() != n )
          throw formic::Exception("When converting to independent variable matrix, function get a matrix of size %i by %i, but expected one of size %i by %i") % mat.rows() % mat.cols() % n % n;

        // get the final matrix dimension
        const int m = 1 + dep_ptr -> n_ind();

        // combine derivatives down columns to get a m by n matrix
        formic::Matrix<S> column_contracted(m, n, 0.0);
        for (int i = 0; i < n; i++) {
          column_contracted.at(0,i) = mat.at(0,i);
          dep_ptr -> compute_ind_derivs(&mat.at(1,i), &column_contracted.at(1,i));
        }

        // transpose this matrix and then combine derivatives down columns to get the final m by m matrix
        column_contracted.tip();
        formic::Matrix<S> final_mat(m, m, zero(S()));
        for (int i = 0; i < m; i++) {
          final_mat.at(0,i) = column_contracted.at(0,i);
          dep_ptr -> compute_ind_derivs(&column_contracted.at(1,i), &final_mat.at(1,i));
        }
        
        // set the input matrix to the final matrix
        mat = final_mat.t();

      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to project ground state wavefunction out of first derivative space
      // 
      //  after this function call, ground state wavefunction will be projected out of first derivative  
      //  space via gram-schmidt, thereby changing hamiltonian and overlap matrix 
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      void proj_ground_out_fder()
      {
        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();
        //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
        //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

        if ( my_rank == 0 ) {
          // get matrix dimension 
          const int m = _hmat.rows();
          const int n = _hmat.cols();

          // get temp matrix 
          formic::Matrix<S> temp_h(m, n);
          formic::Matrix<S> temp_s(m, n);
          formic::Matrix<S> temp_ss(m, n);
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
              if ( !( i == 0 && j == 0) ) {
          temp_s.at(i,j) = _smat.at(i,j) - (_smat.at(j,0) * _smat.at(i,0) / _smat.at(0,0) + _smat.at(i,0) * _smat.at(0,j) / _smat.at(0,0) - _smat.at(0,i) * _smat.at(j,0) / _smat.at(0,0));
                temp_h.at(i,j) = _hmat.at(i,j) - (_smat.at(j,0) * _hmat.at(i,0) / _smat.at(0,0) + _smat.at(0,i) * _hmat.at(0,j) / _smat.at(0,0) - _smat.at(0, i) * _smat.at(j, 0) * _hmat.at(0, 0) / (_smat.at(0, 0) * _smat.at(0, 0))); 
          if ( _ss_build ) 
            temp_ss.at(i, j) = _ssmat.at(i, j) - (_smat.at(j, 0) * _ssmat.at(i, 0) / _smat.at(0, 0) + _smat.at(i, 0) * _ssmat.at(0, j) / _smat.at(0, 0) - _smat.at(0, i) * _smat.at(j, 0) * _ssmat.at(0, 0) / (_smat.at(0, 0) * _smat.at(0, 0)));
        }

        else {
          temp_s.at(i, j) = formic::unity(S());
          temp_h.at(i, j) = _hmat.at(i, j);
          if ( _ss_build ) 
            temp_ss.at(i, j) = _ssmat.at(i, j);
        }
            }
          }

          _hmat = temp_h;
          _smat = temp_s;
          if ( _ss_build ) 
            _ssmat = temp_ss;

        }
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to calculate the pseudo inverse of overlap matrix 
      // 
      //
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      formic::Matrix<S> ovlp_pseudo_inv(const double threshold, std::ostream & fout)
      {
        
        // get rank number and number of ranks 
        int my_rank = formic::mpi::rank();
        int num_rank = formic::mpi::size();
        //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
        //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
       
        // compute SVD on root process 
        formic::Matrix<S> u;
        formic::Matrix<S> vt;
        formic::ColVec<S> sin_vals;
        if (my_rank == 0) {
          _smat.svd(u, sin_vals, vt);

          // print out overlap matrix's singular values 
          fout << boost::format("transformed overlap matrix's singular values are") << std::endl;
          for (int i = 0; i < sin_vals.size(); i++) {
            fout << boost::format("%20.8e") % sin_vals.at(i) << std::endl;
          }
          fout << std::endl;

          // compute and return inverse 
          for (int i = 0; i < vt.rows(); i++) {
            double scaler = (formic::real(sin_vals(i)) > threshold ? formic::unity(S()) / sin_vals(i) : formic::zero(S()));
            vt.scale_row_by(i, scaler);
          }

        }
        return vt.c() * u.c();
      }

      /////////////////////////////////////////////////////////////////////////////////////////////
      // \brief function to analyze each derivative vectors 
      //
      // 
      //
      /////////////////////////////////////////////////////////////////////////////////////////////
      //void derivative_analyze(std::ostream & fout);

      /// \brief returns total weight
      double total_weight() { return _total_weight; }

      /// \brief returns average of |value/guiding|^2 value 
      double vgsa() { return _vgsa; }

      /// \brief returns the approximate derivative vectors
      formic::Matrix<S> & approximate_der_vec() { return _der_rat_appro; }

      /// \brief returns the approximate energy derivatives
      formic::Matrix<S> & approximate_le_der() { return _le_der_appro; }

      /// \brief returns the hamiltonian matrix 
      formic::Matrix<S> & ham() { return _hmat; }

      /// \brief returns the overlap matrix 
      formic::Matrix<S> & ovl() { return _smat; }

      /// \brief return the S^2 matrix 
      formic::Matrix<S> & ssquare() { return _ssmat; }

      /// \brief return the nomral linear method overlap matrix 
      formic::Matrix<S> & lovl() { return _lsmat; }

      ///////////////////////////////////////////////////////////////////////////////////
      // \brief do D^(-1/2) transform on hamiltonian and overlap matrix 
      //
      //
      //
      ////////////////////////////////////////////////////////////////////////////////////

      //void d_inv_half_trans();

    };
  }
}

#endif


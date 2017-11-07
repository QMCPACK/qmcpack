////////////////////////////////////////////////////////////////////////////////////////////////
// \brief matrix build class that builds harmonic davidson matrix
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////

#include<complex>
#include<vector>
#include<numeric>
#include<cassert>
#include<algorithm>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/openmp.h>

#include<formic/utils/exception.h>
#include<formic/utils/lapack_interface.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lmyengine/matrix_builder.h>


//////////////////////////////////////////////////////////////////////////////////////////////
// \brief constructor that initializes the builder and revelent quantities
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
cqmc::engine::HamOvlpBuilderHD::HamOvlpBuilderHD(formic::Matrix<double> & der_rat, 
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
                                                 const bool print_matrix)
:_der_rat(der_rat),
_le_der(le_der),
_ls_der(ls_der),
_vgs(vgs),
_weight(weight),
_hd_shift(hd_shift),
_num_params(num_params),
_appro_degree(appro_degree),
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
  if ( _ss_build ) {
    _ssmat_temp.resize(NumThreads);
  }

  // size the matrix correctly
  int ndim = _num_params + 1;
  for (int ip = 0; ip < NumThreads; ip++) {
    _hmat_temp[ip].reset(ndim, ndim, 0.0);
    _smat_temp[ip].reset(ndim, ndim, 0.0);
    if ( _ss_build ) 
      _ssmat_temp[ip].reset(ndim, ndim, 0.0);
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////
// \brief function that get parameters
// 
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::get_param(const double hd_shift, 
                                               const int num_params,
                                               const int appro_degree, 
                                               const bool spam_use, 
                                               const bool ground_state, 
                                               const bool variance_correct, 
                                               const bool build_lm_matrix, 
                                               const bool ss_build, 
                                               const bool print_matrix) {
  
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

  // number of threads
  int NumThreads = omp_get_max_threads();

  // size the matrix correctly
  int ndim = _num_params + 1;
  for (int ip = 0; ip < NumThreads; ip++) {
    _hmat_temp[ip].reset(ndim, ndim, 0.0);
    _smat_temp[ip].reset(ndim, ndim, 0.0);
    if ( _ss_build ) 
      _ssmat_temp[ip].reset(ndim, ndim, 0.0);
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////
/// \brief constructor used when we do not have derivative vectors stored in memory
///
///
/////////////////////////////////////////////////////////////////////////////////////////////
//cqmc::engine::HamOvlpBuilderHD::HamOvlpBuilderHD(const double hd_shift,
//                                                 const bool ground_state,
//                                                 const bool ss_build,
//                                                 const bool print_matrix) 
//:_hd_shift(hd_shift),
//_ground_state(ground_state),
//_ss_build(ss_build),
//_print_matrix(print_matrix)
//{
//  
//  // intialize derivative vectors
//  formic::Matrix<double> der_rat(1, 1, 0.0) = _der_rat;
//  formic::Matrix<double> le_der(1, 1, 0.0)  = _le_der;
//  formic::Matrix<double> ls_der(1, 1, 0.0)  = _ls_der;
//  std::vector<double> vgs(1, 1.0) = _vgs;
//  std::vector<double> weight(1, 1.0) = _weight;
//}

/////////////////////////////////////////////////////////////////////////////////////////////
// \brief add contribution to matrix from this sample
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::take_sample(std::vector<double> & der_rat_samp,
                                                 std::vector<double> & le_der_samp,
                                                 std::vector<double> & ls_der_samp,
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
    _hmat_temp.at(myThread).reset(ndim, ndim, 0.0);
  if ( _smat_temp.at(myThread).rows() != _smat_temp.at(myThread).cols() || _smat_temp.at(myThread).rows() != der_rat_samp.size() ) 
    _smat_temp.at(myThread).reset(ndim, ndim, 0.0);
  if ( _ss_build ) {
    if ( _ssmat_temp.at(myThread).rows() != _ssmat_temp.at(myThread).cols() || _ssmat_temp.at(myThread).rows() != der_rat_samp.size() ) 
      _ssmat_temp.at(myThread).reset(ndim, ndim, 0.0);
  }

  //std::cout << boost::format("entering take_sample function in matrix build3") << std::endl;

  // include the value to guiding square ratio in the weight 
  const double ww = weight_samp * vgs_samp;

  // for ground state calculations 
  if ( _ground_state ) {
     
    // add contribution to hamiltonian matrix 
    formic::dgemm('N', 'N', ndim, ndim, 1, ww, &der_rat_samp.at(0), ndim, &le_der_samp.at(0), 1, 1.0, &(_hmat_temp.at(myThread).at(0,0)), ndim);

    // add contribution to overlap matrix 
    formic::dgemm('N', 'N', ndim, ndim, 1, ww, &der_rat_samp.at(0), ndim, &der_rat_samp.at(0), 1, 1.0, &(_smat_temp.at(myThread).at(0,0)), ndim);

    // add contribution to spin matrix if requested 
    if ( _ss_build ) 
      formic::dgemm('N', 'N', ndim, ndim, 1, ww, &der_rat_samp.at(0), ndim, &ls_der_samp.at(0), 1, 1.0, &(_ssmat_temp.at(myThread).at(0,0)), ndim);

  }
  
  // for excited state calculations 
  else {
  
    // combine bare and energy derivative ratio with respect to harmonic davidson shift 
    std::vector<double> hle_der_samp(ndim, 0.0);
    for ( int i = 0; i < ndim; i++) 
      hle_der_samp.at(i) = _hd_shift * der_rat_samp.at(i) - le_der_samp.at(i);

    // add contribution to hailtonian matrix 
    formic::dgemm('N', 'N', ndim, ndim, 1, ww, &der_rat_samp.at(0), ndim, &hle_der_samp.at(0), 1, 1.0, &(_hmat_temp.at(myThread).at(0,0)), ndim);

    // add contribution to overlap matrix
    formic::dgemm('N', 'N', ndim, ndim, 1, ww, &hle_der_samp.at(0), ndim, &hle_der_samp.at(0), 1, 1.0, &(_smat_temp.at(myThread).at(0,0)), ndim);

  }
}

/////////////////////////////////////////////////////////////////////////////////////////////
/// \brief finish sample by doing MPI communications 
/// 
///
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::finish_sample(const double total_weight) 
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
    if (_ss_build)
      _ssmat_temp[0] += _ssmat_temp[ip];
  }

  //std::cout << "netering matrix build finish_sample2" << std::endl;

  // temporary matrices used in mpi reduce
  _hmat.reset(_hmat_temp[0].rows(), _hmat_temp[0].cols(), 0.0);
  _smat.reset(_smat_temp[0].rows(), _smat_temp[0].cols(), 0.0);
  if ( _ss_build ) 
    _ssmat.reset(_ssmat_temp[0].rows(), _ssmat_temp[0].cols(), 0.0);
  //std::cout << "netering matrix build finish_sample3" << std::endl;

  // mpi reduce 
  formic::mpi::reduce(&_hmat_temp[0].at(0,0), &_hmat.at(0,0), _hmat_temp[0].size(), MPI_SUM);
  formic::mpi::reduce(&_smat_temp[0].at(0,0), &_smat.at(0,0), _smat_temp[0].size(), MPI_SUM);
  if ( _ss_build ) 
    formic::mpi::reduce(&_ssmat_temp[0].at(0,0), &_ssmat.at(0,0), _ssmat_temp[0].size(), MPI_SUM);

  //std::cout << "netering matrix build finish_sample4" << std::endl;
 
  // compute the average 
  _hmat /= total_weight; 
  _smat /= total_weight;
  //std::cout << "total weight is " << total_weight << std::endl;

  if ( _ss_build ) 
    _ssmat /= total_weight;
 
  //std::cout << "netering matrix build finish_sample4.5" << std::endl;
  // clear temporary matrices 
  for (int ip = 0; ip < NumThreads; ip++) {
    _hmat_temp[ip].reset(_hmat.rows(), _hmat.cols(), 0.0);
    _smat_temp[ip].reset(_smat.rows(), _smat.cols(), 0.0);
      if ( _ss_build ) 
       _ssmat_temp[ip].reset(_ssmat.rows(), _ssmat.cols(), 0.0);
  }

  //std::cout << "netering matrix build finish_sample5" << std::endl;

  // print the matrix if requested
  if ( _print_matrix && my_rank == 0 ) {

    // hamiltonian 
    //output << _hmat.print("%12.6f", "hamiltonian");

    // overlap
    //output << _smat.print("%12.6f", "overlap");
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////
// \brief build harmonic davidson matrix
// 
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::MatrixBuild(std::ostream & output)
{

  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // before computing the matrix, first combine bare derivative ratio and local energy derivative with respect to harmonic davidson shift 
  formic::Matrix<double> _hle_der;
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
  _hmat.reset(M, M, 0.0);
  _smat.reset(M, M, 0.0);
  
  // if doing variance corrected calculation, need to build the normal linear method overlap matrix 
  if ( _variance_correct )
    _lsmat.reset(M, M, 0.0);

  // if requested, also initialize the S^2 matrix 
  if ( _ss_build ) 
    _ssmat.reset(M, M, 0.0);

  formic::Matrix<double> _hmat_temp(M, M, 0.0);
  formic::Matrix<double> _smat_temp(M, M, 0.0);
  formic::Matrix<double> _ssmat_temp(M, M, 0.0);
  formic::Matrix<double> _lsmat_temp(M, M, 0.0);

  // call blas level 3 routine to build the matrix
  if ( !_ground_state ) {
    // build the matrix via dgemm
    formic::dgemm('T', 'N', M, M, N, 1.0, &_der_rat.at(0,0), N, &_hle_der.at(0,0), N, 0.0, &_hmat_temp.at(0,0), M);
    formic::dgemm('T', 'N', M, M, N, 1.0, &_hle_der.at(0,0), N, &_hle_der.at(0,0), N, 0.0, &_smat_temp.at(0,0), M);

    // if doing variance correct calculation, build normal linear method's overlap matrix 
    if ( _variance_correct ) 
      formic::dgemm('T', 'N', M, M, N, 1.0, &_der_rat.at(0,0), N, &_der_rat.at(0,0), N, 0.0, &_lsmat_temp.at(0,0), M);

    // undo the changes to derivative vectors 
    for (int i = 0; i < N; i++) {
      _der_rat.scale_row_by(i, 1.0/std::sqrt(_vgs.at(i) * _weight.at(i)));
    }
  }
  
  else {
    // build the matrix via dgemm
    formic::dgemm('T', 'N', M, M, N, 1.0, &_der_rat.at(0,0), N, &_le_der.at(0,0), N, 0.0, &_hmat_temp.at(0,0), M);
    formic::dgemm('T', 'N', M, M, N, 1.0, &_der_rat.at(0,0), N, &_der_rat.at(0,0), N, 0.0, &_smat_temp.at(0,0), M);

    // if requested, also build the S^2 matrix 
    if ( _ss_build )
      formic::dgemm('T', 'N', M, M, N, 1.0, &_der_rat.at(0,0), N, &_ls_der.at(0,0), N, 0.0, &_ssmat_temp.at(0,0), M);

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
// \brief absorb |value/guiding|^2 value and weight value into derivative(used when 
//        _build_lm_matrix flag to set to be false 
//
// NOTE:After calling this function, derivative vectors have been changed
/////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::HamOvlpBuilderHD::MatrixAbsorb()
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
void cqmc::engine::HamOvlpBuilderHD::MatrixRecover()
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
void cqmc::engine::HamOvlpBuilderHD::reset() 
{
  // clear matrices
  _hmat.reset(0, 0, 0.0);
  _smat.reset(0, 0, 0.0);
  _ssmat.reset(0, 0, 0.0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  Converts the input matrix by combining derivatives for dependent variables into 
//         derivative for independent variables
//
// \param[in]         deps     object decribing the variable dependencies
// \param[in, out]    mat      the matrix to be converted
////////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::convert_to_ind_var_form(const formic::VarDeps * dep_ptr, formic::Matrix<double> & mat)
{
  // get expected initial matrix dimension
  const int n = 1 + dep_ptr -> n_tot();

  // ensure the matirx is the correct size
  if ( mat.rows() != n || mat.cols() != n )
    throw formic::Exception("When converting to independent variable matrix, function get a matrix of size %i by %i, but expected one of size %i by %i") % mat.rows() % mat.cols() % n % n;

  // get the final matrix dimension
  const int m = 1 + dep_ptr -> n_ind();

  // combine derivatives down columns to get a m by n matrix
  formic::Matrix<double> column_contracted(m, n, 0.0);
  for (int i = 0; i < n; i++) {
    column_contracted.at(0,i) = mat.at(0,i);
    dep_ptr -> compute_ind_derivs(&mat.at(1,i), &column_contracted.at(1,i));
  }

  // transpose this matrix and then combine derivatives down columns to get the final m by m matrix
  column_contracted.tip();
  formic::Matrix<double> final_mat(m, m, 0.0);
  for (int i = 0; i < m; i++) {
    final_mat.at(0,i) = column_contracted.at(0,i);
    dep_ptr -> compute_ind_derivs(&column_contracted.at(1,i), &final_mat.at(1,i));
  }
  
  // set the input matrix to the final matrix
  mat = final_mat.t();

}

///////////////////////////////////////////////////////////////////////////////////
// \brief do D^(-1/2) transfrom on hamiltonian and overlap matrix 
//
//
//
////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::d_inv_half_trans()
{  
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  if (my_rank == 0) {
    // first get D^(-1/2) matrix 
    formic::Matrix<double> d_inv_half(_hmat.rows(), _hmat.cols());
    for (int i = 0; i < _hmat.rows(); i++) {
      for (int j = 0; j < _hmat.cols(); j++) {
        if (i != j) 
          d_inv_half.at(i, j) = 0.0;
        else
          d_inv_half.at(i, j) = (_smat.at(i, i) > 1.0e-5 ? (1.0 / std::sqrt(_smat.at(i, i))) : (1.0 / std::sqrt(1.0e-5)));
      }  
    }

    // then do the transformation 
    _hmat = d_inv_half * _hmat * d_inv_half;
    _smat = d_inv_half * _smat * d_inv_half;
  }
}

////////////////////////////////////////////////////////////////////////////////////
// \brief function to project ground state wavefunction out of first derivative space
//
//
//
////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::proj_ground_out_fder()
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
    formic::Matrix<double> temp_h(m, n);
    formic::Matrix<double> temp_s(m, n);
    formic::Matrix<double> temp_ss(m, n);
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        if ( !( i == 0 && j == 0) ) {
	  temp_s.at(i,j) = _smat.at(i,j) - (_smat.at(j,0) * _smat.at(i,0) / _smat.at(0,0) + _smat.at(i,0) * _smat.at(0,j) / _smat.at(0,0) - _smat.at(0,i) * _smat.at(j,0) / _smat.at(0,0));
          temp_h.at(i,j) = _hmat.at(i,j) - (_smat.at(j,0) * _hmat.at(i,0) / _smat.at(0,0) + _smat.at(0,i) * _hmat.at(0,j) / _smat.at(0,0) - _smat.at(0, i) * _smat.at(j, 0) * _hmat.at(0, 0) / (_smat.at(0, 0) * _smat.at(0, 0))); 
	  if ( _ss_build ) 
	    temp_ss.at(i, j) = _ssmat.at(i, j) - (_smat.at(j, 0) * _ssmat.at(i, 0) / _smat.at(0, 0) + _smat.at(i, 0) * _ssmat.at(0, j) / _smat.at(0, 0) - _smat.at(0, i) * _smat.at(j, 0) * _ssmat.at(0, 0) / (_smat.at(0, 0) * _smat.at(0, 0)));
	}

	else {
	  temp_s.at(i, j) = 1.0;
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
 
 

////////////////////////////////////////////////////////////////////////////////////
// \brief calculate the pseudo inverse of overlap matrix 
// 
//
//
//
////////////////////////////////////////////////////////////////////////////////////
formic::Matrix<double> cqmc::engine::HamOvlpBuilderHD::ovlp_pseudo_inv(const double threshold, std::ostream & fout) 
{
  
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
 
  // compute SVD on root process 
  formic::Matrix<double> u;
  formic::Matrix<double> vt;
  formic::ColVec<double> sin_vals;
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
      double scaler = (sin_vals(i) > threshold ? 1.0 / sin_vals(i) : 0.0);
      vt.scale_row_by(i, scaler);
    }

  }
  return vt.t() * u.t();
}

/////////////////////////////////////////////////////////////////////////////////////////////
// \brief function to analyze each derivative vectors 
//
// 
//
/////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::HamOvlpBuilderHD::derivative_analyze(std::ostream & fout)
{
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // analyze derivative vectors on root process
  if ( my_rank == 0 ) {
    
    // get a vector to store the square root of weight 
    std::vector<double> _sqrt_weight;

    // print out weight of each configuration
    fout << boost::format("weight on each configuration is") << std::endl;
    for (int i = 0; i < _weight.size(); i++) {
      _sqrt_weight.push_back( std::sqrt(_weight.at(i)) );
      fout << boost::format("%20.8e") % _sqrt_weight.at(i) << std::endl;
    }

    // for each row of derivative vector matrix, multiply it by weight to get the matrix <n|Psi^x>
    for (int i = 0; i < _der_rat.rows(); i++) {
      _der_rat.scale_row_by(i, _sqrt_weight.at(i));
    }

    // project ground state out of first derivative space 
    for (int i = 1; i < _der_rat.cols(); i++) {
      for (int j = 0; j < _der_rat.rows(); j++) {
        _der_rat.at(j, i) -= _smat.at(i, 0) * _der_rat.at(j, 0);
      }
    }
    this -> proj_ground_out_fder();

    // normalize each derivative vector
    for (int i = 0; i < _der_rat.cols(); i++) {
      _der_rat.scale_col_by(i, 1.0/std::sqrt( _smat(i, i) ));
    }

    // print out the first column
    fout << boost::format("new derivative vector matrix is") << std::endl;
    for (int i = 0; i < _der_rat.cols(); i++) {
      for (int j = 0; j < _der_rat.rows(); j++) {
        fout << boost::format("%20.8e") % _der_rat.at(j, i) << std::endl;
      }
      fout << std::endl;
    }
  }
}

/// \brief returns total weight 
double cqmc::engine::HamOvlpBuilderHD::total_weight() { return _total_weight; }

/// \brief returns average of |value/guiding|^2 value 
double cqmc::engine::HamOvlpBuilderHD::vgsa() { return _vgsa; }

/// \brief returns the approximate derivative vectors
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::approximate_der_vec() { return _der_rat_appro; }

/// \brief returns the approximate energy derivatives
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::approximate_le_der() { return _le_der_appro; }

/// \brief returns the hamiltonian matrix 
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::ham() { return _hmat; }

/// \brief returns the overlap matrix 
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::ovl() { return _smat; }

/// \brief returns the S^2 matrix 
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::ssquare() { return _ssmat; }

/// \brief returns the normal linear method's overlap matrix 
formic::Matrix<double> & cqmc::engine::HamOvlpBuilderHD::lovl() { return _lsmat; }


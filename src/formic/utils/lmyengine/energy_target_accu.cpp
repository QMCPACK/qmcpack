///////////////////////////////////////////////////////////////////////////////////////////////
// \file energy_target_accu.cpp
//
//
// \brief  implementation file for harmonic davidson eigine energy and target function calculate
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include<complex>
#include<vector>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
//#include<mpi.h>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/exception.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lmyengine/engine_numeric.h>
#include<formic/utils/lmyengine/energy_target_accu.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  constructs the computer
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
cqmc::engine::ETCompute::ETCompute(const std::vector<double> & le_history,
                                   const std::vector<double> & vg_history,
                                   const std::vector<double> & w_history,
                                   const bool exact_sampling,
                                   const bool ground_state,
                                   const bool variance_correct,
                                   const double hd_shift,
                                   const double var_weight)
: _le_history(le_history),
  _vg_history(vg_history),
  _w_history(w_history),
  _exact_sampling(exact_sampling),
  _ground_state(ground_state),
  _variance_correct(variance_correct),
  _hd_shift(hd_shift),
  _var_weight(var_weight),
  _energy(0.0),
  _energy_s(0.0),
  _target_fn_val(0.0),
  _variance(0.0),
  _variance_s(0.0),
  _serr(-1.0),
  _tnserr(-1.0),
  _tdserr(-1.0),
  _tserr(-1.0)
{
  // initialize all list based on the input vector 
  this -> history_initialize();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief function that initializes all history list
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::ETCompute::history_initialize()
{
      
  // set constant number 
  const int N = _le_history.size();

  // check to see if input list have the same length
  assert(N == _vg_history.size());
  assert(N == _w_history.size());

  // initialize all history list
  for (int i = 0; i < N; i ++) {
    double current_eng = _le_history.at(i);
    //std::cout << boost::format("%10.8e ") % current_eng;
    double current_engv = _le_history.at(i) * _vg_history.at(i);
    double current_eng_s = current_eng * current_eng;
    double current_eng_sv = current_eng_s * _vg_history.at(i);

    if (!_ground_state) {
      double current_tn = _hd_shift - current_eng;
      double current_tnv = current_tn * _vg_history.at(i);
      double current_td = _hd_shift * _hd_shift - 2 * _hd_shift * current_eng + current_eng_s;
      double current_tdv = current_td * _vg_history.at(i);
      _tn_history.push_back(current_tn);
      _tnv_history.push_back(current_tnv);
      _td_history.push_back(current_td);
      _tdv_history.push_back(current_tdv);

    }

    _lev_history.push_back(current_engv);
    _les_history.push_back(current_eng_s);
    _lesv_history.push_back(current_eng_sv);
  }
  //std::cout << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  returns the variance of the block average for the specific block length
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::ETCompute::bvar(const int nblocks) const 
{
  // check for sane number of blocks
  assert( nblocks > 0);

  // get the number of ranks 
  int rank_num = formic::mpi::size();
  //MPI_Comm_size(MPI_COMM_WORLD, &rank_num);

  // get block length (last block will be longer if division does not result in an integer)
  const int bl = _le_history.size() / nblocks;

  // compute averages of each block on each process 
  std::vector<double> avgs;
  avgs.assign(nblocks, 0.0);
  for (int i = 0; i < nblocks; i++) {
    const int start = i * bl;
    const int end = ( i == nblocks-1 ? _le_history.size() : (i+1) * bl );
    avgs.at(i) = std::accumulate( _le_history.begin() + start, _le_history.begin() + end, 0.0) / double(end - start);
  }

  // compute avergae of each block across all processes
  std::vector<double> full_avgs; 
  full_avgs.assign(nblocks, 0.0);
  formic::mpi::reduce(&avgs.at(0), &full_avgs.at(0), nblocks, MPI_SUM);
  for (int i = 0; i < nblocks; i++) {
    full_avgs.at(i) /= double(rank_num);
  }

  // compute overall average
  const double avg = std::accumulate(full_avgs.begin(), full_avgs.end(), 0.0) / double(nblocks);

  // compute the variance of the different blocks' averages
  double var = 0.0;
  for ( int i = 0; i < nblocks; i++) {
    const double x = avg - full_avgs.at(i);
    var += x * x;
  }
  var /= double(nblocks);

  // return the variance on all processes 
  formic::mpi::bcast(&var, 1);
  return var;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief performs a recursive blocking analysis of the statistical energy error
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::ETCompute::recursive_blocking(std::ostream & fout, const bool print) 
{
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // make sure the sample was not exact sample
  if ( _exact_sampling )
    throw formic::Exception("recursive blocking error analysis should not be performed for exact sampling");

  // get number of processors
  const double nproc = double(num_rank);

  // get copies of weights, numerators and denominators
  std::vector<double> wtv(_w_history);
  std::vector<double> nmv(_lev_history);
  std::vector<double> dnv(_vg_history);

  // get vector to estimate of the errors at different blocking levels
  std::vector<double> estimates;

  // each process gets the sample count and variance for each level of blocking
  int sample_end = nmv.size();
  for (int blk_lvl = 0; sample_end >= 200; blk_lvl++) {

    // compute and print error estimate for this level of blocking
    double eng, var;
    cqmc::mpi_unbiased_ratio_of_means(sample_end, &wtv.at(0), &nmv.at(0), &dnv.at(0), true, eng, var);
    estimates.push_back( std::sqrt( var / sample_end / nproc ) );
    if ( my_rank == 0 && print )
      fout << boost::format("  error estimate for blocks of size 2^(%2i) = %20.12f") % (blk_lvl+1) % (*estimates.rbegin()) << std::endl;

    // create new sample set by averaging adjacent elements of previous sample
    int j = 0;
    //to_reduce.push_back(0.0); // accumulation for variance of new sample
    for (int i = 0; i + 1 < sample_end; i += 2, j += 1) {
      wtv[j] = wtv[i] + wtv[i+1];
      nmv[j] = ( wtv[i] * nmv[i] + wtv[i+1] * nmv[i+1] ) / ( wtv[i] + wtv[i+1] );
      dnv[j] = ( wtv[i] * dnv[i] + wtv[i+1] * dnv[i+1] ) / ( wtv[i] + wtv[i+1] );
    }

    // move the end-of-sample value to the end of the new, smaller sample set
    sample_end = j;

  }
	
  // compute overall error estimate as the average of that of the last few blocks
  if ( estimates.size() < 1 )
    throw formic::Exception("no blocks in energy's statistical error analysis");
  const int n_to_avg = std::min(4, int(estimates.size()));
  _serr = std::accumulate(estimates.rbegin(), estimates.rbegin()+n_to_avg, 0.0) / double(n_to_avg);

  // return the overall error estimate
  return _serr;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  performs a recursive blocking analysis of the target function numerator statistical error
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::ETCompute::target_fn_nuerr(std::ostream & fout, const bool print)
{
    
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // make sure the sample was not exact sample
  if ( _exact_sampling )
    throw formic::Exception("recursive blocking error analysis should not be performed for exact sampling");

  // make sure we are not doing ground state calculation
  if ( _ground_state ) 
    throw formic::Exception("target function numerator error calculation should not be performed when doing ground state calculation");

  // get the number of processes
  const double nproc = double (num_rank);

  // get copies of weights, numerators(omega - E) and denominators(|value/guiding|^2)
  std::vector<double> wtv(_w_history);
  std::vector<double> nmv(_tnv_history);
  std::vector<double> dnv(_vg_history);

  // get vector to estimates of error at different blocking levels 
  std::vector<double> estimates;

  // each process gets the sample count and variance for each level of blocking 
  int sample_end = nmv.size();
  for (int blk_lvl = 0; sample_end >= 200; blk_lvl++) 
  {
    // compute and print error estimate for this level of blocking
    double value, variance;
    cqmc::mpi_unbiased_ratio_of_means(sample_end, &wtv.at(0), &nmv.at(0), &dnv.at(0), true, value, variance);
    estimates.push_back(std::sqrt( variance / sample_end / nproc) );
    if (my_rank == 0 && print)
      fout << boost::format("  tn error estimates for block of size 2^(%2i) = %20.12f") % (blk_lvl+1) % (*estimates.rbegin() ) << std::endl;

    // create new sample by averaging adjacent elements of previous sample 
    int j = 0;
    for (int i = 0; i + 1 < sample_end; i += 2, j += 1)
      {
        wtv[j] = wtv[i] + wtv[i+1];
	nmv[j] = ( wtv[i] * nmv[i] + wtv[i+1] * nmv[i+1] ) / ( wtv[i] + wtv[i+1]);
	dnv[j] = ( wtv[i] * dnv[i] + wtv[i+1] * dnv[i+1] ) / ( wtv[i] + wtv[i+1]);
      }

    sample_end = j;
  }

  // compute the overall error estimate as the average of that of the last few blocks 
  if ( estimates.size() < 1)
    throw formic::Exception("no blocks in target function numerator statistical error analysis");
  const int n_to_avg = std::min(4, int(estimates.size()));
  double serr = std::accumulate(estimates.rbegin(), estimates.rbegin() + n_to_avg, 0.0) / double(n_to_avg);

  // return the overall error estimate 
  return serr;

}
    
//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief performs a recursive blocking analysis of the target function denominator statistical error
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::ETCompute::target_fn_dnerr(std::ostream & fout, const bool print)
{
      
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // make sure the sample was not exact sample 
  if ( _exact_sampling)
    throw formic::Exception("recursive blocking error analysis should not be performed for exact sampling");

  // make sure we are not doing ground state calculation
  if ( _ground_state ) 
   throw formic::Exception("target function denominator error calculation should not be performed when doing ground state calculation");

  // get the number of processes 
  const double nproc = double (num_rank);
	  
  // get copies of weights, numerators(omega - E) and denominators(|value/guiding|^2)
  std::vector<double> wtv(_w_history);
  std::vector<double> nmv(_tdv_history);
  std::vector<double> dnv(_vg_history);

  // get vector to estimates of error at different blocking levels 
  std::vector<double> estimates;

  // each process gets the sample count and variance for each level of blocking 
  int sample_end = nmv.size();
  for (int blk_lvl = 0; sample_end >= 200; blk_lvl++) 
  {
    // compute and print error estimate for this level of blocking
    double value, variance;
    cqmc::mpi_unbiased_ratio_of_means(sample_end, &wtv.at(0), &nmv.at(0), &dnv.at(0), true, value, variance);
    estimates.push_back(std::sqrt( variance / sample_end / nproc) );
    if (my_rank == 0 && print)
      fout << boost::format("  tn error estimates for block of size 2^(%2i) = %20.12f") % (blk_lvl+1) % (*estimates.rbegin() ) << std::endl;

    // create new sample by averaging adjacent elements of previous sample 
    int j = 0;
    for (int i = 0; i + 1 < sample_end; i += 2, j += 1)
    {
      wtv[j] = wtv[i] + wtv[i+1];
      nmv[j] = ( wtv[i] * nmv[i] + wtv[i+1] * nmv[i+1] ) / ( wtv[i] + wtv[i+1]);
      dnv[j] = ( wtv[i] * dnv[i] + wtv[i+1] * dnv[i+1] ) / ( wtv[i] + wtv[i+1]);
    }

    sample_end = j;
  }

  // compute the overall error estimate as the average of that of the last few blocks 
  if ( estimates.size() < 1)
    throw formic::Exception("no blocks in target function numerator statistical error analysis");
  const int n_to_avg = std::min(4, int(estimates.size()));
  double serr = std::accumulate(estimates.rbegin(), estimates.rbegin() + n_to_avg, 0.0) / double(n_to_avg);

  // return the overall error estimate 
  return serr;

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// \brief calculate the statistical error of target function value
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::engine::ETCompute::target_fn_serr(std::ostream & fout, const bool print)
{
  // make sure the sample was not exact sample
  if ( _exact_sampling )
    throw formic::Exception("statistical error analysis should not be performed for exact sampling");

  // make sure we are not doing ground state calculation
  if ( _ground_state )
    throw formic::Exception("target function statistical error should not be performed when doing ground state calculation");

  // calculate the error of numerator
  _tnserr = this -> target_fn_nuerr(fout, print);

  // calculate the error of denominator
  _tdserr = this -> target_fn_dnerr(fout, print);

  // calculate the mean of numerator
  double _tf_numerator = _hd_shift - _energy;

  // calculate the mean of denominator
  double _tf_denominator = _hd_shift * _hd_shift - 2 * _hd_shift * _energy + _energy_s;

  // calculate the statistical error of target function
  _tserr = std::abs(_target_fn_val) * std::sqrt( (_tnserr / _tf_numerator) * (_tnserr / _tf_numerator) + (_tdserr / _tf_denominator) * (_tdserr / _tf_denominator) );

  return _tserr;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief funciton that calculates average energy and target function value
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::ETCompute::en_tar()
{

  // compute the average of energy and its variance 
  cqmc::mpi_unbiased_ratio_of_means(_w_history.size(), &_w_history.at(0), &_lev_history.at(0), &_vg_history.at(0), !_exact_sampling, _energy, _variance);

  // compute the average of energy^2 and its variance 
  cqmc::mpi_unbiased_ratio_of_means(_w_history.size(), &_w_history.at(0), &_lesv_history.at(0), &_vg_history.at(0), !_exact_sampling, _energy_s, _variance_s);

  // if we are doing excted state calculation, compute target function value
  if ( !_ground_state ) {
    cqmc::mpi_unbiased_ratio_of_means(_w_history.size(), &_w_history.at(0), &_tnv_history.at(0), &_tdv_history.at(0), !_exact_sampling, _target_fn_val, _target_fn_var);
    
    // correct for the variance
    //_target_fn_val -= _variance;
  }
  
  // we have not yet calculated the statistical error, so set it to negative to indicate this
  _serr = -1.0;

  // we have not yet calculated the target function statistical error, so set it to negative to indicate this
  _tserr = -1.0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief function that add the "variance correct" term to target function 
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::ETCompute::correct_finite_variance()
{
  // if we are doing ground state calculation, throw out an error
  if ( _ground_state ) 
    throw formic::Exception("variance correct algorithm can only be applied to excited state calculations");

  // evaluate corrections to target function
  double correction = _var_weight * 10.0 * (_energy_s - 2.0 * _hd_shift * _energy + _hd_shift * _hd_shift);

  // add corrections to target function 
  _target_fn_val += correction;

  //_target_fn_val = _energy_s - 2.0 * _hd_shift * _energy + _hd_shift * _hd_shift;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// \brief  prints some statistics about the most recent sample 
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::ETCompute::print_statistics(std::ostream & fout) 
{
    
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank(); 
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

  // compute mean and variance 
  double mean = 0.0;
  if ( !_ground_state )
    mean = _hd_shift - 1.0 / _target_fn_val;

  double le_mean = _energy;
      
  if (my_rank == 0) {
    fout << boost::format("energy and target function accumulation statistics:") << std::endl;
    fout << std::endl;
  }

  // print out recursive blocking information 
  const double err_rec = ( _exact_sampling ? 0.0 : this -> recursive_blocking(fout) );

  double terr_rec = 0.0;
  if ( !_ground_state )
    terr_rec = ( _exact_sampling ? 0.0 : this -> target_fn_serr(fout) );

  if (my_rank == 0) {
    fout << std::endl;
    if ( !_ground_state )
      fout << boost::format("      target function = %20.12f") % _target_fn_val << std::endl;
    else
      fout << boost::format("      target function = %20s") % "N/A" << std::endl;

    fout << boost::format("              le_mean = %20.12f") % le_mean << std::endl;
    fout << boost::format("             les_mean = %20.12f") % _energy_s << std::endl;
    if ( _exact_sampling ) {
      fout << boost::format("             stat err = %20s") % "N/A" << std::endl;
      fout << boost::format("             autocorr = %20s") % "N/A" << std::endl;
    } else {
      fout << boost::format("             stat err = %20.12f") % err_rec << std::endl;
      fout << boost::format("             autocorr = %20.12f") % ( err_rec * err_rec * num_rank * _le_history.size() / _variance ) << std::endl;

      if ( !_ground_state ) {
        fout << boost::format("   target nu stat err = %20.12f") % _tnserr << std::endl;
        fout << boost::format("   target dn stat err = %20.12f") % _tdserr << std::endl;
        fout << boost::format("      target stat err = %20.12f") % terr_rec << std::endl;
      }
    }
    fout << boost::format("              std dev = %20.12f") % std::sqrt(_variance) << std::endl;
    fout << boost::format("             variance = %20.12f") % _variance << std::endl;
    fout << std::endl;
  } 
}

/// \brief function that returns average of local energy 
double cqmc::engine::ETCompute::energy() const { return _energy; }

/// \brief function that returns target function value 
double cqmc::engine::ETCompute::tar_fn_val() const { return _target_fn_val; }

/// \brief function that returns variance 
double cqmc::engine::ETCompute::variance() const { return _variance; }

/// \brief function that returns statistical error of local energy 
double cqmc::engine::ETCompute::eserr(std::ostream & fout) {
  if ( _exact_sampling ) 
    return 0.0;
  if ( _serr < 0.0)
    this -> recursive_blocking(fout, false);
  return _serr;
}

/// \brief function that returns statistical error of target function value 
double cqmc::engine::ETCompute::tserr(std::ostream & fout) {
  if ( _exact_sampling ) 
    return 0.0;
  if ( _tserr < 0.0) 
    this -> target_fn_serr(fout, false);
  return _tserr;
}



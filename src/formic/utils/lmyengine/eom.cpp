///////////////////////////////////////////////////////////////////////////////////////////////
// \brief implementation file for equation of motion JAGP
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include<complex>
#include<vector>
#include<utility>
#include<numeric>
#include<algorithm>
//#include<mpi.h>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/matrix.h>
#include<formic/utils/exception.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lapack_interface.h>
#include<formic/utils/lmyengine/eom.h>

//////////////////////////////////////////////////////////////////////////////////////////////
// \brief constructor that initializes this eom object 
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
cqmc::engine::EOM::EOM(formic::Matrix<double> & hmat, 
                       formic::Matrix<double> & smat,
                       formic::Matrix<double> & ssmat,
                       const bool ssquare,
                       const bool print_matrix,
                       const bool pm_ortho_first,
                       const bool jas_fixed,
                       const int n_sites,
                       const int n_pm,
                       const int n_jas,
                       const double inv_threshold)
:_hmat(hmat),
_smat(smat),
_ssmat(ssmat),
_ssquare(ssquare),
_print_matrix(print_matrix),
_pm_ortho_first(pm_ortho_first),
_jas_fixed(jas_fixed),
_n_sites(n_sites),
_n_pm(n_pm),
_n_pm_zero(0),
_n_jas(n_jas),
_n_jas_zero(0),
_inv_threshold(inv_threshold)
{

  // if the jastrow factor is held fixed, set the number of jastrow factor variable to zero 
  if ( _jas_fixed ) 
    _n_jas = 0;
}

////////////////////////////////////////////////////////////////////////////////////
// \brief a simple way to do eom calculation 
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::EOM::eom_calculation_simple(std::ostream & fout)
{
  
  int my_rank = formic::mpi::rank();

  // solve the eigenvalue problem on root process 
  if (my_rank == 0) {
    
    
    // print out the matrices if requested 
    if ( _print_matrix ) {
     
      // print Hamiltonian matrix 
      fout << boost::format("hamiltonian matrix is") << std::endl;
      for (int i = 0; i < _hmat.rows(); i++) {
        for (int j = 0; j < _hmat.cols(); j++) {
          fout << boost::format("%20.8e  ") % _hmat.at(i, j);
        }
        fout << std::endl;
      }

      fout << std::endl;

      // print overlap matrix 
      fout << boost::format("overlap matrix is") << std::endl;
      for (int i = 0; i < _smat.rows(); i++) {
        for (int j = 0; j < _smat.cols(); j++) {
          fout << boost::format("%20.8e  ") % _smat.at(i, j);
        }
        fout << std::endl;
      }
    }
    
    // set an array to store eigenvalues 
    formic::ColVec<std::complex<double> > evals;


    // solve the eigenvalue problem 
    (ovlp_pseudo_inv(fout) * _hmat).nonsym_eig(evals, _evecs);

    // make a sortable vector of the eigenvalues 
    for (int i = 0; i < evals.size(); i++) 
      _energy_index.push_back(std::make_pair(evals.at(i), i));

    // sort the eigenvalues while preserving their index 
    for (int i = 0; i < _energy_index.size(); i++) {
      for (int j = i + 1; j < _energy_index.size(); j++) {
        
        // get a temp pair to store ith pair 
        std::pair<std::complex<double>, int> temp_pair = _energy_index.at(i);

        // if the value of real part of jth < ith, swap these two pairs 
        if ( (_energy_index.at(j)).first.real() < (_energy_index.at(i)).first.real() ) {
          _energy_index.at(i) = _energy_index.at(j);
          _energy_index.at(j) = temp_pair;
        }
      }
    }

    // calculate S^2 value if requested 
    if ( _ssquare ) { 

      _ss_vals.clear();
      // loop over all states
      for (int i = 0; i < _energy_index.size(); i++) {
        
       	// get state index 
        const int index = _energy_index.at(i).second;
	
        // get denominator and numerator 
        std::complex<double> de(0.0, 0.0);
        std::complex<double> nu(0.0, 0.0);

        // calculate denominator and numerator 
        for (int k = 0; k < _energy_index.size(); k++) {
          for (int l = 0; l < _energy_index.size(); l++) {
            nu += _evecs.at(k, index) * _evecs.at(l, index) * _ssmat.at(k, l);
            de += _evecs.at(k, index) * _evecs.at(l, index) * _smat.at(k, l);
          }
        }
        _ss_vals.push_back(nu / de);
      }   
    }

    // get the real part of the eigenvector
    formic::Matrix<double> _evecs_real(_evecs.rows(), _evecs.cols());
    for (int i = 0; i < _evecs.rows(); i++) {
      for (int j = 0; j < _evecs.cols(); j++) {
        _evecs_real.at(i, j) = _evecs.at(i, j).real();
      }
    }

    // normalize the eigenvectors
    for (int i = 0; i < _evecs.cols(); i++) {
      // calculate the norm of ith eigenvector
      formic::ColVec<double> col_i = _evecs_real.col_as_vec(i);
      double norm = std::sqrt(col_i.norm2());
      // normalize
      _evecs.scale_col_by(i, 1.0/norm);
    }

  }
}


////////////////////////////////////////////////////////////////////////////////////
// \brief calculate the pseudo inverse of overlap matrix(this function is not used  
//        now)
//
//
//
////////////////////////////////////////////////////////////////////////////////////
formic::Matrix<double> cqmc::engine::EOM::ovlp_pseudo_inv(std::ostream & fout) 
{
  
  int my_rank = formic::mpi::rank();
 
  // compute SVD on root process 
  formic::Matrix<double> u;
  formic::Matrix<double> vt;
  formic::Matrix<double> v;
  formic::ColVec<double> sin_vals;
  if (my_rank == 0) {
    
    // performs svd decomposition to the overlap matrix
    _smat.svd(u, sin_vals, vt);
    v = vt.t();

    // print out overlap matrix's singular values 
    fout << boost::format("transformed overlap matrix's singular values are") << std::endl;
    for (int i = 0; i < sin_vals.size(); i++) {
      fout << boost::format("%20.8e") % sin_vals.at(i) << std::endl;
    }
    fout << std::endl;

    // compute and return inverse 
    for (int i = 0; i < v.cols(); i++) {
      double scaler = (sin_vals.at(i) > _inv_threshold ? 1.0 / sin_vals.at(i) : 0.0);
      v.scale_col_by(i, scaler);
    }

  }
  return v * u.t();
}

//////////////////////////////////////////////////////////////////////////////////////////////
// \brief print out eom calculation statistics 
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::EOM::eom_print(std::ostream & fout)
{
  
  int my_rank = formic::mpi::rank();

  // print statistics on root process 
  if (my_rank == 0) {
   
    fout << boost::format("EOM JAGP state statistics:") << std::endl;
    if ( ! _ssquare ) {
      for (int i = 0; i < _energy_index.size(); i++) 
        fout << boost::format("energy =  %20.12f   + %20.12f i") % (_energy_index.at(i).first).real() % (_energy_index.at(i).first).imag() << std::endl;
      fout << std::endl;
    }

    else {
      for (int i = 0; i < _energy_index.size(); i++) 
        fout << boost::format("energy =  %20.12f   + %20.12f i") % (_energy_index.at(i).first).real() % (_energy_index.at(i).first).imag() << boost::format("    S^2 = %20.12f  + %20.12f i") % _ss_vals.at(i).real() % _ss_vals.at(i).imag() << std::endl;
      fout << std::endl;

      return;

    }

  } 
}
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
_n_jas_zero(0),
_n_jas(n_jas),
_inv_threshold(inv_threshold)
{

  // if the jastrow factor is held fixed, set the number of jastrow factor variable to zero 
  if ( _jas_fixed ) 
    _n_jas = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
// \brief do equation of motion calculation 
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_calculation(std::ostream & fout)
//{
//  
//  // get rank number and the number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  // solve the eigenvalue problem on root process 
//  if (my_rank == 0) {
//    
//    
//    // print out the matrices if requested 
//    if ( _print_matrix ) {
//     
//      // print Hamiltonian matrix 
//      fout << boost::format("hamiltonian matrix is") << std::endl;
//      for (int i = 0; i < _hmat.rows(); i++) {
//        for (int j = 0; j < _hmat.cols(); j++) {
//          fout << boost::format("%20.8e  ") % _hmat.at(i, j);
//        }
//        fout << std::endl;
//      }
//
//      fout << std::endl;
//
//      // print overlap matrix 
//      fout << boost::format("overlap matrix is") << std::endl;
//      for (int i = 0; i < _smat.rows(); i++) {
//        for (int j = 0; j < _smat.cols(); j++) {
//          fout << boost::format("%20.8e  ") % _smat.at(i, j);
//        }
//        fout << std::endl;
//      }
//    }
//
//    // transform to a normalized basis 
//    this -> eom_ovl_normalize();
//
//    // project ground state out of first derivative space 
//    this -> eom_ground_proj();
//
//    // normalize the derivative vectors
//    this -> eom_ovl_normalize();
//
//    if ( _pm_ortho_first ) {
//    
//      // orthonormalize pairing matrix space 
//      this -> eom_pm_ortho(pm_transform);
//      _hmat = pm_transform.t() * _hmat * pm_transform;
//      _smat = pm_transform.t() * _smat * pm_transform;
//      _ssmat = pm_transform.t() * _ssmat * pm_transform;
//
//      if ( !_jas_fixed ) { 
//        // project pairing matrix space out of jastrow space 
//        this -> eom_proj_jas_out_of_pm();
//      }
//
//      if ( !_jas_fixed ) { 
//        // orthonormalize jastrow factor space 
//        this -> eom_jas_ortho(jas_transform);
//        _hmat = jas_transform.t() * _hmat * jas_transform;
//        _smat = jas_transform.t() * _smat * jas_transform;
//        _ssmat = jas_transform.t() * _ssmat * jas_transform;
//      }
//    }
//
//
//    else {
//
//      // if jastrow factor variable is held fixed, throw out an error 
//      if ( _jas_fixed ) 
//        throw formic::Exception("Jastrow factor is fixed!");
//      
//      // orthonormalize jastrow factor space 
//      this -> eom_jas_ortho(jas_transform);
//      _hmat = jas_transform.t() * _hmat * jas_transform;
//      _smat = jas_transform.t() * _smat * jas_transform;
//      _ssmat = jas_transform.t() * _ssmat * jas_transform;
//
//      // project pairing matrix component out of jastrow factor space  
//      this -> eom_proj_pm_out_of_jas();
//
//      // orthonormalize pairing matrix space 
//      this -> eom_pm_ortho(pm_transform);
//      _hmat = pm_transform.t() * _hmat * pm_transform;
//      _smat = pm_transform.t() * _smat * pm_transform;
//      _ssmat = pm_transform.t() * _ssmat * pm_transform;
//
//    }
//
//    // solve the eigenvalue problem 
//    Eigen::EigenSolver<Eigen::MatrixXd> es(_hmat);
//
//    // set an array to store eigenvalues 
//    Eigen::ArrayXcd evals;
//    evals = (es.eigenvalues().array());
//
//    // store eigenvectors 
//    _evecs = es.eigenvectors();
//
//    // make a sortable vector of the eigenvalues 
//    for (int i = 0; i < evals.size(); i++) 
//      _energy_index.push_back(std::make_pair(evals(i), i));
//
//    // sort the eigenvalues while preserving their index 
//    for (int i = 0; i < _energy_index.size(); i++) {
//      for (int j = i + 1; j < _energy_index.size(); j++) {
//        
//        // get a temp pair to store ith pair 
//        std::pair<std::complex<double>, int> temp_pair = _energy_index.at(i);
//
//        // if the value of real part of jth < ith, swap these two pairs 
//        if ( (_energy_index.at(j)).first.real() < (_energy_index.at(i)).first.real() ) {
//          _energy_index.at(i) = _energy_index.at(j);
//          _energy_index.at(j) = temp_pair;
//        }
//      }
//    }
//
//    // calculate S^2 value if requested 
//    if ( _ssquare ) { 
//
//      _ss_vals.clear();
//      // loop over all states
//      for (int i = 0; i < _energy_index.size(); i++) {
//        
//       	// get state index 
//        const int index = _energy_index.at(i).second;
//	
//        // get denominator and numerator 
//        std::complex<double> de(0.0, 0.0);
//        std::complex<double> nu(0.0, 0.0);
//
//        // calculate denominator and numerator 
//        for (int k = 0; k < _energy_index.size(); k++) {
//          for (int l = 0; l < _energy_index.size(); l++) {
//            nu += _evecs(k, index) * _evecs(l, index) * _ssmat(k, l);
//            de += _evecs(k, index) * _evecs(l, index) * _smat(k, l);
//          }
//        }
//        _ss_vals.push_back(nu / de);
//      }   
//    } 
//
//    // transfer eigenvectors back to the original normalized and with ground state being projected out space 
//    if ( _pm_ortho_first ) {
//      if ( !_jas_fixed ) { 
//        _evecs = pm_transform * jas_transform * _evecs;
//      }
//      else { 
//        _evecs = pm_transform * _evecs;
//      }
//    }
//    else 
//      _evecs = jas_transform * pm_transform * _evecs;
//
//    // get the real part of the eigenvector
//    Eigen::MatrixXd _evecs_real(_evecs.rows(), _evecs.cols());
//    for (int i = 0; i < _evecs.rows(); i++) {
//      for (int j = 0; j < _evecs.cols(); j++) {
//        _evecs_real(i, j) = _evecs(i, j).real();
//      }
//    }
//
//    // normalize the eigenvectors
//    for (int i = 0; i < _evecs.cols(); i++) {
//      // calculate the norm of ith eigenvector
//      double norm = std::sqrt(_evecs_real.col(i).dot(_evecs_real.col(i)));
//      // normalize
//      _evecs.col(i) /= norm;
//    }
//
//
//  }
//
//  // end of the calculation function 
//}

////////////////////////////////////////////////////////////////////////////////////
// \brief a simple way to do eom calculation 
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::EOM::eom_calculation_simple(std::ostream & fout)
{
  
  // get rank number and the number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

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
  
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
 
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
  
  // get rank number and number of ranks 
  int my_rank = formic::mpi::rank();
  int num_rank = formic::mpi::size();
  //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
  //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);

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

      // print the eigenvectors of the first 10 singlet states 
      int count = 0;
      for (int i = 0; i < 11 && count < _energy_index.size(); ) {
        
        // check to see if it's a singlet state
        if ( _ss_vals.at(count).real() < 2.5 ) {
          
          // print out eigenvectors 
          fout << boost::format("eigenvector of the %d th state(energy = %20.12f   + %20.12f i) is") % i % (_energy_index.at(count).first).real() % (_energy_index.at(count).first).imag() << std::endl;
          for (int j = 0; j < _evecs.rows(); j++) {
            fout << boost::format("%20.12f   + %20.12f i") % _evecs(j, _energy_index.at(count).second).real() % _evecs(j, _energy_index.at(count).second).imag();
            
            // print the corresponding wavefunction variable name 
            if ( j == 0 ) 
              fout << boost::format("  undifferentiated wavefunction") << std::endl;

            else if ( j > 0 && j <= _n_sites * _n_sites ) {
              int row = (j - 1) / _n_sites + 1;
              int col = j - (row - 1) * _n_sites ;
              fout << boost::format("  F %d %d") % row % col << std::endl;
            }
          }
          i++;
          fout << std::endl;
        }
        count++;
      }
    }

  } 
}

//////////////////////////////////////////////////////////////////////////////////////////////
// \brief transform derivative vector to a normalized space by the following transformation 
//        H = D^(-1/2) * H * D^(-1/2) S = D^(-1/2) * S * D^(-1/2)
//        D^(-1/2) is a diagonal matrix which has the 1/ sqrt(norm of derivative vectors) as  
//        diagonal elements
//
//////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_ovl_normalize()
//{
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  if ( my_rank == 0 ) {
//    // calculate D^(-1/2) matrix
//    _d_inv_half.resize(_hmat.rows(), _hmat.cols());
//    for (int i = 0; i < _d_inv_half.rows(); i++) { 
//      for (int j = 0; j < _d_inv_half.cols(); j++) {
//        _d_inv_half(i, j) = ( i == j ? 1.0 / std::sqrt(_smat(i, i)) : 0.0);
//      }
//    }
//
//    // do the transformation 
//    _smat = _d_inv_half * _smat * _d_inv_half;
//    _hmat = _d_inv_half * _hmat * _d_inv_half;
//
//    // if requested, also transform S^2 matrix 
//    if ( _ssquare ) 
//      _ssmat = _d_inv_half * _ssmat * _d_inv_half;
//  }
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//// \brief project ground state out of derivative space 
////
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_ground_proj()
//{
//  
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//
//  if ( my_rank == 0 ) {
//    
//    // get matrix dimension
//    const int m = _hmat.rows();
//    const int n = _smat.cols();
//    
//    // get temp matrix 
//    Eigen::MatrixXd _smatg(m, n);
//    Eigen::MatrixXd _hmatg(m, n);
//    Eigen::MatrixXd _ssmatg(m, n);
//    for (int i = 0; i < m; i++) {
//      for (int j = 0; j < n; j++) {
//        if ( !( i == 0 && j == 0) ) {
//          _smatg(i, j) = _smat(i, j) - (_smat(j, 0) * _smat(i, 0) / _smat(0, 0) + _smat(i, 0) * _smat(0, j) / _smat(0, 0) - _smat(0, i) * _smat(j, 0) / _smat(0, 0));
//          _hmatg(i, j) = _hmat(i, j) - (_smat(j, 0) * _hmat(i, 0) / _smat(0, 0) + _smat(0, i) * _hmat(0, j) / _smat(0, 0) - _smat(0, i) * _smat(j, 0) * _hmat(0, 0) / (_smat(0, 0) * _smat(0, 0)));
//          if ( _ssquare ) 
//            _ssmatg(i, j) = _ssmat(i, j) - (_smat(j, 0) * _ssmat(i, 0) / _smat(0, 0) + _smat(i, 0) * _ssmat(0, j) / _smat(0, 0) - _smat(0, i) * _smat(j, 0) * _ssmat(0, 0) / (_smat(0, 0) * _smat(0, 0)));
//        } else {
//          _smatg(i, j) = 1.0;
//          _hmatg(i, j) = _hmat(i, j);
//          if ( _ssquare ) 
//            _ssmatg(i, j) = _ssmat(i, j);
//        }
//      }
//    }
//
//    _hmat = _hmatg;
//    _smat = _smatg;
//    if ( _ssquare ) 
//      _ssmat = _ssmatg;
//
//  }
//
//}
//	
////////////////////////////////////////////////////////////////////////////////////////////////
//// \brief calculate the transformation matrix that orthonormalizes pairing matrix derivative 
////        subspace 
////   
////
////////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_pm_ortho(Eigen::MatrixXd & pm_transform)
//{
//  
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  if ( my_rank == 0 ) {
//    
//    // get the pairing matrix subspace out of the whole overlap matrix 
//    Eigen::MatrixXd pm_block(_n_pm, _n_pm);
//    pm_block = _smat.block(1, 1, _n_pm, _n_pm);
//
//    // do SVD on this block 
//    Eigen::MatrixXd U;
//    Eigen::MatrixXd V;
//    Eigen::VectorXd sin_vals;
//    Eigen::JacobiSVD<Eigen::MatrixXd> svd(pm_block, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    U = svd.matrixU();
//    V = svd.matrixV();
//    sin_vals = svd.singularValues();
//
//    // compute the inverse of the square root of the singular values and count the number of nonzero singular values 
//    int _n_pm_nonzero = 0;
//    _n_pm_zero = 0;
//    for (int i = 0; i < sin_vals.size(); i++) {
//      if (sin_vals(i) >= _inv_threshold) {
//        sin_vals(i) = 1.0 / std::sqrt(sin_vals(i));
//	_n_pm_nonzero ++;
//      }
//    }
//
//    _n_pm_zero = _n_pm - _n_pm_nonzero;
//
//    // truncate U, V and singular value vector based on the number of non zero singular values 
//    U.conservativeResize(_n_pm, _n_pm_nonzero);
//    V.conservativeResize(_n_pm, _n_pm_nonzero);
//    sin_vals.conservativeResize(_n_pm_nonzero);
//
//    // transfrom singular value vector to diagonal matrix 
//    Eigen::MatrixXd sin_val_matrix = sin_vals.asDiagonal();
//
//    // compute S^(-1/2)
//    Eigen::MatrixXd s_half_inv = V * sin_val_matrix;
//
//    // initialize the transformation matrix 
//    pm_transform.resize(_smat.rows(), _smat.cols() - _n_pm_zero);
//    for (int i = 0; i < pm_transform.rows(); i++) {
//      for (int j = 0; j < pm_transform.cols(); j++) {
//         pm_transform(i, j) = 0.0;
//      }
//    }
//    pm_transform(0, 0) = 1.0;
//
//    // initialize jastrow block as identity matrix 
//    Eigen::MatrixXd jblock(_n_jas, _n_jas);
//    for (int i = 0; i < jblock.rows(); i++) {
//      for (int j = 0; j < jblock.cols(); j++) {
//        jblock(i, j) = (i == j ? 1.0 : 0.0);
//      }
//    }
//    
//    // assign the corresponding block 
//    pm_transform.block(1, 1, _n_pm, _n_pm - _n_pm_zero) = s_half_inv;
//    pm_transform.block(_n_pm + 1, _n_pm + 1 - _n_pm_zero, _n_jas, _n_jas) = jblock;
//
//    // update the number of pairing matrix element 
//    _n_pm -= _n_pm_zero;
//
//  }
//
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//// \brief calculate the transformation matrix that orthonormalizes jastrow factor derivative
////        subspace 
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_jas_ortho(Eigen::MatrixXd & jas_transform)
//{
//  
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  if ( my_rank == 0 ) {
//    
//    // first get the jastrow factor derivative block from the whole overlap matrix 
//    Eigen::MatrixXd jas_block(_n_jas, _n_jas);
//    jas_block = _smat.block(_n_pm + 1, _n_pm + 1, _n_jas, _n_jas);
//
//    // do SVD on this block
//    Eigen::MatrixXd U;
//    Eigen::MatrixXd V;
//    Eigen::VectorXd sin_vals;
//    Eigen::JacobiSVD<Eigen::MatrixXd> svd(jas_block, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    U = svd.matrixU();
//    V = svd.matrixV();
//    sin_vals = svd.singularValues();
//
//    // compute the inverse of the square root of the singular values and count the number of nonzero singular values 
//    int _n_jas_nonzero = 0;
//    for (int i = 0; i < sin_vals.size(); i++) {
//      if (sin_vals(i) >= _inv_threshold) { 
//        sin_vals(i) = 1.0 / std::sqrt(sin_vals(i));
//        _n_jas_nonzero ++;
//      }
//    }
//
//    _n_jas_zero = _n_jas - _n_jas_nonzero;
//
//    U.conservativeResize(_n_jas, _n_jas_nonzero);
//    V.conservativeResize(_n_jas, _n_jas_nonzero);
//    sin_vals.conservativeResize(_n_jas_nonzero);
//
//    // transform to a diagnal matrix 
//    Eigen::MatrixXd sin_val_matrix = sin_vals.asDiagonal();
//
//    // compute S^(-1/2)
//    Eigen::MatrixXd s_half_inv = V * sin_val_matrix;
//
//    // initialize transformation matrix 
//    jas_transform.resize(_smat.rows(), _smat.cols() - _n_jas_zero);
//    for (int i = 0; i < jas_transform.rows(); i++) {
//      for (int j = 0; j < jas_transform.cols(); j++) {
//        jas_transform(i, j) = 0.0;
//      }
//    }
//    jas_transform(0, 0) = 1.0;
//
//    // set the pairing matrix block the identity matrix 
//    Eigen::MatrixXd pblock(_n_pm, _n_pm);
//    for (int i = 0; i < pblock.rows(); i++) {
//      for (int j = 0; j < pblock.cols(); j++) {
//        pblock(i, j) = (i == j ? 1.0 : 0.0);
//      }
//    }
//
//    // assign the corresponding block 
//    jas_transform.block(_n_pm + 1, _n_pm + 1, _n_jas, _n_jas - _n_jas_zero) = s_half_inv;
//    jas_transform.block(1, 1, _n_pm, _n_pm) = pblock;
//   
//    // update the number of jastrow factor elements
//    _n_jas -= _n_jas_zero;
//  }
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//// \brief project pairing matrix derivative vectors out of jastrow factor space 
////
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_proj_pm_out_of_jas()
//{
//  
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  if ( my_rank == 0 ) {
//     
//    // project jastrow factor space out of pairing matrix space 
//
//    // loop over all pairing matrix derivatives 
//    for (int i = 1; i < _n_pm + 1; i++) {
//       
//      // loop over all jastrow factor 
//      for (int j = _n_pm + 1; j < _smat.rows(); j++) {
//         
//	// vector to store new overlap matrix elements 
//	std::vector<double> s_update;
//
//	// vector to store new hamiltonian matrix elements
//	std::vector<double> h_update;
//
//	// vector to store new hamiltonian transpose elements
//	std::vector<double> h_update_t;
//
//	// vector to store new S^2 matrix elements 
//	std::vector<double> ss_update;
//
//	// vector to store new S^2 transpose elements
//	std::vector<double> ss_update_t;
//
//        // loop over all variables 
//        for (int k = 1; k < _smat.rows(); k++) {
//          s_update.push_back( _smat(k, i) - (_smat(i, j) / _smat(j, j)) * _smat(k, j) );
//	  if (i != k ) {
//	    h_update.push_back( _hmat(k, i) - (_smat(i, j) / _smat(j, j)) * _hmat(k, j) );
//	    h_update_t.push_back( _hmat(i, k) - (_smat(j, i) / _smat(j, j)) * _hmat(j, k) );
//
//	    if ( _ssquare ) {
//	      ss_update.push_back( _ssmat(k, i) - (_smat(i, j) / _smat(j, j)) * _ssmat(k, j) );
//	      ss_update_t.push_back( _ssmat(i, k) - (_smat(j, i) / _smat(j, j)) * _ssmat(j, k) );
//	    }
// 
//	  }
//	  else if (i == k) {
//            double temp = _hmat(i, i) - (_smat(i, j) / _smat(j, j)) * _hmat(i, j) - (_smat(j, i) / _smat(j, j)) * _hmat(j, i) + (_smat(j, i) / _smat(j, j)) * (_smat(i, j) / _smat(j, j)) * _hmat(j, j);
//	    h_update.push_back( temp );
//	    h_update_t.push_back( temp );
//
//	    if ( _ssquare ) {
//	      double temp_s = _ssmat(i, i) - (_smat(i, j) / _smat(j, j)) * _ssmat(i, j) - (_smat(j, i) / _smat(j, j)) * _ssmat(j, i) + (_smat(j, i) / _smat(j, j)) * (_smat(i, j) / _smat(j, j)) * _ssmat(j, j);
//	      ss_update.push_back( temp_s );
//	      ss_update_t.push_back( temp_s );
//	    }
//	  }
//        }
//
//        // update corresponding matrix elements
//        for (int k = 1; k < _smat.rows(); k++) {
//          _smat(k, i) = s_update.at(k - 1);
//          _smat(i, k) = s_update.at(k - 1);
//	  _hmat(k, i) = h_update.at(k - 1);
//	  _hmat(i, k) = h_update_t.at(k - 1);
//	  if ( _ssquare ) {
//	    _ssmat(k, i) = ss_update.at(k - 1);
//	    _ssmat(i, k) = ss_update_t.at(k - 1);
//	  }
//        }
//
//        // clear vectors 
//        s_update.clear();
//        h_update.clear();
//        h_update_t.clear();
//	ss_update.clear();
//	ss_update_t.clear();
//      }    
//    }
//  }
//}
//
////////////////////////////////////////////////////////////////////////////////////////////////
//// \brief project jastrow derivative vectors out of pairing matrix space 
////
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
//void cqmc::engine::EOM::eom_proj_jas_out_of_pm()
//{
//  
//  // get rank number and number of ranks 
//  int my_rank;
//  int num_rank;
//  MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
//  MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
//
//  if ( my_rank == 0 ) {
//   
//    // project pairing matrix space out of jastrow factor space 
//
//    // loop over all jastrow factor derivatives 
//    for (int i = _n_pm + 1; i < _smat.rows(); i++) {
//      
//      // loop over all pairing matrix derivatives 
//      for (int j = 1; j < _n_pm + 1; j++) {
//        
//	// vector to store new overlap matrix elements
//	std::vector<double> s_update;
//
//	// vector to store new hamiltonian matrix elements
//	std::vector<double> h_update;
//
//	// vector to store new hamiltonian matrix transpose 
//	std::vector<double> h_update_t;
//
//	// vector to store new S^2 matrix elements
//	std::vector<double> ss_update;
//
//	// vector to store new S^2 matrix transpose 
//	std::vector<double> ss_update_t;
//
//	// loop over all variables 
//	for (int k = 1; k < _smat.rows(); k++) { 
//	  s_update.push_back( _smat(k, i) - (_smat(i, j) / _smat(j, j)) * _smat(k, j) );
//	  if (i != k ) {
//	    h_update.push_back( _hmat(k, i) - (_smat(i, j) / _smat(j, j)) * _hmat(k, j) );
//	    h_update_t.push_back( _hmat(i, k) - (_smat(j, i) / _smat(j, j)) * _hmat(j, k) );
//
//	    if ( _ssquare ) {
//	      ss_update.push_back( _ssmat(k, i) - (_smat(i, j) / _smat(j, j)) * _ssmat(k, j) );
//	      ss_update_t.push_back( _ssmat(i, k) - (_smat(j, i) / _smat(j, j)) * _ssmat(j, k) );
//	    }
// 
//	  }
//	  else if (i == k) {
//            double temp = _hmat(i, i) - (_smat(i, j) / _smat(j, j)) * _hmat(i, j) - (_smat(j, i) / _smat(j, j)) * _hmat(j, i) + (_smat(j, i) / _smat(j, j)) * (_smat(i, j) / _smat(j, j)) * _hmat(j, j);
//	    h_update.push_back( temp );
//	    h_update_t.push_back( temp );
//
//	    if ( _ssquare ) {
//	      double temp_s = _ssmat(i, i) - (_smat(i, j) / _smat(j, j)) * _ssmat(i, j) - (_smat(j, i) / _smat(j, j)) * _ssmat(j, i) + (_smat(j, i) / _smat(j, j)) * (_smat(i, j) / _smat(j, j)) * _ssmat(j, j);
//	      ss_update.push_back( temp_s );
//	      ss_update_t.push_back( temp_s );
//	    }
//	  }
//        }
//
//        // update corresponding matrix elements
//        for (int k = 1; k < _smat.rows(); k++) {
//          _smat(k, i) = s_update.at(k - 1);
//          _smat(i, k) = s_update.at(k - 1);
//	  _hmat(k, i) = h_update.at(k - 1);
//	  _hmat(i, k) = h_update_t.at(k - 1);
//	  if ( _ssquare ) {
//	    _ssmat(k, i) = ss_update.at(k - 1);
//	    _ssmat(i, k) = ss_update_t.at(k - 1);
//	  }
//        }
//
//        // clear vectors 
//        s_update.clear();
//        h_update.clear();
//        h_update_t.clear();
//	ss_update.clear();
//	ss_update_t.clear();
//      }    
//    }
//  }
//}




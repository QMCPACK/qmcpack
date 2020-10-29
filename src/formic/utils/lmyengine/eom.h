////////////////////////////////////////////////////////////////////////////////////////////////
// \brief euqation of motion JAGP class 
// 
//
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_EOM_HEADER
#define ENGINE_EOM_HEADER

#include<complex> 
#include<vector>
#include<numeric>
#include<algorithm>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

//#include<mpi.h>

#include"formic/utils/matrix.h"
#include"formic/utils/exception.h"
#include"formic/utils/mpi_interface.h"
#include"formic/utils/lapack_interface.h"

namespace cqmc {
  
  namespace engine { 
    
    template<typename S> class EOM {

      private: 

        /// \brief [in] linear method hamiltoian matrix 
        formic::Matrix<S> & _hmat;

        /// \brief [in] linear method overlap matrix 
        formic::Matrix<S> & _smat;

        /// \brief [in] linear method S^2 matrix 
        formic::Matrix<S> & _ssmat;

        /// \brief [in] flag to tell whether to do S^2 calculation 
        const bool _ssquare;

        /// \brief [in] flag to tell whether to print out the matrix when constructing this object 
        const bool _print_matrix; 

        /// \brief [in] flag to tell whether to orthonormalize pairing matrix derivative first 
        const bool _pm_ortho_first;

        /// \brief [in] flag to tell whether jastrow factor is held fixed 
        const bool _jas_fixed;

        /// \brief [in] number of sites 
        int _n_sites;

        /// \brief [in] number of pairing matrix elements 
        int _n_pm;

        /// \brief number of linearly dependent pairing matrix element derivative vectors 
        int _n_pm_zero;

        /// \brief number of jastrow factor matrix elements 
        int _n_jas;

        /// \brief number of linearly dependent jastrow factor matrix element derivative vectors
        int _n_jas_zero;

        /// \brief [in] threshold to use when calculating pseudo inverse of overlap matrix 
        const double _inv_threshold; 

        /// \brief [out] a vector storing sorted energies and their index  
        std::vector<std::pair<std::complex<double>, int> > _energy_index;

        /// \brief [out] a vector storing S^2 value 
        std::vector<std::complex<double> > _ss_vals;

        /// \brief [out] a matrix storing eigenvectors 
        formic::Matrix<std::complex<double> > _evecs;

        /// \brief [out] real sorted eigenvectors
        formic::Matrix<S> _evecs_sorted_real;

        /// \brief [optionaly out] D^(-1/2) matrix 
        formic::Matrix<double> _d_inv_half;

        /// \brief [out] pairing matrix space orthonormalization transformation matrix 
        formic::Matrix<double> pm_transform;

        /// \brief [out] jastrow factor space orthonormalization transformation matrix 
        formic::Matrix<double> jas_transform;  

      public:
       
        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief constructor that initializes this eom object 
        //
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        EOM(formic::Matrix<S> & hmat, 
            formic::Matrix<S> & smat,
            formic::Matrix<S> & ssmat,
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

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief do equation of motion calculation 
        //
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_calculation(std::ostream & fout);

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief a simple way to do eom calculation
        //
        //  
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        void eom_calculation_simple(std::ostream & fout)
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

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief calculate pseudo inverse of overlap matrix 
        //
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        formic::Matrix<S> ovlp_pseudo_inv(std::ostream & fout)
        {
          
          // get rank number and number of ranks 
          int my_rank = formic::mpi::rank();
          int num_rank = formic::mpi::size();
          //MPI_Comm_rank(MPI_COMM_WORLD, & my_rank);
          //MPI_Comm_size(MPI_COMM_WORLD, & num_rank);
         
          // compute SVD on root process 
          formic::Matrix<S> u;
          formic::Matrix<S> vt;
          formic::Matrix<S> v;
          formic::ColVec<S> sin_vals;
          if (my_rank == 0) {
            
            // performs svd decomposition to the overlap matrix
            _smat.svd(u, sin_vals, vt);
            v = vt.c();

            // print out overlap matrix's singular values 
            fout << boost::format("transformed overlap matrix's singular values are") << std::endl;
            for (int i = 0; i < sin_vals.size(); i++) {
              fout << boost::format("%20.8e") % sin_vals.at(i) << std::endl;
            }
            fout << std::endl;

            // compute and return inverse 
            for (int i = 0; i < v.cols(); i++) {
              S scaler = (std::abs(sin_vals.at(i)) > _inv_threshold ? 1.0 / sin_vals.at(i) : 0.0);
              v.scale_col_by(i, scaler);
            }

          }
          return v * u.c();
        }

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief transform derivative vector to a normalized space by the following transformation 
        //        H = D^(-1/2) * H * D^(-1/2) S = D^(-1/2) * S * D^(-1/2)
        //        D^(-1/2) is a diagonal matrix which has the 1/ sqrt(norm of derivative vectors) as  
        //        diagonal elements
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_ovl_normalize();

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief project ground state out of derivative space 
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_ground_proj();

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief calculate the transformation matrix that orthonormalizes pairing matrix derivative 
        //        subspace 
        //   
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_pm_ortho(formic::Matrix<double> & pm_transform);

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief calculate the transformation matrix that orthonormalizes jastrow factor derivative
        //        subspace 
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_jas_ortho(formic::Matrix<double> & jas_transform);

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief project jastrow derivative vectors out of pairing matrix space 
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_proj_jas_out_of_pm();

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief project pairing matrix derivative vectors out of jastrow factor space 
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        //void eom_proj_pm_out_of_jas();
        
        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief print out eom calculation statistics 
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        void eom_print(std::ostream & fout)
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
        // \brief output energy list
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<std::pair<std::complex<double>, int> > & energy_list() { return _energy_index; } 

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief output unsorted eigenvectors
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        formic::Matrix<std::complex<double> > & evecs() { return _evecs; } 
    };
  }
}

#endif

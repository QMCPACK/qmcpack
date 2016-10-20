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

#include<formic/utils/matrix.h>

namespace cqmc {
  
  namespace engine { 
    
    class EOM {

      private: 

        /// \brief [in] linear method hamiltoian matrix 
        formic::Matrix<double> & _hmat;

        /// \brief [in] linear method overlap matrix 
        formic::Matrix<double> & _smat;

        /// \brief [in] linear method S^2 matrix 
        formic::Matrix<double> & _ssmat;

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
        formic::Matrix<double> _evecs_sorted_real;

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
        EOM(formic::Matrix<double> & hmat, 
            formic::Matrix<double> & smat,
            formic::Matrix<double> & ssmat,
            const bool ssquare,
            const bool print_matrix,
            const bool pm_ortho_first,
            const bool jas_fixed,
            const int n_sites,
            const int n_pm,
            const int n_jas,
            const double inv_threshold);

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
        void eom_calculation_simple(std::ostream & fout);

        //////////////////////////////////////////////////////////////////////////////////////////////
        // \brief calculate pseudo inverse of overlap matrix 
        //
        //
        //
        //
        //////////////////////////////////////////////////////////////////////////////////////////////
        formic::Matrix<double> ovlp_pseudo_inv(std::ostream & fout);

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
        void eom_print(std::ostream & fout);

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

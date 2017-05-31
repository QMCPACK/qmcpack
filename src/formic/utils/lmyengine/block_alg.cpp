//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file  /utils/lmyengine/block_alg.cpp       
///
/// \brief  Implementation file for the block recursive algorithm class
///
//////////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>
#include<list>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<sstream>
#include<formic/utils/openmp.h>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/matrix.h>
#include<formic/utils/mpi_interface.h>
#include<formic/utils/lmyengine/block_alg.h>
#include<formic/utils/lmyengine/block_detail.h>
#include<formic/utils/lmyengine/eigen_solver.h>
#include<formic/utils/lmyengine/davidson_solver.h>
#include<formic/utils/lmyengine/spam_solver.h>


/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that accumulates data for the object
///
/// \param[in]   d        weight for this sample
/// \param[in]   dr       bare derivative ratio vector
/// \param[in]   er       energy derivative ratio vector
/// \param[in]   ground   whether this is a ground state calculation
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::acc(const double d, const std::vector<double> & dr, const std::vector<double> & er, const bool ground, const double hd_lm_shift) {
  
  // adds up total weight
  # pragma omp critical
  m_tw += d;

  // accumulate hamiltonian data
  m_hdata.acc(d, dr, er, m_ou);

  // accumulate overlap data
  if ( ground ) 
    m_sdata.acc(d, dr, dr, m_ou);
  else { 
    m_sdata.acc(d, er, er, m_ou);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that resets this object
///
/// \param[in]   nv        number of variables
/// \param[in]   nblock    number of blocks
/// \param[in]   ou        vector of vector storing the old update coefficients
///
//////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::reset(const int nv, const int nblock, const std::vector<formic::ColVec<double> > & ou, const bool ground, const bool iterative) {
  
  // reset the total weight to be zero
  m_tw = 0.0;

  // reset the hamiltonian and overlap matrix data
  m_hdata.reset(nv, nblock, ou.size());
  m_sdata.reset(nv, nblock, ou.size());

  // empty the old update matrix
  m_ou.clear();

  // copy the old update coefficients to old update matrix
  if ( ou.size() > 0 ) 
    if ( ou.at(0).size() > 0 )
      m_ou.reset(ou.at(0).size(), ou.size());
  for (int k = 0; k < m_ou.cols(); k++) {
    ou.at(k) >>= m_ou.col_begin(k);
  }

  _ground = ground;  
  _iterative = iterative;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that finalizes data accumulation
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::finalize() {
  
  m_hdata.finalize(m_tw);
  m_sdata.finalize(m_tw);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that finalizes data accumulation across all processors
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::mpi_finalize(const double total_weight) {
  
  // call mpi finalize function for data matrices
  m_hdata.mpi_finalize(total_weight);
  m_sdata.mpi_finalize(total_weight);
}
  

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief 
///
/// \param[in]    b     the block index
/// \param[in]    x     the block index
/// \param[out]   dd    the output matrix
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::prep_lm_block_plus_other_ou_dd_matrix(const int b, const int x, formic::Matrix<double> & dd) {
  
  // begining index of the block
  const int ibeg = 1 + m_hdata.bb(b);

  // block length
  const int len = m_hdata.bl(b);

  // dimension of the matrix
  const int dim = 1 + len + m_ou.cols();

  // size the output matrix correctly
  dd.reset(dim, dim, 0.0);

  // variable part
  for (int i = 0; i < len; i++) 
    dd.at(1+i, 1+i) = 1.0;

  // old updates part
  for (int l = 0; l < m_ou.cols(); l++) {
    for (int k = 0; k < m_ou.cols(); k++) {
      dd.at(1+len+k, 1+len+l) = m_ou_dd.at(x).at(k,l);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that solves each block's eigenvalue problem
///
/// \param[in]   nkps          number of directions to keep
/// \param[in]   shift_i       identity shift
/// \param[in]   shift_s       overlap shift
/// \param[in]   omega         harmonic davidson energy shift
/// \param[in]   shift_scale   scale of each shift
/// \param[out]  block_ups     updates directions
/// \param[out]  output        output stream
/// \param[in]   iterative     whether to use iterative method 
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlocker::solve_for_block_dirs(const formic::VarDeps * dep_ptr,
                                                   const int nkps,
                                                   const double shift_i,
                                                   const double shift_s,
                                                   const std::vector<double> & shift_scale,
                                                   std::vector<std::vector<formic::Matrix<double> > > & block_ups,
                                                   std::ostream & output,
                                                   const double omega) {

  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();
  
  // clear eigenvectors
  m_ou_dd.clear();

  // on root process
  if ( my_rank == 0 ) {

    // loop over blocks
    for (int b = 0; b < m_hdata.nb(); b++) {
      // block begining index
      const int starti = this->bb(b);
      // block end index
      const int stopi  = this->be(b);
    
      m_ou_dd.push_back(formic::Matrix<double>(m_ou.cols(), m_ou.cols(), 0.0));
      for (int l = 0; l < m_ou.cols(); l++) {
        for (int k = 0; k < m_ou.cols(); k++) {
          for (int i = starti; i < stopi; i++) {
            m_ou_dd.rbegin()->at(k,l) += m_ou.at(i,k) * m_ou.at(i,l);
          }
        }
      }
    }
  }

  // initialize update direction vectors
  block_ups.assign(m_hdata.nb(), std::vector<formic::Matrix<double> >());

  // solve for each block on root process
  // loop over blocks
  for (int b = 0; b < m_hdata.nb(); b++) {
    
    // loop over shifts
    for (auto s = shift_scale.begin(); s != shift_scale.end(); s++) {
      
      // begining index of this block
      const int ibeg = 1 + m_hdata.bb(b);

      // block length
      const int len = m_hdata.bl(b);

      // 
      formic::Matrix<double> block_target_vecs(len, m_hdata.nb()-1);

      for (int x = 0, y = 0; x < m_hdata.nb(); x++) {
        if (x == b)
          continue;
        
        // construct block matrix
        formic::Matrix<double> hh, ss, dd, up;

        // build the matrices on root process
        if ( my_rank == 0 ) {

          // hamiltonian
          m_hdata.prep_lm_block_plus_other_ou_matrix(b, x, hh);

          // overlap
          m_sdata.prep_lm_block_plus_other_ou_matrix(b, x, ss);

          // old updates shifts
          this->prep_lm_block_plus_other_ou_dd_matrix(b, x, dd);
        }

        if ( !_iterative ) {
          if ( my_rank == 0 ) 
            up = cqmc::engine::get_important_brlm_dirs(1, this->avg_e(), shift_i * (*s), shift_s * (*s), 1.0e-4, hh, ss, dd, output);
        }
        else  
          up = cqmc::engine::get_important_brlm_dirs_davidson(dep_ptr, 1, omega, this->avg_e(), shift_i * (*s), shift_s * (*s), 1.0e-6, hh, ss, dd, output);

        //copy the resulting update direction to matrix
        if ( my_rank == 0 ) 
          std::copy(up.begin()+1, up.begin()+1+len, block_target_vecs.col_begin(y));
        y++;
      }

      // find the most important direction
      if ( my_rank == 0 ) {
        formic::Matrix<double> u, vt;
        formic::ColVec<double> sig;
        block_target_vecs.svd(u, sig,vt);

        // copy this direction to output matrix
        block_ups.at(b).push_back(formic::Matrix<double>(len, std::min(nkps, int(u.cols()))) <<= u.begin());
      }
    }
  }

  // create a vector containing number of matrix and matrix size for each block
  std::vector<std::vector<std::vector<int> > > mat_size(m_hdata.nb());
  // get matrix size on root process
  for (int b = 0; b < m_hdata.nb(); b++) {
    mat_size.at(b).resize(shift_scale.size());
    for (int s = 0; s < shift_scale.size(); s++) {
      mat_size.at(b).at(s).resize(2);
      if ( my_rank == 0 ) {
        mat_size.at(b).at(s).at(0) = block_ups.at(b).at(s).rows();
        mat_size.at(b).at(s).at(1) = block_ups.at(b).at(s).cols();
      }
    }
  }

  // broadcast thie vector to all process
  for (int b = 0; b < m_hdata.nb(); b++) {
    for (int s = 0; s < shift_scale.size(); s++) {
      formic::mpi::bcast(&mat_size.at(b).at(s).at(0), 2);
      //MPI_Bcast(&mat_size.at(b).at(s).at(0), 2, MPI_INT, 0, MPI_COMM_WORLD);
    }
  }

  // size the matrix on non-root process correctly
  if ( my_rank != 0 ) {
    for (int b = 0; b < m_hdata.nb(); b++) {
      block_ups.at(b).resize(shift_scale.size());
      for ( int s = 0; s < shift_scale.size(); s++) {
        formic::Matrix<double> temp_mat(mat_size.at(b).at(s).at(0), mat_size.at(b).at(s).at(1), 0.0);
        block_ups.at(b).at(s) = temp_mat;
      }
    }
  }

  // broadcast the solve results for each block and each shift to all processes
  for (int b = 0; b < m_hdata.nb(); b++) {
    for (int s = 0; s < shift_scale.size(); s++) {
      formic::mpi::bcast(&block_ups.at(b).at(s).at(0,0), block_ups.at(b).at(s).size());
      //MPI_Bcast(&block_ups.at(b).at(s).at(0,0), block_ups.at(b).at(s).size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }
}


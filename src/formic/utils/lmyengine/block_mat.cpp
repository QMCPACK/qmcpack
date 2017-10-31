/////////////////////////////////////////////////////////////////////////////////////////////////
/// \file  /utils/lmyengine/block_mat.cpp
///
/// \brief  Implementation file for block recursive matrix class
///
/////////////////////////////////////////////////////////////////////////////////////////////////

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
#include<formic/utils/lmyengine/block_mat.h>
#include<formic/utils/lmyengine/block_detail.h>

//////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that contracts each block's previous update components with derivative ratios
///
/// \param[in]   dr        derivative ratio(bare or energy) for this sample
/// \param[in]   ou_mat    matrix holds old updates, size num_old_update * num_variables
/// \param[out]  cont_vec  contracted results
///
//////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::prep_block_ou_contractions(const std::vector<double> & dr, 
                                                                const formic::Matrix<double> & ou_mat, 
                                                                std::vector<formic::ColVec<double> > & cont_vec) {
  
  // loop over blocks
  for (int b = 0; b < this->nb(); b++) {
    
    const int ibeg = 1 + m_block_beg.at(b);
    const int len = m_block_len.at(b);

    formic::ColVec<double> & cont = cont_vec.at(b);
    for (int i = 0; i < cont.size(); i++)
      cont.at(i) = 0.0;

    for (int k = 0; k < m_nou; k++) {
      for (int j = ibeg; j < ibeg+len; j++) {
        cont.at(k) += dr.at(j) * ou_mat.at(j-1,k);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that create a matrix in the basis of a block of variables and older updates 
///         from a different block
///
/// \param[in]   b     block index
/// \param[in]   x     block index
/// \param[out]  mat   output matrix
///
////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::prep_lm_block_plus_other_ou_matrix(const int b, const int x, formic::Matrix<double> & mat) {

  // first check that if x and b is the same block
  if ( x == b ) 
    throw formic::Exception("b and x must be different in prep_lm_block_plus_other_ou_matrix");

  // get the length of this block
  const int len = this->bl(b);

  // get the dimension of this matrix 
  const int dim = 1 + len + m_nou;

  const int y = x - ( x > b ? 1 : 0 ); 
  
  // size the output matrix correctly
  mat.reset(dim, dim);

  // <wfn|wfn> element
  mat.at(0,0) = m_ww[0]; 
  //std::cout << "in prep_lm_block_plus_other_ou_matrix mat(0,0) is " << m_ww << std::endl;

  // <wfn|var>
  for (int i = 0; i < len; i++) 
    mat.at(0, 1+i) = m_wv[0].at(b).at(i);

  // <var|wfn>
  for (int i = 0; i < len; i++)
    mat.at(1+i, 0) = m_vw[0].at(b).at(i);

  // <var|var>
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      mat.at(1+i,1+j) = m_vv[0].at(b).at(i,j);
    }
  }

  // <wfn|old_update>
  for (int k = 0; k < m_nou; k++) 
    mat.at(0, 1+len+k) = m_wo[0].at(x).at(k);

  // <old_update_wfn>
  for (int k = 0; k < m_nou; k++)
    mat.at(1+len+k, 0) = m_ow[0].at(x).at(k);

  // <var|old_update>
  for (int i = 0; i < len; i++) {
    for (int k = 0; k < m_nou; k++) {
      mat.at(1+i, 1+len+k) = m_vo[0].at(b).at(i, y*m_nou+k);
    }
  }

  // <old_update|var>
  for (int k = 0; k < m_nou; k++) {
    for (int i = 0; i < len; i++) {
      mat.at(1+len+k, 1+i) = m_ov[0].at(b).at(y*m_nou+k, i);
    }
  }

  // <old_update|old_update>
  for (int k = 0; k < m_nou; k++) {
    for (int l = 0; l < m_nou; l++) {
      mat.at(1+len+k, 1+len+l) = m_oo[0].at(x).at(k,l);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that accumulates contribution to block matrices for this sample
///
/// \param[in]   d        weight for this sample
/// \param[in]   lr       left vector for outer product
/// \param[in]   rr       right vector for outer product
/// \param[in]   ou_mat   matrix storing the old update coefficients
///
////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::acc(const double d, const std::vector<double> & lr, const std::vector<double> &rr, const formic::Matrix<double> & ou_mat) {
  
  // get the thread number 
  int myThread = omp_get_thread_num();

  // check to see if the input vector is of the correct size
  if ( lr.size() != 1 + this->tot_rows(m_vv[myThread]) )
    throw formic::Exception("lr vector length was %i but should be %i LMBlockerMatData::acc") % lr.size() % ( 1 + this->tot_rows(m_vv[myThread]) );
  if ( rr.size() != 1 + this->tot_rows(m_vv[myThread]) )
    throw formic::Exception("rr vector length was %i but should be %i LMBlockerMatData::acc") % rr.size() % ( 1 + this->tot_rows(m_vv[myThread]) );

  // prep intermediates for old updates
  this->prep_block_ou_contractions(lr, ou_mat, m_boulr[myThread]);
  this->prep_block_ou_contractions(rr, ou_mat, m_bourr[myThread]);

  // <wfn|wfn>
  m_ww[myThread] += d * lr.at(0) * rr.at(0);

  // loop over blocks
  for (int b = 0; b < this->nb(); b++) {
    
    // begining index
    const int ibeg = 1 + m_block_beg.at(b);

    // block length
    const int len = m_block_len.at(b);

    for (int i = 0; i < len; i++) {

      // <wfn|var>
      m_wv[myThread].at(b).at(i) += d * lr.at(0) * rr.at(ibeg + i);

      // <var|wfn>
      m_vw[myThread].at(b).at(i) += d * lr.at(ibeg + i) * rr.at(0);

    }

    for (int k = 0; k < m_nou; k++) {

      // <wfn|old_updates>
      m_wo[myThread].at(b).at(k) += d * lr.at(0) * m_bourr[myThread].at(b).at(k);
      
      // <old_updates|wfn>
      m_ow[myThread].at(b).at(k) += d * m_boulr[myThread].at(b).at(k) * rr.at(0);
    }

    // <var|var>
    { formic::Matrix<double> & vv_mat = m_vv[myThread].at(b);
      for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
          vv_mat.at(i,j) += d * lr.at(ibeg+i) * rr.at(ibeg+j);
        }
      }
    }

    // <var|old_update> note this is for var in the current block and old-update in all other blocks
    { formic::Matrix<double> & vo_mat = m_vo[myThread].at(b);
      for (int x = 0, y = 0; x < this->nb(); x++) {
        if ( x == b )
          continue;
        const formic::ColVec<double> & rr_vec = m_bourr[myThread].at(x);
        for (int o = 0; o < m_nou; o++) {
          for (int i = 0; i < len; i++) {
            vo_mat.at(i,y*m_nou+o) += d * lr.at(ibeg+i) * rr_vec.at(o);
           }
         } 
         y++;
       } 
     } 

     // <old_update|var> note this is for var in the current block and old-update in all other blocks
     { formic::Matrix<double> & ov_mat = m_ov[myThread].at(b);
       for (int i = 0; i < len; i++) {
         for (int x = 0, y = 0; x < this->nb(); x++) {
           if ( x == b ) 
             continue;
           const formic::ColVec<double> & lr_vec = m_boulr[myThread].at(x);
           for (int k = 0; k < m_nou; k++) {
             ov_mat.at(m_nou*y+k, i) += d * lr_vec.at(k) * rr.at(ibeg+i);
           }
           y++;
         }
       }
     }

     // <old_update|old_update>
     { formic::Matrix<double> & oo_mat = m_oo[myThread].at(b);
       for (int k = 0; k < m_nou; k++) {
         for (int l = 0; l < m_nou; l++) {
           oo_mat.at(k,l) += d * m_boulr[myThread].at(b).at(k) * m_bourr[myThread].at(b).at(l);
         }
       }
     }
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that finalizes the accumulation by dividing matrices total weight
///
////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::finalize(const double total_weight) {
  
  // get the number of threads
  int NumThreads = omp_get_max_threads();

  // sum over threads
  for (int ip = 1; ip < NumThreads; ip++) {
    m_ww[0] += m_ww[ip];
    for (int b = 0; b < this->nb(); b++) {
      m_wv[0].at(b) += m_wv[ip].at(b);
      m_vw[0].at(b) += m_vw[ip].at(b);
      m_wo[0].at(b) += m_wo[ip].at(b);
      m_ow[0].at(b) += m_ow[ip].at(b);
      m_ov[0].at(b) += m_ov[ip].at(b);
      m_vo[0].at(b) += m_vo[ip].at(b);
      m_oo[0].at(b) += m_oo[ip].at(b);
    }
  }
  
  // <wfn|wfn>
  m_ww[0] /= total_weight;

  // loop over blocks
  for (int b = 0; b < this->nb(); b++) {
    m_wv[0].at(b) /= total_weight;
    m_vw[0].at(b) /= total_weight;
    m_vv[0].at(b) /= total_weight;
    m_wo[0].at(b) /= total_weight;
    m_ow[0].at(b) /= total_weight;
    m_ov[0].at(b) /= total_weight;
    m_vo[0].at(b) /= total_weight;
    m_oo[0].at(b) /= total_weight;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that finalizes the accumulation by dividing matrices total weight 
///         across the whole processes
///
////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::mpi_finalize(const double total_weight) {
  
  // get rank number and number of ranks
  int my_rank = formic::mpi::rank();

  // get the number of threads 
  int NumThreads = omp_get_max_threads();

  // loop over blocks
  //for (int b = 0; b < this->nb(); b++) {
  //  // print out to debug
  //  for (int i = 0; i < m_ow.at(b).size(); i++) {
  //    std::cout << m_ow.at(b).at(i) << "   ";
  //  }
  //  std::cout << std::endl;
  //}
  //std::cout << std::endl;

  // sum over threads
  for (int ip = 1; ip < NumThreads; ip++) {
    m_ww[0] += m_ww[ip];
    for (int b = 0; b < this->nb(); b++) {
      m_wv[0].at(b) += m_wv[ip].at(b);
      m_vw[0].at(b) += m_vw[ip].at(b);
      m_vv[0].at(b) += m_vv[ip].at(b);
      m_wo[0].at(b) += m_wo[ip].at(b);
      m_ow[0].at(b) += m_ow[ip].at(b);
      m_ov[0].at(b) += m_ov[ip].at(b);
      m_vo[0].at(b) += m_vo[ip].at(b);
      m_oo[0].at(b) += m_oo[ip].at(b);
    }
  }

  // <wfn|wfn>
  double m_ww_tot = 0.0;
  formic::mpi::reduce(&m_ww[0], &m_ww_tot, 1, MPI_SUM);
  //MPI_Reduce(&m_ww, &m_ww_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  m_ww[0] = m_ww_tot / total_weight;

  // get space for MPI reduce 
  std::vector<formic::ColVec<double> > m_wv_tot, m_vw_tot, m_wo_tot, m_ow_tot;
  std::vector<formic::Matrix<double> > m_vv_tot, m_vo_tot, m_ov_tot, m_oo_tot;
  for (int i = 0; i < this->nb(); i++) {
    m_wv_tot.push_back(formic::ColVec<double>(this->bl(i), 0.0));
    m_vw_tot.push_back(formic::ColVec<double>(this->bl(i), 0.0));
    m_vv_tot.push_back(formic::Matrix<double>(this->bl(i), this->bl(i), 0.0));
    m_wo_tot.push_back(formic::ColVec<double>(m_nou, 0.0));
    m_ow_tot.push_back(formic::ColVec<double>(m_nou, 0.0));
    m_ov_tot.push_back(formic::Matrix<double>((this->nb()-1)*m_nou, this->bl(i), 0.0));
    m_vo_tot.push_back(formic::Matrix<double>(this->bl(i), (this->nb()-1)*m_nou, 0.0));
    m_oo_tot.push_back(formic::Matrix<double>(m_nou, m_nou, 0.0));
  }

  // do MPI reduce
  for (int i = 0; i < this->nb(); i++) {
    formic::mpi::reduce(&m_wv[0].at(i).at(0), &m_wv_tot.at(i).at(0), this->bl(i), MPI_SUM);
    formic::mpi::reduce(&m_vw[0].at(i).at(0), &m_vw_tot.at(i).at(0), this->bl(i), MPI_SUM);
    formic::mpi::reduce(&m_vv[0].at(i).at(0,0), &m_vv_tot.at(i).at(0,0), m_vv[0].at(i).size(), MPI_SUM);
    formic::mpi::reduce(&m_wo[0].at(i).at(0), &m_wo_tot.at(i).at(0), m_nou, MPI_SUM);
    formic::mpi::reduce(&m_ow[0].at(i).at(0), &m_ow_tot.at(i).at(0), m_nou, MPI_SUM);
    formic::mpi::reduce(&m_ov[0].at(i).at(0,0), &m_ov_tot.at(i).at(0,0), m_ov[0].at(i).size(), MPI_SUM);
    formic::mpi::reduce(&m_vo[0].at(i).at(0,0), &m_vo_tot.at(i).at(0,0), m_vo[0].at(i).size(), MPI_SUM);
    formic::mpi::reduce(&m_oo[0].at(i).at(0,0), &m_oo_tot.at(i).at(0,0), m_oo[0].at(i).size(), MPI_SUM);
  }

  // evaluate the average across all processors
  if ( my_rank == 0 ) {
    
    // loop over blocks
    for (int b = 0; b < this->nb(); b++) {
      m_vw[0].at(b) = m_vw_tot.at(b) / total_weight;
      m_wv[0].at(b) = m_wv_tot.at(b) / total_weight;
      m_vv[0].at(b) = m_vv_tot.at(b) / total_weight;
      m_ow[0].at(b) = m_ow_tot.at(b) / total_weight;
      m_wo[0].at(b) = m_wo_tot.at(b) / total_weight;
      m_vo[0].at(b) = m_vo_tot.at(b) / total_weight;
      m_ov[0].at(b) = m_ov_tot.at(b) / total_weight;
      m_oo[0].at(b) = m_oo_tot.at(b) / total_weight;
      // print out to debug
      //for (int i = 0; i < m_ow.at(b).size(); i++) {
      //  std::cout << m_ow.at(b).at(i) << "   ";
      //}
      //std::cout << std::endl;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that resets the matrices
///
/////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::engine::LMBlockerMatData::reset(const int nv, const int nblock, const int nou) {
  
  m_nou = nou;

  // get block information
  cqmc::engine::brlm_get_block_info(nv, nblock, m_block_beg, m_block_end, m_block_len);

  // get the maximum number of threads 
  int NumThreads = omp_get_max_threads();
  
  m_ww.assign(NumThreads, 0.0);

  m_wv.clear();
  m_wv.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_wv[ip].push_back(formic::ColVec<double>(this->bl(i), 0.0));
    }
  }

  m_vw.clear();
  m_vw.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_vw[ip].push_back(formic::ColVec<double>(this->bl(i), 0.0));
    }
  }

  m_vv.clear();
  m_vv.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_vv[ip].push_back(formic::Matrix<double>(this->bl(i), this->bl(i), 0.0));
    }
  }

  m_wo.clear();
  m_wo.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_wo[ip].push_back(formic::ColVec<double>(m_nou, 0.0));
    }
  }

  m_ow.clear();
  m_ow.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_ow[ip].push_back(formic::ColVec<double>(m_nou, 0.0));
    }
  }

  m_vo.clear();
  m_vo.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_vo[ip].push_back(formic::Matrix<double>(this->bl(i), (this->nb()-1)*m_nou, 0.0));
    }
  }

  m_ov.clear();
  m_ov.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_ov[ip].push_back(formic::Matrix<double>((this->nb()-1)*m_nou, this->bl(i), 0.0));
    }
  }

  m_oo.clear();
  m_oo.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_oo[ip].push_back(formic::Matrix<double>(m_nou, m_nou, 0.0));
    }
  }

  m_boulr.clear();
  m_bourr.clear();
  m_boulr.resize(NumThreads);
  m_bourr.resize(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++) {
    for (int i = 0; i < this->nb(); i++) {
      m_boulr[ip].push_back(formic::ColVec<double>(m_nou, 0.0));
      m_bourr[ip].push_back(formic::ColVec<double>(m_nou, 0.0));
    }
  }
}


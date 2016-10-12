//////////////////////////////////////////////////////////////////////////////
/// \file  formic/utils/lmyengine/block_detai.h
///
/// \brief  header file for block-recursive algorithm functions
///
//////////////////////////////////////////////////////////////////////////////

#ifndef ENGINE_BLOCK_FUNC_HEADER
#define ENGINE_BLOCK_FUNC_HEADER

#include<iostream>
#include<vector>

#include<formic/utils/matrix.h>

namespace formic { 
  class VarDeps; 
}

namespace cqmc {

  namespace engine {

    void brlm_get_block_info(const int nvar, const int nblock, std::vector<int> & block_beg, std::vector<int> & block_end, std::vector<int> & block_len);

    formic::Matrix<double> get_important_brlm_dirs(const int nkeep,
                                                   const double curr_e,
                                                   const double shift_i,
                                                   const double shift_s,
                                                   const double ovl_thresh,
                                                   const formic::Matrix<double> & hh,
                                                   const formic::Matrix<double> & ss,
                                                   const formic::Matrix<double> & dd,
                                                   std::ostream & output); 

    formic::Matrix<double> get_important_brlm_dirs_davidson(const formic::VarDeps * dep_ptr,
                                                            const int nkeep,
                                                            const double omega,
                                                            const double curre_cost,
                                                            const double shift_i,
                                                            const double shift_s,
                                                            const double ovl_thresh,
                                                            formic::Matrix<double> & hh,
                                                            formic::Matrix<double> & ss,
                                                            formic::Matrix<double> & dd,
                                                            std::ostream & output);

  }
}

#endif

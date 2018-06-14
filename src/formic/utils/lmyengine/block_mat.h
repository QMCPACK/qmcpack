//////////////////////////////////////////////////////////////////////////////////////////////////
/// \file  /utils/lmyengine/block_mat.h       
///
/// \brief  Header file for the block recursive matrix class
///
//////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_RECURSIVE_MAT_HEADER
#define BLOCK_RECURSIVE_MAT_HEADER

#include<vector>
#include<list>
#include<string>
#include<numeric>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<sstream>

#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>

#include<formic/utils/matrix.h>

namespace cqmc {

namespace engine {

class LMBlockerMatData {
  
  protected:
    
    // number of old updates 
    int m_nou;

    // a vector holds the block begin index for each block
    std::vector<int> m_block_beg;

    // a vector holds the block end index for each block
    std::vector<int> m_block_end;

    // a vector holds the block length for each block
    std::vector<int> m_block_len;

    // matrix elemnets of <wfn|wfn>
    std::vector<double> m_ww;

    // matrix elements of <wfn|var>
    std::vector<std::vector<formic::ColVec<double> > > m_wv;

    // matrix elements of <var|wfn>
    std::vector<std::vector<formic::ColVec<double> > > m_vw;

    // matrix elements of <wfn|old_update>
    std::vector<std::vector<formic::ColVec<double> > > m_wo;

    // matrix elements of <old_update|wfn>
    std::vector<std::vector<formic::ColVec<double> > > m_ow;

    // matrix elements of <var|var>
    std::vector<std::vector<formic::Matrix<double> > > m_vv;

    // matrix elements of <var|old_update>
    std::vector<std::vector<formic::Matrix<double> > > m_vo;

    // matrix elements of <old_update|var>
    std::vector<std::vector<formic::Matrix<double> > > m_ov;

    // matrix elements of <old_update|old_update>
    std::vector<std::vector<formic::Matrix<double> > > m_oo;

    // used to hold contractions of each block's component of old updates with derivative ratios
    std::vector<std::vector<formic::ColVec<double> > > m_boulr;

    // used to hold contractions of each block's component of old updates with derivative ratios
    std::vector<std::vector<formic::ColVec<double> > > m_bourr;
    
  protected:

    int tot_rows(const std::vector<formic::Matrix<double> > & vom) {
      int retval = 0;
      for (int i = 0; i < vom.size(); i++)
        retval += vom.at(i).rows();
      return retval;
    }

  public:
    
    // function that returns the number of blocks
    int nb() const { return m_block_beg.size(); }

    // function that returns the begining index of ith block
    int bb(const int i) const { return m_block_beg.at(i); }

    // function that returns the end index of ith block
    int be(const int i) const { return m_block_end.at(i); }

    // function that returns the length of ith block
    int bl(const int i) const { return m_block_len.at(i); }

    // function that returns <wfn|wfn>
    double ww_element() const { return m_ww[0]; }

    // function that reset the object by clear all data arrays
    void reset(const int nv, const int nblock, const int nou);

    // functon that contracts the each block's previous update components with derivatives
    void prep_block_ou_contractions(const std::vector<double> & dr, const formic::Matrix<double> & ou_mat, std::vector<formic::ColVec<double> > & cont_vec);

    // create a matrix in the basis of a block of variables and older updates from a different block
    void prep_lm_block_plus_other_ou_matrix(const int b, const int x, formic::Matrix<double> & mat);

    // accumulate function 
    void acc(const double d, const std::vector<double> & lr, const std::vector<double> & rr, const formic::Matrix<double> & ou_mat);

    // finalize accumulation by dividing each matrix the total weight
    void finalize(const double total_weight);

    // finalize accumulaion by doing mpi reduce
    void mpi_finalize(const double total_weight);
};

}

}

#endif 

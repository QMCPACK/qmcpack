/////////////////////////////////////////////////////////////////////////////////////////////////
/// \file  /utils/lmyengine/block_alg.h
///
/// \brief  Header file for the block recursive algorithm class
///
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_RECURSIVE_ALG_HEADER
#define BLOCK_RECURSIVE_ALG_HEADER

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
#include<formic/utils/lmyengine/block_mat.h>
#include<formic/utils/lmyengine/var_dependencies.h>

namespace cqmc {

namespace engine {

class LMBlocker {
  
  protected:
    
    // Hamiltonian data
    cqmc::engine::LMBlockerMatData m_hdata;

    // overlap data
    cqmc::engine::LMBlockerMatData m_sdata;

    // matrix storing the old updates coefficients
    formic::Matrix<double> m_ou;

    // directions for each block
    std::vector<formic::Matrix<double> >  m_ou_dd;

    // total weight
    double m_tw;

    // whether this calculation is ground or excited 
    bool _ground;

    // whether to solve the eigenvalue problem iteratively
    bool _iterative;

  public:
    
    // default constructor
    LMBlocker() {}

    // function that returns the number of blocks
    int nb() const { return m_hdata.nb(); }

    // function that returns the begining index of ith block
    int bb(const int i) const { return m_hdata.bb(i); }

    // function that returns the end index of ith block
    int be(const int i) const { return m_hdata.be(i); }

    // function that returns the length of ith block
    int bl(const int i) const { return m_hdata.bl(i); }

    // function that returns the old updates matrix
    formic::Matrix<double> & ou_mat() { return m_ou; }

    // function that returns whether we use iterative method
    bool iterative() const { return _iterative; }

    // function that resets the object
    void reset(const int nv, const int nblock, const std::vector<formic::ColVec<double> > & ou, const bool ground = true, const bool iterative = false); 

    // function that returns the average of local energy
    double avg_e() const { return m_hdata.ww_element() / m_sdata.ww_element(); }

    // function that returns the total weight
    double total_weight() const { return m_tw; }

    // function that accumlates data for Hamiltonian and overlap
    void acc(const double d, const std::vector<double> & dr, const std::vector<double> & er, const bool ground, const double hd_lm_shift=0.0);

    // function that finalizes the data accumulation
    void finalize();

    // function that finalize the data accumulation across all processors
    void mpi_finalize(const double total_weight);

    // function that
    void prep_lm_block_plus_other_ou_dd_matrix(const int b, const int x, formic::Matrix<double> & dd);

    // function that solves the update direction for each block
    void solve_for_block_dirs(const formic::VarDeps * dep_ptr,
                              const int nkps, 
                              const double shift_i,
                              const double shift_s,
                              const std::vector<double> & shift_scale,
                              std::vector<std::vector<formic::Matrix<double> > > & block_ups,
                              std::ostream & output,
                              const double omega=0.0);
  };
}
}

#endif

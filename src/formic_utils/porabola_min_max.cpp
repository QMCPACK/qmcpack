///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/porabola_min_max.cpp
///
/// \brief   implementation for the porabola_min_max function
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>

#include<formic/utils/numeric.h>
#include<formic/utils/lapack_interface.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Finds the minima or maxima of the porabola defined by the three supplied points
///
/// \param[in]      xvec     length 3.  the x coordinates of the points
/// \param[in]      yvec     length 3.  the y coordinates of the points
///
///////////////////////////////////////////////////////////////////////////////////////////////////
double formic::porabola_min_max(const double * const xvec, const double * const yvec) {

        const int n = 3;
        std::vector<double> lin_reg_mat(n*n, 0.0);
        std::vector<double> lin_reg_vec(n, 0.0);
        std::vector<double> lin_reg_sol(n, 0.0);
        for (int i = 0; i < n; i++) {
          lin_reg_mat.at(i*n+0) = xvec[i] * xvec[i];
          lin_reg_mat.at(i*n+1) = xvec[i];
          lin_reg_mat.at(i*n+2) = 1.0;
          lin_reg_vec.at(i) = yvec[i];
        }
        double d = 0.0;
        std::vector<double> work(n*n, 0.0);
        std::vector<int> iwork(2*n, 0);
        formic::matrix_inverse_lu(n, d, &lin_reg_mat.at(0), &work.at(0), &iwork.at(0));
        formic::xgemm('T', 'N', n, 1, n, 1.0, &lin_reg_mat.at(0), n, &lin_reg_vec.at(0), n, 0.0, &lin_reg_sol.at(0), n);
        return ( -0.5 * lin_reg_sol.at(1) / lin_reg_sol.at(0) );

}

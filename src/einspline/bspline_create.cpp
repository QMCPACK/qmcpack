/////////////////////////////////////////////////////////////////////////////
//  einspline:  a library for creating and evaluating B-splines            //
//  Copyright (C) 2007 Kenneth P. Esler, Jr.                               //
//  Released under the BSD-3-clause license                                //
/////////////////////////////////////////////////////////////////////////////

#include "bspline_create.h"
#include <array>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <stdexcept>
#include <vector>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////       Helper functions for spline creation         ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
void init_sse_data();

void find_coefs_1d_d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////     Double-Precision, Real Creation Routines       ////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void solve_deriv_interp_1d_d(double bands[], double coefs[], int M, int cstride)
{
  // Solve interpolating equations
  // First and last rows are different
  bands[4 * (0) + 1] /= bands[4 * (0) + 0];
  bands[4 * (0) + 2] /= bands[4 * (0) + 0];
  bands[4 * (0) + 3] /= bands[4 * (0) + 0];
  bands[4 * (0) + 0] = 1.0;
  bands[4 * (1) + 1] -= bands[4 * (1) + 0] * bands[4 * (0) + 1];
  bands[4 * (1) + 2] -= bands[4 * (1) + 0] * bands[4 * (0) + 2];
  bands[4 * (1) + 3] -= bands[4 * (1) + 0] * bands[4 * (0) + 3];
  bands[4 * (0) + 0] = 0.0;
  bands[4 * (1) + 2] /= bands[4 * (1) + 1];
  bands[4 * (1) + 3] /= bands[4 * (1) + 1];
  bands[4 * (1) + 1] = 1.0;

  // Now do rows 2 through M+1
  for (int row = 2; row < (M + 1); row++)
  {
    bands[4 * (row) + 1] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 2];
    bands[4 * (row) + 3] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 3];
    bands[4 * (row) + 2] /= bands[4 * (row) + 1];
    bands[4 * (row) + 3] /= bands[4 * (row) + 1];
    bands[4 * (row) + 0] = 0.0;
    bands[4 * (row) + 1] = 1.0;
  }

  // Do last row
  bands[4 * (M + 1) + 1] -= bands[4 * (M + 1) + 0] * bands[4 * (M - 1) + 2];
  bands[4 * (M + 1) + 3] -= bands[4 * (M + 1) + 0] * bands[4 * (M - 1) + 3];
  bands[4 * (M + 1) + 2] -= bands[4 * (M + 1) + 1] * bands[4 * (M) + 2];
  bands[4 * (M + 1) + 3] -= bands[4 * (M + 1) + 1] * bands[4 * (M) + 3];
  bands[4 * (M + 1) + 3] /= bands[4 * (M + 1) + 2];
  bands[4 * (M + 1) + 2] = 1.0;

  coefs[(M + 1) * cstride] = bands[4 * (M + 1) + 3];
  // Now back substitute up
  for (int row = M; row > 0; row--)
    coefs[row * cstride] = bands[4 * (row) + 3] - bands[4 * (row) + 2] * coefs[cstride * (row + 1)];

  // Finish with first row
  coefs[0] = bands[4 * (0) + 3] - bands[4 * (0) + 1] * coefs[1 * cstride] - bands[4 * (0) + 2] * coefs[2 * cstride];
}

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void solve_periodic_interp_1d_d(double bands[], double coefs[], double lastCol[], int M, size_t cstride)
{
  // Now solve:
  // First and last rows are different
  bands[4 * (0) + 2] /= bands[4 * (0) + 1];
  bands[4 * (0) + 0] /= bands[4 * (0) + 1];
  bands[4 * (0) + 3] /= bands[4 * (0) + 1];
  bands[4 * (0) + 1] = 1.0;
  bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 2] * bands[4 * (0) + 0];
  bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 2] * bands[4 * (0) + 3];
  bands[4 * (M - 1) + 2] = -bands[4 * (M - 1) + 2] * bands[4 * (0) + 2];
  lastCol[0]             = bands[4 * (0) + 0];

  for (int row = 1; row < (M - 1); row++)
  {
    bands[4 * (row) + 1] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 2];
    bands[4 * (row) + 3] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 3];
    lastCol[row]         = -bands[4 * (row) + 0] * lastCol[row - 1];
    bands[4 * (row) + 0] = 0.0;
    bands[4 * (row) + 2] /= bands[4 * (row) + 1];
    bands[4 * (row) + 3] /= bands[4 * (row) + 1];
    lastCol[row] /= bands[4 * (row) + 1];
    bands[4 * (row) + 1] = 1.0;
    if (row < (M - 2))
    {
      bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 2] * bands[4 * (row) + 3];
      bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 2] * lastCol[row];
      bands[4 * (M - 1) + 2] = -bands[4 * (M - 1) + 2] * bands[4 * (row) + 2];
    }
  }

  // Now do last row
  // The [2] element and [0] element are now on top of each other
  bands[4 * (M - 1) + 0] += bands[4 * (M - 1) + 2];
  bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 0] * (bands[4 * (M - 2) + 2] + lastCol[M - 2]);
  bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 0] * bands[4 * (M - 2) + 3];
  bands[4 * (M - 1) + 3] /= bands[4 * (M - 1) + 1];
  coefs[M * cstride] = bands[4 * (M - 1) + 3];
  for (int row = M - 2; row >= 0; row--)
    coefs[(row + 1) * cstride] =
        bands[4 * (row) + 3] - bands[4 * (row) + 2] * coefs[(row + 2) * cstride] - lastCol[row] * coefs[M * cstride];

  coefs[0 * cstride]       = coefs[M * cstride];
  coefs[(M + 1) * cstride] = coefs[1 * cstride];
  coefs[(M + 2) * cstride] = coefs[2 * cstride];
}

// On input, bands should be filled with:
// row 0   :  abcdInitial from boundary conditions
// rows 1:M:  basis functions in first 3 cols, data in last
// row M+1 :  abcdFinal   from boundary conditions
// cstride gives the stride between values in coefs.
// On exit, coefs with contain interpolating B-spline coefs
void solve_antiperiodic_interp_1d_d(double bands[], double coefs[], double lastCol[], int M, size_t cstride)
{
  bands[4 * 0 + 0] *= -1.0;
  bands[4 * (M - 1) + 2] *= -1.0;
  // Now solve:
  // First and last rows are different
  bands[4 * (0) + 2] /= bands[4 * (0) + 1];
  bands[4 * (0) + 0] /= bands[4 * (0) + 1];
  bands[4 * (0) + 3] /= bands[4 * (0) + 1];
  bands[4 * (0) + 1] = 1.0;
  bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 2] * bands[4 * (0) + 0];
  bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 2] * bands[4 * (0) + 3];
  bands[4 * (M - 1) + 2] = -bands[4 * (M - 1) + 2] * bands[4 * (0) + 2];
  lastCol[0]             = bands[4 * (0) + 0];

  for (int row = 1; row < (M - 1); row++)
  {
    bands[4 * (row) + 1] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 2];
    bands[4 * (row) + 3] -= bands[4 * (row) + 0] * bands[4 * (row - 1) + 3];
    lastCol[row]         = -bands[4 * (row) + 0] * lastCol[row - 1];
    bands[4 * (row) + 0] = 0.0;
    bands[4 * (row) + 2] /= bands[4 * (row) + 1];
    bands[4 * (row) + 3] /= bands[4 * (row) + 1];
    lastCol[row] /= bands[4 * (row) + 1];
    bands[4 * (row) + 1] = 1.0;
    if (row < (M - 2))
    {
      bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 2] * bands[4 * (row) + 3];
      bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 2] * lastCol[row];
      bands[4 * (M - 1) + 2] = -bands[4 * (M - 1) + 2] * bands[4 * (row) + 2];
    }
  }

  // Now do last row
  // The [2] element and [0] element are now on top of each other
  bands[4 * (M - 1) + 0] += bands[4 * (M - 1) + 2];
  bands[4 * (M - 1) + 1] -= bands[4 * (M - 1) + 0] * (bands[4 * (M - 2) + 2] + lastCol[M - 2]);
  bands[4 * (M - 1) + 3] -= bands[4 * (M - 1) + 0] * bands[4 * (M - 2) + 3];
  bands[4 * (M - 1) + 3] /= bands[4 * (M - 1) + 1];
  coefs[M * cstride] = bands[4 * (M - 1) + 3];
  for (int row = M - 2; row >= 0; row--)
    coefs[(row + 1) * cstride] =
        bands[4 * (row) + 3] - bands[4 * (row) + 2] * coefs[(row + 2) * cstride] - lastCol[row] * coefs[M * cstride];

  coefs[0 * cstride]       = -coefs[M * cstride];
  coefs[(M + 1) * cstride] = -coefs[1 * cstride];
  coefs[(M + 2) * cstride] = -coefs[2 * cstride];
}


namespace
{
// The first three entries multiply the endpoint's three neighboring
// coefficients; the fourth is the prescribed derivative.
std::array<double, 4> derivative_boundary_row(bc_code code, double value, double delta_inv)
{
  switch (code)
  {
  case FLAT:
  case DERIV1:
    return {-0.5 * delta_inv, 0.0, 0.5 * delta_inv, code == FLAT ? 0.0 : value};
  case NATURAL:
  case DERIV2:
    return {delta_inv * delta_inv, -2.0 * delta_inv * delta_inv, delta_inv * delta_inv, code == NATURAL ? 0.0 : value};
  default:
    throw std::invalid_argument("Expected a derivative B-spline boundary condition");
  }
}
} // namespace

void find_coefs_1d_d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride)
{
  int M           = grid.num;
  double basis[4] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 0.0};
  if (bc.lCode == PERIODIC || bc.lCode == ANTIPERIODIC)
  {
    double* bands   = static_cast<double*>(std::malloc(4 * M * sizeof(double)));
    double* lastCol = static_cast<double*>(std::malloc(M * sizeof(double)));

    for (int i = 0; i < M; i++)
    {
      bands[4 * i + 0] = basis[0];
      bands[4 * i + 1] = basis[1];
      bands[4 * i + 2] = basis[2];
      bands[4 * i + 3] = data[i * dstride];
    }

    if (bc.lCode == ANTIPERIODIC)
      solve_antiperiodic_interp_1d_d(bands, coefs, lastCol, M, cstride);
    else
      solve_periodic_interp_1d_d(bands, coefs, lastCol, M, cstride);

    std::free(bands);
    std::free(lastCol);
  }
  else
  {
    const auto abcd_left  = derivative_boundary_row(bc.lCode, bc.lVal, grid.delta_inv);
    const auto abcd_right = derivative_boundary_row(bc.rCode, bc.rVal, grid.delta_inv);

    double* bands = static_cast<double*>(std::malloc((M + 2) * 4 * sizeof(double)));

    for (int i = 0; i < 4; i++)
    {
      bands[4 * (0) + i]     = abcd_left[i];
      bands[4 * (M + 1) + i] = abcd_right[i];
    }
    for (int i = 0; i < M; i++)
    {
      for (int j = 0; j < 3; j++)
        bands[4 * (i + 1) + j] = basis[j];
      bands[4 * (i + 1) + 3] = data[i * dstride];
    }
    // Now, solve for coefficients
    solve_deriv_interp_1d_d(bands, coefs, M, cstride);

    std::free(bands);
  }
}


UBspline_1d_d* create_UBspline_1d_d(Ugrid x_grid, BCtype_d xBC, double* data)
{
  // Create new spline
  UBspline_1d_d* restrict spline = static_cast<UBspline_1d_d*>(std::malloc(sizeof(UBspline_1d_d)));
  spline->spcode                 = U1D;
  spline->tcode                  = DOUBLE_REAL;
  spline->xBC                    = xBC;

  // Setup internal variables
  int M = x_grid.num;
  int N;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
  {
    x_grid.delta = (x_grid.end - x_grid.start) / (double)(x_grid.num);
    N            = M + 3;
  }
  else
  {
    x_grid.delta = (x_grid.end - x_grid.start) / (double)(x_grid.num - 1);
    N            = M + 2;
  }

  x_grid.delta_inv = 1.0 / x_grid.delta;
  spline->x_grid   = x_grid;

#ifndef HAVE_SSE2
  spline->coefs = static_cast<double*>(std::malloc(sizeof(double) * N));
#else
  void* coefs_ptr = nullptr;
  if (posix_memalign(&coefs_ptr, 16, sizeof(double) * N) != 0)
    coefs_ptr = nullptr;
  spline->coefs = static_cast<double*>(coefs_ptr);
#endif
  if (data != nullptr) // only data is provided
    find_coefs_1d_d(spline->x_grid, xBC, data, 1, spline->coefs, 1);

  init_sse_data();
  return spline;
}


UBspline_3d_d* create_UBspline_3d_d(Ugrid x_grid,
                                    Ugrid y_grid,
                                    Ugrid z_grid,
                                    BCtype_d xBC,
                                    BCtype_d yBC,
                                    BCtype_d zBC,
                                    double* data)
{
  // Create new spline
  UBspline_3d_d* restrict spline = static_cast<UBspline_3d_d*>(std::malloc(sizeof(UBspline_3d_d)));
  spline->spcode                 = U3D;
  spline->tcode                  = DOUBLE_REAL;
  spline->xBC                    = xBC;
  spline->yBC                    = yBC;
  spline->zBC                    = zBC;

  // Setup internal variables
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx + 3;
  else
    Nx = Mx + 2;
  x_grid.delta     = (x_grid.end - x_grid.start) / (double)(Nx - 3);
  x_grid.delta_inv = 1.0 / x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My + 3;
  else
    Ny = My + 2;
  y_grid.delta     = (y_grid.end - y_grid.start) / (double)(Ny - 3);
  y_grid.delta_inv = 1.0 / y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
    Nz = Mz + 3;
  else
    Nz = Mz + 2;
  z_grid.delta     = (z_grid.end - z_grid.start) / (double)(Nz - 3);
  z_grid.delta_inv = 1.0 / z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny * Nz;
  spline->y_stride = Nz;

  spline->coefs_size = (size_t)Nx * (size_t)Ny * (size_t)Nz;

#ifndef HAVE_SSE2
  spline->coefs = static_cast<double*>(std::malloc(sizeof(double) * spline->coefs_size));
#else
  void* coefs_ptr = nullptr;
  if (posix_memalign(&coefs_ptr, 16, (sizeof(double) * spline->coefs_size)) != 0)
    coefs_ptr = nullptr;
  spline->coefs = static_cast<double*>(coefs_ptr);
#endif

  if (data != nullptr) // only data is provided
  {
    // First, solve in the X-direction
#pragma omp parallel for
    for (int iy = 0; iy < My; iy++)
      for (int iz = 0; iz < Mz; iz++)
      {
        intptr_t doffset = iy * Mz + iz;
        intptr_t coffset = iy * Nz + iz;
        find_coefs_1d_d(spline->x_grid, xBC, data + doffset, My * Mz, spline->coefs + coffset, Ny * Nz);
      }

    // Now, solve in the Y-direction
#pragma omp parallel for
    for (int ix = 0; ix < Nx; ix++)
      for (int iz = 0; iz < Nz; iz++)
      {
        intptr_t doffset = ix * Ny * Nz + iz;
        intptr_t coffset = ix * Ny * Nz + iz;
        find_coefs_1d_d(spline->y_grid, yBC, spline->coefs + doffset, Nz, spline->coefs + coffset, Nz);
      }

    // Now, solve in the Z-direction
#pragma omp parallel for
    for (int ix = 0; ix < Nx; ix++)
      for (int iy = 0; iy < Ny; iy++)
      {
        intptr_t doffset = (ix * Ny + iy) * Nz;
        intptr_t coffset = (ix * Ny + iy) * Nz;
        find_coefs_1d_d(spline->z_grid, zBC, spline->coefs + doffset, 1, spline->coefs + coffset, 1);
      }
  }
  init_sse_data();
  return spline;
}

namespace
{
// Helper class to factor the cubic B-spline interpolation system once per axis.
// Reuse is possible because all lines share the interior weights (1/6, 2/3,
// 1/6) and the same boundary equations. Periodic and antiperiodic systems have
// corner couplings. Derivative conditions instead add one equation at each end
// without wrapping.
class SplineLineSolver
{
  // Factors for each row except the last, handled separately.
  struct Row
  {
    double inv_diag;          // Reciprocal pivot before normalizing the row.
    double upper;             // Normalized coefficient of the next unknown.
    double wrap        = 0.0; // Cyclic coupling to the final unknown.
    double last_weight = 0.0; // Cyclic elimination weight in the final row.
  };

  static constexpr double a = 1.0 / 6.0; // Off-diagonal interpolation weight.
  static constexpr double b = 2.0 / 3.0; // Main-diagonal interpolation weight.
  const int num_points_;
  const double sign_; // +1 for periodic, -1 for antiperiodic, zero for derivative boundaries.
  std::vector<Row> rows_;
  double last_lower_    = 0.0; // Final row's coupling to the penultimate unknown before elimination.
  double last_inv_diag_ = 0.0; // Reciprocal pivot of the final row after elimination.
  // Additional factors and RHS values from the derivative endpoint equations.
  double left_second_upper_  = 0.0;
  double left_rhs_           = 0.0;
  double right_second_lower_ = 0.0;
  double right_rhs_          = 0.0;

  // Compute the elimination pivot of a tridiagonal interior row with weights
  // (a, b, a) after eliminating the sub-diagonal entry `a` against the
  // normalized upper-diagonal entry `previous_upper` of the preceding row.
  static double interior_pivot(double previous_upper) { return b - a * previous_upper; }

  // Precompute the Gauss elimination factors for a 1D system of size (m + 2)
  // with derivative boundary conditions (FLAT, NATURAL, DERIV1, or DERIV2) at
  // both endpoints. Populates `rows_` (size m + 1) along with the endpoint
  // coupling scalars `left_second_upper_`, `left_rhs_`, `right_second_lower_`,
  // `last_lower_`, `last_inv_diag_`, and `right_rhs_`.
  void factor_derivative(Ugrid grid, BCtype_d bc)
  {
    const auto left  = derivative_boundary_row(bc.lCode, bc.lVal, grid.delta_inv);
    const auto right = derivative_boundary_row(bc.rCode, bc.rVal, grid.delta_inv);
    const int m      = num_points_;
    rows_.resize(m + 1);

    // The left boundary row couples the first three coefficients. Eliminate
    // the first coefficient from the first sample row before the regular sweep.
    rows_[0]                 = {1.0 / left[0], left[1] / left[0]};
    left_second_upper_       = left[2] / left[0];
    left_rhs_                = left[3] / left[0];
    const double first_pivot = interior_pivot(rows_[0].upper);
    rows_[1]                 = {1.0 / first_pivot, (a - a * left_second_upper_) / first_pivot};
    for (int i = 2; i <= m; ++i)
    {
      const double pivot = interior_pivot(rows_[i - 1].upper);
      rows_[i]           = {1.0 / pivot, a / pivot};
    }

    // The right boundary row couples the last three coefficients. Eliminate the
    // first two using their normalized sample rows.
    right_second_lower_ = right[0];
    last_lower_         = right[1] - right_second_lower_ * rows_[m - 1].upper;
    last_inv_diag_      = 1.0 / (right[2] - last_lower_ * rows_[m].upper);
    right_rhs_          = right[3];
  }

  // Solve a 1D derivative-boundary system for a single grid line using the
  // precomputed factors from `factor_derivative`. Reads `m` samples from `data`
  // with stride `dstride` into `rhs` during forward elimination, then writes
  // `m + 2` B-spline coefficients to `coefs` with stride `cstride` via back
  // substitution. Safe when `data` and `coefs` alias.
  void solve_derivative(const double* data,
                        intptr_t dstride,
                        double* coefs,
                        intptr_t cstride,
                        std::vector<double>& rhs) const
  {
    const int m = num_points_;
    rhs[0]      = left_rhs_;
    for (int i = 1; i <= m; ++i)
      rhs[i] = (data[(i - 1) * dstride] - a * rhs[i - 1]) * rows_[i].inv_diag;

    // All samples have been read so back substitute in place.
    double next              = (right_rhs_ - right_second_lower_ * rhs[m - 1] - last_lower_ * rhs[m]) * last_inv_diag_;
    coefs[(m + 1) * cstride] = next;
    for (int i = m; i > 0; --i)
    {
      next               = rhs[i] - rows_[i].upper * next;
      coefs[i * cstride] = next;
    }
    coefs[0] = left_rhs_ - rows_[0].upper * next - left_second_upper_ * coefs[2 * cstride];
  }

  // Precompute the Gauss elimination factors for a cyclic (PERIODIC or
  // ANTIPERIODIC) 1D system of `m` independent equations. Tracks the cyclic
  // corner couplings in each row's `wrap` (coupling to the final unknown) and
  // `last_weight` (fill-in weight accumulated into the final row), populating
  // `rows_` (size m - 1), `last_lower_`, and `last_inv_diag_`.
  void factor_cyclic()
  {
    const int m = num_points_;
    rows_.resize(m - 1);
    // Normalize the first row and eliminate the first unknown from the final row.
    rows_[0]         = {1.0 / b, a / b, sign_ * a / b, sign_ * a};
    double last_diag = b - rows_[0].last_weight * rows_[0].wrap; // Running pivot of the final row.
    double weight    = -rows_[0].last_weight * rows_[0].upper; // Induced coupling in the final row to the next unknown.

    // Walk the remaining rows index i except the last.
    for (int i = 1; i < m - 1; ++i)
    {
      const double diag = interior_pivot(rows_[i - 1].upper);
      rows_[i]          = {1.0 / diag, a / diag, -a * rows_[i - 1].wrap / diag, weight};
      if (i < m - 2)
      {
        last_diag -= weight * rows_[i].wrap;
        weight = -weight * rows_[i].upper;
      }
    }
    // At the penultimate row, the induced coupling joins the original lower
    // diagonal. That row's upper and wrap terms both refer to the final unknown.
    last_lower_    = a + weight;
    last_inv_diag_ = 1.0 / (last_diag - last_lower_ * (rows_.back().upper + rows_.back().wrap));
  }

public:
  // Construct and pre-factor the 1D cubic B-spline interpolation solver for the
  // given axis grid `grid` and boundary conditions `bc`. Dispatches to
  // `factor_cyclic()` for PERIODIC/ANTIPERIODIC or `factor_derivative()` otherwise.
  SplineLineSolver(Ugrid grid, BCtype_d bc)
      : num_points_(grid.num),
        sign_(bc.lCode == PERIODIC           ? 1.0
                  : bc.lCode == ANTIPERIODIC ? -1.0
                                             : 0.0)
  {
    if (sign_ == 0.0)
      factor_derivative(grid, bc);
    else
      factor_cyclic();
  }

  // Allocate a per-thread scratch vector sized to `rows_.size()` (`m + 1` for
  // derivative BCs, `m - 1` for cyclic BCs) to hold the transformed right-hand
  // side during `solve()`; reuse across lines to avoid per-line allocations.
  std::vector<double> make_workspace() const { return std::vector<double>(rows_.size()); }

  // Solve/interpolate one line of `m` samples (`data` with stride `dstride`)
  // into B-spline coefficients (`coefs` with stride `cstride`), using `rhs`
  // from `make_workspace()` as scratch space. Dispatches to `solve_derivative`
  // for derivative BCs, or executes the cyclic forward elimination, back
  // substitution, and periodic/antiperiodic boundary padding (`m + 3` total
  // coefficients). Safe for in-place use (`data == coefs`) because all input
  // samples are consumed during forward elimination before any coefficients are
  // written.
  void solve(const double* data, intptr_t dstride, double* coefs, intptr_t cstride, std::vector<double>& rhs) const
  {
    if (sign_ == 0.0)
    {
      solve_derivative(data, dstride, coefs, cstride, rhs);
      return;
    }

    const int m = num_points_;

    // Apply the stored elimination to the right-hand side. rhs holds the
    // normalized entries for rows 0 through m - 2; last holds the final entry.
    rhs[0]      = data[0] * rows_[0].inv_diag;
    double last = data[(m - 1) * dstride] - rows_[0].last_weight * rhs[0];
    for (int i = 1; i < m - 2; ++i)
    {
      rhs[i] = (data[i * dstride] - a * rhs[i - 1]) * rows_[i].inv_diag;
      last -= rows_[i].last_weight * rhs[i];
    }
    if (m > 2) // For two samples, the first entry is already the penultimate one.
      rhs[m - 2] = (data[(m - 2) * dstride] - a * rhs[m - 3]) * rows_.back().inv_diag;

    // Solve for the final unknown, then walk backward. Keep that value in
    // last for the wrap terms; next holds the most recently solved coefficient.
    last               = (last - last_lower_ * rhs.back()) * last_inv_diag_;
    double next        = last;
    coefs[m * cstride] = last;
    for (int i = m - 2; i >= 0; --i)
    {
      next                     = rhs[i] - rows_[i].upper * next - rows_[i].wrap * last;
      coefs[(i + 1) * cstride] = next;
    }
    // Independent coefficients occupy slots 1 through m. Fill the three
    // boundary slots with signed copies for periodic or antiperiodic evaluation.
    coefs[0]                 = sign_ * last;
    coefs[(m + 1) * cstride] = sign_ * coefs[cstride];
    coefs[(m + 2) * cstride] = sign_ * coefs[2 * cstride];
  }
};
} // namespace

// Recompute coefficients using one factorization per axis for either cyclic
// or derivative boundary conditions.
void recompute_UBspline_3d_d(UBspline_3d_d* spline, double* data)
{
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)
    Nx = Mx + 3;
  else
    Nx = Mx + 2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)
    Ny = My + 3;
  else
    Ny = My + 2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)
    Nz = Mz + 3;
  else
    Nz = Mz + 2;

  // Solve the spline system one axis at a time. The axis solvers share
  // read-only factors across OpenMP iterations; each thread owns its
  // rhs workspace and reuses it across its assigned grid lines. X reads the
  // input samples; Y and Z solve in place in the coefficient table.
  {
    const SplineLineSolver x_solver(spline->x_grid, spline->xBC);
#pragma omp parallel
    {
      auto rhs = x_solver.make_workspace();
#pragma omp for
      for (int iy = 0; iy < My; iy++)
        for (int iz = 0; iz < Mz; iz++)
        {
          intptr_t doffset = iy * Mz + iz;
          intptr_t coffset = iy * Nz + iz;
          x_solver.solve(data + doffset, My * Mz, spline->coefs + coffset, Ny * Nz, rhs);
        }
    }
  }

  { // In-place Y-direction solve
    const SplineLineSolver y_solver(spline->y_grid, spline->yBC);
#pragma omp parallel
    {
      auto rhs = y_solver.make_workspace();
#pragma omp for
      for (int ix = 0; ix < Nx; ix++)
        for (int iz = 0; iz < Nz; iz++)
        {
          intptr_t doffset = ix * Ny * Nz + iz;
          intptr_t coffset = ix * Ny * Nz + iz;
          y_solver.solve(spline->coefs + doffset, Nz, spline->coefs + coffset, Nz, rhs);
        }
    }
  }


  { // In-place Z-direction solve
    const SplineLineSolver z_solver(spline->z_grid, spline->zBC);
#pragma omp parallel
    {
      auto rhs = z_solver.make_workspace();
#pragma omp for
      for (int ix = 0; ix < Nx; ix++)
        for (int iy = 0; iy < Ny; iy++)
        {
          intptr_t doffset = (ix * Ny + iy) * Nz;
          intptr_t coffset = (ix * Ny + iy) * Nz;
          z_solver.solve(spline->coefs + doffset, 1, spline->coefs + coffset, 1, rhs);
        }
    }
  }
}

void destroy_UBspline(Bspline* spline)
{
  std::free(spline->coefs);
  std::free(spline);
}

void destroy_multi_UBspline(Bspline* spline);

void destroy_Bspline(void* spline)
{
  Bspline* sp = static_cast<Bspline*>(spline);
  if (sp->sp_code <= U3D)
    destroy_UBspline(sp);
  else if (sp->sp_code <= MULTI_U3D)
    destroy_multi_UBspline(sp);
  else
  {
    std::fprintf(stderr, "Error in destroy_Bspline:  invalide spline code %d.\n", sp->sp_code);
    std::abort();
  }
}

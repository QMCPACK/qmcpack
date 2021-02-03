#include <CL/sycl.hpp>
#include "QMCWaveFunctions/detail/SYCL/delayed_update_helper.h"
//#include "oneapi/mkl.hpp"
#include "oneapi/mkl/lapack.hpp"
#include "oneapi/mkl/blas.hpp"
namespace sycl = cl::sycl;

int main()
{
  std::vector<std::complex<float>> mat_float(100);
  std::vector<std::complex<double>> mat_double(100);
  std::vector<int> delay_list_gpu(100);

  sycl::queue q;

  copy_matrix_sycl(0, 0, mat_float.data(), 1, mat_double.data(), 0, q);
  extract_matrix_diagonal_sycl(0, mat_double.data(), 0, mat_double.data(), q);
  make_identity_matrix_sycl(0, mat_double.data(), 0, q);
  applyW_stageV_sycl(delay_list_gpu.data(), 0, mat_float.data(), 0, 0, mat_float.data(), mat_float.data(), q);
}

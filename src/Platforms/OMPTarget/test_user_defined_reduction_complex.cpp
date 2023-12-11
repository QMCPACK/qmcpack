#include <vector>
#include <complex>
#include <cassert>

#if !defined(OPENMP_NO_UDR)
#pragma omp declare reduction(+: std::complex<float>: omp_out += omp_in)
#endif

int main()
{
  const int N = 100;
  std::vector<std::complex<float>> array(N);

  auto array_ptr = array.data();
  for (int i = 0; i < N; i++)
    array_ptr[i] = std::complex<float>(i);

  std::complex<float> sum;
  #pragma omp parallel for reduction(+: sum)
  for (int i = 0; i < N; i++)
    sum += array_ptr[i];

  assert(std::real(sum) == 4950);
  assert(std::imag(sum) == 0);
}

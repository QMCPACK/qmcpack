/* boost histogram.cpp graphical verification of distribution functions
 *
 * Copyright Jens Maurer 2000
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * $Id: histogram.cpp,v 1.8 2004/07/27 03:43:34 dgregor Exp $
 *
 * This test program allows to visibly examine the results of the
 * distribution functions.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>
#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsApp/RandomNumberControl.h"
using namespace qmcplusplus;

void plot_histogram(const std::vector<int>& slots, int samples,
                    double from, double to)
{
  int m = *std::max_element(slots.begin(), slots.end());
  const int nRows = 20;
  std::cout.setf(std::ios::fixed|std::ios::left);
  std::cout.precision(5);
  for(int r = 0; r < nRows; r++)
  {
    double y = ((nRows - r) * double(m))/(nRows * samples);
    std::cout << std::setw(10) << y << "  ";
    for(unsigned int col = 0; col < slots.size(); col++)
    {
      char out = ' ';
      if(slots[col]/double(samples) >= y)
        out = 'x';
      std::cout << out;
    }
    std::cout << std::endl;
  }
  std::cout << std::setw(12) << " "
            << std::setw(10) << from;
  std::cout.setf(std::ios::right, std::ios::adjustfield);
  std::cout << std::setw(slots.size()-10) << to << std::endl;
}

// I am not sure whether these two should be in the library as well

// maintain sum of NumberGenerator results
template<class NumberGenerator,
         class Sum = typename NumberGenerator::result_type>
class sum_result
{
public:
  typedef NumberGenerator base_type;
  typedef typename base_type::result_type result_type;
  explicit sum_result(const base_type & g) : gen(g), _sum(0) { }
  result_type operator()()
  {
    result_type r = gen();
    _sum += r;
    return r;
  }
  base_type & base()
  {
    return gen;
  }
  Sum sum() const
  {
    return _sum;
  }
  void reset()
  {
    _sum = 0;
  }
private:
  base_type gen;
  Sum _sum;
};


// maintain square sum of NumberGenerator results
template<class NumberGenerator,
         class Sum = typename NumberGenerator::result_type>
class squaresum_result
{
public:
  typedef NumberGenerator base_type;
  typedef typename base_type::result_type result_type;
  explicit squaresum_result(const base_type & g) : gen(g), _sum(0) { }
  result_type operator()()
  {
    result_type r = gen();
    _sum += r*r;
    return r;
  }
  base_type & base()
  {
    return gen;
  }
  Sum squaresum() const
  {
    return _sum;
  }
  void reset()
  {
    _sum = 0;
  }
private:
  base_type gen;
  Sum _sum;
};


template<class RNG>
void histogram(RNG base, int samples, double from, double to,
               const std::string & name)
{
  typedef squaresum_result<sum_result<RNG, double>, double > SRNG;
  SRNG gen((sum_result<RNG, double>(base)));
  const int nSlots = 60;
  std::vector<int> slots(nSlots,0);
  for(int i = 0; i < samples; i++)
  {
    double val = gen();
    if(val < from || val >= to)    // early check avoids overflow
      continue;
    int slot = int((val-from)/(to-from) * nSlots);
    if(slot < 0 || slot > (int)slots.size())
      continue;
    slots[slot]++;
  }
  std::cout << name << std::endl;
  plot_histogram(slots, samples, from, to);
  double mean = gen.base().sum() / samples;
  std::cout << "mean: " << mean
            << " sigma: " << std::sqrt(gen.squaresum()/samples-mean*mean)
            << "\n" << std::endl;
}

template<class RNG, class SEEDARRAY>
void histogram_OMP(RNG base, int samples, double from, double to,
                   SEEDARRAY& seeds, const std::string & name)
{
  typedef squaresum_result<sum_result<RNG, double>, double > SRNG;
  const int nSlots = 60;
  std::vector<int> slots_tot(nSlots,0);
  double sum1=0.0;
  double sum2=0.0;
  #pragma omp parallel
  {
    std::vector<int> slots(nSlots,0);
    RNG base_copy(base);
    //base_copy.reset();
    base_copy.seed(seeds[omp_get_thread_num()]);
    SRNG gen((sum_result<RNG, double>(base_copy)));
    for(int i = 0; i < samples; i++)
    {
      double val = gen();
      if(val < from || val >= to)    // early check avoids overflow
        continue;
      int slot = int((val-from)/(to-from) * nSlots);
      if(slot < 0 || slot > (int)slots.size())
        continue;
      slots[slot]++;
    }
    #pragma omp critical
    {
      for(int i=0; i<nSlots; i++)
        slots_tot[i]+=slots[i];
      sum1 += gen.base().sum();
      sum2 += gen.squaresum();
      double mean = gen.base().sum() / samples;
      std::cout << "random seed " << seeds[omp_get_thread_num()]
                << " mean: " << mean
                << " sigma: " << std::sqrt(gen.squaresum()/samples-mean*mean)
                << std::endl;
    }
  }
  std::cout << name << std::endl;
  plot_histogram(slots_tot, samples, from, to);
  int samples_tot=samples*omp_get_max_threads();
  double mean = sum1/samples_tot;
  //double mean = gen.base().sum() / samples;
  std::cout << "mean: " << mean
            << " sigma: " << std::sqrt(sum2/samples_tot-mean*mean)
            << "\n" << std::endl;
}

template<class PRNG, class Dist>
inline boost::variate_generator<PRNG&, Dist> make_gen(PRNG & rng, Dist d)
{
  return boost::variate_generator<PRNG&, Dist>(rng, d);
}

template<class PRNG>
void histograms()
{
  PRNG rng;
  using namespace boost;
  //histogram(make_gen(rng, uniform_smallint<>(0, 5)), 100000, -1, 6,
  //          "uniform_smallint(0,5)");
  //histogram(make_gen(rng, uniform_int<>(0, 5)), 100000, -1, 6,
  //          "uniform_int(0,5)");
  histogram(make_gen(rng, uniform_real<>(0,1)), 1000000, -0.5, 1.5,
            "uniform_real(0,1)");
  //histogram(make_gen(rng, bernoulli_distribution<>(0.2)), 100000, -0.5, 1.5,
  //          "bernoulli(0.2)");
  //histogram(make_gen(rng, binomial_distribution<>(4, 0.2)), 100000, -1, 5,
  //          "binomial(4, 0.2)");
  //histogram(make_gen(rng, triangle_distribution<>(1, 2, 8)), 100000, 0, 10,
  //          "triangle(1,2,8)");
  //histogram(make_gen(rng, geometric_distribution<>(5.0/6.0)), 100000, 0, 10,
  //          "geometric(5/6)");
  //histogram(make_gen(rng, exponential_distribution<>(0.3)), 100000, 0, 10,
  //          "exponential(0.3)");
  //histogram(make_gen(rng, cauchy_distribution<>()), 100000, -5, 5,
  //          "cauchy");
  //histogram(make_gen(rng, lognormal_distribution<>(3, 2)), 100000, 0, 10,
  //          "lognormal");
  histogram(make_gen(rng, normal_distribution<>()), 1000000, -3, 3,
            "normal");
  histogram(make_gen(rng, normal_distribution<>(0.5, 0.5)), 1000000, -3, 3,
            "normal(0.5, 0.5)");
  //histogram(make_gen(rng, poisson_distribution<>(1.5)), 100000, 0, 5,
  //          "poisson(1.5)");
  //histogram(make_gen(rng, poisson_distribution<>(10)), 100000, 0, 20,
  //          "poisson(10)");
  //histogram(make_gen(rng, gamma_distribution<>(0.5)), 100000, 0, 0.5,
  //          "gamma(0.5)");
  //histogram(make_gen(rng, gamma_distribution<>(1)), 100000, 0, 3,
  //          "gamma(1)");
  //histogram(make_gen(rng, gamma_distribution<>(2)), 100000, 0, 6,
  //          "gamma(2)");
}

void simple_test()
{
  typedef qmcplusplus::RandomGenerator_t generator_type;
  typedef generator_type::uint_type uint_type;
  generator_type uni;
  uni.seed(static_cast<uint_type>(std::time(0))%1024);
  //copy a generator and seed it by random number from the original generator
  generator_type uni_copy(uni);
  uni_copy.seed(static_cast<uint_type>(uni()*1000));
  for(int i=0; i<20; i++)
    std::cout << uni() << " " << uni_copy() << std::endl;
  //check serial histogram
  histogram(uni, 1000000, -0.5, 1.5, "uniform_real(0,1)");
  //test PrimeNumberSet class
  PrimeNumberSet<uint_type> primes;
  std::vector<uint_type> myprimes;
  //use a random number as an offset
  int n=primes.get(static_cast<uint_type>(std::time(0))%1024,omp_get_max_threads(),myprimes);
  //int n=primes.get(10000,myprimes,5);
  histogram_OMP(uni, 100000, -0.5, 1.5, myprimes, "uniform_real(0,1)");
  for(int i=0; i<myprimes.size(); i++)
    std::cout << i << std::setw(12) << myprimes[i] << std::endl;
}

int main(int argc, char** argv)
{
  //histograms<boost::mt19937>();
  //histograms<boost::lagged_fibonacci607>();
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo Welcome(argc,argv,OHMMS::Controller->rank());
  RandomNumberControl rng;
  rng.put(NULL);
  //check serial histogram
  histogram(Random, 10000, -0.5, 1.5, "uniform_real(0,1)");
  std::vector<uint32_t> myprimes;
  //rng.PrimeNumbers.get(Random.offset()+omp_get_max_threads(),omp_get_max_threads(),myprimes);
  rng.PrimeNumbers.get(1024,512,myprimes);
  std::cout <<"============================  " << std::endl;
  for(int i=0; i<myprimes.size();)
  {
    for(int j=0; j<8; j++, i++)
      std::cout << std::setw(12) << myprimes[i];
    std::cout << std::endl;
  }
  histogram_OMP(Random, 10000, -0.5, 1.5, myprimes, "uniform_real(0,1)");
  std::cout <<"============================  " << std::endl;
  int imax=8*(rng.PrimeNumbers.size()/8);
  for(int i=0; i<imax;)
  {
    for(int j=0; j<8; j++, i++)
      std::cout << std::setw(12) << rng.PrimeNumbers[i];
    std::cout << std::endl;
  }
  //std::stringstream a;
  //Random.write(a);
  //cout << "Size of std::string " << a.str().size() << std::endl;
  //vector<uint32_t> v;
  //uint32_t vt;
  //while(!a.eof())
  //{
  //  if(a>>vt) v.push_back(vt);
  //}
  //for(int i=0; i<v.size(); i++) std::cout << v[i] << std::endl;
  //Random.write(std::cout);
  //cout << std::endl;
  //cout << " size of data " << v.size() << std::endl;
  OHMMS::Controller->finalize();
  return 0;
}


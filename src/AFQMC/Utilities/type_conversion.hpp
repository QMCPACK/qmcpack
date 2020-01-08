////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_TYPE_CONVERSION_HPP 
#define AFQMC_TYPE_CONVERSION_HPP 

namespace qmcplusplus
{
namespace afqmc
{

// convert to single precision
template<typename T>
struct to_single_precision { using value_type = T; };
 
template<>
struct to_single_precision<double> { using value_type = float; };

template<>
struct to_single_precision<std::complex<double>> { using value_type = std::complex<float>; };

template<>
struct to_single_precision<long> { using value_type = int; };

template<>
struct to_single_precision<unsigned long> { using value_type = unsigned int; };

// convert to double precision
template<typename T>
struct to_double_precision { using value_type = T; };

template<>
struct to_double_precision<float> { using value_type = double; };

template<>
struct to_double_precision<std::complex<float>> { using value_type = std::complex<double>; };

template<>
struct to_double_precision<int> { using value_type = long; };

template<>
struct to_double_precision<unsigned int> { using value_type = unsigned long; };


// remove complex
template<typename T>
struct remove_complex { using value_type = T; };

template<>
struct remove_complex<std::complex<float>> { using value_type = float; };

template<>
struct remove_complex<std::complex<double>> { using value_type = double; };

// convert a set of types into strings
template<class T>
inline std::string type_to_string() { return std::string(""); } 

template<>
inline std::string type_to_string<double>() { return std::string("double"); }

template<>
inline std::string type_to_string<int>() { return std::string("int"); }

template<>
inline std::string type_to_string<std::string>() { return std::string("std::string"); }

// create std::vector<T> from std::vector<Q>
template<class Q, class T>
std::vector<Q> make_vector(std::vector<T> const& other)
{
  std::vector<Q> res;
  res.reserve(other.size());
  for(auto& v:other) res.emplace_back(v);
  return res;
}

//   can use: return std::vector<Q>{std::make_move_iterator(other.begin()),
//                                  std::make_move_iterator(other.end()));
// and can replace move_vector with this outside
template<class Q, class T>
std::vector<Q> move_vector(std::vector<T> && other)
{
  std::vector<Q> res;
  res.reserve(other.size());
  for(auto& v:other) res.emplace_back(std::move(v));
  return res;
}

template<class Q, class T, class Aux>
std::vector<Q> move_vector(std::vector<T> && other, Aux param)
{
  std::vector<Q> res;
  res.reserve(other.size());
  for(auto& v:other) res.emplace_back(std::move(v),param);
  return res;
}

}

}

#endif

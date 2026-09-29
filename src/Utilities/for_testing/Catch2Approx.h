//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CATCH2_APPROX_H
#define QMCPLUSPLUS_CATCH2_APPROX_H

#include <catch2/catch_approx.hpp>

#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

#include <config.h>

class Approx : public Catch::Approx
{
public:
  template<typename T>
  explicit Approx(const T& value) : Catch::Approx(value)
  {
    scale(1.0);
  }

  friend std::ostream& operator<<(std::ostream& os, const Approx& value)
  {
    return os << value.toString();
  }
};

class ComplexApprox
{
public:
  explicit ComplexApprox(const std::complex<double>& value) : value_(value) {}
  explicit ComplexApprox(const std::complex<float>& value) : value_(value) {}
  explicit ComplexApprox(double value) : value_(value) {}
  explicit ComplexApprox(float value) : value_(value) {}

  ComplexApprox& epsilon(double new_epsilon)
  {
    epsilon_ = new_epsilon;
    return *this;
  }

  double epsilon() const { return epsilon_; }

  friend bool operator==(const std::complex<double>& lhs, const ComplexApprox& rhs) { return rhs.compare(lhs); }

  friend bool operator==(const std::complex<float>& lhs, const ComplexApprox& rhs)
  {
    return rhs.compare(std::complex<double>(lhs));
  }

  friend bool operator==(const ComplexApprox& lhs, const std::complex<double>& rhs) { return rhs == lhs; }
  friend bool operator==(const ComplexApprox& lhs, const std::complex<float>& rhs) { return rhs == lhs; }

  friend std::ostream& operator<<(std::ostream& os, const ComplexApprox& value)
  {
    return os << "ComplexApprox( " << value.value_ << " )";
  }

private:
  bool compare(const std::complex<double>& other) const
  {
    const double distance  = std::abs(other - value_);
    const double magnitude = std::abs(value_);
    return distance == Catch::Approx(0.0).epsilon(0.0).margin(epsilon_ * (1.0 + magnitude));
  }

  std::complex<double> value_;
  double epsilon_ = std::numeric_limits<float>::epsilon() * 100;
};

class LogComplexApprox
{
public:
  explicit LogComplexApprox(const std::complex<double>& value) : value_(value) {}
  explicit LogComplexApprox(const std::complex<float>& value) : value_(value) {}
  explicit LogComplexApprox(double value) : value_(value) {}
  explicit LogComplexApprox(float value) : value_(value) {}

  LogComplexApprox& epsilon(double new_epsilon)
  {
    epsilon_ = new_epsilon;
    return *this;
  }

  double epsilon() const { return epsilon_; }

  friend bool operator==(const std::complex<double>& lhs, const LogComplexApprox& rhs)
  {
    return std::exp(lhs) == ComplexApprox(std::exp(rhs.value_)).epsilon(rhs.epsilon_);
  }

  friend bool operator==(const std::complex<float>& lhs, const LogComplexApprox& rhs)
  {
    return std::exp(lhs) == ComplexApprox(std::exp(rhs.value_)).epsilon(rhs.epsilon_);
  }

  friend bool operator==(const LogComplexApprox& lhs, const std::complex<double>& rhs) { return rhs == lhs; }
  friend bool operator==(const LogComplexApprox& lhs, const std::complex<float>& rhs) { return rhs == lhs; }

  friend std::ostream& operator<<(std::ostream& os, const LogComplexApprox& value)
  {
    return os << "LogComplexApprox( " << value.value_ << " )";
  }

private:
  std::complex<double> value_;
  double epsilon_ = std::numeric_limits<float>::epsilon() * 100;
};

#ifdef QMC_COMPLEX
using ValueApprox = ::ComplexApprox;
#else
using ValueApprox = ::Approx;
#endif

#endif

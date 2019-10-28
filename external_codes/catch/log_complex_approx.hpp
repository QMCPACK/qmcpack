#ifndef CATCH_LOG_COMPLEX_APPROX
#define CATCH_LOG_COMPLEX_APPROX

#include <complex>
#include <cmath>

// Copy and modify the ComplexApprox class to handle complex numbers for log(complex)

namespace Catch {
class LogComplexApprox
{
    std::complex<double> m_value;
    double m_epsilon;

    bool approx_compare(const double lhs, const double rhs) const
    {
        return Approx(lhs).epsilon(m_epsilon) == rhs;
    }

    bool approx_compare_polar(const double lhs, const double rhs) const
    {
        return Approx(std::sin(lhs)).epsilon(m_epsilon) == std::sin(rhs) &&
               Approx(std::cos(lhs)).epsilon(m_epsilon) == std::cos(rhs);
    }

public:

    LogComplexApprox(const std::complex<double> &value) : m_value(value) {
      init_epsilon();
    }

    LogComplexApprox(const std::complex<float> &value) : m_value(value) {
      init_epsilon();
    }

    LogComplexApprox(const double &value) : m_value(value) {
      init_epsilon();
    }

    LogComplexApprox(const float &value) : m_value(value) {
      init_epsilon();
    }

   void init_epsilon() {
      // Copied from catch.hpp - would be better to copy it from Approx object
      m_epsilon = std::numeric_limits<float>::epsilon()*100;
    }

    LogComplexApprox& epsilon(double new_epsilon)
    {
      m_epsilon = new_epsilon;
      return *this;
    }

    double epsilon() const
    {
      return m_epsilon;
    }

    friend bool operator == (std::complex<double> const& lhs, LogComplexApprox const& rhs)
    {
        return rhs.approx_compare(lhs.real(), rhs.m_value.real()) &&
               rhs.approx_compare_polar(lhs.imag(), rhs.m_value.imag());
    }

    friend bool operator == (std::complex<float> const& lhs, LogComplexApprox const& rhs)
    {
        return rhs.approx_compare(lhs.real(), rhs.m_value.real()) &&
               rhs.approx_compare_polar(lhs.imag(), rhs.m_value.imag());
    }

    friend bool operator == (LogComplexApprox const &lhs, std::complex<double> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (LogComplexApprox const &lhs, std::complex<float> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    std::string toString() const {
        std::ostringstream oss;
        oss <<"LogComplexApprox( " << ::Catch::Detail::stringify(m_value) << " )";
        return oss.str();
    }

    friend std::ostream& operator << ( std::ostream& os, LogComplexApprox const& ca )
    {
       os << ca.toString();
       return os;
    }
};

template<>
struct StringMaker<LogComplexApprox> {
  static std::string convert(LogComplexApprox const &value);
};

#ifdef CATCH_IMPL
std::string StringMaker<LogComplexApprox>::convert(LogComplexApprox const& value)
{
  return value.toString();
}
#endif
}

using Catch::LogComplexApprox;

#endif

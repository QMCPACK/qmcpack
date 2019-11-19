#ifndef CATCH_COMPLEX_APPROX
#define CATCH_COMPLEX_APPROX

#include <complex>

// Copy and modify the Approx class to handle complex numbers

namespace Catch {
class ComplexApprox
{
    std::complex<double> m_value;
    double m_epsilon;

    bool approx_compare(const double lhs, const double rhs) const
    {
        return Approx(lhs).epsilon(m_epsilon) == rhs;
    }

public:
    ComplexApprox(const std::complex<double> &value) : m_value(value) {
      init_epsilon();
    }

    ComplexApprox(const std::complex<float> &value) : m_value(value) {
      init_epsilon();
    }

    ComplexApprox(const double &value) : m_value(value) {
      init_epsilon();
    }

    ComplexApprox(const float &value) : m_value(value) {
      init_epsilon();
    }

   void init_epsilon() {
      // Copied from catch.hpp - would be better to copy it from Approx object
      m_epsilon = std::numeric_limits<float>::epsilon()*100;
    }

    ComplexApprox& epsilon(double new_epsilon)
    {
      m_epsilon = new_epsilon;
      return *this;
    }

    double epsilon() const
    {
      return m_epsilon;
    }

    friend bool operator == (std::complex<double> const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = rhs.approx_compare(lhs.real(), rhs.m_value.real());
        is_equal &= rhs.approx_compare(lhs.imag(), rhs.m_value.imag());
        return is_equal;
    }

    friend bool operator == (std::complex<float> const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = rhs.approx_compare(lhs.real(), rhs.m_value.real());
        is_equal &= rhs.approx_compare(lhs.imag(), rhs.m_value.imag());
        return is_equal;
    }

    friend bool operator == (ComplexApprox const &lhs, std::complex<double> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (ComplexApprox const &lhs, std::complex<float> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    std::string toString() const {
        std::ostringstream oss;
        oss <<"ComplexApprox( " << ::Catch::Detail::stringify(m_value) << " )";
        return oss.str();
    }

    friend std::ostream& operator << ( std::ostream& os, ComplexApprox const& ca )
    {
       os << ca.toString();
       return os;
    }
};

template<>
struct StringMaker<ComplexApprox> {
  static std::string convert(ComplexApprox const &value);
};

#ifdef CATCH_IMPL
std::string StringMaker<ComplexApprox>::convert(ComplexApprox const& value)
{
  return value.toString();
}
#endif
}

using Catch::ComplexApprox;

#endif

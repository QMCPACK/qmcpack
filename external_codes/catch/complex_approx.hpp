#ifndef CATCH_COMPLEX_APPROX
#define CATCH_COMPLEX_APPROX

#include <complex>

// Copy and modify the Approx class to handle complex numbers

namespace Catch {
class ComplexApprox
{
public:
    ComplexApprox(std::complex<double> value) : m_value(value), m_compare_real_only(false) {
      // Copied from catch.hpp - would be better to copy it from Approx object
      m_epsilon = std::numeric_limits<float>::epsilon()*100;
    }

    std::complex<double> m_value;
    bool m_compare_real_only;
    double m_epsilon;

    bool approx_compare(const double lhs, const double rhs) const
    {
        return Approx(lhs).epsilon(m_epsilon) == rhs;
    }

    friend bool operator == (double const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = rhs.approx_compare(lhs, rhs.m_value.real());
        if (!rhs.m_compare_real_only)
        {
          is_equal &= rhs.approx_compare(0.0, rhs.m_value.imag());
        }
        return is_equal;
    }

    friend bool operator == (ComplexApprox const& lhs, double const &rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (std::complex<double> const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = rhs.approx_compare(lhs.real(), rhs.m_value.real());
        if (!rhs.m_compare_real_only)
        {
          is_equal &= rhs.approx_compare(lhs.imag(), rhs.m_value.imag());
        }
        return is_equal;
    }

    friend bool operator == (std::complex<float> const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = rhs.approx_compare(lhs.real(), rhs.m_value.real());
        if (!rhs.m_compare_real_only)
        {
          is_equal &= rhs.approx_compare(lhs.imag(), rhs.m_value.imag());
        }
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

    ComplexApprox &compare_real_only()
    {
      m_compare_real_only = true;
      return *this;
    }

    ComplexApprox &epsilon(double new_epsilon)
    {
      m_epsilon = new_epsilon;
      return *this;
    }

    double epsilon() const
    {
      return m_epsilon;
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

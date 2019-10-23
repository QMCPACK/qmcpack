#ifndef CATCH_POLAR_APPROX
#define CATCH_POLAR_APPROX

#include <complex>
#include <cmath>

// Copy and modify the ComplexApprox class to handle complex numbers for polar coordinates

namespace Catch {
class PolarApprox : public ComplexApprox
{
public:
    using ComplexApprox::ComplexApprox;

    bool approx_compare_polar(const double lhs, const double rhs) const
    {
        return Approx(std::sin(lhs)).epsilon(m_epsilon) == std::sin(rhs) &&
               Approx(std::cos(lhs)).epsilon(m_epsilon) == std::cos(rhs);
    }

    friend bool operator == (std::complex<double> const& lhs, PolarApprox const& rhs)
    {
        return rhs.approx_compare(lhs.real(), rhs.m_value.real()) &&
               rhs.approx_compare_polar(lhs.imag(), rhs.m_value.imag());
    }

    friend bool operator == (std::complex<float> const& lhs, PolarApprox const& rhs)
    {
        return rhs.approx_compare(lhs.real(), rhs.m_value.real()) &&
               rhs.approx_compare_polar(lhs.imag(), rhs.m_value.imag());
    }

    friend bool operator == (PolarApprox const &lhs, std::complex<double> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (PolarApprox const &lhs, std::complex<float> const& rhs)
    {
        return operator==( rhs, lhs );
    }

    std::string toString() const {
        std::ostringstream oss;
        oss <<"PolarApprox( " << ::Catch::Detail::stringify(m_value) << " )";
        return oss.str();
    }

    friend std::ostream& operator << ( std::ostream& os, PolarApprox const& ca )
    {
       os << ca.toString();
       return os;
    }
};

template<>
struct StringMaker<PolarApprox> {
  static std::string convert(PolarApprox const &value);
};

#ifdef CATCH_IMPL
std::string StringMaker<PolarApprox>::convert(PolarApprox const& value)
{
  return value.toString();
}
#endif
}

using Catch::PolarApprox;

#endif

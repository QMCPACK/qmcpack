#ifndef CATCH_COMPLEX_APPROX
#define CATCH_COMPLEX_APPROX

#include <complex>

// Copy and modify the Approx class to handle complex numbers

namespace Catch {
class ComplexApprox
{
public:
    ComplexApprox(std::complex<double> value) : m_value(value), m_compare_real_only(false) {}
    std::complex<double> m_value;
    bool m_compare_real_only;

    friend bool operator == (double const& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = Approx(lhs) == rhs.m_value.real();
        if (!rhs.m_compare_real_only)
        {
          is_equal &= Approx(0.0) == rhs.m_value.imag();
        }
        return is_equal;
    }

    friend bool operator == (ComplexApprox const& lhs, double const &rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (std::complex<double>& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = Approx(lhs.real()) == rhs.m_value.real();
        if (!rhs.m_compare_real_only)
        {
          is_equal &= Approx(lhs.imag()) == rhs.m_value.imag();
        }
        return is_equal;
    }

    friend bool operator == (std::complex<float>& lhs, ComplexApprox const& rhs)
    {
        bool is_equal = Approx(lhs.real()) == rhs.m_value.real();
        if (!rhs.m_compare_real_only)
        {
          is_equal &= Approx(lhs.imag()) == rhs.m_value.imag();
        }
        return is_equal;
    }

    friend bool operator == (ComplexApprox const &lhs, std::complex<double>& rhs)
    {
        return operator==( rhs, lhs );
    }

    friend bool operator == (ComplexApprox const &lhs, std::complex<float>& rhs)
    {
        return operator==( rhs, lhs );
    }

    ComplexApprox &compare_real_only()
    {
      m_compare_real_only = true;
      return *this;
    }


    std::string toString() const {
        std::ostringstream oss;
        oss <<"ComplexApprox( " << m_value << " )";
        return oss.str();
    }

};

template<>
inline std::string toString<ComplexApprox>( ComplexApprox const& value ) {
    return value.toString();
}
}

using Catch::ComplexApprox;

#endif

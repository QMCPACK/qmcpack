#ifndef CATCH_COMPLEX_APPROX
#define CATCH_COMPLEX_APPROX

#include <complex>
#include <limits>

// Copy and modify the Approx class to handle complex numbers

namespace Catch {
namespace Detail {
class ComplexApprox
{
    std::complex<double> m_value;
    double m_epsilon;

    bool approx_compare_zero(const double dist, const double magnitude) const
    {
      /* need to workaround catch to provide consistent behavior for a real value.
        bool Approx::equalityComparisonImpl(const double other) const {
        return marginComparison(m_value, other, m_margin) || marginComparison(m_value, other, m_epsilon * (m_scale + std::fabs(m_value)));
        }
      */
      return Approx(dist).epsilon(0.0).margin(m_epsilon * (1.0 + magnitude)) == 0.0;
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
        return rhs.approx_compare_zero(std::abs(lhs - rhs.m_value), std::abs(rhs.m_value));
    }

    friend bool operator == (std::complex<float> const& lhs, ComplexApprox const& rhs)
    {
        return rhs.approx_compare_zero(std::abs(std::complex<double>(lhs) - rhs.m_value), std::abs(rhs.m_value));
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
}
  
template<>
struct StringMaker<Catch::Detail::ComplexApprox> {
  static std::string convert(Catch::Detail::ComplexApprox const &value);
};

#ifdef CATCH_IMPL
std::string StringMaker<Catch::Detail::ComplexApprox>::convert(Catch::Detail::ComplexApprox const& value)
{
  return value.toString();
}
#endif

}

using Catch::Detail::ComplexApprox;

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/mantissa_exponent.h
///
/// \brief   header file for a class to store numbers in the mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_UTILS_MANTISSA_EXPONENT_HEADER
#define FORMIC_UTILS_MANTISSA_EXPONENT_HEADER

#include<formic/utils/exception.h>
#include<formic/utils/zero_one.h>

namespace formic {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   A class to store numbers in the format \f$ x = y * 10^z \f$ in order to avoid
  ///          problems from very large numbers and very small numbers that have exponents
  ///          outside the range that can be represented by a double.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class T> class MantExp {

    private:

      T _y; ///< the mantissa of the number, \f$ 1 \le |y| < 10 \f$
      int _z; ///< the exponent of the number

    public:

      T mantissa() const { return _y; } ///< get the mantissa
      int exponent() const { return _z; } ///< get the exponent

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   default constructor creates the number zero
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp() : _y(formic::zero(T())), _z(0) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   constructor that converts a regular number into mantissa/exponent format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      explicit MantExp(const T x) {
        const double absolute_value = std::abs(x);
        if (absolute_value == 0.0) {
          _y = formic::zero(T());
          _z = 0;
        } else {
          const T phase = x / absolute_value;
          const double logarithm = std::log10(absolute_value);
          _z = int(logarithm);
          _y = phase * std::pow(double(10.0), logarithm - _z);
          while(std::abs(_y) >= 10.0) {
            _y /= 10.0;
            _z += 1;
          }
          while(std::abs(_y) < 1.0) {
            _y *= 10.0;
            _z -= 1;
          }
        }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   constructor that creates the number \f$ x * 10^n \f$
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp(const T x, const int n) {
        *this = MantExp<T>(x);
        if (std::abs(_y) > 0.0)
          _z += n;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   function to return the stored number in standard format
      ///
      /// \return the stored number in standard format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      T get() const {
        if (_z < -100)
          return formic::zero(T());
        if (_z > 100)
          throw formic::Exception("large exponent (%i) detected in formic::MantExp::get") % _z;
        return _y * std::pow(double(10.0), _z);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to multiply the stored number by a regular number
      ///
      /// \return reference to the stored number after performing the multiplication
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> & operator*=(const T x) {
        *this = MantExp<T>(_y * x, _z);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to divide the stored number by a regular number
      ///
      /// \return reference to the stored number after performing the division
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> & operator/=(const T x) {
        *this = MantExp<T>(_y / x, _z);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to multiply the stored number by another mantissa/exponent number
      ///
      /// \return reference to the stored number after performing the multiplication
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> & operator*=(const MantExp<T> & other) {
        *this = MantExp<T>(_y * other._y, _z + other._z);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to divide the stored number by another mantissa/exponent number
      ///
      /// \return reference to the stored number after performing the division
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> & operator/=(const MantExp<T> & other) {
        *this = MantExp<T>(_y / other._y, _z - other._z);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to add a mantissa/exponent number to the stored number
      ///
      /// \return reference to the stored number after performing the addition
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> & operator+=(const MantExp<T> & other) {
        *this = MantExp<T>(_y + MantExp<T>(other._y, other._z - _z).get(), _z);
        return *this;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to multiply two numbers in mantissa/exponent format
      ///
      /// \return the result of the multiplication in mantissa/exponent format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> operator*(const MantExp<T> & other) const { return MantExp<T>(_y * other._y, _z + other._z); }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to compute the quotient of two numbers in mantissa/exponent format
      ///
      /// \return the result of the division in mantissa/exponent format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> operator/(const MantExp<T> & other) const { return MantExp<T>(_y / other._y, _z - other._z); }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to add two numbers in mantissa/exponent format
      ///
      /// \return the result of the addition in mantissa/exponent format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> operator+(const MantExp<T> & other) const {
        if ( _z > other._z ) return (*this) * ( formic::unity(T()) + ( (other) / (*this) ).get() );
                             return (other) * ( formic::unity(T()) + ( (*this) / (other) ).get() );
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   operator to subtract two numbers in mantissa/exponent format
      ///
      /// \return the result of the subtraction in mantissa/exponent format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      MantExp<T> operator-(const MantExp<T> & other) const {
        if ( _z > other._z ) return (*this) * (  formic::unity(T()) - ( (other) / (*this) ).get() );
                             return (other) * ( -formic::unity(T()) + ( (*this) / (other) ).get() );
      }

  };

} // end namespace formic

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   operator to compute the product of a mantissa/exponent number and a regular number
///
/// \return the result of the multiplication in mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> inline formic::MantExp<T> operator*(const formic::MantExp<T> & m, const T x) {
  return formic::MantExp<T>(m.mantissa() * x, m.exponent());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   operator to compute the product of a regular number and a mantissa/exponent number
///
/// \return the result of the multiplication in mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> inline formic::MantExp<T> operator*(const T x, const formic::MantExp<T> & m) {
  return formic::MantExp<T>(m.mantissa() * x, m.exponent());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   operator to compute the division of a mantissa/exponent number by a regular number
///
/// \return the result of the division in mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> inline formic::MantExp<T> operator/(const formic::MantExp<T> & m, const T x) {
  return formic::MantExp<T>(m.mantissa() / x, m.exponent());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   operator to compute the division of a regular number by a mantissa/exponent number
///
/// \return the result of the division in mantissa/exponent format
///
///////////////////////////////////////////////////////////////////////////////////////////////////
template<class T> inline formic::MantExp<T> operator/(const T x, const formic::MantExp<T> & m) {
  return formic::MantExp<T>(x / m.mantissa(), -m.exponent());
}

namespace formic {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   computes the square root of a mantissa/exponent number
  ///
  /// \return  the square root in mantissa/exponent format
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template<class T> inline formic::MantExp<T> sqrt(const formic::MantExp<T> & m) {
    if ( m.exponent() % 2 == 0 )
      return formic::MantExp<T>(std::sqrt(m.mantissa()), m.exponent() / 2 );
    return formic::MantExp<T>(std::sqrt(10.0 * m.mantissa()), ( m.exponent() - 1 ) / 2 );
  }

}

#endif

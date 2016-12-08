///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/random.h
///
/// \brief   header file for functions and objects associated with random numbers
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_RANDOM_HEADER
#define FORMIC_RANDOM_HEADER

#include <cmath>
#include <typeinfo>
#include <complex>

#include <formic/utils/exception.h>

namespace formic {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   The integer type used in random number generators.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  typedef unsigned int lcg_int_t;

  void set_seed(unsigned int seed);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Class for a linear congruential random number generator.
  ///
  /// Next number x = ( a * x + c ) mod m
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  class LinearCongruentialGenerator {

    private:

      const lcg_int_t _m; ///< integer used to define generator
      const lcg_int_t _a; ///< integer used to define generator
      const lcg_int_t _c; ///< integer used to define generator
      lcg_int_t _x;       ///< the current value of the random number

    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   constructor
      ///
      /// \param[in]     m        value for the modulus parameter
      /// \param[in]     a        value for the multiplicative parameter
      /// \param[in]     c        value for the addative parameter
      /// \param[in]     seed     the random seed
      ///
      /// The use of m = 2^31 - 1, a = 16807, c = 0 is recommended by S. K. Park and K. W. Miller,
      /// ``Random Number Generators: Good Ones are Hard to Find,'' Transactions of the ACM, Nov. 1988.
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      LinearCongruentialGenerator(const lcg_int_t m = 2147483647, //pcps::pow<pcps::lcg_int_t>(2,31)-1,
                                  const lcg_int_t a = 16807,
                                  const lcg_int_t c = 0,
                                  const lcg_int_t seed = 1)
        :  _m(m),
           _a(a),
           _c(c),
           _x(seed)
      {
        if (sizeof(lcg_int_t) != 4)
          throw formic::Exception("formic::LinearCongruentialGenerator requires lcg_int_t to be 32-bit");
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns the maximum possible value of the generator
      ///
      /// \return the maximum possible value of the generator
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t max() const { return _m - 1; }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns the minimum possible value of the generator
      ///
      /// \return the minimum possible value of the generator
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t min() const { return 0; }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns a random number and increments the generator's state
      ///
      /// \return a random number
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t operator()() {
        const lcg_int_t retval = ( _a * _x + _c ) % _m;
        if (retval == _x)
          throw formic::Exception("formic::LinearCongruentialGenerator got stuck.  Try a different seed.");
        _x = retval;
        return retval;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   sets the random seed
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      void set_seed(const lcg_int_t seed) { _x = seed; }

  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   A global linear congruential random number generator.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  extern formic::LinearCongruentialGenerator global_lcg;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Class for a lagged fibonacci random number generator.
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <lcg_int_t p = 1279, lcg_int_t q = 418> class LaggedFibonacci {

    private:

      lcg_int_t _word_length;        ///< range of the possible random numbers
      lcg_int_t _data_length;        ///< length of the array of previous values
      lcg_int_t _current_element;    ///< index of the current element
      std::vector<lcg_int_t> _data;  ///< the array of previous values

    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   constructs the generator
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      LaggedFibonacci()
        : _word_length(1073741824), // 2^30 // (pcps::pow<lcg_int_t>(2,30)),
          _data_length(p+1),
          _current_element(p),
          _data(p+1, 0)
      {

        // check that values are sane
        if (sizeof(lcg_int_t) != 4)
          throw formic::Exception("formic::LaggedFibonacci requires lcg_int_t to be 32-bit");
        if (q >= p)
          throw formic::Exception("formic::LaggedFibonacci requires p > q");
        if (p <= 0)
          throw formic::Exception("formic::LaggedFibonacci requires p > 0");
        if (q <= 0)
          throw formic::Exception("formic::LaggedFibonacci requires q > 0");

        // initialize the data vector
        bool have_odd = false;
        for (lcg_int_t i = 0; i < p; i++) {
          _data[i] = (lcg_int_t)( 0.1 + std::fabs(   double(formic::global_lcg())
                                                   * double(_word_length) / double(formic::global_lcg.max()) ) ) % _word_length;
          if (_data[i] % 2 == 1)
            have_odd = true;
        }
        if (!have_odd) {
          const lcg_int_t to_be_odd = formic::global_lcg() % p;
          if (_data[to_be_odd] == 0)
            _data[to_be_odd] += 1;
          else
            _data[to_be_odd] -= 1;
        }

      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns the maximum possible value of the generator
      ///
      /// \return the maximum possible value of the generator
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t max() const { return _word_length - 1; }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns the minimum possible value of the generator
      ///
      /// \return the minimum possible value of the generator
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t min() const { return 0; }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   returns a random number and increments the generator's state
      ///
      /// \return a random number
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      lcg_int_t operator()() {
        const lcg_int_t retval = (   _data[ ( (_current_element + _data_length) - p) % _data_length ]
                                   + _data[ ( (_current_element + _data_length) - q) % _data_length ] ) % _word_length;
        _data[_current_element] = retval;
        _current_element++;
        if (_current_element == _data_length)
          _current_element = 0;
        return retval;
      }

  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns a random number in the appropriate type
  ///
  /// \return a random number in the appropriate type
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // no specialization
  template <class S> inline S random_number() {
    throw formic::Exception("no specialization of formic::random_number() for type \"%s\"") % typeid(S).name();
  }
  // specialization for double
  template <> inline double random_number<double>() {
    return double(formic::global_lcg()) / double(formic::global_lcg.max());
  }
  // specialization for complex
  template <> inline std::complex<double> random_number< std::complex<double> >() {
    return std::complex<double>(formic::random_number<double>(), formic::random_number<double>());
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   returns a random number chosen uniformly from the interval 0 to 1
  ///
  /// \return a random number chosen uniformly from the interval 0 to 1
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class GENERATOR> inline double uni_01(GENERATOR & gen) {
    const double retval = std::fabs( double(gen()) / double(gen.max()) );
    assert( retval >= 0.0 && retval <= 1.0 );
    return retval;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   Given a range containing a discreet probability distribution, whose elements sum to 1,
  ///          pick at random a position in the range according to the probability distribution.
  ///
  /// \param[in]      start    beginning of the probability distribution range
  /// \param[in]      end      end of the probability distribution range
  /// \param[in,out]  gen      the random number generator to use
  ///
  /// \return  the randomly chosen position in the range [0, 1, 2, ..., end-start)
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  template <class ITERATOR, class GENERATOR> inline int rand_pos_in_prob_vec(ITERATOR start, ITERATOR end, GENERATOR & gen) {

    // if the range is empty return 0
    if ( start == end )
      return 0;

    // pick a random number between 0 and 1
    double rn = formic::uni_01(gen);

    // loop over the range until the cumulative probability exceeds the randomly drawn number
    // or the range is exhausted (in which case the last element with nonzero probability is selected)
    int retval = -1;
    do {
      rn -= *start;
      retval++;
      start++;
    } while ( rn >= 0.0 && start != end );
    if ( rn >= 1.0e-8 )
      throw formic::Exception("unexpectedly high residual probability in formic::rand_pos_in_prob_vec");
    start--;
    while ( retval >= 0 ) {
      if ( *start != 0.0 )
       break;
      retval--;
      start--;
    }
    if ( retval < 0 )
      throw formic::Exception("failed to choose an element in formic::rand_pos_in_prob_vec");

    // return the randomly chosen position
    return retval;

  }

} // end namespace formic

#endif

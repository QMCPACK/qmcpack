///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/utils/exception.h
///
/// \brief   defines a general exception class
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_EXCEPTION_HEADER
#define FORMIC_EXCEPTION_HEADER

#include <string>
#include <stdexcept>

#include <boost/format.hpp>

namespace formic {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   a general exception class
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  class Exception : public std::runtime_error {

    private:

      boost::format _f; ///< object used for formatting the error message

    public:

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   constructs the exception from a string
      ///
      /// \param[in]      s    the error message, which may contain formatting used by boost::format
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      Exception(const std::string & s = "unspecified formic exception") : std::runtime_error(s), _f(s) {}

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   replaces one of the formatting specifiers with the supplied object
      ///
      /// \param[in]      arg  the object to use for the first of the remaining formatting specifiers
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      template<class T> formic::Exception & operator%(const T & arg) {

        // feed the argument to the formatting object
        _f % arg;

        // try to complete the formatting
        bool complete = false;
        std::string formatted_message;
        try {
          formatted_message = _f.str(); // boost will throw an error here if the formatting is incomplete
          complete = true;
        } catch (...) {
        }

        // if we switch to boost 1.46.1, we can replace this try/catch statement with a test on _f.remaining_args()

        // if the formatting is complete, reset the exception using the formatted error message
        if (complete) *this = formic::Exception(formatted_message);

        // return a reference to the exception
        return *this;

      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      /// \brief   destructor, which is declared with throw() in order to match the virtual destructor
      ///          of the parent class std::runtime_error
      ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////
      ~Exception() throw() {}

  };

}

#endif

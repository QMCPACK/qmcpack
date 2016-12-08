///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/timing/timing.h
///
/// \brief   declarations of functions and classes used in timing the code
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef FORMIC_TIMING_HEADER
#define FORMIC_TIMING_HEADER

#include <string>

#include <boost/scoped_ptr.hpp>

#include <formic/utils/exception.h>

namespace cqmc {

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   A class for keeping track of elapsed time
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  class Stopwatch {
    public:
      Stopwatch();
      void reset(const std::string & name = std::string(""));
      void start(const std::string & name = std::string(""));
      void stop(const std::string & name = std::string(""));
      double elapsed_seconds() const;
    private:
      bool _running; ///< whether the stopwatch is running
      double _start_time; ///< when the clock was started
      double _elapsed_time; ///< elapsed time after most recent stop
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief   A class for a timer that starts on construction and ends on destruction
  ///
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  class ScopedTimer {
    public:
      ScopedTimer(std::pair<std::string,cqmc::Stopwatch> * pp) : _pp(pp) {
        if ( !_pp )
          throw formic::Exception("illegal construction of an empty ScopedTimer");
        _pp->second.start(_pp->first); // start the stopwatch upon construction
      }
      void clear() {
        if ( !_pp )
          throw formic::Exception("empty ScopedTimer cannot be cleared");
        _pp->second.stop(_pp->first); // stop the stopwatch
        _pp = 0; // empty the object
      }
      ~ScopedTimer() {
        if ( _pp )
          this->clear(); // if the timer has not already been stopped, stop it upon destruction
      }
    private:
      std::pair<std::string,cqmc::Stopwatch> * _pp; ///< pointer to name and underlying stopwatch
      ScopedTimer(const ScopedTimer &); // disable copy construction
      ScopedTimer & operator=(const ScopedTimer &); // disable assignment
  };

  std::string print_timers();
  void start_timer(const std::string & name);
  void stop_timer(const std::string & name);
  void reset_timer(const std::string & name);

  namespace timers {
    extern boost::scoped_ptr<std::pair<std::string,cqmc::Stopwatch> > propose_1e;
    extern boost::scoped_ptr<std::pair<std::string,cqmc::Stopwatch> > propose_2e;
  }

} // end namespace cqmc

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file formic/timing/timing.cpp
///
/// \brief   implementation of functions and classes used in timing the code
///
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <list>

#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <formic/utils/timing.h>
#include <formic/utils/mpi_interface.h>
#include <formic/utils/exception.h>

boost::scoped_ptr<std::pair<std::string,formic::Stopwatch> > formic::timers::propose_1e( new std::pair<std::string,formic::Stopwatch>(std::string("propose_1e"), formic::Stopwatch()) );
boost::scoped_ptr<std::pair<std::string,formic::Stopwatch> > formic::timers::propose_2e( new std::pair<std::string,formic::Stopwatch>(std::string("propose_2e"), formic::Stopwatch()) );

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   constructs a new stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
formic::Stopwatch::Stopwatch() : _running(false), _start_time(0), _elapsed_time(0) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   resets the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::Stopwatch::reset(const std::string & name) {
  if (_running)
    throw formic::Exception("cannot reset Stopwatch \"%s\" because it is running") % name;
  _elapsed_time = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   starts the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::Stopwatch::start(const std::string & name) {
  if (_running)
    throw formic::Exception("cannot start Stopwatch \"%s\" because it is already running") % name;
  _start_time = MPI_Wtime();
  _running = true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   stops the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::Stopwatch::stop(const std::string & name) {
  if (!_running)
    throw formic::Exception("cannot stop Stopwatch \"%s\" because it is not running") % name;
  _elapsed_time = _elapsed_time + ( MPI_Wtime() - _start_time );
  _running = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   returns the elapsed time in seconds
///
///////////////////////////////////////////////////////////////////////////////////////////////////
double formic::Stopwatch::elapsed_seconds() const {
  double total = _elapsed_time;
  if (_running)
    total = total + ( MPI_Wtime() - _start_time );
  return total;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   container to hold timers in use by the formic code
///
///////////////////////////////////////////////////////////////////////////////////////////////////
static std::map<std::string, formic::Stopwatch> formic_timer_container;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   prints all timers in formic_timer_container
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string formic::print_timers() {

  struct local_funcs {

    static void print_timer(std::stringstream & ss,
                            const formic::Stopwatch & sw,
                            const std::string & name) {

      // do not print the timer if it recorded no time
      if ( sw.elapsed_seconds() == 0.0 )
        return;

      // print the timer
      ss << boost::format("  %-30s      %15.6f seconds") % name % sw.elapsed_seconds() << std::endl;

    }

  };

  std::stringstream ss;

  ss << std::endl;
  ss << boost::format("Printing all timers used during this execution:") << std::endl;

  local_funcs::print_timer(ss, formic::timers::propose_1e->second, formic::timers::propose_1e->first);
  local_funcs::print_timer(ss, formic::timers::propose_2e->second, formic::timers::propose_2e->first);

  for (std::map<std::string, formic::Stopwatch>::const_iterator it = formic_timer_container.begin();
       it != formic_timer_container.end(); it++)
    local_funcs::print_timer(ss, it->second, it->first);
    //ss << boost::format("  %-30s      %15.6f seconds") % it->first % it->second.elapsed_seconds() << std::endl;

  ss << std::endl;

  return ss.str();

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   starts a timer
///
/// \param[in]     name     a string identifying the timer
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::start_timer(const std::string & name) {
  formic_timer_container[name].start(name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   stops a timer
///
/// \param[in]     name     a string identifying the timer
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::stop_timer(const std::string & name) {
  formic_timer_container[name].stop(name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   resets a timer
///
/// \param[in]     name     a string identifying the timer
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void formic::reset_timer(const std::string & name) {
  formic_timer_container[name].reset(name);
}


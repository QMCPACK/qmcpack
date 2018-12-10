///////////////////////////////////////////////////////////////////////////////////////////////////
/// \file cqmc/timing/timing.cpp
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

#include<formic/utils/mpi_interface.h>

#include<formic/utils/lmyengine/engine_timing.h>
#include<formic/utils/exception.h>

boost::scoped_ptr<std::pair<std::string,cqmc::Stopwatch> > cqmc::timers::propose_1e( new std::pair<std::string,cqmc::Stopwatch>(std::string("propose_1e"), cqmc::Stopwatch()) );
boost::scoped_ptr<std::pair<std::string,cqmc::Stopwatch> > cqmc::timers::propose_2e( new std::pair<std::string,cqmc::Stopwatch>(std::string("propose_2e"), cqmc::Stopwatch()) );

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   constructs a new stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
cqmc::Stopwatch::Stopwatch() : _running(false), _start_time(0), _elapsed_time(0) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   resets the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::Stopwatch::reset(const std::string & name) {
  if (_running)
    throw formic::Exception("cannot reset Stopwatch \"%s\" because it is running") % name;
  _elapsed_time = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   starts the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::Stopwatch::start(const std::string & name) {
  if (_running)
    throw formic::Exception("cannot start Stopwatch \"%s\" because it is already running") % name;
  _start_time = MPI_Wtime();
  _running = true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   stops the stopwatch
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::Stopwatch::stop(const std::string & name) {
  if (!_running)
    throw formic::Exception("cannot stop Stopwatch \"%s\" because it is not running") % name;
  _elapsed_time = _elapsed_time + ( MPI_Wtime() - _start_time );
  _running = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   returns the elapsed time in seconds
///
///////////////////////////////////////////////////////////////////////////////////////////////////
double cqmc::Stopwatch::elapsed_seconds() const {
  double total = _elapsed_time;
  if (_running)
    total = total + ( MPI_Wtime() - _start_time );
  return total;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   container to hold timers in use by the cqmc code
///
///////////////////////////////////////////////////////////////////////////////////////////////////
static std::map<std::string, cqmc::Stopwatch> cqmc_timer_container;

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   prints all timers in cqmc_timer_container
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string cqmc::print_timers() {

  struct local_funcs {

    static void print_timer(std::stringstream & ss,
                            const cqmc::Stopwatch & sw,
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

  local_funcs::print_timer(ss, cqmc::timers::propose_1e->second, cqmc::timers::propose_1e->first);
  local_funcs::print_timer(ss, cqmc::timers::propose_2e->second, cqmc::timers::propose_2e->first);

  for (std::map<std::string, cqmc::Stopwatch>::const_iterator it = cqmc_timer_container.begin();
       it != cqmc_timer_container.end(); it++)
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
void cqmc::start_timer(const std::string & name) {
  cqmc_timer_container[name].start(name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   stops a timer
///
/// \param[in]     name     a string identifying the timer
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::stop_timer(const std::string & name) {
  cqmc_timer_container[name].stop(name);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief   resets a timer
///
/// \param[in]     name     a string identifying the timer
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void cqmc::reset_timer(const std::string & name) {
  cqmc_timer_container[name].reset(name);
}


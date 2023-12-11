//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef MPIEXECPTIONWRAPPER_H
#define MPIEXECPTIONWRAPPER_H

#include "Message/Communicate.h"

namespace qmcplusplus
{

class MPIExceptionWrapper
{
public:
  MPIExceptionWrapper() {};

  /** Call an arbitrary function that can throw an exception catch it, extract message,
   *  and call APP_ABORT
   *
   */
  template<typename F, typename... Args>
  void operator()(F&& f, Args&&... args);

};

template<typename F, typename... Args>
void MPIExceptionWrapper::operator()(F&& f, Args&&... args)
{
  try
  {
    f(std::forward<Args>(args)...);
  }
  catch (const std::exception& re)
  {
    APP_ABORT(re.what());
  }
  catch (...)
  {
    APP_ABORT("MPIExceptionWrapper caught an unknown exception!");
  }

}

}


#endif /* MPIEXECPTIONWRAPPER_H */

////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by:
// Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by:
// Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_TASKSOPENMP_HPP
#define QMCPLUSPLUS_TASKSOPENMP_HPP

#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{


template<>
template<typename F, typename... Args>
void TasksOneToOne<Threading::OPENMP>::operator()(F&& f, Args&&... args)
{
#pragma omp parallel num_threads(num_threads_)
  {
      f(omp_get_thread_num(), std::forward<Args>(args)...);
  }
}

}

#endif

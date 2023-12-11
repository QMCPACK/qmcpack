//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_LISTENER_HPP
#define QMCPLUSPLUS_LISTENER_HPP
/** \file
 *  Listener types that allow Estimators to register with QMCHamiltonian to have "trace" values from
 *  operators reported.
 *  This is aims to be lightweight and minimal but still use much of the trace manager based
 *  implementation of observables that rely on per particle Hamiltonian values.
 */

#include <functional>

#include "type_traits/template_types.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

/** An object of this type is a listener expecting a callback to the report function  with a vector of values,
 *  the convention being 1 per particle, called once per walker per component.
 *  The name is primarily for debugging purposes, if have decided against having the QMCHamiltonian use it to estable
 *  routing. Instead the register functions in the QMCHamiltonian are specfic to the sort of value that the listener wants
 *  to listen to.
 */
template<typename T>
class ListenerVector
{
public:
  /** "Callback" function type for an operator to report a vector of values to a listener
   *  \param[in] walker_index    a numeric walker id, could be important to listener
   *  \param[in] name            of operator reporting, could be important to listener
   *  \param[in] values          vector of values, per particle by convention. Also by convention
   *                             the receiver should not assume the reference values have any persistence
   *                             beyond the scope of the callback.
   */
  using ReportingFunction =
      std::function<void(const int walker_index, const std::string& name, const Vector<T>& values)>;
  /** constructor that should appear in the code
   *  \param[in]  name           intended to be strictly information
   *  \param[in]  report_func    intended to be a lambda callback function that will be used to report.
   */
  ListenerVector(const std::string& name, ReportingFunction report_func) : report(report_func), name_(name) {}
  ReportingFunction report;
  const std::string& get_name() const { return name_; }

private:
  const std::string name_;
};

/** Convenience container for common optional element to mw_eval.._impl.
 *  This allows the per_particle and reduced mw_eval_... to share the same
 *  implementation method.
 *
 *  member naming is such that on usage:
 *       ListenerOption listeners
 *       ...
 *       if (listeners)
 *         for (const ListenerVector<Real>& listener : listeners->electron_values)
 *           listener.report(iw, O_leader.getName(), ve_sample);
 *
 *  see NonLocalECPotential for example of usage.
 */
template<typename T>
struct ListenerOption
{
  ListenerOption(const std::vector<ListenerVector<T>>& le, const std::vector<ListenerVector<T>>& li)
      : electron_values(le), ion_values(li)
  {}
  const std::vector<ListenerVector<T>>& electron_values;
  const std::vector<ListenerVector<T>>& ion_values;
};

} // namespace qmcplusplus


#endif

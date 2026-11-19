//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_DRIFTMODIFIER_BUILDER_H
#define QMCPLUSPLUS_DRIFTMODIFIER_BUILDER_H

/** @file
 *  @brief Factory functions for constructing DriftModifierBase instances.
 */

#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBase.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
/** create a DriftModifierBase instance from an XML input node (used by legacy drivers)
 * @param cur XML node pointer to the QMC driver section
 * @param myComm MPI communicator
 * @return pointer to a newly allocated DriftModifierBase object
 */
DriftModifierBase* createDriftModifier(xmlNodePtr cur, const Communicate& myComm);

/** create a DriftModifierBase instance from parsed input parameters (used by batched drivers)
 * @param drift_modifier_str name of the drift modifier scheme (e.g., "UNR", case-insensitive)
 * @param unr_a scaling parameter 'a' for the UNR drift modifier
 * @return pointer to a newly allocated DriftModifierBase object
 */
DriftModifierBase* createDriftModifier(const std::string& drift_modifier_str, QMCTraits::RealType unr_a);
} // namespace qmcplusplus

#endif
